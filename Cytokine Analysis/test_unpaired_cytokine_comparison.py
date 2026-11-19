import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


class UnpairedCytokineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        notebook_path = Path(__file__).with_name("260527_Cytokine Analysis.ipynb")
        notebook = json.loads(notebook_path.read_text())
        sources = [
            "".join(cell["source"])
            for cell in notebook["cells"]
            if "".join(cell["source"]).startswith("def Plot_Unpaired_Cytokine_Comparison(")
        ]
        if len(sources) != 1:
            raise AssertionError("Expected one unpaired comparison function.")
        namespace = {}
        exec(compile(sources[0], str(notebook_path), "exec"), namespace)
        cls.compare = staticmethod(namespace["Plot_Unpaired_Cytokine_Comparison"])

    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.source = Path(self.directory.name) / "input.xlsx"

    def run_comparison(self, metrics, *, repeats=False, save=False, exclude=None,
                       sheet="Data", end_row=None, output_dir=None):
        count = len(next(iter(metrics.values()))[0])
        raw = pd.DataFrame(
            np.full((1 + count * 2 * (2 if repeats else 1), 6 + 3 * len(metrics)), None)
        )
        for index, name in enumerate(metrics):
            raw.iat[0, 5 + 3 * index] = name
        row = 1
        for group_index, group in enumerate(("6hr Ctrl", "48hr 5ug")):
            for rep in range(count):
                for offset in ((-0.5, 0.5) if repeats else (0,)):
                    raw.iat[row, 0] = f"Rep {rep + 1}"
                    raw.iat[row, 1] = group
                    for index, values in enumerate(metrics.values()):
                        raw.iat[row, 6 + 3 * index] = values[group_index][rep] + offset
                    row += 1
        raw.to_excel(self.source, sheet_name=sheet, header=False, index=False)
        with contextlib.redirect_stdout(io.StringIO()):
            return self.compare(
                self.source, sheet_name=sheet, header_row=1, data_start_row=2,
                data_end_row=end_row, output_dir=output_dir,
                expected_replicates=count, save_results=save, display_results=False,
                exclude_replicates=exclude,
            )

    def test_excluding_one_treatment_replicate_preserves_controls_and_audit(self):
        points, results, paths = self.run_comparison(
            {"IL-6": ([1, 2, 3], [4, 1000, 6])}, repeats=True, save=True,
            exclude=[("48hr 5ug", "Rep 2")],
        )
        self.assertEqual(set(points.loc[points.Sample.eq("6hr Ctrl"), "Rep"]),
                         {"Rep 1", "Rep 2", "Rep 3"})
        self.assertEqual(set(points.loc[points.Sample.eq("48hr 5ug"), "Rep"]),
                         {"Rep 1", "Rep 3"})
        self.assertEqual(results.loc[0, "N_Control"], 3)
        self.assertEqual(results.loc[0, "N_Treatment"], 2)
        self.assertEqual(results.loc[0, "Mean_Treatment"], 5)
        self.assertAlmostEqual(results.loc[0, "P_Value"], 0.2)
        self.assertAlmostEqual(results.loc[0, "Q_BH"], 0.2)
        self.assertEqual(results.loc[0, "Status"], "tested")
        self.assertEqual(Path(paths[-1]).parent.name,
                         "Unpaired_Cytokine_Comparison_Excluded_Replicates__Data__rows2-13")
        with pd.ExcelFile(paths[-1]) as book:
            excluded = pd.read_excel(book, sheet_name="Excluded_Observations")
            self.assertEqual(len(excluded), 2)
            self.assertTrue(excluded.Sample.eq("48hr 5ug").all())
            self.assertTrue(excluded.Rep.eq("Rep 2").all())
            retained = pd.read_excel(book, sheet_name="Input_Observations")
            self.assertEqual(len(retained), 10)
            notes = pd.read_excel(book, sheet_name="Analysis_Notes").set_index("Setting")
            self.assertEqual(notes.loc["Excluded replicates", "Value"], "48hr 5ug: Rep 2")

    def test_exclusion_does_not_relax_missing_retained_value_requirement(self):
        _, results, _ = self.run_comparison(
            {"IL-6": ([1, 2, 3], [4, 1000, np.nan])},
            exclude=[("48hr 5ug", "Rep 2")],
        )
        self.assertTrue(results.loc[0, "Status"].startswith("not tested:"))
        self.assertTrue(np.isnan(results.loc[0, "P_Value"]))
        self.assertEqual(results.loc[0, "Missing_Treatment_IDs"], "Rep 3")

    def test_invalid_exclusions_raise(self):
        for exclusions in (
            [("48hr 5ug", "Rep 9")],
            [("48hr 5ug", "Rep 1"), ("48hr 5ug", "Rep 2")],
            [("48hr 5ug", "Rep 2"), ("48hr 5ug", "Rep 2")],
            [("unknown", "Rep 2")],
        ):
            with self.subTest(exclusions=exclusions), self.assertRaises(ValueError):
                self.run_comparison(
                    {"IL-6": ([1, 2, 3], [4, 5, 6])}, exclude=exclusions,
                )

    def test_exact_ties_constants_missing_and_fdr(self):
        _, results, paths = self.run_comparison({
            "Separated": ([-3, -2, -1], [0, 1, 2]),
            "Ties": ([1, 1, 2], [2, 3, 3]),
            "Identical": ([1, 1, 1], [1, 1, 1]),
            "ConstantSeparated": ([0, 0, 0], [1, 1, 1]),
            "Missing": ([np.nan, 1, 2], [3, 4, 5]),
        })
        table = results.set_index("Cytokine")
        np.testing.assert_allclose(table.P_Value.iloc[:4], [0.1, 0.2, 1, 0.1])
        np.testing.assert_allclose(table.U_Statistic.iloc[:4], [9, 8.5, 4.5, 9])
        np.testing.assert_allclose(
            table.Q_BH.iloc[:4], stats.false_discovery_control([0.1, 0.2, 1, 0.1])
        )
        self.assertTrue(table.FDR_Family_Size.eq(4).all())
        self.assertTrue(table.P_Method.iloc[:4].eq("exact permutation (ties included)").all())
        self.assertTrue(table.loc["Missing", "Status"].startswith("not tested:"))
        self.assertEqual(table.loc["Missing", "Missing_Control_IDs"], "Rep 1")
        self.assertTrue(np.isnan(table.loc["Missing", "P_Value"]))
        self.assertTrue(np.isnan(table.loc["Missing", "Q_BH"]))
        self.assertFalse({"T_Statistic", "DF", "CI95_Low", "CI95_High"} & set(results))
        self.assertEqual(paths, [])
        self.assertEqual(plt.get_fignums(), [])

    def test_technical_repeats_and_exports(self):
        points, results, paths = self.run_comparison(
            {"IL-6": ([-3, -2, -1], [0, 1, 2])}, repeats=True, save=True
        )
        self.assertEqual(len(points), 6)
        self.assertTrue(points.Number_Values.eq(2).all())
        self.assertTrue(points.Number_Rows.eq(2).all())
        np.testing.assert_allclose(points.Value, [-3, -2, -1, 0, 1, 2])
        self.assertEqual(results.loc[0, "N_Control"], 3)
        self.assertEqual(results.loc[0, "Mean_Difference"], 3)
        self.assertEqual(results.loc[0, "Median_Control"], -2)
        self.assertEqual(results.loc[0, "Median_Treatment"], 1)
        self.assertEqual(len(paths), 2)
        self.assertTrue(all(Path(path).stat().st_size > 0 for path in paths))
        self.assertTrue(all(
            Path(path).parent.name == "Unpaired_Cytokine_Comparison__Data__rows2-13"
            for path in paths
        ))
        with pd.ExcelFile(paths[-1]) as book:
            self.assertEqual(book.sheet_names, [
                "Unpaired_Statistics", "Biological_Replicates",
                "Input_Observations", "Analysis_Notes",
            ])
            exported = pd.read_excel(book, sheet_name="Unpaired_Statistics")
            self.assertAlmostEqual(exported.loc[0, "P_Value"], 0.1)
            self.assertEqual(exported.loc[0, "U_Statistic"], 9)
            notes = pd.read_excel(book, sheet_name="Analysis_Notes")
            self.assertIn("Mann-Whitney U", notes.set_index("Setting").loc["Test", "Value"])
            self.assertEqual(
                notes.set_index("Setting").loc["Scale", "Value"],
                "Values as supplied; no transformation, normalization, or imputation",
            )

    def test_explicit_bounds_and_sanitized_sheet_in_default_folder(self):
        _, _, paths = self.run_comparison(
            {"IL-6": ([1, 2, 3], [4, 5, 6])},
            save=True, sheet="261018 Updated Analysis", end_row=7,
        )
        self.assertTrue(all(
            Path(path).parent.name ==
            "Unpaired_Cytokine_Comparison__261018_Updated_Analysis__rows2-7"
            for path in paths
        ))

    def test_explicit_output_directory_is_unchanged(self):
        output = Path(self.directory.name) / "custom"
        _, _, paths = self.run_comparison(
            {"IL-6": ([1, 2, 3], [4, 5, 6])}, save=True, output_dir=output,
        )
        self.assertTrue(all(Path(path).parent == output.resolve() for path in paths))

    def test_raw_values_are_not_transformed(self):
        points, results, _ = self.run_comparison({
            "IL-6": ([1, 10, 100], [2, 20, 200]),
        })
        np.testing.assert_allclose(points.Value, [1, 10, 100, 2, 20, 200])
        self.assertEqual(results.loc[0, "Mean_Control"], 37)
        self.assertEqual(results.loc[0, "Mean_Treatment"], 74)

    def test_larger_samples_use_tie_corrected_asymptotic_test(self):
        control = [0] * 5 + [1] * 5
        treatment = [1] * 5 + [2] * 5
        _, results, _ = self.run_comparison({"IL-6": (control, treatment)})
        expected = stats.mannwhitneyu(
            treatment, control, alternative="two-sided",
            method="asymptotic", use_continuity=True,
        )
        self.assertAlmostEqual(results.loc[0, "P_Value"], expected.pvalue)
        self.assertEqual(results.loc[0, "U_Statistic"], expected.statistic)
        self.assertEqual(
            results.loc[0, "P_Method"], "asymptotic (tie and continuity corrected)"
        )


class StandaloneUnpairedCytokineTests(UnpairedCytokineTests):
    @classmethod
    def setUpClass(cls):
        from Cytokine_48hr_5ug_vs_Control import Plot_Unpaired_Cytokine_Comparison

        cls.compare = staticmethod(Plot_Unpaired_Cytokine_Comparison)

    def test_display_is_disabled_by_default(self):
        import inspect

        self.assertIs(
            inspect.signature(self.compare).parameters["display_results"].default,
            False,
        )

    def test_clean_copies_preserve_original_style_and_statistics(self):
        from functools import partial
        from unittest.mock import patch
        from matplotlib.collections import PathCollection
        from matplotlib.figure import Figure
        from matplotlib.colors import to_rgba

        _, baseline_statistics, baseline_paths = self.run_comparison(
            {"IL-6": ([1, 2, 3], [4, 5, 6])}, save=True,
        )
        original_bytes = Path(baseline_paths[0]).read_bytes()
        snapshots = []
        savefig = Figure.savefig

        def capture(figure, path, **kwargs):
            ax = figure.axes[0]
            snapshots.append((
                ax.get_legend() is not None,
                [collection.get_facecolors().copy() for collection in ax.collections
                 if isinstance(collection, PathCollection)],
            ))
            return savefig(figure, path, **kwargs)

        self.compare = partial(self.compare, save_clean_plots=True)
        with patch.object(Figure, "savefig", capture):
            _, statistics, paths = self.run_comparison(
                {"IL-6": ([1, 2, 3], [4, 5, 6])}, save=True,
            )
        pd.testing.assert_frame_equal(statistics, baseline_statistics)
        self.assertEqual(len(paths), 3)
        self.assertEqual(Path(paths[0]).read_bytes(), original_bytes)
        self.assertEqual(Path(paths[1]).parent.name, "No_Legend_Black_Dots")
        self.assertTrue(Path(paths[1]).stat().st_size > 0)
        self.assertTrue(snapshots[0][0])
        self.assertFalse(snapshots[1][0])
        self.assertEqual(len(snapshots[1][1]), 6)
        for facecolors in snapshots[1][1]:
            np.testing.assert_allclose(facecolors, [to_rgba("black")])

    def test_saving_without_display_does_not_create_gui_figures(self):
        from unittest.mock import patch

        with patch.object(plt, "subplots", side_effect=AssertionError("GUI figure created")), \
                patch.object(plt, "show", side_effect=AssertionError("Plot window opened")):
            _, _, paths = self.run_comparison(
                {"IL-6": ([1, 2, 3], [4, 5, 6])}, save=True,
            )
        self.assertEqual(len(paths), 2)
        self.assertTrue(all(Path(path).stat().st_size > 0 for path in paths))

    def test_significance_bracket_uses_adjusted_q_value_in_both_plot_styles(self):
        from functools import partial
        from unittest.mock import patch
        from matplotlib.figure import Figure

        self.compare = partial(self.compare, save_clean_plots=True)
        separated = {"IL-6": ([1, 2, 3, 4], [5, 6, 7, 8])}
        with patch.object(Figure, "savefig", autospec=True) as save:
            _, significant, _ = self.run_comparison(separated, save=True)
        self.assertLess(significant.iloc[0].Q_BH, 0.05)
        self.assertEqual(len(save.call_args_list), 2)
        for call in save.call_args_list:
            ax = call.args[0].axes[0]
            self.assertEqual([text.get_text() for text in ax.texts], ["*"])
            np.testing.assert_array_equal(ax.lines[-1].get_xdata(), [0, 0, 1, 1])
            self.assertLess(ax.texts[0].get_position()[1], ax.get_ylim()[1])

        with patch.object(Figure, "savefig", autospec=True) as save:
            _, adjusted, _ = self.run_comparison(
                dict(separated, Other=([1, 2, 3, 4], [1, 2, 3, 4])), save=True,
            )
        self.assertLess(adjusted.iloc[0].P_Value, 0.05)
        self.assertGreaterEqual(adjusted.iloc[0].Q_BH, 0.05)
        self.assertTrue(all(not call.args[0].axes[0].texts for call in save.call_args_list))


if __name__ == "__main__":
    unittest.main()
