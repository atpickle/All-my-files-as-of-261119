import ast
import contextlib
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from scipy import stats

from Cytokine_Heatmap_Statistics import analyze_heatmap_groups, annotate_heatmap_bar
from Cytokine_Heatmap_Bar_Plots import plot_cytokine_heatmap


class HeatmapBarTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = Path(__file__).with_name("260527_Cytokine Analysis.ipynb")
        notebook = json.loads(path.read_text())
        cls.functions = []
        for cell in notebook["cells"]:
            source = "".join(cell["source"])
            if source.startswith("def plot_cytokine_heatmap("):
                node = ast.parse(source).body[0]
                namespace = {}
                exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), "exec"), namespace)
                cls.functions.append(namespace["plot_cytokine_heatmap"])
        if len(cls.functions) != 2:
            raise AssertionError("Expected both heatmap definitions.")
        cls.functions.append(plot_cytokine_heatmap)

    def test_standalone_file_runs_without_notebook_or_helper_module(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            shutil.copyfile(
                Path(__file__).with_name("Cytokine_Heatmap_Bar_Plots.py"),
                folder / "Cytokine_Heatmap_Bar_Plots.py",
            )
            raw = pd.DataFrame([[None] * 7 for _ in range(6)])
            raw.iat[0, 5] = "IL-6"
            for row, (sample, rep, value) in enumerate([
                ("Control", "1", 1), ("Control", "2", 2),
                ("Treatment", "1", 3), ("Treatment", "2", 4), ("Treatment", "3", 5),
            ], start=1):
                raw.iat[row, 0], raw.iat[row, 1], raw.iat[row, 6] = rep, sample, value
            raw.to_excel(folder / "input.xlsx", header=False, index=False)
            code = """
from pathlib import Path
from Cytokine_Heatmap_Bar_Plots import plot_cytokine_heatmap
assert not Path("output").exists()
plot_cytokine_heatmap(
    "input.xlsx", header_row=1, samples=["Control", "Treatment"],
    control_sample="Control", scale_mode="zscore_row", output_dir="output",
)
assert len(list(Path("output").glob("*.png"))) == 1
assert len(list(Path("output/Raw Bar Plots").glob("*.png"))) == 1
assert len(list(Path("output/Raw Bar Plots").glob("*.xlsx"))) == 1
"""
            result = subprocess.run(
                [sys.executable, "-c", code], cwd=folder, capture_output=True,
                text=True, timeout=60, env={**os.environ, "MPLBACKEND": "Agg"},
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_raw_bars_match_filtered_means_in_every_scaling_mode(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            source = folder / "input.xlsx"
            raw = pd.DataFrame([[None] * 13 for _ in range(7)])
            for col, name in [(5, "IL-6"), (8, "Flat"), (11, "Excluded")]:
                raw.iat[0, col] = name
            for row, (rep, sample, value) in enumerate([
                ("Rep 1", "A", 2),
                ("Rep 1", "A", 4),
                ("Rep 2", "A", 6),
                ("Rep 1", "B", 8),
                ("Rep 2", "B", None),
                ("Rep 1", "Excluded sample", 1000),
            ], start=1):
                raw.iat[row, 0], raw.iat[row, 1] = rep, sample
                raw.iat[row, 6], raw.iat[row, 9], raw.iat[row, 12] = value, 5, 999
            raw.to_excel(source, header=False, index=False)
            original_savefig = Figure.savefig

            for function in self.functions:
                for mode in ("none", "log2fc", "zscore_row"):
                    for preserve_order in (True, False):
                        with self.subTest(function=function, mode=mode, preserve_order=preserve_order):
                            recorded = {}

                            def capture(figure, path, *args, **kwargs):
                                path = Path(path)
                                if path.parent.name == "Raw Bar Plots":
                                    ax = figure.axes[0]
                                    labels = [tick.get_text() for tick in ax.get_xticklabels()]
                                    recorded[ax.get_title()] = (
                                        labels, [bar.get_height() for bar in ax.patches]
                                    )
                                if mode == "none" and preserve_order:
                                    original_savefig(figure, path, *args, **kwargs)

                            with patch.object(Figure, "savefig", capture), contextlib.redirect_stdout(io.StringIO()):
                                heatmap, long_df = function(
                                    source, header_row=1, samples=["B", "A"],
                                    cytokines=["IL-6", "Flat"], scale_mode=mode,
                                    ref_sample="A", preserve_sample_order=preserve_order,
                                    cluster_rows=False, output_dir=folder,
                                    run_statistics=False,
                                )
                            self.assertEqual(set(long_df.Sample), {"A", "B"})
                            self.assertEqual(set(long_df.Cytokine), {"IL-6", "Flat"})
                            self.assertEqual(set(recorded), {"IL-6", "Flat"})
                            labels, heights = recorded["IL-6"]
                            expected = {"B\n(n=1)": 8, "A\n(n=2)": 4}
                            self.assertEqual(dict(zip(labels, heights)), expected)
                            if preserve_order:
                                self.assertEqual(labels, list(expected))
                            np.testing.assert_allclose(recorded["Flat"][1], [5, 5])
                            if mode == "zscore_row":
                                np.testing.assert_allclose(
                                    heatmap.loc["IL-6", ["B", "A"]].to_numpy(),
                                    np.array([1, -1]) / np.sqrt(2),
                                )
                                self.assertNotIn("Flat", heatmap.index)
                            self.assertEqual(plt.get_fignums(), [])
            images = list((folder / "Raw Bar Plots").glob("*.png"))
            self.assertEqual(len(images), 2)
            self.assertTrue(all(image.stat().st_size > 10000 for image in images))

    def test_statistics_export_accepts_two_and_three_replicates_in_both_definitions(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            source = folder / "input.xlsx"
            raw = pd.DataFrame([[None] * 7 for _ in range(6)])
            raw.iat[0, 5] = "IL-6"
            for row, (rep, sample, value) in enumerate([
                ("1", "Control", 1), ("2", "Control", 2),
                ("1", "Treatment", 3), ("2", "Treatment", 4), ("3", "Treatment", 5),
            ], start=1):
                raw.iat[row, 0], raw.iat[row, 1], raw.iat[row, 6] = rep, sample, value
            raw.to_excel(source, header=False, index=False)
            for function in self.functions:
                for mode in ("none", "zscore_row"):
                    with patch.object(Figure, "savefig"), contextlib.redirect_stdout(io.StringIO()):
                        function(
                            source, header_row=1, samples=["Control", "Treatment"],
                            control_sample="Control", scale_mode=mode, output_dir=folder,
                        )
                    workbooks = list((folder / "Raw Bar Plots").glob("*Results.xlsx"))
                    self.assertEqual(len(workbooks), 1)
                    overall = pd.read_excel(workbooks[0], sheet_name="Kruskal_Wallis")
                    pairs = pd.read_excel(workbooks[0], sheet_name="Dunn_vs_Control")
                    self.assertEqual(overall.loc[0, "Status"], "tested")
                    self.assertEqual(pairs.loc[0, "N_Control"], 2)
                    self.assertEqual(pairs.loc[0, "N_Group"], 3)
                    self.assertAlmostEqual(
                        overall.loc[0, "P_Value"], stats.kruskal([1, 2], [3, 4, 5]).pvalue,
                    )
                    self.assertEqual(plt.get_fignums(), [])

    def test_unknown_control_raises_explicit_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            source = folder / "input.xlsx"
            raw = pd.DataFrame([[None] * 7 for _ in range(3)])
            raw.iat[0, 5] = "IL-6"
            for row, sample in enumerate(["A", "B"], start=1):
                raw.iat[row, 0], raw.iat[row, 1], raw.iat[row, 6] = "1", sample, row
            raw.to_excel(source, header=False, index=False)
            for function in self.functions:
                with self.assertRaisesRegex(ValueError, "must be included"), contextlib.redirect_stdout(io.StringIO()):
                    function(source, header_row=1, output_dir=folder)

    def test_censored_observations_do_not_enter_inference_as_zero(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            source = folder / "input.xlsx"
            raw = pd.DataFrame([[None] * 7 for _ in range(5)])
            raw.iat[0, 5] = "IL-6"
            for row, (sample, rep, value) in enumerate([
                ("Control", "1", 1), ("Control", "2", 2),
                ("Treatment", "1", "OOR <"), ("Treatment", "2", 4),
            ], start=1):
                raw.iat[row, 0], raw.iat[row, 1], raw.iat[row, 6] = rep, sample, value
            raw.to_excel(source, header=False, index=False)
            for function in self.functions:
                with patch.object(Figure, "savefig"), contextlib.redirect_stdout(io.StringIO()):
                    _, parsed = function(
                        source, header_row=1, samples=["Control", "Treatment"],
                        control_sample="Control", output_dir=folder,
                    )
                self.assertIn(0, parsed.Obs_Conc.to_list())
                workbook = next((folder / "Raw Bar Plots").glob("*Results.xlsx"))
                overall = pd.read_excel(workbook, sheet_name="Kruskal_Wallis")
                pairs = pd.read_excel(workbook, sheet_name="Dunn_vs_Control")
                self.assertIn("censored", overall.loc[0, "Status"])
                self.assertTrue(pairs.P_Holm.isna().all())
                self.assertFalse(pairs.Significant.any())


class HeatmapStatisticsTests(unittest.TestCase):
    def analyze(self, data, censored=None):
        with contextlib.redirect_stdout(io.StringIO()):
            return analyze_heatmap_groups(
                data, ["Control", "A", "B"], list(data.Cytokine.unique()),
                "Control", 0.05, set() if censored is None else censored,
            )

    def data(self):
        rows = []
        for cytokine in ["IL-6", "Constant"]:
            for sample, values in [("Control", [1, 2]), ("A", [3, 4, 5]), ("B", [6, 7, 8])]:
                for rep, value in enumerate(values, start=1):
                    rows.append((sample, str(rep), cytokine, value if cytokine == "IL-6" else 1))
        return pd.DataFrame(rows, columns=["Sample", "Rep", "Cytokine", "Obs_Conc"])

    def test_unequal_n_dunn_holm_bh_and_constants(self):
        _, overall, pairs = self.analyze(self.data())
        expected = stats.kruskal([1, 2], [3, 4, 5], [6, 7, 8])
        self.assertAlmostEqual(overall.loc[0, "H_Statistic"], expected.statistic)
        self.assertAlmostEqual(overall.loc[0, "P_Value"], expected.pvalue)
        self.assertAlmostEqual(overall.loc[0, "Q_BH"], min(1, 2 * expected.pvalue))
        first, second = 2 * stats.norm.sf(2.5 / np.sqrt(5)), 2 * stats.norm.sf(5.5 / np.sqrt(5))
        np.testing.assert_allclose(pairs.P_Unadjusted.iloc[:2], [first, second])
        np.testing.assert_allclose(pairs.P_Holm.iloc[:2], [max(first, 2 * second), 2 * second])
        self.assertTrue(overall.Status.eq("tested").all())
        self.assertEqual(overall.loc[1, "P_Value"], 1)
        self.assertTrue(pairs.P_Holm.iloc[2:].eq(1).all())

    def test_technical_repeats_are_not_extra_biological_replicates(self):
        data = self.data()
        repeated = pd.concat([data, data.iloc[[0]]], ignore_index=True)
        points, actual, pairs = self.analyze(repeated)
        _, expected, expected_pairs = self.analyze(data)
        pd.testing.assert_frame_equal(actual, expected)
        pd.testing.assert_frame_equal(pairs, expected_pairs)
        self.assertEqual(points.iloc[0].Number_Rows, 2)

    def test_tied_values_match_existing_kruskal_wallis_implementation(self):
        data = self.data()
        data.Obs_Conc = data.Obs_Conc // 2
        _, overall, pairs = self.analyze(data)
        groups = [data.loc[(data.Cytokine == "IL-6") & (data.Sample == sample), "Obs_Conc"]
                  for sample in ["Control", "A", "B"]]
        self.assertAlmostEqual(overall.loc[0, "P_Value"], stats.kruskal(*groups).pvalue)
        variance = 6 - 18 / 84
        expected_z = np.array([23 / 6 - 7 / 4, 7 - 7 / 4]) / np.sqrt(variance * (1 / 2 + 1 / 3))
        np.testing.assert_allclose(pairs.Z_Statistic.iloc[:2], expected_z)
        self.assertTrue(pairs.P_Holm.notna().all())

    def test_bh_family_excludes_skipped_cytokines(self):
        data = self.data()
        data.loc[(data.Cytokine == "IL-6") & (data.Sample == "Control") & (data.Rep == "2"), "Obs_Conc"] = np.nan
        _, overall, _ = self.analyze(data)
        self.assertTrue(overall.FDR_Family_Size.eq(1).all())
        self.assertTrue(np.isnan(overall.loc[0, "Q_BH"]))
        self.assertEqual(overall.loc[1, "Q_BH"], 1)

    def test_small_groups_censoring_and_missing_ids_have_explicit_skips(self):
        for reason in ["small", "censored", "missing ID"]:
            data = self.data()
            censored = set()
            if reason == "small":
                data.loc[(data.Sample == "Control") & (data.Rep == "2"), "Obs_Conc"] = np.nan
            elif reason == "censored":
                censored = {("A", "1", "IL-6"), ("A", "1", "Constant")}
            else:
                data.loc[data.Sample == "Control", "Rep"] = None
            _, overall, pairs = self.analyze(data, censored)
            self.assertTrue(overall.Status.str.startswith("not tested:").all())
            self.assertTrue(pairs.Status.str.startswith("not tested:").all())
            self.assertTrue(pairs.P_Holm.isna().all())
            self.assertFalse(pairs.Significant.any())

    def test_nonfinite_data_raises(self):
        data = self.data()
        data["Obs_Conc"] = data.Obs_Conc.astype(float)
        data.loc[0, "Obs_Conc"] = np.inf
        with self.assertRaisesRegex(ValueError, "Nonfinite"):
            self.analyze(data)

    def test_significance_annotation_and_skip_message(self):
        _, overall, pairs = self.analyze(self.data())
        pairs.loc[0, ["P_Holm", "Significant"]] = [0.001, True]
        fig, ax = plt.subplots()
        try:
            annotate_heatmap_bar(ax, "IL-6", ["Control", "A", "B"], overall, pairs, "Control")
            self.assertTrue(any(text.get_text() == "**" for text in ax.texts))
            overall.loc[0, "Status"] = "not tested: insufficient replicates"
            annotate_heatmap_bar(ax, "IL-6", ["Control", "A", "B"], overall, pairs, "Control")
            self.assertIn("not tested", ax.get_title())
        finally:
            plt.close(fig)


if __name__ == "__main__":
    unittest.main()
