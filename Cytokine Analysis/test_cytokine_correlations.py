import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from Cytokine_Morphology_Correlation import (
    Plot_Cytokine_Morphology_Correlation_Plots, read_cell_area_points,
    read_paired_morphology_points,
)


class CorrelationTests(unittest.TestCase):
    def setUp(self):
        cyto, morph = [], []
        for dose, slope in [(5, 2), (100, -3)]:
            for time, value in [(6, 1), (24, 2), (72, 4)]:
                sample = f"{time}hr {dose}ng"
                for rep, offset in [("A", -0.5), ("B", 0.5)]:
                    # Different technical-repeat counts must not weight replicates.
                    for _ in range(3 if rep == "A" else 1):
                        cyto.append(dict(
                            Cytokine="IL-6", Sample=sample, Rep=rep, Time_hr=time,
                            Obs_Conc=20 + slope * value + offset,
                        ))
                for rep, offset in [("unrelated-1", -0.25), ("unrelated-2", 0.25)]:
                    for scale, result in [("ln", value + offset), ("raw", np.exp(value + offset))]:
                        morph.append(dict(
                            Metric="Area", Sample=sample, Biological_Rep=rep,
                            Time_hr=time, Results=result, Morphology_Scale=scale, SEM=0.1,
                        ))
        self.cyto = pd.DataFrame(cyto)
        self.morph = pd.DataFrame(morph)

    def call(self, **kwargs):
        return Plot_Cytokine_Morphology_Correlation_Plots(
            self.cyto, self.morph, display_results=False, **kwargs,
        )

    def test_separate_panels_save_independent_figures_by_dose(self):
        self.cyto = pd.concat([self.cyto, self.cyto.assign(Cytokine="IL-10")], ignore_index=True)
        self.morph = self.morph.assign(Metric="Area (µm²)")
        with tempfile.TemporaryDirectory() as tmp:
            with patch.object(plt.Figure, "savefig", autospec=True) as save:
                _, stats, paths = self.call(
                    separate_panels=True, descriptive_only=True,
                    concentration_unit="normalized",
                    label_inset_points=False, output_dir=tmp,
                )
            self.assertEqual(len(paths), 2 * len(stats))
            for index, row in enumerate(stats.itertuples()):
                (main, main_path), (correlation, correlation_path) = [
                    call.args for call in save.call_args_list[2 * index:2 * index + 2]
                ]
                folder = Path(tmp) / "Area" / f"{row.LPS_ng_mL:g} ng"
                self.assertEqual(main_path.parent, folder / "Main")
                self.assertEqual(correlation_path.parent, folder / "Correlation")
                self.assertEqual(len(main.axes), 2)
                self.assertFalse(main.axes[0].child_axes)
                self.assertEqual(len(correlation.axes), 1)
                self.assertFalse(main.legends)
                self.assertFalse(correlation.legends)
                self.assertFalse(correlation.texts)
                self.assertEqual(main.axes[0].get_title(), f"{row.Cytokine} and Area")
                self.assertEqual(main.axes[0].get_ylabel(), "ln(Area (µm²))")
                self.assertEqual(correlation.axes[0].get_title(), main.axes[0].get_title())
                self.assertEqual(main.axes[1].get_ylabel(), f"{row.Cytokine} Concentration, normalized")
                self.assertIn("Unadjusted R²", correlation.axes[0].texts[0].get_text())
            legends = {path.parent: figure for figure, path in (
                call.args for call in save.call_args_list if call.args[1].name == "Legend.png"
            )}
            self.assertEqual(len(legends), 4)
            for dose in (5, 100):
                folder = Path(tmp) / "Area" / f"{dose} ng"
                main_key = legends[folder / "Main"]
                correlation_key = legends[folder / "Correlation"]
                self.assertFalse(main_key.axes or correlation_key.axes)
                self.assertEqual(len(main_key.legends[0].get_texts()), 5)
                self.assertEqual(
                    {text.get_text().splitlines()[0] for text in correlation_key.legends[0].get_texts()},
                    {"IL-6", "IL-10"},
                )

    def test_paired_model_uses_complete_pairs_and_reports_missing_cytokine(self):
        cyto, morph = [], []
        for i, time in enumerate([6, 24, 72]):
            for j, rep in enumerate(["1", "2", "3"]):
                interaction = i * j
                morph.append(dict(Metric="Area", Sample=f"{time}hr_5ng", Time_hr=time,
                                  Biological_Rep=rep, Morphology_Scale="ln",
                                  Results=100, Area_Log_Adjusted=10 + i + 3 * j + interaction))
                cyto.append(dict(Cytokine="IL-6", Sample=f"{time}hr 5ng", Time_hr=time,
                                 Rep=f"Rep {rep}", Obs_Conc=(np.nan if (i, j) == (1, 2)
                                                             else 2.0 ** (20 + 4 * i - j + 2 * interaction))))
        self.cyto, self.morph = pd.DataFrame(cyto), pd.DataFrame(morph)
        mapping = {f"Rep {rep}": rep for rep in ["1", "2", "3"]}
        with tempfile.TemporaryDirectory() as tmp:
            with patch.object(plt.Figure, "savefig", autospec=True) as save:
                matched, stats, _ = self.call(control_for_time=True, replicate_map=mapping, output_dir=tmp)
            self.assertEqual(len(matched), 8)
            row = stats.iloc[0]
            self.assertEqual((row.n_pairs, row.n_biological_reps, row.Residual_df), (8, 3, 2))
            self.assertAlmostEqual(row.Partial_R_squared, 1)
            self.assertAlmostEqual(row.Residual_Slope, 0.5)
            figure = save.call_args.args[0]
            inset = figure.axes[0].child_axes[0]
            self.assertIn("Partial R² = 1.000", [text.get_text() for text in inset.texts])
            self.assertEqual(inset.get_xlabel(), "log₂(IL-6 concentration) residual")
            self.assertEqual(inset.get_ylabel(), "ln(area + 0.5) residual")
            np.testing.assert_allclose(matched.Cytokine_Log2, np.log2(matched.Obs_Conc))
            np.testing.assert_allclose(matched.Morphology_Residual, 0.5 * matched.Cytokine_Residual, atol=1e-12)
            self.assertIn("Missing cytokine: 24 h, replicate 3",
                          [text.get_text() for text in figure.legends[0].get_texts()])
            self.assertIn("6h / R1", [text.get_text() for text in inset.texts])
            # The assays spell the same condition differently, but share one legend key.
            legend = figure.legends[0]
            keys = [handle for handle in legend.get_lines() if handle.get_label() == "6hr 5ng"]
            self.assertEqual(len(keys), 1)
            from matplotlib.colors import to_rgba
            for axis in figure.axes:
                np.testing.assert_allclose(axis.collections[0].get_facecolors()[0],
                                           to_rgba(keys[0].get_color()))
            with tempfile.TemporaryDirectory() as unlabeled_dir:
                with patch.object(plt.Figure, "savefig", autospec=True) as save:
                    _, unlabeled_stats, _ = self.call(
                        control_for_time=True, replicate_map=mapping,
                        output_dir=unlabeled_dir, label_inset_points=False,
                    )
                unlabeled = save.call_args.args[0].axes[0].child_axes[0]
                self.assertEqual([text.get_text() for text in unlabeled.texts], ["Partial R² = 1.000"])
                self.assertEqual(len(unlabeled.collections[0].get_offsets()), 8)
                self.assertEqual(len(unlabeled.lines), 1)
                pd.testing.assert_frame_equal(stats, unlabeled_stats)
        with self.assertRaisesRegex(ValueError, "one-to-one replicate_map"):
            self.call(control_for_time=True, replicate_map={"Rep 1": "1", "Rep 2": "1"})
        self.cyto = pd.concat([self.cyto, self.cyto.iloc[[0]]], ignore_index=True)
        with self.assertRaisesRegex(ValueError, "one measured cytokine value"):
            self.call(control_for_time=True, replicate_map=mapping, save_results=False)

    def test_paired_area_logs_each_cell_and_weights_wells_equally(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "Area.csv"
            pd.DataFrame({
                "Treatment": ["6hr_5ng"] * 4, "Biological_Rep": [1] * 4,
                "Well": [1, 1, 1, 2], "Results": [0.5, 1.5, 3.5, 15.5],
            }).to_csv(path, index=False)
            row = read_cell_area_points(path).iloc[0]
        self.assertAlmostEqual(row.Area_Log_Adjusted, 2.5 * np.log(2))
        self.assertAlmostEqual(row.Results, 5.25)
        self.assertEqual((row.Area_Measurements, row.Area_Wells), (4, 2))

    def test_paired_metric_transforms_with_unequal_well_sizes(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "Metric.csv"
            for metric, values, expected in (
                ("Perimeter (µm)", [2, 8, 32], 3.5 * np.log(2)),
                ("AUC", [2, 8, 32], 3.5 * np.log(2)),
                ("Circularity", [0.2, 0.4, 0.9], 0.6),
            ):
                with self.subTest(metric=metric):
                    pd.DataFrame({
                        "Treatment": ["6hr_5ng"] * 3, "Biological_Rep": [1] * 3,
                        "Well": [1, 1, 2], "Results": values,
                    }).to_csv(path, index=False)
                    row = read_paired_morphology_points(path, metric).iloc[0]
                    self.assertAlmostEqual(row.Morphology_Transformed, expected)
                    self.assertAlmostEqual(row.Results, np.mean(values))

    def test_descriptive_inset_shows_the_unadjusted_fit_for_each_dose(self):
        with tempfile.TemporaryDirectory() as tmp:
            with patch.object(plt.Figure, "savefig", autospec=True) as save:
                _, stats, _ = self.call(descriptive_only=True, output_dir=tmp)
            for call, row in zip(save.call_args_list, stats.itertuples()):
                figure = call.args[0]
                inset = figure.axes[0].child_axes[0]
                self.assertEqual(len(inset.lines), 1)
                self.assertTrue({"6h", "24h", "72h", "Unadjusted R² = 1.000"}
                                .issubset(text.get_text() for text in inset.texts))
                self.assertTrue(any("Unadjusted linear fit" in text.get_text()
                                    for text in figure.legends[0].get_texts()))
                slope = {5: 2, 100: -3}[row.LPS_ng_mL]
                self.assertAlmostEqual(row.Slope, slope)
                line = inset.lines[0]
                np.testing.assert_allclose(line.get_xdata()[[0, -1]], [1, 4])
                np.testing.assert_allclose(line.get_ydata(), 20 + slope * line.get_xdata())

    def test_descriptive_inset_suppresses_undefined_correlation(self):
        self.cyto["Obs_Conc"] = 1
        _, stats, _ = self.call(descriptive_only=True, save_results=False)
        self.assertTrue(stats.R_squared.isna().all())
        self.assertTrue(stats.Status.str.contains("Constant condition means").all())

    def test_exact_coefficients_separate_doses_and_equal_rep_weight(self):
        before_cyto, before_morph = self.cyto.copy(), self.morph.copy()
        matched, stats, paths = self.call(save_results=False)
        self.assertEqual(len(matched), 6)
        self.assertEqual(paths, [])
        self.assertEqual(set(matched.Morphology_Scale), {"ln"})
        self.assertTrue(matched.Cytokine_n.eq(2).all())
        self.assertTrue(matched.Morphology_n.eq(2).all())
        stats = stats.set_index("LPS_ng_mL")
        for dose, slope, r in [(5, 2, 1), (100, -3, -1)]:
            self.assertAlmostEqual(stats.loc[dose, "Slope"], slope)
            self.assertAlmostEqual(stats.loc[dose, "Intercept"], 20)
            self.assertAlmostEqual(stats.loc[dose, "Pearson_r"], r)
            self.assertAlmostEqual(stats.loc[dose, "R_squared"], 1)
            self.assertEqual(stats.loc[dose, "n_conditions"], 3)
        pd.testing.assert_frame_equal(self.cyto, before_cyto)
        pd.testing.assert_frame_equal(self.morph, before_morph)
        self.assertEqual(plt.get_fignums(), [])

    def test_missing_conditions_not_imputed_and_sparse_fit_reported(self):
        self.morph = self.morph.loc[self.morph.Sample.ne("72hr 5ng")]
        matched, stats, _ = self.call(save_results=False)
        sparse = stats.loc[stats.LPS_ng_mL.eq(5)].iloc[0]
        self.assertEqual(sparse.n_conditions, 2)
        self.assertIn("3 matched conditions", sparse.Status)
        self.assertTrue(pd.isna(sparse.Pearson_r))
        self.assertFalse(matched.Sample.eq("72hr 5ng").any())

    def test_constant_data_not_given_a_success_shaped_fit(self):
        self.cyto["Obs_Conc"] = 1
        _, stats, _ = self.call(save_results=False)
        self.assertTrue(stats.Status.str.startswith("Constant").all())
        self.assertTrue(stats.R_squared.isna().all())

    def test_partial_correlation_with_three_conditions_is_not_reported(self):
        matched, stats, _ = self.call(control_for_time=True, save_results=False)
        self.assertTrue(stats.Residual_df.eq(0).all())
        self.assertTrue(stats.Partial_R_squared.isna().all())
        self.assertTrue(stats.Status.str.contains("4 matched conditions").all())
        self.assertTrue(matched.Morphology_Residual.isna().all())

    def test_partial_correlation_removes_time_and_labels_inset(self):
        times = np.array([6, 24, 42, 60])
        x_residual = np.array([1, -1, -1, 1])
        y_residual = -2 * x_residual + np.array([1, -3, 3, -1])
        shared = {"Sample": [f"{time}hr 5ng" for time in times], "Time_hr": times}
        self.cyto = pd.DataFrame(shared | {
            "Cytokine": "IL-6", "Rep": "A", "Obs_Conc": 40 + 3 * times + y_residual,
        })
        self.morph = pd.DataFrame(shared | {
            "Metric": "Area", "Biological_Rep": "B", "Morphology_Scale": "ln",
            "Results": 10 + times + x_residual,
        })
        with tempfile.TemporaryDirectory() as tmp:
            with patch.object(plt.Figure, "savefig", autospec=True) as save:
                matched, stats, paths = self.call(control_for_time=True, output_dir=tmp)
            row = stats.iloc[0]
            self.assertEqual(row.Residual_df, 1)
            self.assertAlmostEqual(row.Partial_r, -2 / 3)
            self.assertAlmostEqual(row.Partial_R_squared, 4 / 9)
            ordered = matched.sort_values("Time_hr")
            np.testing.assert_allclose(ordered.Morphology_Residual, x_residual)
            np.testing.assert_allclose(ordered.Cytokine_Residual, y_residual)
            self.assertTrue(paths[0].endswith("_partial_correlation.png"))
            figure = save.call_args.args[0]
            inset = figure.axes[0].child_axes[0]
            self.assertIn("Partial R² = 0.444", [text.get_text() for text in inset.texts])
            self.assertTrue(any("Partial r = -0.667" in text.get_text()
                                for text in figure.legends[0].get_texts()))
            workbook = Path(tmp) / "Cytokine_Morphology_Partial_Correlations.xlsx"
            self.assertAlmostEqual(pd.read_excel(workbook, sheet_name="Partial Correlations")
                                   .iloc[0].Partial_R_squared, 4 / 9)

        self.cyto["Obs_Conc"] = 40 + 3 * times
        _, stats, _ = self.call(control_for_time=True, save_results=False)
        self.assertTrue(stats.Partial_R_squared.isna().all())
        self.assertIn("No residual variation", stats.iloc[0].Status)

    def test_invalid_inputs_and_no_match_raise(self):
        with self.assertRaisesRegex(ValueError, "Unknown Metric"):
            self.call(metrics=["missing"], save_results=False)
        with self.assertRaisesRegex(ValueError, "morphology_scale"):
            self.call(morphology_scale="log10", save_results=False)
        self.morph.loc[0, "Results"] = np.inf
        with self.assertRaisesRegex(ValueError, "invalid/nonfinite Results"):
            self.call(save_results=False)
        self.morph.loc[0, "Results"] = 1
        self.morph["Sample"] = "48hr 5ng"
        self.morph["Time_hr"] = 48
        self.morph["Biological_Rep"] = self.morph.index.astype(str)
        with self.assertRaisesRegex(ValueError, "No numeric conditions match"):
            self.call(save_results=False)

    def test_rejects_duplicate_morphology_replicates_and_inconsistent_times(self):
        self.morph = pd.concat([self.morph, self.morph.iloc[[0]]], ignore_index=True)
        with self.assertRaisesRegex(ValueError, "one row per"):
            self.call(save_results=False)
        self.morph = self.morph.iloc[:-1].copy()
        self.morph.loc[0, "Time_hr"] = 24
        with self.assertRaisesRegex(ValueError, "disagree"):
            self.call(save_results=False)

    def test_saves_inset_images_and_workbook(self):
        with tempfile.TemporaryDirectory() as tmp:
            with patch.object(plt.Figure, "savefig", autospec=True, wraps=None) as save:
                self.call(output_dir=tmp, save_results=True)
                self.assertEqual(len(save.call_args_list), 2)
                figure = save.call_args_list[0].args[0]
                self.assertEqual(len(figure.axes), 2)
                left, right = figure.axes
                self.assertEqual(left.get_ylabel(), "ln(Area)")
                self.assertEqual(right.get_ylabel(), "IL-6 observed concentration")
                np.testing.assert_allclose(left.collections[0].get_offsets()[:, 1], [0.75, 1.25])
                np.testing.assert_allclose(right.collections[0].get_offsets()[:, 1], [16.5, 17.5])
                self.assertEqual(len(left.containers), 3)
                self.assertTrue(any(line.get_linestyle() == "--" for line in left.lines))
                self.assertTrue(any(line.get_linestyle() == "-" for line in right.lines))
                from matplotlib.collections import PathCollection

                for axis in (left, right):
                    scatter = [item for item in axis.collections if isinstance(item, PathCollection)]
                    self.assertEqual(len(scatter), 3)
                    self.assertEqual(sum(len(item.get_offsets()) for item in scatter), 6)
                self.assertEqual(len(figure.axes[0].child_axes), 1)
                inset = figure.axes[0].child_axes[0]
                self.assertEqual(len(inset.collections[0].get_offsets()), 3)
                self.assertTrue(any("slope =" in text.get_text() for text in inset.texts))
                self.assertIsNone(left.get_legend())
                self.assertEqual(len(figure.legends), 1)
                self.assertTrue(any(
                    "dashed trend" in text.get_text().lower()
                    for text in figure.legends[0].get_texts()
                ))
                figure.canvas.draw()
                renderer = figure.canvas.get_renderer()
                legend_bounds = figure.legends[0].get_window_extent(renderer)
                self.assertGreater(legend_bounds.x0, inset.get_tightbbox(renderer).x1)
                self.assertLessEqual(legend_bounds.x1, figure.bbox.x1)
            matched, stats, paths = self.call(output_dir=tmp, save_results=True)
            self.assertEqual(len(paths), 2)
            self.assertTrue(all(Path(path).stat().st_size > 10000 for path in paths))
            workbook = Path(tmp) / "Cytokine_Morphology_Correlations.xlsx"
            pd.testing.assert_frame_equal(
                pd.read_excel(workbook, sheet_name="Exploratory Regressions"),
                stats, check_dtype=False, check_exact=False, rtol=1e-12, atol=1e-12,
            )
            self.assertEqual(
                len(pd.read_excel(workbook, sheet_name="Matched Condition Means")), len(matched),
            )
        self.assertEqual(plt.get_fignums(), [])

    def test_unmatched_conditions_stay_in_timecourse_not_regression(self):
        self.cyto = self.cyto.loc[self.cyto.Sample.ne("72hr 100ng")]
        with tempfile.TemporaryDirectory() as tmp:
            with patch.object(plt.Figure, "savefig", autospec=True) as save:
                matched, statistics, _ = self.call(output_dir=tmp, save_results=True)
                figure = save.call_args_list[0].args[0]
                left, right = figure.axes
                self.assertEqual(len(left.containers), 3)
                self.assertEqual(len(right.collections), 2)
                inset = left.child_axes[0]
                self.assertEqual(len(inset.collections[0].get_offsets()), 2)
        self.assertFalse(matched.Sample.eq("72hr 100ng").any())
        self.assertEqual(statistics.iloc[0].n_conditions, 2)

    def test_invalid_sem_raises(self):
        self.morph.loc[0, "SEM"] = -0.1
        with self.assertRaisesRegex(ValueError, "SEM must be nonnegative"):
            self.call(save_results=False)


if __name__ == "__main__":
    unittest.main()
