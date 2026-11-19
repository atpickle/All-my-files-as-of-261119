"""Standalone 48hr 5ug vs 6hr Ctrl analysis from 260527_Cytokine Analysis.ipynb.

Run this script to analyze "261018 Updated 48hr", rows 194-235, in the
workbook configured under __main__. Requires pandas, numpy, scipy, matplotlib,
and openpyxl. Obs Conc. values are used as supplied, without transformation;
technical repeats are averaged within each biological replicate.

Plots show mean +/- SEM and replicate dots. Brackets mark BH-adjusted q-values
below 0.05, 0.01, 0.001, and 0.0001 with one to four asterisks, respectively.
Independent, two-sided Mann-Whitney
U tests use exact label permutations (including ties) up to 100,000 allocations,
then tie/continuity-corrected asymptotic tests. BH FDR spans valid cytokine tests.
With three replicates per group, the smallest exact two-sided p-value is 0.10.
Do not use this independent-groups test on values normalized to a shared control.

PNGs and a workbook containing statistics, biological replicates, input values,
and analysis notes are saved beside the input workbook in a sheet/row-labeled
Unpaired_Cytokine_Comparison subfolder. The original notebook is unchanged.
Plots are saved without opening windows by default; display_results=True opts
into interactive display.
Running the script also saves legend-free copies with uniform black dots in
No_Legend_Black_Dots. Original plots retain their legends and replicate colors.
"""


def Plot_Unpaired_Cytokine_Comparison(
    file_path,
    sheet_name="Morphology Paper Cytokines",
    *,
    control="6hr Ctrl",
    treatment="48hr 5ug",
    header_row=3,
    data_start_row=6,
    data_end_row=None,
    sample_col="B",
    rep_col="A",
    start_col="F",
    col_step=3,
    obs_conc_offset=1,
    cytokines=None,
    expected_replicates=3,
    exclude_replicates=None,
    output_dir=None,
    save_results=True,
    display_results=False,
    save_clean_plots=False,
    plot_group_labels=None,
    plot_value_label="Value (as supplied)",
):
    """Compare independent groups of supplied values with Mann-Whitney U.

    Returns (bio_points, statistics, saved_paths). Tests require the expected
    replicates per group; missing/censored values are never substituted.
    No transformation or normalization is applied. Technical repeats are
    averaged on the supplied scale. Mean differences are descriptive, in
    supplied units, treatment minus control. U uses treatment
    first. Exact permutation p-values handle ties for up to 100,000 allocations;
    larger samples use the tie- and continuity-corrected asymptotic test.
    exclude_replicates is a list of (group, replicate ID) pairs. Expected
    counts are validated before exclusion, then reduced only in that group.
    save_clean_plots saves additional legend-free, black-dot PNGs when
    save_results=True, without changing the original plot style.
    plot_group_labels overrides the two tick labels and uses cytokine-only
    titles. plot_value_label changes only the y-axis text.
    """
    from pathlib import Path
    from math import comb
    import re
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    from matplotlib.collections import PathCollection
    import numpy as np
    import pandas as pd
    from openpyxl.utils import column_index_from_string
    from scipy import stats

    def mann_whitney_u(x, y):
        ranks = stats.rankdata(np.concatenate((x, y)))
        return ranks[:len(x)].sum() - len(x) * (len(x) + 1) / 2

    for label, value in (("header_row", header_row), ("data_start_row", data_start_row),
                         ("col_step", col_step), ("expected_replicates", expected_replicates)):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"{label} must be a positive integer.")
    if expected_replicates < 2:
        raise ValueError("expected_replicates must be at least 2.")
    if isinstance(obs_conc_offset, bool) or not isinstance(obs_conc_offset, int) or not 0 <= obs_conc_offset < col_step:
        raise ValueError("obs_conc_offset must be an integer within each cytokine block.")
    if any(not isinstance(name, str) or not name.strip() for name in (control, treatment)):
        raise ValueError("Treatment names must be nonempty strings.")
    control, treatment = control.strip(), treatment.strip()
    if control == treatment:
        raise ValueError("Control and treatment must be different.")
    exclusions = [] if exclude_replicates is None else exclude_replicates
    if (not isinstance(exclusions, list)
            or any(not isinstance(pair, (tuple, list)) or len(pair) != 2
                   or any(not isinstance(item, str) or not item.strip() for item in pair)
                   for pair in exclusions)):
        raise ValueError("exclude_replicates must be a list of (group, replicate ID) pairs.")
    exclusions = [(group.strip(), rep.strip()) for group, rep in exclusions]
    if len(set(exclusions)) != len(exclusions) or any(group not in (control, treatment) for group, rep in exclusions):
        raise ValueError("Exclusions must be unique and refer to the selected groups.")
    if cytokines is not None:
        if (not isinstance(cytokines, list) or not cytokines
                or any(not isinstance(name, str) or not name.strip() for name in cytokines)):
            raise ValueError("cytokines must be a nonempty list of names, or None.")
        cytokines = [name.strip() for name in cytokines]
        if len(set(cytokines)) != len(cytokines):
            raise ValueError("cytokines contains duplicate names.")
    indices = []
    for column in (sample_col, rep_col, start_col):
        if not isinstance(column, str) or not re.fullmatch(r"[A-Za-z]+", column):
            raise ValueError("Column settings must be Excel column letters.")
        indices.append(column_index_from_string(column) - 1)
    sample_idx, rep_idx, first_idx = indices
    source = Path(file_path).expanduser().resolve()
    with pd.ExcelFile(source) as book:
        if not isinstance(sheet_name, str) or sheet_name not in book.sheet_names:
            raise ValueError(f"Unknown sheet {sheet_name!r}. Available: {book.sheet_names}")
        raw = pd.read_excel(book, sheet_name=sheet_name, header=None)
    last_row = len(raw) if data_end_row is None else data_end_row
    if isinstance(last_row, bool) or not isinstance(last_row, int) or not 1 <= header_row < data_start_row <= last_row <= len(raw):
        raise ValueError("Choose valid inclusive Excel data rows after header_row.")
    if max(indices) >= raw.shape[1]:
        raise ValueError("A selected column lies outside the sheet.")
    headers = {}
    for idx in range(first_idx, raw.shape[1], col_step):
        value = raw.iat[header_row - 1, idx]
        if pd.isna(value) or not str(value).strip():
            continue
        name = str(value).strip()
        if name in headers or idx + obs_conc_offset >= raw.shape[1]:
            raise ValueError(f"Duplicate header or missing value column for {name!r}.")
        headers[name] = idx + obs_conc_offset
    names = list(headers) if cytokines is None else cytokines
    if not names or set(names) - set(headers):
        raise ValueError(f"Choose cytokines from these headers: {list(headers)}")
    rows = raw.iloc[data_start_row - 1:last_row].copy()
    labels = rows.iloc[:, sample_idx].fillna("").astype(str).str.strip()
    absent = [name for name in (control, treatment) if not labels.eq(name).any()]
    if absent:
        raise ValueError(f"Missing groups {absent}. Available: {sorted(set(labels) - {''})}")
    records = []
    for row_idx in rows.index[labels.isin([control, treatment])]:
        rep = raw.iat[row_idx, rep_idx]
        if pd.isna(rep) or not str(rep).strip():
            raise ValueError(f"Missing biological replicate ID at Excel row {row_idx + 1}.")
        for name in names:
            value = raw.iat[row_idx, headers[name]]
            if pd.isna(value) or (isinstance(value, str) and not value.strip()):
                value = np.nan
            else:
                try:
                    value = float(value)
                except (TypeError, ValueError) as exc:
                    raise ValueError(f"Nonnumeric/censored {name} at Excel row {row_idx + 1}: {value!r}. Resolve on the supplied scale first.") from exc
                if not np.isfinite(value):
                    raise ValueError(f"Nonfinite {name} at Excel row {row_idx + 1}.")
            records.append((labels.loc[row_idx], str(rep).strip(), name, value, row_idx + 1))
    observations = pd.DataFrame(records, columns=["Sample", "Rep", "Cytokine", "Value", "Excel_Row"])
    rep_sets = [set(observations.loc[observations.Sample.eq(group), "Rep"]) for group in (control, treatment)]
    if any(len(ids) != expected_replicates for ids in rep_sets):
        raise ValueError(f"Expected {expected_replicates} biological IDs per group; found {control}: {sorted(rep_sets[0])}; {treatment}: {sorted(rep_sets[1])}.")
    excluded_observations = observations.iloc[:0].copy()
    expected_counts = {control: expected_replicates, treatment: expected_replicates}
    excluded_mask = pd.Series(False, index=observations.index)
    for group, rep in exclusions:
        matches = observations.Sample.eq(group) & observations.Rep.eq(rep)
        if not matches.any():
            raise ValueError(f"Excluded replicate not found: {group}: {rep}.")
        excluded_mask |= matches
        expected_counts[group] -= 1
    if min(expected_counts.values()) < 2:
        raise ValueError("Exclusions must leave at least two biological replicates in each group.")
    exclusion_note = "; ".join(f"{group}: {rep}" for group, rep in exclusions) or "none"
    if exclusions:
        excluded_observations = observations.loc[excluded_mask].copy()
        observations = observations.loc[~excluded_mask].copy()
        print(f"Sensitivity analysis; excluded {exclusion_note}.")
    bio_points = observations.groupby(["Cytokine", "Sample", "Rep"], sort=False, as_index=False).agg(
        Value=("Value", "mean"), Number_Values=("Value", "count"), Number_Rows=("Value", "size"))
    missing_count = int(observations.Value.isna().sum())
    print(f"Independent groups: {expected_counts}. Missing input values: {missing_count}. Values used as supplied, without transformation; technical repeats averaged on the supplied scale.")
    results, group_values = [], {}
    for name in names:
        metric = bio_points.loc[bio_points.Cytokine.eq(name)]
        groups = {group: metric.loc[metric.Sample.eq(group)].set_index("Rep").Value for group in (control, treatment)}
        group_values[name] = groups
        control_values, treatment_values = (groups[group].dropna() for group in (control, treatment))
        nc, nt = len(control_values), len(treatment_values)
        difference = treatment_values.mean() - control_values.mean()
        row = dict(Cytokine=name, Control=control, Treatment=treatment,
                   N_Control=nc, N_Treatment=nt,
                   Mean_Control=control_values.mean(), SEM_Control=control_values.sem(),
                   Mean_Treatment=treatment_values.mean(), SEM_Treatment=treatment_values.sem(),
                   Median_Control=control_values.median(), Median_Treatment=treatment_values.median(),
                   Mean_Difference=difference, U_Statistic=np.nan,
                   P_Value=np.nan, Q_BH=np.nan, P_Method="not tested",
                   Status="tested",
                   Missing_Control_IDs=", ".join(groups[control].index[groups[control].isna()]),
                   Missing_Treatment_IDs=", ".join(groups[treatment].index[groups[treatment].isna()]))
        if nc != expected_counts[control] or nt != expected_counts[treatment]:
            row["Status"] = f"not tested: requires {expected_counts[control]} numeric control and {expected_counts[treatment]} numeric treatment replicates"
        else:
            if comb(nc + nt, nt) <= 100_000:
                test = stats.permutation_test(
                    (treatment_values.to_numpy(), control_values.to_numpy()),
                    mann_whitney_u, permutation_type="independent",
                    alternative="two-sided", n_resamples=np.inf, vectorized=False)
                row["P_Method"] = "exact permutation (ties included)"
            else:
                test = stats.mannwhitneyu(treatment_values, control_values,
                                          alternative="two-sided", method="asymptotic",
                                          use_continuity=True)
                row["P_Method"] = "asymptotic (tie and continuity corrected)"
            if not np.isfinite(test.statistic) or not np.isfinite(test.pvalue):
                raise ValueError(f"Nonfinite Mann-Whitney U result for {name!r}.")
            row.update(U_Statistic=float(test.statistic), P_Value=float(test.pvalue))
        results.append(row)
    statistics = pd.DataFrame(results)
    valid = statistics.P_Value.dropna().sort_values()
    if len(valid):
        adjusted = valid.to_numpy() * len(valid) / np.arange(1, len(valid) + 1)
        statistics.loc[valid.index, "Q_BH"] = np.minimum.accumulate(adjusted[::-1])[::-1].clip(0, 1)
    statistics["FDR_Family_Size"] = len(valid)
    print(f"Two-sided Mann-Whitney U tests: {len(valid)}/{len(names)} cytokines; BH FDR across {len(valid)} valid tests. Exploratory; planned retained counts: {expected_counts}.")
    if list(expected_counts.values()) == [3, 3]:
        print("At n=3 per group, the smallest possible two-sided exact p-value is 0.10.")
    elif sorted(expected_counts.values()) == [2, 3]:
        print("At n=3 versus n=2, the smallest possible two-sided exact p-value is 0.20.")
    for row in statistics.itertuples():
        if row.Status != "tested":
            print(f"{row.Cytokine}: {row.Status}; missing control IDs: {row.Missing_Control_IDs or 'none'}; missing treatment IDs: {row.Missing_Treatment_IDs or 'none'}.")
    default_folder = "Unpaired_Cytokine_Comparison_Excluded_Replicates" if exclusions else "Unpaired_Cytokine_Comparison"
    sheet_tag = re.sub(r"[^\w-]+", "_", sheet_name).strip("_") or "sheet"
    default_folder += f"__{sheet_tag}__rows{data_start_row}-{last_row}"
    output = Path(output_dir).expanduser().resolve() if output_dir is not None else source.parent / default_folder
    saved_paths = []
    if save_results:
        output.mkdir(parents=True, exist_ok=True)
    colors = plt.get_cmap("tab10")
    group_colors = ("#BDBDBD", "#E69F00")
    for index, row in enumerate(statistics.itertuples(), start=1):
        if display_results:
            fig, ax = plt.subplots(figsize=(6, 5))
        else:
            fig = Figure(figsize=(6, 5))
            ax = fig.subplots()
        try:
            groups = group_values[row.Cytokine]
            for x, group in enumerate((control, treatment)):
                values = groups[group].dropna()
                if len(values):
                    ax.bar(x, values.mean(), width=0.6, color=group_colors[x],
                           edgecolor="black", linewidth=1, alpha=0.65, zorder=1)
                    if len(values) > 1:
                        ax.errorbar(x, values.mean(), yerr=values.sem(), fmt="none",
                                    color="black", capsize=4, zorder=3)
                offsets = np.linspace(-0.16, 0.16, len(groups[group]))
                for rep_index, (rep, value) in enumerate(groups[group].sort_index().items()):
                    ax.scatter(x + offsets[rep_index], value, color=colors((x * expected_replicates + rep_index) % 10),
                               edgecolors="black", linewidths=0.5, s=45, zorder=4, label=f"{group}: {rep}")
            ax.set_xticks([0, 1])
            if plot_group_labels is None:
                annotation = (f"Mann-Whitney U: p={row.P_Value:.4g}; BH q={row.Q_BH:.4g}"
                              if row.Status == "tested" else row.Status)
                title = f"{row.Cytokine}\n{annotation}"
                if exclusions:
                    title += f"\nExcluded: {exclusion_note}"
                ax.set_title(title, fontsize=10)
                ax.set_xticklabels([f"{control}\n(n={row.N_Control})", f"{treatment}\n(n={row.N_Treatment})"])
            else:
                ax.set_title(row.Cytokine, fontsize=10)
                ax.set_xticklabels(plot_group_labels)
            ax.set_ylabel(plot_value_label)
            ax.set_xlim(-0.6, 1.6)
            ax.set_xlabel("Treatment group (bars: mean +/- SEM)")
            ax.axhline(0, color="0.8", linewidth=0.7, zorder=0)
            if row.Status == "tested" and row.Q_BH < 0.05:
                stars = "*" * sum(row.Q_BH < cutoff for cutoff in (0.05, 0.01, 0.001, 0.0001))
                lower, upper = ax.get_ylim()
                span = upper - lower
                top = upper + 0.04 * span
                ax.plot([0, 0, 1, 1], [upper, top, top, upper], color="black", linewidth=1)
                ax.text(0.5, top + 0.02 * span, stars, ha="center", va="bottom")
                ax.set_ylim(lower, upper + 0.19 * span)
                ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1, 1))
            else:
                ax.legend(fontsize=8)
            fig.tight_layout()
            if save_results:
                safe_name = re.sub(r"[^\w-]+", "_", row.Cytokine).strip("_")
                path = output / f"Unpaired_{index:02d}_{safe_name}.png"
                fig.savefig(path, dpi=300, bbox_inches="tight")
                saved_paths.append(str(path))
            if display_results:
                plt.show()
            if save_results and save_clean_plots:
                ax.get_legend().remove()
                for collection in ax.collections:
                    if isinstance(collection, PathCollection):
                        collection.set_facecolor("black")
                fig.tight_layout()
                clean_output = output / "No_Legend_Black_Dots"
                clean_output.mkdir(parents=True, exist_ok=True)
                clean_path = clean_output / f"Unpaired_{index:02d}_{safe_name}_NoLegend_BlackDots.png"
                fig.savefig(clean_path, dpi=300, bbox_inches="tight")
                saved_paths.append(str(clean_path))
        finally:
            plt.close(fig)
    if save_results:
        workbook_path = output / "Unpaired_Cytokine_Results.xlsx"
        if workbook_path.resolve() == source:
            raise ValueError("Output workbook must not overwrite the source workbook.")
        with pd.ExcelWriter(workbook_path) as writer:
            statistics.to_excel(writer, sheet_name="Unpaired_Statistics", index=False)
            bio_points.to_excel(writer, sheet_name="Biological_Replicates", index=False)
            observations.to_excel(writer, sheet_name="Input_Observations", index=False)
            if exclusions:
                excluded_observations.to_excel(writer, sheet_name="Excluded_Observations", index=False)
            pd.DataFrame({"Setting": ["Source", "Sheet", "Test", "Scale", "FDR family", "Expected replicates per group", "Retained replicate counts", "Excluded replicates", "Excel data rows"],
                          "Value": [str(source), sheet_name, "Two-sided Mann-Whitney U; treatment first; exact permutations up to 100,000 allocations, otherwise tie/continuity-corrected asymptotic; see P_Method. Mean differences descriptive only; no mean-difference CI.",
                                    "Values as supplied; no transformation, normalization, or imputation",
                                    f"BH across {len(valid)} valid tests", str(expected_replicates), str(expected_counts), exclusion_note, f"{data_start_row}-{last_row}"]}).to_excel(writer, sheet_name="Analysis_Notes", index=False)
        saved_paths.append(str(workbook_path))
        print(f"Saved {len(saved_paths) - 1} PNGs and results workbook to {output}.")
    return bio_points, statistics, saved_paths


if __name__ == "__main__":
    unpaired_bio_points, unpaired_statistics, unpaired_saved_paths = Plot_Unpaired_Cytokine_Comparison(
        file_path="/Users/allisonpickle/Desktop/261007_Cytokine Reanalysis/261118_Compiled_Cytokines_3BioReps.xlsx",
        sheet_name="261018 Updated 48hr",
        control="6hr Ctrl",
        treatment="48hr 5ug",
        header_row=1,
        data_start_row=194,
        data_end_row=235,
        expected_replicates=3,
        cytokines=None,
        save_results=True,
        display_results=False,
        save_clean_plots=True,
        plot_group_labels=("6-hr Control", "48-hr 5 µg/mL LPS"),
        plot_value_label="log2 normalized concentration",
    )
    print(unpaired_statistics.to_string(index=False))
