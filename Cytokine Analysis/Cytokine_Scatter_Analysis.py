"""Cytokine scatter plots and exploratory independent-group time analyses."""


def Plot_Cytokine_Scatter_Plots(
    file_path,
    sheet_name=0,
    sample_col="B",
    rep_col="A",
    header_row=16,
    start_col="F",
    col_step=3,
    obs_conc_offset=1,
    data_start_row=None,
    data_end_row=None,
    output_dir=None,
    sample_order=None,
    show_mean=True,
    save_results=True,
    display_results=True,
    include_48hr=True,
    run_statistics=True,
    show_significance=True,
    lower_quantification_limits=None,
    log2_offset=0.5,
    samples=None,
    sheet_selections=None,
    include_all_samples=False,
):
    """Plot 5/100 ng/mL LPS time courses and optionally test time effects.

    48 hr is used only when present and can be explicitly skipped with
    ``include_48hr=False``. Ordinary one-way ANOVA uses all available biological
    replicate measurements independently at each time, with at least two per
    analyzed time. No replicate blocking or cross-time matching is performed.
    This assumes independent biological units across times and equal variances.
    OOR < values remain flagged and retain the legacy
    zero value for descriptive plots. For inference, censored values require
    an assay-specific LLOQ and use the explicit LLOQ/sqrt(2) substitution.
    Positive, skewed panels are transformed as log2(value + log2_offset), with
    one shared offset (default 0.5 concentration units, half an assumed
    1-unit reporting increment), never a cytokine-specific offset. Override
    that global assumption if the assay concentration units/resolution differ.
    Inference is exploratory at n=3.

    Returns the legacy tuple: parsed observations, summary, saved plot paths.
    PNG plots and the analysis workbook default to the input Excel file's
    containing folder; output_dir overrides this location.

    ``samples`` optionally selects a nonempty list of exact sample names before
    replicate averaging, summaries, statistics, and plotting. None includes all.
    Unknown names raise ValueError with available names in the selected rows.
    The numeric time-course eligibility rules still apply: controls and NP
    treatments are not plotted and excluded names are reported explicitly.
    This selects groups, not x-axis order; numeric time determines plot order.

    ``sheet_selections`` optionally supplies a nonempty list of dictionaries
    with sheet_name and overrides for header_row, data_start_row, data_end_row,
    start_col, and obs_conc_offset. Omitted settings inherit function defaults.
    Row bounds are inclusive Excel rows. All other layout settings are shared.
    Cytokines align by name; matching sample/replicate IDs across sheets are
    technical repeats, averaged before statistics and plotting. Combined parsed
    observations include Sheet and Excel_Row provenance. Overlapping ranges on
    the same sheet are rejected. None preserves the single-sheet interface.
    Returned observations also retain original sheet/row selections in
    attrs["source_selections"] for downstream output folder naming, including
    single-sheet input and rows excluded by sample/time eligibility filters.

    include_all_samples=True plots every selected sample with its own legend
    entry. Supported explicit hour labels determine placement; the last supported
    hour label wins when a pretreatment duration also appears. Names without any
    hour label are placed at 0 hr. Unsupported hours raise ValueError (48 hr is
    excluded when include_48hr=False). Only plain-LPS sample names contribute to
    concentration trends and statistical tests; NP treatments/controls remain
    separate descriptive points. The default retains time-course-only behavior.
    """
    from pathlib import Path
    import re

    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from matplotlib.lines import Line2D
    from openpyxl.utils import column_index_from_string
    from scipy import stats
    from scipy.interpolate import PchipInterpolator

    if not np.isfinite(log2_offset) or log2_offset <= 0:
        raise ValueError("log2_offset must be a finite positive concentration value.")
    if samples is not None:
        if not isinstance(samples, list) or not samples:
            raise ValueError("samples must be a nonempty list of sample names, or None.")
        if any(not isinstance(name, str) or not name.strip() for name in samples):
            raise ValueError("Each sample name must be a nonempty string.")
        samples = [name.strip() for name in samples]
        if len(set(samples)) != len(samples):
            raise ValueError("samples must not contain duplicate names.")
    lloq = {} if lower_quantification_limits is None else dict(lower_quantification_limits)
    for name, limit in lloq.items():
        if not np.isfinite(limit) or limit <= 0:
            raise ValueError(f"LLOQ for {name!r} must be finite and positive.")

    sample_idx = column_index_from_string(sample_col) - 1
    rep_idx = column_index_from_string(rep_col) - 1
    if not isinstance(col_step, int) or isinstance(col_step, bool) or col_step < 1:
        raise ValueError("col_step must be a positive integer.")
    multi_sheet = sheet_selections is not None
    selections = sheet_selections if multi_sheet else [{"sheet_name": sheet_name}]
    if not isinstance(selections, list) or not selections:
        raise ValueError("sheet_selections must be a nonempty list of dictionaries.")
    allowed_keys = {
        "sheet_name", "header_row", "data_start_row", "data_end_row",
        "start_col", "obs_conc_offset",
    }
    for selection in selections:
        if not isinstance(selection, dict) or "sheet_name" not in selection:
            raise ValueError("Each sheet selection must contain sheet_name.")
        if set(selection) - allowed_keys:
            raise ValueError(f"Unknown sheet selection keys: {set(selection) - allowed_keys}")

    def parse_sheet(raw, selected_sheet, selected_header, first_row, last_row,
                    selected_start, selected_offset):
        for label, value in (
            ("header_row", selected_header), ("data_start_row", first_row),
            ("data_end_row", last_row), ("obs_conc_offset", selected_offset),
        ):
            if not isinstance(value, int) or isinstance(value, bool):
                raise ValueError(f"{label} must be an integer in sheet {selected_sheet!r}.")
        if not (1 <= selected_header < first_row <= last_row <= len(raw)):
            raise ValueError(f"Select valid inclusive Excel data rows after the header in sheet {selected_sheet!r}.")
        if selected_offset < 0:
            raise ValueError("obs_conc_offset must be nonnegative.")
        if not isinstance(selected_start, str) or not re.fullmatch(r"[A-Za-z]+", selected_start):
            raise ValueError(f"start_col must be an Excel column letter in sheet {selected_sheet!r}.")
        start_idx = column_index_from_string(selected_start) - 1
        if max(sample_idx, rep_idx, start_idx) >= raw.shape[1]:
            raise ValueError(f"Selected columns are outside sheet {selected_sheet!r}.")
        headers = []
        for idx in range(start_idx, raw.shape[1], col_step):
            name = raw.iat[selected_header - 1, idx]
            if pd.notna(name) and str(name).strip():
                if multi_sheet and (
                    isinstance(name, (int, float, np.number))
                    or re.fullmatch(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", str(name).strip())
                ):
                    raise ValueError(f"Numeric cytokine header {name!r} in sheet {selected_sheet!r}; check header_row.")
                if idx + selected_offset >= raw.shape[1]:
                    raise ValueError(f"Missing concentration column for {name} in sheet {selected_sheet!r}.")
                headers.append((str(name).strip(), idx + selected_offset))
        if not headers or len({name for name, _ in headers}) != len(headers):
            raise ValueError(f"Expected unique cytokine headers in sheet {selected_sheet!r}.")
        sheet_records = []
        for row in range(first_row - 1, last_row):
            sample, rep = raw.iat[row, sample_idx], raw.iat[row, rep_idx]
            if pd.isna(sample) or not str(sample).strip():
                continue
            if pd.isna(rep) or not str(rep).strip():
                raise ValueError(f"Missing biological-replicate ID in sheet {selected_sheet!r} at Excel row {row + 1}.")
            for cytokine, col in headers:
                original = raw.iat[row, col]
                below = isinstance(original, str) and bool(
                    re.match(r"^OOR\s*<", original.strip(), re.IGNORECASE)
                )
                if below:
                    observed = 0.0
                    statistical = lloq[cytokine] / np.sqrt(2) if cytokine in lloq else np.nan
                elif pd.isna(original) or not str(original).strip():
                    observed = statistical = np.nan
                else:
                    try:
                        observed = float(original)
                    except (TypeError, ValueError) as exc:
                        raise ValueError(
                            f"Invalid {cytokine} value in sheet {selected_sheet!r} at Excel row {row + 1}: {original!r}"
                        ) from exc
                    if not np.isfinite(observed):
                        raise ValueError(f"Nonfinite {cytokine} concentration in sheet {selected_sheet!r} at Excel row {row + 1}.")
                    statistical = observed
                record = {
                    "Sample": str(sample).strip(), "Rep": str(rep).strip(),
                    "Cytokine": cytokine, "Obs_Conc": observed,
                    "BelowDetection": below, "Statistical_Conc": statistical,
                }
                if multi_sheet:
                    record.update(Sheet=selected_sheet, Excel_Row=row + 1)
                sheet_records.append(record)
        if not sheet_records:
            raise ValueError(f"No observations parsed from sheet {selected_sheet!r}, rows {first_row}-{last_row}.")
        return sheet_records, headers

    cytokines = []
    names = set()
    records = []
    selected_rows = set()
    source_ranges = []
    source_selections = []
    with pd.ExcelFile(file_path) as workbook:
        sheet_cache = {}
        for selection in selections:
            selected_sheet = selection["sheet_name"]
            if isinstance(selected_sheet, int) and not isinstance(selected_sheet, bool):
                if not 0 <= selected_sheet < len(workbook.sheet_names):
                    raise ValueError(f"Invalid sheet index: {selected_sheet}")
                selected_sheet = workbook.sheet_names[selected_sheet]
            if not isinstance(selected_sheet, str) or selected_sheet not in workbook.sheet_names:
                raise ValueError(f"Unknown sheet {selected_sheet!r}; available: {workbook.sheet_names}")
            if selected_sheet not in sheet_cache:
                sheet_cache[selected_sheet] = pd.read_excel(workbook, sheet_name=selected_sheet, header=None)
            raw = sheet_cache[selected_sheet]
            selected_header = selection.get("header_row", header_row)
            first_row = selection.get("data_start_row", data_start_row)
            last_row = selection.get("data_end_row", data_end_row)
            if not isinstance(selected_header, int) or isinstance(selected_header, bool):
                raise ValueError(f"header_row must be an integer in sheet {selected_sheet!r}.")
            first_row = selected_header + 1 if first_row is None else first_row
            last_row = len(raw) if last_row is None else last_row
            sheet_records, headers = parse_sheet(
                raw, selected_sheet, selected_header, first_row, last_row,
                selection.get("start_col", start_col),
                selection.get("obs_conc_offset", obs_conc_offset),
            )
            row_keys = {(selected_sheet, row) for row in range(first_row, last_row + 1)}
            if selected_rows & row_keys:
                raise ValueError(f"Overlapping row selections in sheet {selected_sheet!r}.")
            selected_rows.update(row_keys)
            records.extend(sheet_records)
            for name, col in headers:
                if name not in names:
                    cytokines.append((name, col))
                    names.add(name)
            source_ranges.append(f"{selected_sheet!r} rows {first_row}-{last_row}")
            source_selections.append({
                "sheet_name": selected_sheet,
                "data_start_row": first_row,
                "data_end_row": last_row,
            })
    unknown_lloq = set(lloq) - names
    if unknown_lloq:
        raise ValueError(f"LLOQ values provided for unknown cytokines: {sorted(unknown_lloq)}")
    long_df = pd.DataFrame(records)
    long_df.attrs["source_selections"] = source_selections
    if samples is not None:
        available = sorted(long_df["Sample"].unique())
        unknown = [name for name in samples if name not in available]
        if unknown:
            raise ValueError(f"Unknown samples: {unknown}. Available samples: {available}")
        long_df = long_df.loc[long_df["Sample"].isin(samples)].copy()
        print(f"Analyzing selected samples: {samples}")
    labels = long_df["Sample"].str.extract(
        r"^\s*(\d+)\s*(?:hr|h)\s*[_ ]+([\d.]+)\s*ng(?:/mL)?(?:\s*LPS)?\s*$",
        flags=re.IGNORECASE,
    )
    long_df["Time_hr"] = pd.to_numeric(labels[0], errors="coerce")
    long_df["LPS_ng_mL"] = pd.to_numeric(labels[1], errors="coerce")
    if include_all_samples:
        long_df["Plain_LPS"] = labels[0].notna() & labels[1].notna()
        def sample_time(name):
            hours = [int(value) for value in re.findall(r"(?<!\d)(\d+)\s*(?:hrs?|hours?|h)(?![A-Za-z])", name, re.I)]
            supported = [value for value in hours if value in {6, 24, 48, 72}]
            if supported:
                return supported[-1]
            if hours:
                raise ValueError(f"No supported 6/24/48/72 hr label in sample {name!r}.")
            return 0
        time_by_sample = {name: sample_time(name) for name in long_df["Sample"].unique()}
        long_df["Time_hr"] = long_df["Sample"].map(time_by_sample)
        doses = long_df["Sample"].str.extract(r"(?<![\d.])(5|100)\s*ng\b", flags=re.I)[0]
        long_df["LPS_ng_mL"] = pd.to_numeric(doses, errors="coerce")
        no_time = [name for name, time in time_by_sample.items() if time == 0]
        if no_time:
            print(f"Samples without hour labels placed at 0 hr: {no_time}")
    requested_times = [6, 24, 48, 72] if include_48hr else [6, 24, 72]
    concentrations = [5, 100]
    selected = long_df["Time_hr"].isin(requested_times) & long_df["LPS_ng_mL"].isin(concentrations)
    if include_all_samples:
        selected = long_df["Time_hr"].isin([0, *requested_times])
    if samples is not None:
        excluded_names = list(dict.fromkeys(long_df.loc[~selected, "Sample"]))
        if excluded_names:
            print(f"Selected samples excluded by time-course eligibility rules: {excluded_names}")
    print(f"Excluded {int((~selected).sum())} rows by {'hour' if include_all_samples else '5/100 ng/mL time-course'} eligibility.")
    long_df = long_df.loc[selected].copy()
    if long_df.empty:
        raise ValueError("No samples remain after time eligibility filtering.")
    times = [time for time in ([0, *requested_times] if include_all_samples else requested_times) if long_df["Time_hr"].eq(time).any()]
    absent = [time for time in requested_times if time not in times]
    if absent:
        print(f"Time point(s) {absent} hr not provided; skipped.")
    valid = long_df.dropna(subset=["Obs_Conc"])
    if valid.empty:
        raise ValueError("No numeric concentrations in the selected rows.")
    if samples is not None:
        missing = sorted(set(long_df["Sample"]) - set(valid["Sample"]))
        if missing:
            raise ValueError(f"Selected samples have no numeric concentrations: {missing}")
    print(
        f"Reading {'; '.join(source_ranges)}; "
        f"{int(long_df['BelowDetection'].sum())} below-detection observations retained and flagged."
    )
    print(f"Missing concentrations omitted from descriptive summaries: {int(long_df['Obs_Conc'].isna().sum())}.")

    group_keys = ["Cytokine", "Time_hr", "LPS_ng_mL", "Rep"]
    if multi_sheet or include_all_samples:
        group_keys.insert(1, "Sample")
    if include_all_samples:
        group_keys.append("Plain_LPS")
    points = valid.groupby(group_keys, sort=False, as_index=False, dropna=False).agg(
        Obs_Conc=("Obs_Conc", "mean"),
        Statistical_Conc=(
            "Statistical_Conc",
            lambda values: values.mean() if values.notna().all() else np.nan,
        ),
        TechnicalRows=("Obs_Conc", "count"),
        BelowDetection=("BelowDetection", "any"),
    )
    if (points["TechnicalRows"] > 1).any():
        print("Technical repeats averaged within each biological replicate.")
    summary_keys = ["Cytokine", "Time_hr", "LPS_ng_mL"]
    if include_all_samples:
        summary_keys.insert(1, "Sample")
    summary = points.groupby(
        summary_keys, sort=False, as_index=False, dropna=False
    ).agg(mean=("Obs_Conc", "mean"), std=("Obs_Conc", "std"), n=("Obs_Conc", "count"))
    detection = long_df.groupby(
        summary_keys if include_all_samples else ["Cytokine", "LPS_ng_mL", "Time_hr"],
        sort=False, as_index=False, dropna=False
    ).agg(
        Number_of_Observations=("Obs_Conc", "count"),
        Number_Below_Detection=("BelowDetection", "sum"),
    )
    detection["Percent_Below_Detection"] = np.where(
        detection["Number_of_Observations"] > 0,
        100 * detection["Number_Below_Detection"] / detection["Number_of_Observations"],
        np.nan,
    )
    if sample_order is not None:
        raise ValueError("sample_order is not used for numeric time-trend plots; leave it as None.")

    alpha = 0.05
    comparisons = [(6, 24), (24, 48), (48, 72), (6, 48), (6, 72), (24, 72)]

    def p_label(p_value):
        if pd.isna(p_value):
            return "not tested"
        if p_value < 0.0001:
            return "****"
        if p_value < 0.001:
            return "***"
        if p_value < 0.01:
            return "**"
        if p_value < alpha:
            return "*"
        return "n.s."

    overall_rows, pair_rows = [], []
    for cytokine, _ in cytokines:
        for concentration in concentrations:
            if not run_statistics:
                overall_rows.append({
                    "Cytokine": cytokine, "LPS_ng_mL": concentration,
                    "Biological_Replicates": 0, "F_statistic": np.nan,
                    "Total_Observations": 0, "Replicates_Per_Time": "",
                    "Numerator_df": np.nan, "Denominator_df": np.nan,
                    "Raw_p_value": np.nan, "Significance": "not tested",
                    "Model_Status_or_Warnings": "Statistical testing disabled.",
                })
                for a, b in comparisons:
                    pair_rows.append({
                        "Cytokine": cytokine, "LPS_ng_mL": concentration,
                        "Time_Comparison": f"{a} vs {b} hr",
                        "Earlier_Time_n": 0, "Later_Time_n": 0,
                        "Estimated_Difference": np.nan, "Standard_Error": np.nan,
                        "CI_Lower": np.nan, "CI_Upper": np.nan,
                        "Test_Statistic": np.nan, "Degrees_of_Freedom": np.nan,
                        "Raw_p_value": np.nan, "Holm_adjusted_p_value": np.nan,
                        "Significance": "not tested", "Comparison_Status": "Statistical testing disabled.",
                    })
                continue

            observed = long_df[
                long_df["Cytokine"].eq(cytokine) & long_df["LPS_ng_mL"].eq(concentration)
            ]
            subset = points[
                points["Cytokine"].eq(cytokine)
                & points["LPS_ng_mL"].eq(concentration)
                & points["Time_hr"].isin(times)
            ]
            if include_all_samples:
                observed = observed.loc[observed["Plain_LPS"]]
                subset = subset.loc[subset["Plain_LPS"]]
            groups = {
                time: subset.loc[
                    subset["Time_hr"].eq(time), "Statistical_Conc"
                ].dropna().to_numpy(dtype=float)
                for time in times
            }
            counts = {time: len(values) for time, values in groups.items()}
            count_description = "; ".join(f"{time} hr: n={counts[time]}" for time in times)
            total_n = sum(counts.values())
            has_censored = bool(observed["BelowDetection"].any())
            if len(times) < 2:
                status = "Not run: fewer than two requested time points are available."
            elif has_censored and cytokine not in lloq:
                status = "Not run: censored values; supply this cytokine's LLOQ. No values imputed or excluded."
            elif any(n < 2 for n in counts.values()):
                status = "Not run: fewer than two biological replicate measurements at an analyzed time."
            else:
                status = "Exploratory ordinary one-way ANOVA; independent time groups, no replicate block."
            status += f" {count_description}."
            print(f"{cytokine}, {concentration} ng/mL: {status}")
            if absent:
                status += f" Times analyzed: {times} hr; unavailable/skipped: {absent} hr."

            f_value = df_num = df_den = raw_p = np.nan
            transformed = False
            mse = np.nan
            if status.startswith("Exploratory"):
                pooled = np.concatenate(list(groups.values()))
                skew = 0.0
                if np.all(pooled >= 0) and np.ptp(pooled) > 0:
                    skew = stats.skew(pooled, bias=False) if pooled.size >= 3 else 0.0
                    if np.isfinite(skew) and skew > 1:
                        groups = {
                            time: np.log2(values + log2_offset)
                            for time, values in groups.items()
                        }
                        transformed = True
                df_num, df_den = len(times) - 1, total_n - len(times)
                ss_error = sum(
                    np.sum((values - values.mean()) ** 2) for values in groups.values()
                )
                mse = ss_error / df_den
                if mse > 0:
                    result = stats.f_oneway(*groups.values())
                    f_value, raw_p = float(result.statistic), float(result.pvalue)
                else:
                    status = f"Model fitting issue: within-group variance is zero. {count_description}."
                if transformed:
                    status += f" log2(value + {log2_offset:g}) used (shared offset; skew={skew:.3g})."
                if any(n < 3 for n in counts.values()):
                    status += " At least one time has fewer than the planned three biological replicates."
            overall_rows.append({
                "Cytokine": cytokine, "LPS_ng_mL": concentration,
                "Biological_Replicates": subset.loc[
                    subset["Statistical_Conc"].notna(), "Rep"
                ].nunique(),
                "Total_Observations": total_n, "Replicates_Per_Time": count_description,
                "F_statistic": f_value,
                "Numerator_df": df_num, "Denominator_df": df_den,
                "Raw_p_value": raw_p, "Significance": p_label(raw_p),
                "Model_Status_or_Warnings": status,
            })

            contrast_p = {}
            for a, b in comparisons:
                row = {
                    "Cytokine": cytokine, "LPS_ng_mL": concentration,
                    "Time_Comparison": f"{a} vs {b} hr",
                    "Earlier_Time_n": counts.get(a, 0),
                    "Later_Time_n": counts.get(b, 0),
                    "Estimated_Difference": np.nan, "Standard_Error": np.nan,
                    "CI_Lower": np.nan, "CI_Upper": np.nan,
                    "Test_Statistic": np.nan, "Degrees_of_Freedom": np.nan,
                    "Raw_p_value": np.nan, "Holm_adjusted_p_value": np.nan,
                    "Significance": "not tested", "Comparison_Status": "",
                }
                if a not in times or b not in times:
                    row["Comparison_Status"] = "Not run: a time point was not provided or was skipped."
                elif has_censored and cytokine not in lloq:
                    row["Comparison_Status"] = "Not run: censored observations require a supplied LLOQ."
                elif not status.startswith("Exploratory"):
                    row["Comparison_Status"] = f"Not run: one-way ANOVA unavailable. {status}"
                else:
                    difference = groups[b].mean() - groups[a].mean()
                    se = np.sqrt(mse * (1 / counts[a] + 1 / counts[b]))
                    df = df_den
                    statistic = difference / se if se > 0 else (np.inf if difference else 0.0)
                    p_value = 2 * stats.t.sf(abs(statistic), df) if se > 0 else (0.0 if difference else 1.0)
                    critical = stats.t.ppf(1 - alpha / 2, df)
                    row.update({
                        "Estimated_Difference": difference, "Standard_Error": se,
                        "CI_Lower": difference - critical * se, "CI_Upper": difference + critical * se,
                        "Test_Statistic": statistic, "Degrees_of_Freedom": df,
                        "Raw_p_value": p_value,
                        "Comparison_Status": (
                            f"Model-based independent-group contrast on "
                            f"{'log2-transformed concentration' if transformed else 'raw concentration'} scale; "
                            "later time minus earlier time; pooled one-way ANOVA residual variance."
                        ),
                    })
                    contrast_p[(a, b)] = p_value
                pair_rows.append(row)

            # Holm correction is scoped to the comparisons that could actually be tested.
            ordered = sorted(contrast_p, key=contrast_p.get)
            running = 0.0
            for rank, pair in enumerate(ordered):
                adjusted = min(1.0, (len(ordered) - rank) * contrast_p[pair])
                running = max(running, adjusted)
                target = f"{pair[0]} vs {pair[1]} hr"
                result = next(
                    row for row in reversed(pair_rows)
                    if row["Cytokine"] == cytokine
                    and row["LPS_ng_mL"] == concentration
                    and row["Time_Comparison"] == target
                )
                result["Holm_adjusted_p_value"] = running
                result["Significance"] = p_label(running)

    overall = pd.DataFrame(overall_rows)
    pairwise = pd.DataFrame(pair_rows)
    colors = {5: "#E65400", 100: "#0072B2"}
    rep_labels = sorted(points["Rep"].unique())
    rep_colors = dict(zip(
        rep_labels, plt.get_cmap("hsv")(np.linspace(0, 1, len(rep_labels), endpoint=False))
    ))
    rep_colors.update({"Rep 1": "#E74C3C", "Rep 2": "#2ECC71", "Rep 3": "#9B59B6"})
    sample_labels = list(dict.fromkeys(long_df["Sample"]))
    sample_colors = dict(zip(
        sample_labels, plt.get_cmap("tab20")(np.linspace(0, 1, len(sample_labels)))
    ))
    output = Path(output_dir) if output_dir else Path(file_path).resolve().parent
    if save_results:
        output.mkdir(parents=True, exist_ok=True)
    saved_paths = []
    rng = np.random.default_rng(42)
    for index, (cytokine, _) in enumerate(cytokines, start=1):
        cyto_points = points[points["Cytokine"].eq(cytokine)]
        if cyto_points.empty:
            print(f"No numeric data for {cytokine}; no plot generated.")
            continue
        fig, ax = plt.subplots(figsize=(8, 6))
        if include_all_samples:
            for sample in sample_labels:
                sample_points = cyto_points.loc[cyto_points["Sample"].eq(sample)]
                if sample_points.empty:
                    continue
                ax.scatter(
                    sample_points["Time_hr"] + rng.uniform(-0.6, 0.6, len(sample_points)),
                    sample_points["Obs_Conc"], color=sample_colors[sample],
                    s=55, edgecolors="white", linewidths=0.6, zorder=3,
                )
                if show_mean:
                    ax.scatter(
                        [sample_points["Time_hr"].iloc[0]], [sample_points["Obs_Conc"].mean()],
                        color=sample_colors[sample], marker="D", s=65,
                        edgecolors="black", zorder=4,
                    )
        for concentration in concentrations:
            subset = cyto_points[cyto_points["LPS_ng_mL"].eq(concentration)]
            if include_all_samples:
                subset = subset.loc[subset["Plain_LPS"]]
            if subset.empty:
                continue
            jitter = rng.uniform(-0.6, 0.6, len(subset))
            if not include_all_samples:
                ax.scatter(
                    subset["Time_hr"] + jitter, subset["Obs_Conc"],
                    color=subset["Rep"].map(rep_colors).tolist(), s=55,
                    edgecolors="white", linewidths=0.6, zorder=3,
                )
            means = subset.groupby("Time_hr")["Obs_Conc"].mean().sort_index()
            color = colors[concentration]
            if show_mean and not include_all_samples:
                ax.scatter(means.index, means.values, color=color, marker="D",
                           s=65, edgecolors="black", zorder=4)
            if len(means) >= 3:
                grid = np.linspace(means.index.min(), means.index.max(), 200)
                ax.plot(grid, PchipInterpolator(means.index, means.values, extrapolate=False)(grid),
                        color=color, linewidth=2)
            elif len(means) == 2:
                ax.plot(means.index, means.values, color=color, linewidth=2)

        if run_statistics and show_significance:
            values = cyto_points["Obs_Conc"].dropna()
            ymin, ymax = (float(values.min()), float(values.max())) if not values.empty else (0.0, 1.0)
            span = max(ymax - ymin, abs(ymax) * 0.1, 1e-6)
            level = 0
            for concentration in concentrations:
                for a, b in zip(times, times[1:]):
                    rows = pairwise[
                        pairwise["Cytokine"].eq(cytokine)
                        & pairwise["LPS_ng_mL"].eq(concentration)
                        & pairwise["Time_Comparison"].eq(f"{a} vs {b} hr")
                    ]
                    if rows.empty or pd.isna(rows.iloc[0]["Holm_adjusted_p_value"]):
                        continue
                    y = ymax + span * (0.08 + 0.09 * level)
                    h = span * 0.025
                    color = colors[concentration]
                    ax.plot([a, a, b, b], [y, y + h, y + h, y], color=color, linewidth=1.3)
                    ax.text((a + b) / 2, y + h, rows.iloc[0]["Significance"],
                            color=color, ha="center", va="bottom", fontsize=9)
                    level += 1
            if level:
                ax.set_ylim(top=ymax + span * (0.15 + 0.11 * level))

        if include_all_samples:
            sample_handles = [
                Line2D([0], [0], linestyle="none", marker="o",
                       markerfacecolor=sample_colors[sample], markeredgecolor="white",
                       label=sample + (" (no numeric data)" if not cyto_points["Sample"].eq(sample).any() else ""))
                for sample in sample_labels
            ]
            ax.legend(handles=sample_handles, title="Sample", loc="upper left",
                      bbox_to_anchor=(1.02, 1), fontsize=9)
        concentration_handles = [
            Line2D([0], [0], color=colors[c], linewidth=2,
                   marker="D" if show_mean else None, label=f"{c} ng/mL")
            for c in concentrations if c in set(cyto_points["LPS_ng_mL"])
        ]
        concentration_legend = None if include_all_samples else ax.legend(
            handles=concentration_handles, title="LPS concentration",
            loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=9,
        )
        if concentration_legend is not None:
            ax.add_artist(concentration_legend)
        rep_handles = [
            Line2D([0], [0], linestyle="none", marker="o",
                   markerfacecolor=rep_colors[rep], markeredgecolor="white",
                   markersize=8, label=rep)
            for rep in rep_labels if rep in set(cyto_points["Rep"])
        ]
        if not include_all_samples:
            ax.legend(handles=rep_handles, title="Bio Rep", loc="upper left",
                      bbox_to_anchor=(1.02, 0.65), fontsize=9)
        ax.set_xticks(times)
        ax.set_xticklabels([f"{time} hr" for time in times])
        ax.set_xlim(min(times) - 3, max(times) + 3)
        ax.set(xlabel="Exposure time", ylabel="Observed concentration", title=cytokine)
        ax.grid(axis="y", alpha=0.3, linestyle="--")
        fig.tight_layout()
        if save_results:
            safe_name = re.sub(r"[^\w-]+", "_", cytokine, flags=re.UNICODE).strip("_")
            path = output / f"{index:02d}_{safe_name}.png"
            fig.savefig(path, dpi=300)
            saved_paths.append(str(path))
        if display_results:
            plt.show()
        plt.close(fig)

    if save_results:
        with pd.ExcelWriter(output / "Cytokine_Scatter_Data.xlsx") as writer:
            long_df.to_excel(writer, sheet_name="Parsed Observations", index=False)
            points.to_excel(writer, sheet_name="Biological Replicate Points", index=False)
            summary.to_excel(writer, sheet_name="Summary", index=False)
            overall.to_excel(writer, sheet_name="Overall Time ANOVA", index=False)
            pairwise.to_excel(writer, sheet_name="Pairwise Time Comparisons", index=False)
            detection.to_excel(writer, sheet_name="Detection Limit Summary", index=False)
        print(f"Saved {len(saved_paths)} scatter plots and the data workbook to {output}.")
    return long_df, summary, saved_paths
