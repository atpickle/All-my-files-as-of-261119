"""Paired qPCR comparison after 6 hours at 50 ng/mL LPS."""

from argparse import ArgumentParser
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import permutation_test, ttest_ind, ttest_rel

from Cytokine_Heatmap_Statistics import _holm


def mRNA_50ngLPS_comparison(file_path, output_dir=None):
    """Export paired qPCR plots, mean ± SD bars, and statistics; return summary and plate pairs.

    P2 and P4 are treated as independent experimental units, as assumed in
    the workbook. Measurements within each plate and condition are averaged.
    The two-sided exact sign-flip and paired t-tests compare LPS minus control
    on the supplied normalized scale. Holm adjustment is applied separately
    to each test across the three targets.
    """
    file_path = Path(file_path)
    output = Path(output_dir) if output_dir is not None else file_path.parent / "mRNA_50ngLPS_Comparison"
    targets = ("TNF-α", "IL-1β", "IL-6")
    conditions = ["Ctrl", "LPS 50ng/mL"]
    values = ["Value 1", "Value 2", "Value 3"]
    data = pd.read_excel(file_path, sheet_name="Prism Data", header=2)
    data = data.loc[data["Target / measure"].isin(targets)].copy()
    if data["Plate"].isna().any() or not data["Condition"].isin(conditions).all():
        raise ValueError("Selected qPCR rows must have plate IDs and valid control or LPS conditions.")
    measurements = data[values].apply(pd.to_numeric, errors="raise")
    if not np.isfinite(measurements.to_numpy()).all() or measurements.le(0).any().any():
        raise ValueError("Selected qPCR rows must have three positive, finite measurements.")
    data["Plate mean"] = measurements.mean(axis=1)
    paired = data.pivot(index=["Target / measure", "Plate"], columns="Condition", values="Plate mean")
    paired = paired.reindex(columns=conditions).rename(columns=dict(zip(conditions, ["Control", "LPS"])))
    if paired.isna().any().any() or set(paired.index.get_level_values(0)) != set(targets):
        raise ValueError("Each selected target and plate must have both control and LPS measurements.")
    paired["LPS minus control"] = paired["LPS"] - paired["Control"]
    paired["LPS / control"] = paired["LPS"] / paired["Control"]

    summaries = []
    for target in targets:
        group = paired.loc[target]
        differences = group["LPS minus control"].to_numpy()
        if len(differences) < 2:
            raise ValueError(f"{target}: at least two independent plate pairs are required.")
        test = permutation_test((differences,), np.mean, permutation_type="samples",
                                n_resamples=np.inf, alternative="two-sided")
        t_test = ttest_rel(group["LPS"], group["Control"]) if np.ptp(differences) > 0 else None
        summaries.append({
            "Target": target,
            "Plate pairs": len(group),
            "Pairs with increase": int(np.count_nonzero(differences > 0)),
            "Mean paired difference": differences.mean(),
            "Geometric mean LPS / control": np.exp(np.log(group["LPS / control"]).mean()),
            "Exact paired p": float(test.pvalue),
            "Paired t": t_test.statistic if t_test is not None else np.nan,
            "Paired t df": len(differences) - 1,
            "Paired t p": t_test.pvalue if t_test is not None else np.nan,
        })
    summary = pd.DataFrame(summaries)
    summary["Holm-adjusted p"] = _holm(summary["Exact paired p"].to_numpy())
    t_p = summary["Paired t p"]
    summary["Paired t Holm-adjusted p"] = pd.Series(
        _holm(t_p.fillna(1).to_numpy()), index=summary.index,
    ).where(t_p.notna())
    output.mkdir(parents=True, exist_ok=True)
    notes = [
        "Exposure: 50 ng/mL LPS for 6 h. Values are normalized qPCR measurements, not protein concentrations.",
        "Plate independence is assumed. The three values per plate and condition were averaged.",
        "Two-sided exact sign-flip test on within-plate differences; Holm adjustment across three targets.",
        "With two independent pairs, the smallest attainable two-sided exact p-value is 0.50.",
        "Two-sided paired t-tests use the same plate differences. With two pairs (df = 1), "
        "normality cannot be meaningfully assessed; t-test results are exploratory.",
        "A t-test is undefined when paired differences have zero variance. Undefined tests "
        "are left blank and counted in the Holm family.",
    ]
    with pd.ExcelWriter(output / "mRNA_50ngLPS_Results.xlsx") as writer:
        summary.to_excel(writer, sheet_name="Summary", index=False)
        paired.reset_index().to_excel(writer, sheet_name="Plate Pairs", index=False)
        pd.DataFrame({"Analysis notes": notes}).to_excel(writer, sheet_name="Notes", index=False)

    fig, axes = plt.subplots(1, len(targets), figsize=(10, 4), layout="constrained")
    try:
        colors = {plate: f"C{i % 10}" for i, plate in enumerate(paired.index.get_level_values("Plate").unique())}
        for ax, result in zip(axes, summaries):
            target = result["Target"]
            for plate, row in paired.loc[target].iterrows():
                ax.plot([0, 1], [row["Control"], row["LPS"]], "o-", color=colors[plate], label=plate)
            ax.set(xticks=[0, 1], xticklabels=["Control", "LPS"],
                   title=f"{target}\nHolm-adjusted p = {summary.loc[summary['Target'].eq(target), 'Holm-adjusted p'].iloc[0]:.2f}")
            ax.legend(title="Plate", frameon=False, loc="upper left")
        axes[0].set_ylabel("Normalized mRNA expression")
        fig.suptitle("Control vs. 50 ng/mL LPS at 6 h")
        fig.savefig(output / "mRNA_50ngLPS_Paired_Plots.png", dpi=300)
    finally:
        plt.close(fig)

    fig, axes = plt.subplots(1, len(targets), figsize=(10, 4), layout="constrained")
    try:
        for ax, result in zip(axes, summary.to_dict("records")):
            target = result["Target"]
            group = paired.loc[target, ["Control", "LPS"]]
            x = np.arange(2)
            ax.bar(x, group.mean(), yerr=group.std(ddof=1), width=0.58,
                   color=["#dddddd", "#7fa7cd"], edgecolor="black", capsize=4)
            for offset, (plate, row) in zip(np.linspace(-0.09, 0.09, len(group)), group.iterrows()):
                ax.scatter(x + offset, row, color=colors[plate], edgecolors="black",
                           s=32, zorder=3, label=plate)
            ax.set(xticks=x, xticklabels=["Control", "50 ng/mL LPS"], ylim=(0, None),
                   title=f"{target}\nHolm-adjusted p = {result['Holm-adjusted p']:.2f}")
            ax.spines[["top", "right"]].set_visible(False)
            ax.legend(title="Plate", frameon=False, loc="upper left")
        axes[0].set_ylabel("Normalized mRNA expression")
        fig.suptitle("6 h qPCR: mean ± SD across plate means")
        fig.savefig(output / "mRNA_50ngLPS_Bar_Plots.png", dpi=300)
    finally:
        plt.close(fig)
    return summary, paired.reset_index()


def plot_p2_normalized(file_path, output_dir=None, *, annotate_welch=False):
    """Plot P2 measurements, optionally marking technical-replicate Welch tests with stars."""
    file_path = Path(file_path)
    output = Path(output_dir) if output_dir is not None else file_path.parent / "mRNA_50ngLPS_Comparison"
    sheet = pd.read_excel(file_path, sheet_name="P2 Normalized", header=None)
    targets = ("TNF-α", "IL-1β", "IL-6")
    groups = {}
    for row, target, header in ((10, "TNF-α", "tnf-a"), (5, "IL-1β", "il-1b"), (0, "IL-6", "il-6")):
        labels = [str(value).strip().lower() for value in sheet.iloc[row, :3]]
        if labels[1:] != ["ctrl", "50 nglps"] or (
            labels[0] != header and not (target == "TNF-α" and labels[0] == "il-1b")
        ):
            raise ValueError(f"Unexpected P2 Normalized header in Excel row {row + 1}: {labels}")
        values = sheet.iloc[row + 1:row + 4, 1:3].to_numpy(dtype=float)
        if values.shape != (3, 2) or not np.isfinite(values).all():
            raise ValueError(f"{target}: expected three finite P2 values per condition.")
        if target == "TNF-α" and labels[0] == "il-1b":
            # Confirm the duplicated header is a mislabel before plotting it as TNF-α.
            source = pd.read_excel(file_path, sheet_name="Prism Data", header=2)
            reference = source.loc[source["Plate"].eq("P2") & source["Target / measure"].eq(target)]
            reference = reference.set_index("Condition").reindex(["Ctrl", "LPS 50ng/mL"])
            expected = reference[["Value 1", "Value 2", "Value 3"]].to_numpy(dtype=float).T
            if not np.allclose(values, expected, rtol=0, atol=0.001):
                raise ValueError("The second il-1b block does not match the P2 TNF-α source values.")
        groups[target] = values

    output.mkdir(parents=True, exist_ok=True)
    suffix = "_Welch" if annotate_welch else ""
    path = output / f"mRNA_50ngLPS_P2_Normalized_Bar_Plots{suffix}.png"
    fig, axes = plt.subplots(1, len(targets), figsize=(10, 4), layout="constrained")
    try:
        x = np.arange(2)
        for ax, target in zip(axes, targets):
            values = groups[target]
            ax.bar(x, values.mean(axis=0), yerr=values.std(axis=0, ddof=1), width=0.58,
                   color=["#dddddd", "#7fa7cd"], edgecolor="black", capsize=4)
            ax.scatter((x + np.linspace(-0.09, 0.09, len(values))[:, None]).ravel(),
                       values.ravel(), color="black", s=18, zorder=3)
            ax.set(xticks=x, xticklabels=["Control", "50 ng/mL LPS"], ylim=(0, None), title=target)
            ax.spines[["top", "right"]].set_visible(False)
            if annotate_welch:
                test = ttest_ind(values[:, 1], values[:, 0], equal_var=False, alternative="two-sided")
                if not np.isfinite(test.pvalue):
                    raise ValueError(f"{target}: Welch's t-test did not return a finite p-value.")
                ymax = ax.get_ylim()[1]
                ax.set_ylim(top=1.28 * ymax)
                ax.plot([0, 0, 1, 1], [ymax, 1.04 * ymax, 1.04 * ymax, ymax],
                        color="black", linewidth=1)
                stars = "*" * sum(test.pvalue < cutoff for cutoff in (0.05, 0.01, 0.001, 0.0001))
                ax.text(0.5, 1.07 * ymax, stars or "ns", ha="center", va="bottom")
        axes[0].set_ylabel("Normalized mRNA expression")
        title = "P2 normalized qPCR at 6 h: mean ± SD within one plate"
        if annotate_welch:
            title += "\nTwo-sided Welch tests, unadjusted p-values; 3 technical repeats per group"
            title += "\n* p < 0.05; ** p < 0.01; *** p < 0.001; **** p < 0.0001"
        fig.suptitle(title)
        fig.savefig(path, dpi=300)
    finally:
        plt.close(fig)
    return path


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("workbook", nargs="?", type=Path, default=Path.home() / "Desktop" / "LPS50_vs_Ctrl_qPCR.xlsx")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--sheet", choices=["P2 Normalized", "Prism Data"], default="P2 Normalized",
                        help="Plot P2 only (default), or run the original paired-plate analysis.")
    parser.add_argument("--annotate-welch", action="store_true",
                        help="Export a P2-only bar plot with stars for technical-replicate Welch tests.")
    args = parser.parse_args()
    if args.annotate_welch and args.sheet != "P2 Normalized":
        parser.error("--annotate-welch requires the P2 Normalized sheet.")
    if args.sheet == "P2 Normalized":
        print(plot_p2_normalized(args.workbook, args.output_dir, annotate_welch=args.annotate_welch))
    else:
        results, _ = mRNA_50ngLPS_comparison(args.workbook, args.output_dir)
        print(results.to_string(index=False))
