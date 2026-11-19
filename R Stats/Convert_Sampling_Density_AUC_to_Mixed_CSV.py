"""Add mixed-model identifiers to long sampling-density AUC worksheets.

Input columns: Key, Step Size, AUC (raw AUC, not log-transformed AUC).
Output: Key, Treatment, Biological_Rep, Well, Image, StepSize, Results.
StepSize values: Continuous, S2.5, S5.

Well mapping follows the original converter: filename positions 1-4 -> well 1,
5-8 -> well 2, 9-12 -> well 3, repeating every 12 positions.
Biological replicate is supplied explicitly, not inferred from Key.
"""

import os
from collections.abc import Mapping
from numbers import Integral

import numpy as np
import pandas as pd


def _unique_output_path(path):
    if not os.path.exists(path):
        return path
    base, ext = os.path.splitext(path)
    counter = 1
    while os.path.exists(f"{base} ({counter}){ext}"):
        counter += 1
    return f"{base} ({counter}){ext}"


def _sanitize_filename_part(value):
    value = str(value).strip()
    for char in '/\\:*?"<>|':
        value = value.replace(char, "_")
    return " ".join(value.split())


def split_sampling_density_sheet_to_mixed_csv(
    df,
    sheet_name,
    input_excel_path,
    key_col="Key",
    stepsize_col="Step Size",
    auc_col="AUC",
    overwrite=False,
    biological_rep=None,
    output_dir=None,
    include_sheet_name_in_output=True,
):
    """Export one row per Key/density; accept Continuous, 2.5, 5 or S2.5/S5.

    Every cell must have all three densities, with no duplicate observations.
    """
    if (
        isinstance(biological_rep, bool)
        or not isinstance(biological_rep, Integral)
        or biological_rep < 1
    ):
        raise ValueError(f"Sheet '{sheet_name}': biological_rep must be a positive integer.")

    data = df.copy()
    data.columns = [str(col).strip() for col in data.columns]
    if data.columns.duplicated().any():
        raise ValueError(f"Sheet '{sheet_name}': duplicate column headers.")
    required = [key_col, stepsize_col, auc_col]
    missing = [col for col in required if col not in data.columns]
    if missing:
        raise ValueError(f"Sheet '{sheet_name}': missing columns {missing}. Found: {list(data.columns)}")
    if data.empty:
        raise ValueError(f"Sheet '{sheet_name}': no rows to export.")
    if data[key_col].isna().any():
        raise ValueError(f"Sheet '{sheet_name}': missing Key.")
    data[key_col] = data[key_col].astype(str).str.strip()
    if data[key_col].eq("").any():
        raise ValueError(f"Sheet '{sheet_name}': Keys must be nonempty.")
    labels = {
        "continuous": "Continuous", "2.5": "S2.5", "s2.5": "S2.5",
        "5": "S5", "5.0": "S5", "s5": "S5",
    }
    data[stepsize_col] = data[stepsize_col].astype(str).str.strip().str.lower().map(labels)
    if data[stepsize_col].isna().any():
        raise ValueError(f"Sheet '{sheet_name}': Step Size must be Continuous, 2.5 or 5.")
    if data.duplicated([key_col, stepsize_col]).any():
        raise ValueError(f"Sheet '{sheet_name}': duplicate Key/Step Size observations.")
    for key, group in data.groupby(key_col):
        if set(group[stepsize_col]) != {"Continuous", "S2.5", "S5"}:
            raise ValueError(f"Sheet '{sheet_name}': '{key}' needs all three step sizes.")
    data[auc_col] = pd.to_numeric(data[auc_col], errors="raise")
    if not np.isfinite(data[auc_col].to_numpy(dtype=float)).all() or (data[auc_col] < 0).any():
        raise ValueError(f"Sheet '{sheet_name}': AUC must be finite and nonnegative.")

    rows = []
    for _, row in data.iterrows():
        key = row[key_col]
        # Do not splitext the cell stem: it may contain a decimal-valued token.
        stem = os.path.basename(key)
        if stem.lower().endswith((".tif", ".tiff", ".png")):
            stem = os.path.splitext(stem)[0]
        parts = stem.split("_")
        if len(parts) < 4 or parts[0] != "Tri" or parts[1] not in {"Ctrl", "5ugLPS"}:
            raise ValueError(f"Sheet '{sheet_name}': invalid Key '{key}'; expected Tri_Ctrl_01_1 or Tri_5ugLPS_01_1.")
        try:
            original_well = int(parts[2])
            image = int(parts[3])
        except ValueError as error:
            raise ValueError(f"Sheet '{sheet_name}': noninteger well/image tokens in '{key}'.") from error
        if original_well < 1 or image < 1:
            raise ValueError(f"Sheet '{sheet_name}': well/image tokens must be positive in '{key}'.")
        well = ((original_well - 1) // 4) % 3 + 1
        rows.append({
            "Key": key,
            "Treatment": f"{parts[0]}_{parts[1]}",
            "Biological_Rep": int(biological_rep),
            "Well": well,
            "Image": image,
            "StepSize": row[stepsize_col],
            "Results": float(row[auc_col]),
        })
    result = pd.DataFrame(rows)
    root = os.path.abspath(output_dir or os.path.dirname(os.path.abspath(input_excel_path)))
    os.makedirs(root, exist_ok=True)
    stem = os.path.splitext(os.path.basename(input_excel_path))[0]
    sheet_tag = f"_{_sanitize_filename_part(sheet_name)}" if include_sheet_name_in_output else ""
    path = os.path.join(root, f"{stem}{sheet_tag}_SamplingDensity_AUC_Long.csv")
    if not overwrite:
        path = _unique_output_path(path)
    result.to_csv(path, index=False)
    print(f"Saved: {path} ({data[key_col].nunique()} cells; {len(result)} cell/density rows)")
    return result, path


def convert_sampling_density_excel_all_sheets_to_mixed_csv(
    input_excel_path,
    key_col="Key",
    stepsize_col="Step Size",
    auc_col="AUC",
    overwrite=False,
    biological_reps=None,
    output_dir=None,
    include_sheet_name_in_output=True,
):
    """biological_reps: positive integer, sheet-name mapping, or None to prompt."""
    results = {}
    with pd.ExcelFile(input_excel_path) as workbook:
        for sheet in workbook.sheet_names:
            if isinstance(biological_reps, Mapping):
                if sheet not in biological_reps:
                    raise ValueError(f"No biological replicate supplied for sheet '{sheet}'.")
                rep = biological_reps[sheet]
            elif biological_reps is not None:
                rep = biological_reps
            else:
                while True:
                    try:
                        rep = int(input(f"Biological replicate for sheet '{sheet}': "))
                    except ValueError:
                        print("Enter a positive integer.")
                        continue
                    if rep > 0:
                        break
                    print("Enter a positive integer.")
            results[sheet] = split_sampling_density_sheet_to_mixed_csv(
                pd.read_excel(workbook, sheet_name=sheet),
                sheet,
                input_excel_path,
                key_col=key_col,
                stepsize_col=stepsize_col,
                auc_col=auc_col,
                overwrite=overwrite,
                biological_rep=rep,
                output_dir=output_dir,
                include_sheet_name_in_output=include_sheet_name_in_output,
            )
    return results


if __name__ == "__main__":
    EXCEL_FILE_PATH = '/Users/allisonpickle/Desktop/Sampling Density/Compiled Step Sizes.xlsx'
    results = convert_sampling_density_excel_all_sheets_to_mixed_csv(
        EXCEL_FILE_PATH,
        biological_reps=1,  # Or {"BioRep1": 1, "BioRep2": 2}; use actual sheet names.
        output_dir='/Users/allisonpickle/Desktop/Sampling Density',
        overwrite=False,
    )
