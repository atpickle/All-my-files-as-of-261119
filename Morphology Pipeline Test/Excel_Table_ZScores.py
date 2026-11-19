"""Column-wise descriptive z-scores and cell coloring for Excel tables."""

from copy import copy
import math
from pathlib import Path
from statistics import mean, stdev
import warnings

from openpyxl import load_workbook
from openpyxl.comments import Comment
from openpyxl.styles import Font, PatternFill
from openpyxl.utils import column_index_from_string, get_column_letter


def color_excel_tables_by_zscore(
    file_path,
    output_path=None,
    *,
    sheet_names=None,
    header_row=1,
    label_col="A",
    numeric_cols=None,
    color_limit=2.0,
    decimal_places=None,
):
    """Save a new workbook with within-column z-score cell colors.

    Each sheet contains one table: a header row, one label column, and numeric
    readout columns. By default all sheets and all named columns except the
    label column are processed. numeric_cols optionally lists Excel letters
    (e.g. ["B", "C"]); sheet_names optionally lists sheet names.

    Rows below the header with nonblank labels are included; entirely blank
    rows are ignored. Blank measurements are excluded column by column.
    Numeric strings, including whitespace, are accepted. Invalid values,
    booleans, infinities and Excel errors raise ValueError. Formula values
    require a numeric cache saved by Excel, otherwise an error is raised.

    z = (value - column mean) / sample SD, computed separately in each sheet,
    using ddof=1. Positive means above average, NOT necessarily better.
    Fewer than two numeric values or zero SD leaves scores undefined/blank,
    clears cell fills, and emits a warning. These are descriptive scores,
    not hypothesis tests, effect sizes against a control, or p-values.

    Original table values/layouts are retained with blue (negative), white
    (zero), and red (positive) fills. Colors saturate at +/-color_limit;
    actual scores are not clipped. Companion sheets contain numeric z-scores,
    and a key sheet records means, SDs, counts and undefined-score reasons.
    Existing conditional formatting on scored sheets is removed so it cannot
    override the new colors; a warning reports this change.

    decimal_places=3 displays readout values, z-scores, means and SDs to three
    decimal places without rounding stored numbers or changing calculations.
    Numeric text in scored cells is converted to numbers for Excel formatting;
    formulas are retained. None preserves original readout formatting and
    displays companion z-scores to three decimals as before.

    Defaults to <input stem>_ZScores.xlsx beside the source. Neither the source
    nor an existing output is overwritten. Returns the saved absolute Path.
    With decimal_places set, the default name ends in _ZScores_<n>dp.xlsx.

    Example:
        result = color_excel_tables_by_zscore(
            "/Users/allisonpickle/Desktop/Morphology Paper Tables.xlsx"
        )
    """
    source = Path(file_path).expanduser().resolve()
    if decimal_places is not None and (
        not isinstance(decimal_places, int)
        or isinstance(decimal_places, bool)
        or not 0 <= decimal_places <= 15
    ):
        raise ValueError("decimal_places must be an integer from 0 to 15, or None.")
    suffix = "" if decimal_places is None else f"_{decimal_places}dp"
    number_format = (
        "0" if decimal_places == 0
        else "0." + "0" * (3 if decimal_places is None else decimal_places)
    )
    destination = (
        Path(output_path).expanduser().resolve()
        if output_path is not None
        else source.with_name(f"{source.stem}_ZScores{suffix}.xlsx")
    )
    if source.suffix.lower() != ".xlsx" or destination.suffix.lower() != ".xlsx":
        raise ValueError("Input and output must be .xlsx workbooks.")
    if destination == source:
        raise ValueError("Output must differ from the source workbook.")
    if destination.exists():
        raise FileExistsError(f"Output already exists; choose another output_path: {destination}")
    if not isinstance(header_row, int) or isinstance(header_row, bool) or header_row < 1:
        raise ValueError("header_row must be a positive Excel row number.")
    if isinstance(color_limit, bool) or not math.isfinite(color_limit) or color_limit <= 0:
        raise ValueError("color_limit must be finite and positive.")
    label_index = column_index_from_string(label_col)
    if numeric_cols is not None:
        if not isinstance(numeric_cols, list) or not numeric_cols:
            raise ValueError("numeric_cols must be a nonempty list of Excel column letters.")
        column_indices = [column_index_from_string(col) for col in numeric_cols]
        if len(set(column_indices)) != len(column_indices) or label_index in column_indices:
            raise ValueError("numeric_cols must be unique and exclude label_col.")
    else:
        column_indices = None

    def unique_title(workbook, base):
        title = base[:31]
        index = 2
        while title.casefold() in {name.casefold() for name in workbook.sheetnames}:
            suffix = f" {index}"
            title = base[:31 - len(suffix)] + suffix
            index += 1
        return title

    def score_fill(score):
        if score is None:
            return PatternFill(fill_type=None)
        fraction = min(abs(score) / color_limit, 1.0)
        endpoint = (49, 130, 189) if score < 0 else (203, 24, 29)
        rgb = tuple(round(255 + fraction * (channel - 255)) for channel in endpoint)
        return PatternFill("solid", fgColor="FF" + "".join(f"{channel:02X}" for channel in rgb))

    workbook = load_workbook(source)
    cached = load_workbook(source, data_only=True)
    try:
        selected = workbook.sheetnames[:] if sheet_names is None else sheet_names
        if not isinstance(selected, list) or not selected or len(set(selected)) != len(selected):
            raise ValueError("sheet_names must be a nonempty list of unique sheet names.")
        unknown = set(selected) - set(workbook.sheetnames)
        if unknown:
            raise ValueError(f"Unknown sheet names: {sorted(unknown)}")
        key = workbook.create_sheet(unique_title(workbook, "Z-score Key"))
        key.append(["Within-column descriptive z-scores"])
        key.append(["Formula", "(value - column mean) / sample SD (ddof=1)"])
        key.append(["Interpretation", "Positive = above average, not necessarily better; not a statistical test."])
        key.append(["Colors", f"Blue <= -{color_limit:g}; white = 0; red >= +{color_limit:g}; scores not clipped."])
        key.append(["Missing / undefined", "Blank scores; no fill. Missing values excluded from mean and SD."])
        key.append(["Source", str(source)])
        key.append([])
        key.append(["Sheet", "Column", "Readout", "Numeric n", "Mean", "Sample SD", "Status", "Z-score sheet"])
        for sheet_name in selected:
            sheet = workbook[sheet_name]
            if header_row > sheet.max_row:
                raise ValueError(f"{sheet_name}: header_row is outside the sheet.")
            if sheet.cell(header_row, label_index).value is None:
                raise ValueError(f"{sheet_name}: label-column header is blank.")
            cols = column_indices if column_indices is not None else [
                col for col in range(1, sheet.max_column + 1)
                if col != label_index and sheet.cell(header_row, col).value is not None
            ]
            if not cols:
                raise ValueError(f"{sheet_name}: no numeric readout columns selected.")
            rows = []
            for row in range(header_row + 1, sheet.max_row + 1):
                label = sheet.cell(row, label_index).value
                if label is None or not str(label).strip():
                    if any(sheet.cell(row, col).value is not None for col in cols):
                        raise ValueError(f"{sheet_name}: data at row {row} has no row label.")
                else:
                    rows.append(row)
            if not rows:
                raise ValueError(f"{sheet_name}: no labeled data rows.")
            if len(sheet.conditional_formatting):
                warnings.warn(f"{sheet_name}: removing existing conditional formatting to show z-score colors.", stacklevel=2)
                sheet.conditional_formatting = type(sheet.conditional_formatting)()
            scores = workbook.create_sheet(unique_title(workbook, f"{sheet_name} ZScores"))
            scores.cell(1, 1, sheet.cell(header_row, label_index).value)
            for output_row, row in enumerate(rows, start=2):
                scores.cell(output_row, 1, sheet.cell(row, label_index).value)
            for output_col, col in enumerate(cols, start=2):
                header = sheet.cell(header_row, col).value
                if header is None:
                    raise ValueError(f"{sheet_name}: column {get_column_letter(col)} has no header.")
                scores.cell(1, output_col, str(header).strip())
                values = {}
                for row in rows:
                    cell = sheet.cell(row, col)
                    value = cached[sheet_name].cell(row, col).value if cell.data_type == "f" else cell.value
                    if cell.data_type == "f" and value is None:
                        raise ValueError(f"{sheet_name}!{cell.coordinate}: formula has no cached result; recalculate and save in Excel.")
                    if value is None or isinstance(value, str) and not value.strip():
                        continue
                    try:
                        number = float(value)
                    except (TypeError, ValueError) as exc:
                        raise ValueError(f"{sheet_name}!{cell.coordinate}: expected numeric value, got {value!r}.") from exc
                    if isinstance(value, bool) or cell.data_type == "e" or not math.isfinite(number):
                        raise ValueError(f"{sheet_name}!{cell.coordinate}: invalid/nonfinite numeric value.")
                    values[row] = number
                average = mean(values.values()) if values else None
                sd = stdev(values.values()) if len(values) >= 2 else None
                status = "Scored"
                if sd is None or sd == 0:
                    status = "Undefined: fewer than two numeric values" if sd is None else "Undefined: zero standard deviation"
                    warnings.warn(f"{sheet_name}, {header}: {status}.", stacklevel=2)
                key.append([sheet_name, get_column_letter(col), str(header).strip(), len(values), average, sd, status, scores.title])
                if decimal_places is not None:
                    for key_col in (5, 6):
                        key.cell(key.max_row, key_col).number_format = number_format
                for output_row, row in enumerate(rows, start=2):
                    score = (values[row] - average) / sd if row in values and sd is not None and sd > 0 else None
                    original = sheet.cell(row, col)
                    if decimal_places is not None:
                        if row in values and original.data_type != "f":
                            original.value = values[row]
                        original.number_format = number_format
                    fill = score_fill(score)
                    original.fill = copy(fill)
                    text = f"Within-column z-score: {score:.4f} SD from mean." if score is not None else f"Z-score undefined/missing. {status}."
                    previous = original.comment.text + "\n\n" if original.comment else ""
                    original.comment = Comment(previous + text, "Z-score analysis")
                    target = scores.cell(output_row, output_col, score)
                    target.number_format = number_format
                    target.fill = copy(fill)
            scores.freeze_panes = "B2"
            scores.column_dimensions["A"].width = 28
            for col in range(2, len(cols) + 2):
                scores.column_dimensions[get_column_letter(col)].width = 22
            for cell in scores[1]:
                cell.font = Font(bold=True)
        key.column_dimensions["A"].width = 28
        for col in range(2, 9):
            key.column_dimensions[get_column_letter(col)].width = 24
        key.column_dimensions["B"].width = 72
        key.freeze_panes = "A9"
        workbook.save(destination)
    finally:
        workbook.close()
        cached.close()
    print(f"Saved colored tables and z-scores to: {destination}")
    return destination

if __name__ == "__main__":
    color_excel_tables_by_zscore(
        "/Users/allisonpickle/Desktop/Morphology Paper Tables.xlsx",
        output_path="/Users/allisonpickle/Desktop/Morphology Paper Tables_ZScores_3dp_v2.xlsx",
        decimal_places=3,
    )