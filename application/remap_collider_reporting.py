"""Apply the archived Collider reporting revision without re-estimating effects."""
from collections import Counter
from pathlib import Path
import csv


ROOT = Path(__file__).resolve().parent
COLLIDER_IDS = {
    "T042_BMI_SBP_HTN",
    "T071_BMI_HDL_C_CAD",
    "T083_BMI_HDL_C_AF",
    "T114_BMI_Glucose_HF",
    "T117_BMI_Glucose_AAA",
    "T145_BMI_Smoking_AAA",
    "T153_BMI_Smoking_AF",
    "T156_BMI_Alcohol_HF",
    "T167_BMI_Alcohol_AF",
}
EXPECTED_COUNTS = {
    "Necessary": 5,
    "Unnecessary": 129,
    "Overadjustment": 22,
    "Unclassified": 12,
}
DECISION_REASON = (
    "Collider role alone does not determine the consequence of candidate-trait "
    "inclusion; additional retained-instrument evidence was not established."
)


def main():
    source = ROOT / "provenance/final_application_audit_before_collider_20261006.csv"
    destination = ROOT / "reporting_input/final_application_audit_168.csv"
    with source.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or [])
        rows = list(reader)

    if len(rows) != 168 or len({row["triad_id"] for row in rows}) != 168:
        raise ValueError("Expected 168 distinct application triads.")
    observed_ids = {
        row["triad_id"] for row in rows if row["display_role"] == "Collider"
    }
    if observed_ids != COLLIDER_IDS:
        raise ValueError("Collider IDs differ from the archived nine-triad revision.")

    for row in rows:
        if row["triad_id"] not in COLLIDER_IDS:
            continue
        if row["adjustment_decision"] != "Overadjustment":
            raise ValueError(f"Unexpected initial decision for {row['triad_id']}.")
        row["adjustment_decision"] = "Unclassified"
        if "display_adjustment_decision" in row:
            row["display_adjustment_decision"] = "Unclassified"
        for field in ("decision_basis", "display_decision_reason"):
            if field in row:
                row[field] = DECISION_REASON
        if "decision_status" in row:
            row["decision_status"] = "unclassified_collider_instrument_structure"

    if Counter(row["adjustment_decision"] for row in rows) != EXPECTED_COUNTS:
        raise ValueError("Adjustment counts do not match 5/129/22/12.")
    destination.parent.mkdir(exist_ok=True)
    with destination.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print("PASS: nine Collider decisions remapped; estimates and diagnostics retained.")


if __name__ == "__main__":
    main()
