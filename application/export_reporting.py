"""Export the eight application fields from the archived reporting register."""
from collections import Counter
from pathlib import Path
import csv


ROOT = Path(__file__).resolve().parent
FIELD_MAP = {
    "triad_id": "triad_id",
    "candidate_trait": "candidate",
    "outcome": "outcome",
    "working_structural_class": "display_role",
    "mediator_like_flag": "mediator_like_flag",
    "target_specific_adjustment_label": "adjustment_decision",
    "observed_total_marginal_logOR_response": "observed_response_same_set",
    "minimum_conditional_F": "F_min",
}
EXPECTED_COUNTS = {
    "Necessary": 5,
    "Unnecessary": 129,
    "Overadjustment": 22,
    "Unclassified": 12,
}


def main():
    source = ROOT / "reporting_input/final_application_audit_168.csv"
    destination = ROOT / "check_output/Supplementary_Data_1.csv"
    expected = ROOT / "expected/Supplementary_Data_1.csv"

    with source.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = set(FIELD_MAP.values()) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing reporting columns: {', '.join(sorted(missing))}")
        rows = list(reader)

    if len(rows) != 168 or len({row["triad_id"] for row in rows}) != 168:
        raise ValueError("Expected 168 distinct application triads.")
    if Counter(row["adjustment_decision"] for row in rows) != EXPECTED_COUNTS:
        raise ValueError("Adjustment counts do not match 5/129/22/12.")

    from io import StringIO
    buffer = StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=list(FIELD_MAP))
    writer.writeheader()
    writer.writerows(
        {target: row[source_field] for target, source_field in FIELD_MAP.items()}
        for row in rows
    )
    exported = buffer.getvalue().encode("utf-8-sig")
    if exported != expected.read_bytes():
        raise ValueError("Reporting export differs from the archived expected CSV.")

    destination.parent.mkdir(exist_ok=True)
    destination.write_bytes(exported)
    print("PASS: 168-row export matches the archived CSV; counts 5/129/22/12.")


if __name__ == "__main__":
    main()
