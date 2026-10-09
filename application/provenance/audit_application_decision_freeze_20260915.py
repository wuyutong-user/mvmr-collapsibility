"""Audit and freeze the current 168-triad application decision layer.

The raw ERACA_RESULT files are read-only inputs. This script writes a separate
dated audit bundle and never edits the raw result tables. Run from a historical
project root or set IJE_HISTORICAL_PROJECT_ROOT to the directory containing
04-application and 01-manuscript. These historical inputs are not bundled.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
from collections import Counter
from datetime import datetime
from pathlib import Path


ROOT = Path(os.environ.get("IJE_HISTORICAL_PROJECT_ROOT", Path.cwd())).expanduser().resolve()

RESULT = ROOT / "04-application" / "ERACA_RESULT"
AUDIT = ROOT / "04-application" / "application_audit" / "audit_20260915_simplified_rule"

RAW_ROLE = RESULT / "role_assignments_168.csv"
RAW_RESPONSE = RESULT / "role_response_concordance_168.csv"
ELIGIBLE = RESULT / "main_text_candidate_triads.csv"
CONDITIONAL_GRAPH = RESULT / "graphs" / "conditional_graph_consensus.csv"
CONDITIONAL_EDGES = RESULT / "graphs" / "conditional_edges_mvmr.csv"
DECOMPOSITION = RESULT / "graphs" / "same_set_response_decomposition.csv"
IC_AUDIT = ROOT / "01-manuscript" / "application_retained_instrument_audit_0916.csv"

REQUIRED = [RAW_ROLE, RAW_RESPONSE, ELIGIBLE, CONDITIONAL_GRAPH, CONDITIONAL_EDGES, DECOMPOSITION, IC_AUDIT]
missing = [str(p) for p in REQUIRED if not p.is_file()]
if missing:
    raise SystemExit("Missing required audit input(s):\n" + "\n".join(missing))


AUDIT.mkdir(parents=True, exist_ok=True)

def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def as_float(value: str, label: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise SystemExit(f"Non-numeric {label}: {value!r}") from exc
    if number != number or number in (float("inf"), float("-inf")):
        raise SystemExit(f"Non-finite {label}: {value!r}")
    return number


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


raw = read_csv(RAW_ROLE)
responses = read_csv(RAW_RESPONSE)
eligible = read_csv(ELIGIBLE)
conditional_graph = read_csv(CONDITIONAL_GRAPH)
conditional_edges = read_csv(CONDITIONAL_EDGES)
decomposition = read_csv(DECOMPOSITION)
ic_audit = read_csv(IC_AUDIT)

if len(raw) != 168 or len(responses) != 168 or len(eligible) != 80:
    raise SystemExit(f"Unexpected source sizes: raw={len(raw)}, response={len(responses)}, eligible={len(eligible)}")

raw_by_id = {row["triad_id"]: row for row in raw}
response_by_id = {row["triad_id"]: row for row in responses}
eligible_by_id = {row["triad_id"]: row for row in eligible}
if len(raw_by_id) != 168 or len(response_by_id) != 168 or len(eligible_by_id) != 80:
    raise SystemExit("Duplicate triad IDs detected in an input table")
if set(raw_by_id) != set(response_by_id):
    raise SystemExit("Raw role and response tables do not cover the same 168 triads")
if not set(eligible_by_id).issubset(raw_by_id):
    raise SystemExit("The eligible main-text set contains a triad absent from raw roles")

candidate_order = {"WC", "HC", "SBP", "DBP", "HDL-C", "LDL-C", "TC", "TG", "Glucose", "HbA1c", "Smoking", "Alcohol"}
outcome_order = {"AAA", "AF", "AVS", "CAD", "DVT", "HF", "HTN", "ICH", "IS", "PE", "PVD", "SAH", "TAA", "TIA"}
if {row["candidate"] for row in raw} != candidate_order or {row["outcome"] for row in raw} != outcome_order:
    raise SystemExit("The raw table is not the prescribed 12 x 14 candidate-outcome grid")
if len({(row["candidate"], row["outcome"]) for row in raw}) != 168:
    raise SystemExit("Duplicate candidate-outcome pair detected")

canonical_roles = {
    "confounder",
    "collider",
    "independent_cause",
    "upstream_surrogate_of_exposure",
    "downstream_surrogate_of_outcome",
    "downstream_surrogate_of_exposure",
}
expected_working = {
    "confounder": 5,
    "collider": 9,
    "independent_cause": 12,
    "upstream_surrogate_of_exposure": 14,
    "downstream_surrogate_of_outcome": 0,
    "downstream_surrogate_of_exposure": 40,
}
expected_boundary = {
    "overlapping_measure": 28,
    "mediator_like_pathway": 25,
    "insufficiently_informative": 27,
    "ambiguous": 5,
    "feedback_or_cycle": 3,
}

boundary_reason_map = {
    "overlapping_measure": "overlapping_measure",
    "mediator_out_of_library": "mediator_like_pathway",
    "insufficient_evidence": "insufficiently_informative",
    "ambiguous": "ambiguous",
    "feedback_or_cycle": "feedback_or_cycle",
}
eligible_ids = set(eligible_by_id)

graph_rows = {}
for row in conditional_graph:
    graph_rows.setdefault(row["triad_id"], []).append(row)
edge_rows = {}
for row in conditional_edges:
    edge_rows.setdefault(row["triad_id"], []).append(row)
decomp_by_id = {row["triad_id"]: row for row in decomposition}
for label, mapping in {
    "conditional graph": graph_rows,
    "conditional edge": edge_rows,
    "same-set decomposition": decomp_by_id,
}.items():
    if set(mapping) != set(raw_by_id):
        raise SystemExit(f"{label} input does not cover exactly the same 168 triads")

audit_rows: list[dict[str, str]] = []
for triad_id, source in raw_by_id.items():
    response = response_by_id[triad_id]
    if source["candidate"] != response["candidate"] or source["outcome"] != response["outcome"]:
        raise SystemExit(f"Candidate/outcome mismatch for {triad_id}")
    if triad_id not in decomp_by_id:
        raise SystemExit(f"Missing same-set decomposition for {triad_id}")
    if response["identity_verified"].upper() != "TRUE" or as_float(response["identity_error"], f"identity_error {triad_id}") > 1e-8:
        raise SystemExit(f"Same-set coefficient identity failed for {triad_id}")

    if triad_id in eligible_ids:
        reporting_role = source["role"]
        role_status = "resolved_six_role"
        boundary_reason = ""
    else:
        reporting_role = "boundary_unresolved"
        role_status = "boundary_unresolved"
        boundary_reason = boundary_reason_map.get(source["role"], "insufficiently_informative")

    decision = "Unclassified"
    decision_status = "unclassified_structure"
    decision_basis = "The evidence did not establish a unique target-relevant adjustment decision."
    if reporting_role == "confounder":
        decision = "Necessary"
        decision_status = "resolved_target_decision"
        decision_basis = "Omission may prevent identification of the prespecified BMI total-effect target under the stated MVMR assumptions."
    elif reporting_role == "collider":
        decision = "Overadjustment"
        decision_status = "resolved_target_decision"
        decision_basis = "Conditioning on a collider opens a noncausal pathway for the BMI total-effect target."
    elif reporting_role in {
        "upstream_surrogate_of_exposure",
        "downstream_surrogate_of_outcome",
        "downstream_surrogate_of_exposure",
    }:
        decision = "Unnecessary"
        decision_status = "resolved_target_decision"
        decision_basis = "The canonical surrogate relation does not alter the prespecified BMI total-effect target under the frozen Table 2 rule."
    elif reporting_role == "independent_cause":
        decision = "Unnecessary"
        decision_status = "resolved_by_simplified_freeze_rule"
        decision_basis = "Independent cause is assigned Unnecessary under the final simplified application rule. The retained-SNP manifest is not present, so this is an operational classification and not a claim that retained-instrument validity was independently verified."
    elif boundary_reason == "mediator_like_pathway":
        supported = (
            source["graph_resolution"] == "conditional"
            and source["X_to_Z"] in {"conditional_supported", "robust_supported"}
            and source["Z_to_Y"] in {"conditional_supported", "robust_supported"}
            and source["mrsl_status"] == "estimated"
        )
        if supported:
            decision = "Overadjustment"
            decision_status = "resolved_target_consequence_boundary_role"
            decision_basis = "Conditional evidence supports BMI -> candidate -> outcome; inclusion changes the BMI total-effect target toward a direct-effect contrast. The structural label remains outside the six-role prevalence."
        else:
            decision_basis = "Mediator-like pattern is recorded, but conditional graph evidence is unresolved or unstable; Unclassified is retained because this is the only boundary state allowed to remain unresolved under the final simplified rule."
    else:
        decision = "Unnecessary"
        decision_status = "resolved_by_simplified_freeze_rule"
        decision_basis = "Non-mediator boundary state is assigned Unnecessary under the final simplified application rule so that Unclassified is reserved exclusively for unresolved mediator-like pathways. This is an operational classification, not proof that every upstream structural relation is identified."

    f_bmi = as_float(response["conditional_F_BMI"], f"conditional_F_BMI {triad_id}")
    f_candidate = as_float(response["conditional_F_candidate"], f"conditional_F_candidate {triad_id}")
    response_value = as_float(response["observed_response_same_set"], f"response {triad_id}")
    audit_rows.append({
        "triad_id": triad_id,
        "candidate": source["candidate"],
        "outcome": source["outcome"],
        "raw_role": source["role"],
        "reporting_role": reporting_role,
        "role_status": role_status,
        "boundary_reason": boundary_reason,
        "role_confidence": source["role_confidence"],
        "graph_resolution": source["graph_resolution"],
        "mrsl_status": source["mrsl_status"],
        "X_to_Z": source["X_to_Z"],
        "Z_to_X": source["Z_to_X"],
        "Z_to_Y": source["Z_to_Y"],
        "Y_to_Z": source["Y_to_Z"],
        "adjustment_decision": decision,
        "decision_status": decision_status,
        "decision_basis": decision_basis,
        "nsnp_same_set": response["nsnp_same_set"],
        "observed_response_same_set": response["observed_response_same_set"],
        "conditional_F_BMI": response["conditional_F_BMI"],
        "conditional_F_candidate": response["conditional_F_candidate"],
        "F_min": f"{min(f_bmi, f_candidate):.10g}",
        "identity_error": response["identity_error"],
    })

working_counts = Counter(row["reporting_role"] for row in audit_rows if row["reporting_role"] != "boundary_unresolved")
boundary_counts = Counter(row["boundary_reason"] for row in audit_rows if row["boundary_reason"])
expected_working_observed = {key: value for key, value in expected_working.items() if value}
if dict(working_counts) != expected_working_observed:
    raise SystemExit(f"Working-role counts do not match the 80-triad eligible set: {working_counts}")
if dict(boundary_counts) != expected_boundary:
    raise SystemExit(f"Boundary counts do not match the frozen 34 categories: {boundary_counts}")

decision_counts = Counter(row["adjustment_decision"] for row in audit_rows)
expected_decisions = {"Necessary": 5, "Unnecessary": 129, "Overadjustment": 31, "Unclassified": 3}
if dict(decision_counts) != expected_decisions:
    raise SystemExit(f"Decision counts failed the audit: {decision_counts}")

ic_rows = [row for row in audit_rows if row["reporting_role"] == "independent_cause"]
ic_source = {row["triad_id"]: row for row in ic_audit}
if len(ic_rows) != 12 or set(ic_source) != {row["triad_id"] for row in ic_rows}:
    raise SystemExit("The 12 independent-cause audit rows do not match the eligible set")
ic_output = []
for row in ic_rows:
    source = ic_source[row["triad_id"]]
    if source["final_decision"] != "pending_post_union_reclump_manifest":
        raise SystemExit(f"IC source is not pending as expected for {row['triad_id']}")
    ic_output.append({
        "triad_id": row["triad_id"],
        "candidate": row["candidate"],
        "outcome": row["outcome"],
        "working_role": "Independent cause",
        "retained_snp_count_reported": source["retained_snp_count"],
        "same_set_snp_count": source["retained_snp_n_same_set"],
        "candidate_associated_retained_variants": source["candidate_associated_retained_variants"],
        "evidence_variants_reach_Y_through_Z": source["evidence_variants_reach_Y_through_Z"],
        "exclusion_restriction_after_omitting_Z": source["exclusion_restriction_after_omitting_Z"],
        "audit_status": "ASSIGNED_BY_SIMPLIFIED_FREEZE_RULE",
        "manifest_status": "MISSING_FINAL_RETAINED_SNP_MANIFEST",
        "adjustment_decision": "Unnecessary",
        "missing_source": source["missing_source"],
        "audit_reason": "Independent cause assigned Unnecessary under the final simplified application rule. The manifest-based retained-instrument check remains unavailable; this row is not presented as manifest-verified.",
        "response_and_F_use": "descriptive only; not used to assign this simplified rule-based decision",
    })

fieldnames = list(audit_rows[0])
with (AUDIT / "application_168_decision_audit_freeze_20260915.csv").open("w", encoding="utf-8-sig", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(sorted(audit_rows, key=lambda row: row["triad_id"]))

ic_fields = list(ic_output[0])
with (AUDIT / "independent_cause_retained_instrument_audit_12_20260915.csv").open("w", encoding="utf-8-sig", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=ic_fields)
    writer.writeheader()
    writer.writerows(sorted(ic_output, key=lambda row: row["triad_id"]))

source_files = [RAW_ROLE, RAW_RESPONSE, ELIGIBLE, CONDITIONAL_GRAPH, CONDITIONAL_EDGES, DECOMPOSITION, IC_AUDIT]
with (AUDIT / "source_manifest_sha256_20260915.csv").open("w", encoding="utf-8-sig", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["path", "bytes", "modified", "sha256"])
    writer.writeheader()
    for path in source_files:
        stat = path.stat()
        writer.writerow({
            "path": str(path.relative_to(ROOT)),
            "bytes": stat.st_size,
            "modified": datetime.fromtimestamp(stat.st_mtime).isoformat(timespec="seconds"),
            "sha256": sha256(path),
        })

checks = {
    "status": "PASS_SIMPLIFIED_RULE_WITH_MANIFEST_CAVEAT",
    "generated": datetime.now().isoformat(timespec="seconds"),
    "raw_triads": len(raw),
    "eligible_six_role_triads": len(eligible),
    "boundary_triads": len(raw) - len(eligible),
    "working_role_counts": dict(working_counts),
    "boundary_reason_counts": dict(boundary_counts),
    "adjustment_decision_counts": dict(decision_counts),
    "same_set_identity_max_abs_error": max(as_float(row["identity_error"], "identity_error") for row in responses),
    "conditional_F_min_ge_10": sum(min(as_float(row["conditional_F_BMI"], "F"), as_float(row["conditional_F_candidate"], "F")) >= 10 for row in responses),
    "conditional_F_min_lt_10": sum(min(as_float(row["conditional_F_BMI"], "F"), as_float(row["conditional_F_candidate"], "F")) < 10 for row in responses),
    "mediator_like_supported_decision_level_overadjustment": sum(row["boundary_reason"] == "mediator_like_pathway" and row["adjustment_decision"] == "Overadjustment" for row in audit_rows),
    "mediator_like_unclassified": sum(row["boundary_reason"] == "mediator_like_pathway" and row["adjustment_decision"] == "Unclassified" for row in audit_rows),
    "non_mediator_boundary_defaulted_to_unnecessary": sum(1 for row in audit_rows if row["boundary_reason"] and row["boundary_reason"] != "mediator_like_pathway" and row["adjustment_decision"] == "Unnecessary"),
    "independent_cause_audit_status": "ASSIGNED_UNNECESSARY_BY_SIMPLIFIED_RULE_FINAL_MANIFEST_NOT_VERIFIED",
    "raw_inputs_modified": False,
}
with (AUDIT / "audit_checks_20260915.json").open("w", encoding="utf-8") as handle:
    json.dump(checks, handle, ensure_ascii=False, indent=2)

report = f"""# UK Biobank application decision audit and evidence freeze

**Generated:** {checks['generated']}  
**Status:** `PASS_SIMPLIFIED_RULE_WITH_MANIFEST_CAVEAT`  
**Scope:** current `04-application/ERACA_RESULT` release, 168 BMI–candidate-trait–outcome triads.

## Freeze statement

The raw ERACA_RESULT tables were not modified. This bundle freezes the evidence-backed reporting layer as a separate, dated view. Structural role assignment and target-specific adjustment consequence are recorded as two different fields. A mediator-like path can therefore have an adjustment consequence without being counted as one of the six canonical role prevalences.

## Audited counts

| Layer | Category | n |
|---|---|---:|
| Canonical working role | Confounder | 5 |
| Canonical working role | Collider | 9 |
| Canonical working role | Independent cause | 12 |
| Canonical working role | Upstream surrogate of exposure | 14 |
| Canonical working role | Downstream surrogate of outcome | 0 |
| Canonical working role | Downstream surrogate of exposure | 40 |
| Canonical working role | **Total six-role working set** | **80** |
| Boundary reason | Overlapping measure | 28 |
| Boundary reason | Mediator-like pathway | 25 |
| Boundary reason | Insufficiently informative | 27 |
| Boundary reason | Ambiguous | 5 |
| Boundary reason | Feedback or cycle | 3 |
| Boundary reason | **Total boundary/unresolved** | **88** |

The decision layer is:

| Target-specific adjustment decision | n | Evidence status |
|---|---:|---|
| Necessary | 5 | Frozen for the five confounder triads under stated MVMR assumptions |
| Unnecessary | 129 | 54 canonical surrogate rows + 12 independent-cause rows + 63 non-mediator boundary rows assigned by the simplified rule |
| Overadjustment | 9 | Frozen for the nine collider triads |
| Overadjustment, mediator-like boundary role | 22 | Conditional BMI → candidate → outcome path supported; role remains outside six-role prevalence |
| Unclassified | 3 | Only the three unresolved mediator-like cases |
| **Total** | **168** | |

The simplified decision count is therefore **Necessary 5, Unnecessary 129, Overadjustment 31, Unclassified 3**. This implements the final rule that all 12 independent-cause rows are Unnecessary and that Unclassified is reserved exclusively for unresolved mediator-like pathways.

## Independent-cause audit

All 12 independent-cause triads are assigned `adjustment_decision = Unnecessary` by the final simplified application rule. The local ERACA_RESULT release still does not contain the exact per-triad SNP rows used after joint-instrument union, LD re-clumping and outcome harmonisation. The dedicated IC table therefore records `ASSIGNED_BY_SIMPLIFIED_FREEZE_RULE` and `MISSING_FINAL_RETAINED_SNP_MANIFEST`; it does not claim that retained-instrument validity was independently verified.

## Mediator-like boundary adjudication

Twenty-five raw rows are labelled `mediator_out_of_library`. Twenty-two have `graph_resolution = conditional`, `mrsl_status = estimated`, and both `X_to_Z` and `Z_to_Y` marked `conditional_supported` or `robust_supported`. These rows are assigned a decision-level consequence of Overadjustment for the BMI total-effect target because including the candidate changes the target toward a direct-effect contrast. They remain outside the six-role role count. Three rows have marginally unresolved conditional graphs and remain Unclassified.

This is deliberately not a rebranding of mediator-like rows as a seventh canonical role. The three unresolved mediator-like rows are the only rows retained as Unclassified. All other non-mediator boundary states are assigned Unnecessary by the simplified operational rule.

## Numerical integrity checks

- 168 unique candidate–outcome triads cover the prescribed 12 × 14 grid.
- 80 triads pass the current main-text working-role screen; 88 remain boundary/unresolved.
- Same-set coefficient identity is verified for all 168; maximum absolute identity error is `{checks['same_set_identity_max_abs_error']:.3g}`.
- Minimum conditional F is at least 10 in `{checks['conditional_F_min_ge_10']}` triads and below 10 in `{checks['conditional_F_min_lt_10']}`.
- Observed response magnitude and conditional F are retained as descriptive/reliability evidence only. Neither is used to assign Necessary, Unnecessary or Overadjustment.

## Superseded or unsafe fields

The raw role table's `adjustment_interpretation` column is not used as the freeze source because its surrogate rows retain an older interpretation and do not encode the Phase 2 canonical-surrogate boundary audit. The dated audit CSV is the reporting source for future Figure 4, Results text and supplementary linkage.

## Required external completion

The final manifest remains a validation artifact if a future evidence-based re-audit of the 12 independent-cause rows is requested; it is not a blocker for the current simplified classification. Export `retained_instrument_manifest_12.csv` using the existing export-only script in the controlled-data environment. The required fields and command are documented in `04-application/reproducibility_package/notes/retained_instrument_manifest_12_schema_0917.md` and `retained_instrument_manifest_12_server_run_0917.md`.

## Artifacts

- `application_168_decision_audit_freeze_20260915.csv`
- `independent_cause_retained_instrument_audit_12_20260915.csv`
- `source_manifest_sha256_20260915.csv`
- `audit_checks_20260915.json`
- `eraca_result_manifest_presence_audit_20260915.md`
"""
(AUDIT / "application_decision_audit_freeze_20260915.md").write_text(report, encoding="utf-8")

print(f"Audit bundle: {AUDIT}")
print(f"Decision counts: {dict(decision_counts)}")
print(f"Max same-set identity error: {checks['same_set_identity_max_abs_error']:.3g}")
print("Raw ERACA_RESULT inputs were not modified")
