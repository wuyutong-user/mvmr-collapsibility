import csv
import re
from copy import deepcopy
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

from lxml import etree


W_NS = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
M_NS = "http://schemas.openxmlformats.org/officeDocument/2006/math"
XML_NS = "http://www.w3.org/XML/1998/namespace"
NS = {"w": W_NS, "m": M_NS}
W = "{" + W_NS + "}"


def text_content(element):
    return "".join(element.itertext())


def set_plain_paragraph(paragraph, text):
    if paragraph.xpath(".//m:oMath", namespaces=NS):
        raise ValueError("The paragraph contains an OMML formula")
    nodes = paragraph.xpath(".//w:t", namespaces=NS)
    if not nodes:
        raise ValueError("The paragraph has no text nodes")
    nodes[0].text = text
    for node in nodes[1:]:
        node.text = ""


def text_run(text, run_properties=None):
    run = etree.Element(W + "r")
    if run_properties is not None:
        run.append(deepcopy(run_properties))
    node = etree.SubElement(run, W + "t")
    if text[:1].isspace() or text[-1:].isspace():
        node.set("{" + XML_NS + "}space", "preserve")
    node.text = text
    return run


def make_paragraph(text):
    paragraph = etree.Element(W + "p")
    paragraph.append(text_run(text))
    return paragraph


def insert_after(anchor, element):
    parent = anchor.getparent()
    parent.insert(parent.index(anchor) + 1, element)
    return element


def insert_before(anchor, element):
    parent = anchor.getparent()
    parent.insert(parent.index(anchor), element)
    return element


def make_table(headers, rows):
    table = etree.Element(W + "tbl")
    properties = etree.SubElement(table, W + "tblPr")
    style = etree.SubElement(properties, W + "tblStyle")
    style.set(W + "val", "TableGrid")
    width = etree.SubElement(properties, W + "tblW")
    width.set(W + "w", "0")
    width.set(W + "type", "auto")
    grid = etree.SubElement(table, W + "tblGrid")
    for _ in headers:
        etree.SubElement(grid, W + "gridCol")

    for values in [headers] + rows:
        row = etree.SubElement(table, W + "tr")
        for value in values:
            cell = etree.SubElement(row, W + "tc")
            cell_properties = etree.SubElement(cell, W + "tcPr")
            cell_width = etree.SubElement(cell_properties, W + "tcW")
            cell_width.set(W + "w", "0")
            cell_width.set(W + "type", "auto")
            cell.append(make_paragraph(str(value)))
    return table


def set_cell_text(cell, text):
    cell_properties = cell.find(W + "tcPr")
    for child in list(cell):
        if child is not cell_properties:
            cell.remove(child)
    cell.append(make_paragraph(str(text)))


def replace_caption(root, number, text):
    prefix = f"Figure S{number}."
    for paragraph in root.xpath(".//w:p", namespaces=NS):
        if text_content(paragraph).startswith(prefix):
            set_plain_paragraph(paragraph, text)
            return
    raise ValueError(f"Caption not found: {prefix}")


def compact_audit_table(table, display_map):
    rows = table.xpath("./w:tr", namespaces=NS)
    headers = [text_content(cell).strip() for cell in rows[0].xpath("./w:tc", namespaces=NS)]
    required = [
        "triad_id",
        "candidate",
        "outcome",
        "role_status",
        "reporting_role",
        "adjustment_decision",
        "observed_response_same_set",
        "F_min",
    ]
    if any(item not in headers for item in required):
        raise ValueError(f"Unexpected Table S27 header: {headers}")
    index = {name: headers.index(name) for name in headers}

    compact_headers = [
        "triad_id",
        "candidate",
        "outcome",
        "evidence status",
        "display role",
        "mediator-like flag",
        "adjustment decision",
        "observed total marginal logOR response",
        "F_min",
    ]
    header_cells = rows[0].xpath("./w:tc", namespaces=NS)
    for cell in header_cells[9:]:
        rows[0].remove(cell)
    grid = table.find(W + "tblGrid")
    if grid is not None:
        for column in list(grid)[9:]:
            grid.remove(column)
    set_row_values(rows[0], compact_headers)

    for row in rows[1:]:
        cells = row.xpath("./w:tc", namespaces=NS)
        values = [text_content(cell).strip() for cell in cells]
        triad_id = values[index["triad_id"]]
        if triad_id not in display_map:
            raise ValueError(f"No display metadata for {triad_id}")
        compact_values = [
            values[index["triad_id"]],
            values[index["candidate"]],
            values[index["outcome"]],
            values[index["role_status"]],
            display_map[triad_id]["figure4_display_role"],
            display_map[triad_id]["mediator_like_display_marker"],
            values[index["adjustment_decision"]],
            values[index["observed_response_same_set"]],
            values[index["F_min"]],
        ]
        for cell in cells[9:]:
            row.remove(cell)
        set_row_values(row, compact_values)


def set_row_values(row, values):
    cells = row.xpath("./w:tc", namespaces=NS)
    if len(cells) != len(values):
        raise ValueError(f"Expected {len(values)} cells but found {len(cells)}")
    for cell, value in zip(cells, values):
        set_cell_text(cell, value)


def export_full_audit(table, output_path):
    rows = table.xpath("./w:tr", namespaces=NS)
    headers = [text_content(cell).strip() for cell in rows[0].xpath("./w:tc", namespaces=NS)]
    with output_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)
        for row in rows[1:]:
            writer.writerow([text_content(cell).strip() for cell in row.xpath("./w:tc", namespaces=NS)])


def read_display_map(csv_path):
    with csv_path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    required = {"triad_id", "figure4_display_role", "mediator_like_display_marker"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError("The Figure 4 display CSV lacks required fields")
    return {row["triad_id"]: row for row in rows}


def remove_paragraph_range(root, start_prefix, end_prefix):
    paragraphs = root.xpath(".//w:p", namespaces=NS)
    start = next(i for i, p in enumerate(paragraphs) if text_content(p).startswith(start_prefix))
    end = next(i for i, p in enumerate(paragraphs) if text_content(p).startswith(end_prefix))
    for paragraph in paragraphs[start:end]:
        parent = paragraph.getparent()
        if parent is not None:
            parent.remove(paragraph)


def remove_paragraph(paragraph):
    parent = paragraph.getparent()
    if parent is not None:
        parent.remove(paragraph)


def revise(input_path, output_path, full_audit_csv, display_csv):
    display_map = read_display_map(display_csv)
    with ZipFile(input_path, "r") as source_zip:
        root = etree.fromstring(source_zip.read("word/document.xml"))

        tables = root.xpath(".//w:tbl", namespaces=NS)
        full_audit_table = tables[26]
        export_full_audit(full_audit_table, full_audit_csv)
        compact_audit_table(full_audit_table, display_map)

        for paragraph in list(root.xpath(".//w:p", namespaces=NS)):
            text = text_content(paragraph)

            if text.startswith("BMI, candidate trait and outcome summary associations came from the same resource"):
                set_plain_paragraph(
                    paragraph,
                    "BMI, candidate-trait and outcome summary associations were derived from the same UK Biobank resource, so participant overlap across the GWAS could not be excluded and may have been substantial. Exact pairwise overlap could not be reconstructed from the available summary files. The covariance between SNP–BMI and SNP–candidate-trait estimates was unavailable and was therefore set to zero for the conditional-F calculations. No correction for sample overlap was applied; overlap-related weak-instrument bias may therefore remain.",
                )
            elif text.startswith("Cardiovascular outcomes were ascertained using the following codes"):
                set_plain_paragraph(
                    paragraph,
                    "Definitions, coding sources, case counts and source-cohort prevalences for the 14 cardiovascular outcomes are provided in Table S10.",
                )
            elif text.startswith("Coronary artery disease was defined by ICD-9"):
                remove_paragraph(paragraph)
            elif text.startswith("S18. MR analysis and reliability diagnostics"):
                set_plain_paragraph(paragraph, "S18. Directional MR, working-role assignment and reliability diagnostics")
            elif text.startswith("For BMI and each candidate trait, variants associated at"):
                set_plain_paragraph(
                    paragraph,
                    "Directional UVMR evidence. The protocol-defined directional-analysis family comprised 388 exposure–outcome tasks: 24 BMI–candidate-trait directions, 336 candidate-trait–cardiovascular-outcome directions and 28 BMI–cardiovascular-outcome directions. Instruments were selected at P<5×10−8, LD-clumped within 10,000 kb at r²<0.001, and harmonised with effect-allele alignment using action=2. A Wald ratio was used for single-variant analyses, IVW was the primary estimator when at least two variants were available, and weighted-median and MR-Egger estimates were additionally obtained when at least three variants were available. Steiger directionality, heterogeneity, the MR-Egger intercept and leave-one-out analyses were retained as directionality and stability diagnostics.",
                )
            elif text.startswith("We aligned effect alleles across BMI, the candidate trait and the outcome"):
                set_plain_paragraph(
                    paragraph,
                    "Directional evidence tiers. The prespecified multiplicity rule was Benjamini–Hochberg adjustment across the complete family of 388 primary directional tests, rather than separate adjustment within candidate traits or outcomes. The protocol-defined family therefore comprises 388 tasks. The currently available primary-results archive contains 336 records, leaving 52 protocol-defined tasks without an estimate or explicit not-estimable or failure status. Complete 388-task classification is not claimed until those statuses are reconciled.",
                )
            elif text.startswith("Conditional F statistics were calculated with MVMR::strength_mvmr"):
                set_plain_paragraph(
                    paragraph,
                    "Local conditional analysis. For each BMI–candidate-trait–outcome triad, marginal directional evidence was assembled into a three-node graph. The joint instrument set was formed by union, LD re-clumping and multivariable harmonisation. Conditional MVMR was fitted with TwoSampleMR::mv_multiple. A conditional edge was treated as supported at the prespecified operational threshold when P<0.05 and the corresponding conditional F≥10. MRSL was then used to prune the marginal graph with cutoff 0.05, adjustment method 1, use_eggers_step2=0 and vary_mvmr_adj=0. When the conditional graph was not estimable or remained cyclic, the marginal evidence state was retained and the protocol-defined reporting rule was applied. A direction labelled not_supported indicates lack of support at the stated threshold and should not be interpreted as proof that the corresponding causal edge is absent.",
                )
            elif text.startswith("For each BMI-outcome-candidate trait combination"):
                remove_paragraph(paragraph)
            elif not text and paragraph.xpath(".//m:oMath", namespaces=NS):
                math_text = text_content(paragraph)
                if "S32" in math_text or "ΔOR" in math_text:
                    remove_paragraph(paragraph)
            elif text.startswith("The factor 100 expresses the difference"):
                remove_paragraph(paragraph)
            elif text.startswith("The six prespecified local structures and their adjustment interpretations"):
                set_plain_paragraph(
                    paragraph,
                    "Working-role mapping. The six prespecified local structures were mapped from the directional and local-conditional evidence according to the canonical patterns shown below. Observed coefficient response was not used as an input to working-role assignment. A direction labelled not_supported indicates lack of support at the stated threshold and should not be interpreted as proof that the corresponding causal edge is absent.",
                )
            elif text.startswith("S20. Historical response summaries and application diagnostics"):
                set_plain_paragraph(paragraph, "S20. Application diagnostics")
            elif text.startswith("Table S21 summarizes historical descriptive responses"):
                set_plain_paragraph(
                    paragraph,
                    "Table S7 reports descriptive phenotype correlations. Table S8 reports conditional instrument-strength diagnostics, and Table S9 reports pre-harmonisation SNP overlap. These quantities are descriptive or reliability diagnostics; none was used as a standalone causal role classifier.",
                )
            elif text.startswith("In the historical estimates, SBP and DBP produced"):
                set_plain_paragraph(
                    paragraph,
                    "Pre-harmonisation SNP overlap records exact variant overlap between the clumped genome-wide-significant instrument sets. It is distinct from the population summary-association set G_XZ and from the weighted projection δ_XZ. Phenotypic correlations are descriptive and were not used for role classification. Conditional F was evaluated after target specification and was used only as a reliability diagnostic.",
                )
            elif text.startswith("S21. Frozen role mapping and target-specific decision audit"):
                set_plain_paragraph(paragraph, "S21. Application role mapping and target-specific decision audit")
            elif text.startswith("The 12 Independent cause triads and the other non-mediator boundary states"):
                set_plain_paragraph(
                    paragraph,
                    "The 12 Independent cause triads and the other non-mediator boundary states were reported as Unnecessary adjustment under the simplified reporting rule. These operational decision assignments should be distinguished from independent verification of all retained-instrument identification conditions.",
                )
            elif text.startswith("The observed same-set total marginal logOR response was used"):
                set_plain_paragraph(
                    paragraph,
                    "The observed same-set total marginal logOR response was used as descriptive compatibility evidence, and the minimum conditional F statistic was used only as a reliability diagnostic. Neither quantity defined the adjustment class. The complete field-level audit is provided in Supplementary Data 1, whereas the compact table below reports the submission-level audit fields. The identity-error field assesses numerical implementation of the same-set regression identity and is not an independent validation of the working causal role. The C codes remain source-state labels because the current literature-edge source still contains an Example row; they should not be interpreted as a completed confirmatory external-evidence layer until that source is cleaned.",
                )
            elif text.startswith("Each row is one BMI–candidate-trait–outcome triad"):
                set_plain_paragraph(
                    paragraph,
                    "Each row is one BMI–candidate-trait–outcome triad. Evidence status records the audit-layer structural state, and display role records the complete-grid visual category. The mediator-like flag identifies the mediator-like display states. The adjustment decision records the target-specific reporting decision and uses Necessary, Unnecessary, Overadjustment or Unclassified. The observed total marginal logOR response is descriptive compatibility evidence, and F_min is a reliability diagnostic. The compact table does not replace the field-level Supplementary Data 1 audit.",
                )
            elif text.startswith("Table S21.") or text.startswith("Table S23."):
                remove_paragraph(paragraph)

        remove_paragraph_range(root, "S19. Historical outcome-specific MR and MVMR estimates", "S20. Application diagnostics")

        # Renumber the retained tables after removing the historical tables and Table S23.
        for paragraph in root.xpath(".//w:p", namespaces=NS):
            text = text_content(paragraph)
            replacements = {
                "Table S22.": "Table S7.",
                "Table S24.": "Table S8.",
                "Table S25.": "Table S9.",
                "Table S26.": "Table S10.",
                "Table S27.": "Table S11.",
            }
            for old, new in replacements.items():
                if text.startswith(old):
                    set_plain_paragraph(paragraph, new + text[len(old):])
                    break

        inline_replacements = {
            "Table S22": "Table S7",
            "Table S24": "Table S8",
            "Table S25": "Table S9",
            "Table S26": "Table S10",
            "Table S27": "Table S11",
        }
        for paragraph in root.xpath(".//w:p", namespaces=NS):
            text = text_content(paragraph)
            updated = text
            for old, new in inline_replacements.items():
                updated = updated.replace(old, new)
            if updated != text:
                set_plain_paragraph(paragraph, updated)

        # Remove historical Tables S7–S21 and redundant Table S23 from the XML.
        tables = root.xpath(".//w:tbl", namespaces=NS)
        remove_indices = list(range(6, 21)) + [22]
        for index in sorted(remove_indices, reverse=True):
            table = tables[index]
            table.getparent().remove(table)

        # Insert the evidence-tier table after its explanatory paragraph.
        paragraphs = root.xpath(".//w:p", namespaces=NS)
        evidence_anchor = next(p for p in paragraphs if text_content(p).startswith("Directional evidence tiers."))
        evidence_heading = make_paragraph("Directional evidence tiers")
        evidence_table = make_table(
            ["Directional state", "Operational criterion"],
            [
                ["robust_supported", "BH-adjusted P<0.05, at least 3 SNPs, mean F≥10, at least one concordant secondary estimator with P<0.05 and no instability flag"],
                ["sensitivity_qualified", "Primary P<0.05, at least 2 SNPs and mean F≥10, without meeting the robust tier"],
                ["nominal", "Primary P<0.05, fewer than 2 usable SNPs and no instability flag"],
                ["not_supported", "Primary P≥0.05"],
                ["weak_instrument", "Mean F<10"],
                ["unstable", "Conflicting significant estimators, MR-Egger intercept P<0.05 or prespecified evidence conflict"],
                ["not_estimable", "No valid estimate"],
            ],
        )
        insert_after(evidence_anchor, evidence_heading)
        insert_after(evidence_heading, evidence_table)

        local_anchor = next(p for p in root.xpath(".//w:p", namespaces=NS) if text_content(p).startswith("Local conditional analysis."))
        reliability_paragraph = make_paragraph(
            "Reliability and software. Conditional F statistics were calculated with MVMR::strength_mvmr after MVMR::format_mvmr. The cross-trait covariance term was set to zero because it was unavailable. Conditional F was used only to assess reliability, not to define a causal role or adjustment class. The current application response was calculated as the total marginal logOR difference ΔlogOR = β_unadj − β_adj from a fixed retained-SNP set. Analyses were performed in R version 4.3.2 on Ubuntu 22.04.3 LTS. Package versions for TwoSampleMR, MVMR and MRSL were not recoverable from the current application archive because a package-version manifest was unavailable.",
        )
        insert_after(local_anchor, reliability_paragraph)

        role_anchor = next(p for p in root.xpath(".//w:p", namespaces=NS) if text_content(p).startswith("Working-role mapping."))
        role_heading = make_paragraph("Canonical working-role patterns")
        role_table = make_table(
            ["Working role", "Canonical supported pattern"],
            [
                ["Confounder", "Z→X and Z→Y"],
                ["Collider", "X→Z and Y→Z"],
                ["Independent cause", "Z→Y, with no supported X↔Z direction"],
                ["Upstream surrogate of exposure", "Z→X, with no supported Z→Y direction"],
                ["Downstream surrogate of outcome", "Y→Z, with no supported X↔Z or Z→Y direction"],
                ["Downstream surrogate of exposure", "X→Z, with no supported Z→Y or Y→Z direction"],
            ],
        )
        insert_after(role_anchor, role_heading)
        insert_after(role_heading, role_table)

        # Add the V/R/S/C/Q codebook before the compact audit table.
        compact_caption = next(p for p in root.xpath(".//w:p", namespaces=NS) if text_content(p).startswith("Table S11."))
        codebook_heading = make_paragraph("State-code definitions")
        codebook_table = make_table(
            ["Code", "Meaning"],
            [
                ["V0", "Historical/current response lineage not fully reconciled"],
                ["R2", "Binary outcome reported as observed total marginal logOR response"],
                ["S1", "Local conditional graph resolved"],
                ["S3", "Marginal, unresolved or cyclic conditional graph"],
                ["C1", "External evidence compatible with the working classification"],
                ["C3", "External evidence in conflict or tension"],
                ["C0", "Insufficient or boundary external evidence"],
                ["Q1", "Both conditional F values met the prespecified reliability threshold"],
                ["Q2", "Minimum conditional F<10 or estimation instability"],
            ],
        )
        data_note = make_paragraph("The full 26-field audit is provided in Supplementary Data 1. The compact table retains the fields needed to interpret evidence status, display role, target-specific decision and reliability.")
        insert_before(compact_caption, codebook_heading)
        insert_before(compact_caption, codebook_table)
        insert_before(compact_caption, data_note)
        set_plain_paragraph(
            compact_caption,
            "Table S11. Compact application audit for the 168 BMI–candidate-trait–cardiovascular-outcome triads.",
        )

        # Remove the historical crosswalk formula paragraph if it survived the first pass.
        for paragraph in list(root.xpath(".//w:p", namespaces=NS)):
            if text_content(paragraph).startswith("The factor 100 expresses the difference"):
                remove_paragraph(paragraph)

        document_xml = etree.tostring(root, xml_declaration=True, encoding="UTF-8", standalone=True)
        with ZipFile(output_path, "w", ZIP_DEFLATED) as target_zip:
            for item in source_zip.infolist():
                if item.filename == "word/document.xml":
                    target_zip.writestr(item, document_xml)
                else:
                    target_zip.writestr(item, source_zip.read(item.filename))


def replace_table(root, table, headers, rows):
    replacement = make_table(headers, rows)
    parent = table.getparent()
    parent.replace(table, replacement)
    return replacement


def clean_formal_audit_csv(input_path, output_path, display_csv):
    """Remove provisional version/evidence fields from the formal audit export."""
    with input_path.open(newline="", encoding="utf-8-sig") as source:
        reader = csv.DictReader(source)
        if reader.fieldnames is None:
            raise ValueError("The application audit CSV has no header")
        drop = {"V", "C"}
        fieldnames = [name for name in reader.fieldnames if name not in drop]
        if len(fieldnames) != len(reader.fieldnames) - len(drop):
            raise ValueError("Unexpected application audit fields")
        rows = list(reader)

    display_map = read_display_map(display_csv)
    if not rows or any(row.get("triad_id") not in display_map for row in rows):
        raise ValueError("The display-role register does not cover the application audit rows")
    display_fields = ["display_role", "mediator_like_flag"]
    insertion_index = fieldnames.index("role_status")
    fieldnames[insertion_index:insertion_index] = display_fields

    with output_path.open("w", newline="", encoding="utf-8-sig") as target:
        writer = csv.DictWriter(target, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            display_row = display_map[row["triad_id"]]
            row["display_role"] = display_row["figure4_display_role"]
            row["mediator_like_flag"] = display_row["mediator_like_display_marker"]
            writer.writerow({name: row.get(name, "") for name in fieldnames})

    return len(rows), fieldnames


def revise_v7(input_path, output_path, input_audit_csv, output_audit_csv, display_csv):
    """Apply the final S17-S20 structural and field-level corrections to v6."""
    row_count, cleaned_fields = clean_formal_audit_csv(input_audit_csv, output_audit_csv, display_csv)
    if row_count != 168:
        raise ValueError(f"Expected 168 application audit rows, found {row_count}")
    if "V" in cleaned_fields or "C" in cleaned_fields:
        raise ValueError("Provisional V/C fields remain in the formal audit export")

    with ZipFile(input_path, "r") as source_zip:
        root = etree.fromstring(source_zip.read("word/document.xml"))

        # S17: Table S10 supplies counts and prevalences; phenotype definitions remain
        # anchored to the source study and its coding systems.
        for paragraph in root.xpath(".//w:p", namespaces=NS):
            text = text_content(paragraph)
            if text.startswith("Definitions, coding sources, case counts and source-cohort prevalences for the 14 cardiovascular outcomes"):
                set_plain_paragraph(
                    paragraph,
                    "Case counts and source-cohort prevalences for the 14 cardiovascular outcomes are provided in Table S10; phenotype definitions followed the source study and referenced coding systems.",
                )
            elif text.startswith("Directional evidence tiers. The prespecified multiplicity rule"):
                set_plain_paragraph(
                    paragraph,
                    "Directional evidence tiers. The prespecified multiplicity rule was Benjamini–Hochberg adjustment across the complete family of 388 primary directional tests, rather than separate adjustment within candidate traits or outcomes. The protocol-defined family therefore comprises 388 tasks. The currently available primary-results archive contains 336 records, leaving 52 protocol-defined tasks without an estimate or explicit not-estimable or failure status. Complete 388-task classification is not claimed until those statuses are reconciled. The current archive therefore supports an incomplete application audit, not a complete confirmatory rerun or a fully reproducible 388-task result set.",
                )

        # S18: keep the 388-task limitation explicit and separate mutually exclusive
        # evidence tiers from orthogonal reliability and stability flags.
        evidence_table = None
        for table in root.xpath(".//w:tbl", namespaces=NS):
            rows = table.xpath("./w:tr", namespaces=NS)
            if not rows:
                continue
            header = [text_content(cell).strip() for cell in rows[0].xpath("./w:tc", namespaces=NS)]
            if header == ["Directional state", "Operational criterion"]:
                evidence_table = table
                break
        if evidence_table is None:
            raise ValueError("Directional evidence-tier table not found")
        replace_table(
            root,
            evidence_table,
            ["Field", "State", "Criterion"],
            [
                ["Evidence tier", "robust_supported", "BH-adjusted P<0.05, at least 3 SNPs, mean F≥10, at least one concordant secondary estimator with P<0.05 and no instability flag"],
                ["Evidence tier", "sensitivity_qualified", "Primary P<0.05, at least 2 SNPs and mean F≥10, without meeting the robust tier"],
                ["Evidence tier", "nominal", "Exactly one usable SNP and primary P<0.05"],
                ["Evidence tier", "not_supported", "Primary P≥0.05"],
                ["Evidence tier", "not_estimable", "No valid estimate"],
                ["Reliability flag", "weak_instrument", "Mean F<10"],
                ["Stability flag", "unstable", "Estimator conflict, MR-Egger intercept P<0.05 or prespecified evidence conflict"],
            ],
        )

        # Final section numbering after removal of historical S19.
        for paragraph in root.xpath(".//w:p", namespaces=NS):
            text = text_content(paragraph)
            if text.startswith("S20. Application diagnostics"):
                set_plain_paragraph(paragraph, "S19. Application diagnostics")
            elif text.startswith("S21. Application role mapping and target-specific decision audit"):
                set_plain_paragraph(paragraph, "S20. Application role mapping and target-specific decision audit")

        # S18: make complete-grid display precedence explicit and keep trait-specific
        # reasons available in the display register.
        role_table = None
        for table in root.xpath(".//w:tbl", namespaces=NS):
            rows = table.xpath("./w:tr", namespaces=NS)
            if not rows:
                continue
            header = [text_content(cell).strip() for cell in rows[0].xpath("./w:tc", namespaces=NS)]
            if header == ["Working role", "Canonical supported pattern"]:
                role_table = table
                break
        if role_table is None:
            raise ValueError("Canonical working-role table not found")
        insert_after(
            role_table,
            make_paragraph(
                "Complete-grid display precedence. Canonical edge patterns were assigned first. Mediator-like X→Z→Y states were displayed within the downstream-surrogate-of-exposure category and marked M. Boundary states were then mapped according to the protocol-defined precedence rules used for the complete-grid display. Response magnitude was not used in this mapping. Trait-specific display reasons are retained in the display-role register.",
            ),
        )

        # S20: retain only scientifically interpretable audit states in the formal
        # data file. Provisional version and external-evidence codes are not submitted.
        for paragraph in root.xpath(".//w:p", namespaces=NS):
            text = text_content(paragraph)
            if text.startswith("The observed same-set total marginal logOR response was used"):
                set_plain_paragraph(
                    paragraph,
                    "The observed same-set total marginal logOR response was used as descriptive compatibility evidence, and the minimum conditional F statistic was used only as a reliability diagnostic. Neither quantity defined the adjustment class. The field-level audit is provided in Supplementary Data 1, whereas the compact table below reports the submission-level audit fields. The identity-error field assesses numerical implementation of the same-set regression identity and is not an independent validation of the working causal role. Display role, mediator-like flag, reliability and conditional-graph state fields are retained in the formal audit file. External-evidence codes were excluded because the current literature-edge source still contains an Example row.",
                )
            elif text.startswith("State-code definitions"):
                set_plain_paragraph(paragraph, "Formal audit-state definitions")
            elif text.startswith("The full 26-field audit is provided in Supplementary Data 1"):
                set_plain_paragraph(
                    paragraph,
                    "The formal audit file retains the scientific state, role, decision, response and reliability fields. The compact table reports the fields needed for the submission-level audit.",
                )
            elif text.startswith("Figure S46. Application workflow"):
                set_plain_paragraph(
                    paragraph,
                    "Figure S46. Application workflow from structural evidence to target-specific reporting. The workflow separates structure source, role hypothesis, role-implied response, observed response and diagnostics, compatibility assessment and target-specific reporting. Final adjustment classification is reported for the prespecified BMI total-effect target; the working role remains a structural interpretation rather than proof of a true DAG.",
                )

        codebook = None
        for table in root.xpath(".//w:tbl", namespaces=NS):
            rows = table.xpath("./w:tr", namespaces=NS)
            if not rows:
                continue
            header = [text_content(cell).strip() for cell in rows[0].xpath("./w:tc", namespaces=NS)]
            if header == ["Code", "Meaning"]:
                codebook = table
                break
        if codebook is None:
            raise ValueError("Application audit codebook not found")
        replace_table(
            root,
            codebook,
            ["Code", "Meaning"],
            [
                ["R2", "Binary outcome reported as observed total marginal logOR response"],
                ["S1", "Local conditional graph resolved"],
                ["S3", "Marginal, unresolved or cyclic conditional graph"],
                ["Q1", "Both conditional F values met the prespecified reliability threshold"],
                ["Q2", "Minimum conditional F<10 or estimation instability"],
            ],
        )

        # Word applies automatic numbering to the reference list. Remove only the
        # typed duplicate numbers from the paragraph text.
        reference_heading = None
        for paragraph in root.xpath(".//w:p", namespaces=NS):
            if text_content(paragraph).strip() == "Supplementary References":
                reference_heading = paragraph
                break
        if reference_heading is not None:
            paragraphs = root.xpath(".//w:p", namespaces=NS)
            start = paragraphs.index(reference_heading)
            for paragraph in paragraphs[start + 1:]:
                text = text_content(paragraph)
                if re.match(r"^\d+\.\s+", text):
                    set_plain_paragraph(paragraph, re.sub(r"^\d+\.\s+", "", text))

        document_xml = etree.tostring(root, xml_declaration=True, encoding="UTF-8", standalone=True)
        with ZipFile(output_path, "w", ZIP_DEFLATED) as target_zip:
            for item in source_zip.infolist():
                if item.filename == "word/document.xml":
                    target_zip.writestr(item, document_xml)
                else:
                    target_zip.writestr(item, source_zip.read(item.filename))


if __name__ == "__main__":
    import os
    base = Path(os.environ.get("IJE_HISTORICAL_PROJECT_ROOT", ".")).resolve()
    supplementary_dir = base / "02-supplementary"
    revise_v7(
        supplementary_dir / "Supplementary_0919_WYT_S1-S21_revision_v6_20260922.docx",
        supplementary_dir / "Supplementary_0919_WYT_S1-S20_revision_v7_20260922.docx",
        supplementary_dir / "Supplementary_Data_1_application_audit_168.csv",
        supplementary_dir / "Supplementary_Data_1_application_audit_168_v7_20260922.csv",
        supplementary_dir / "00-current" / "Supplementary_Data_Figure4_role_decision_168.csv",
    )
