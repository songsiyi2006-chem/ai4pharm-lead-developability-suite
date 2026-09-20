#!/usr/bin/env python3
"""Audit archived versus fresh omnibus CSV inputs without modifying either.

Rows are aligned by position, not silently sorted or joined. Numeric differences
use |fresh - archived| <= atol + rtol * |archived|. Missing cells are never zero.
Expected scientific/numerical differences are reported and return exit code 0;
an invalid provenance file or unreadable/malformed input is an audit failure (2).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
MISSING = {"", "na", "n/a", "nan", "null", "none"}
STATUSES = ("identical", "tolerant", "differ", "missing", "audit_fail")


def digest(path):
    if not path.is_file():
        return None
    checksum = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024*1024), b""):
            checksum.update(chunk)
    return checksum.hexdigest()


def parse_cell(raw):
    text = raw.strip()
    if text.lower() in MISSING:
        return "missing", None
    try:
        value = float(text)
    except ValueError:
        return "text", raw
    return "numeric", value


def finite_json(value):
    """Keep valid JSON while explicitly representing infinite numeric results."""
    if math.isfinite(value):
        return value
    return "Infinity" if value > 0 else "-Infinity"


def read_table(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        table = list(csv.reader(handle))
    if not table or not table[0]:
        raise ValueError(f"Empty CSV: {path}")
    # ViennaRNA BPP inputs are headerless square matrices. Numeric header names
    # are absent in the other source-backed tables selected by the atlas.
    headerless = all(parse_cell(cell)[0] == "numeric" for cell in table[0])
    columns = [f"matrix_col_{i+1:03d}" for i in range(len(table[0]))] if headerless else table.pop(0)
    if len(set(columns)) != len(columns):
        raise ValueError(f"Duplicate CSV headers: {path}")
    bad = [i+1+int(not headerless) for i, row in enumerate(table) if len(row) != len(columns)]
    if bad:
        raise ValueError(f"Ragged CSV rows at lines {bad[:5]}: {path}")
    return columns, table, headerless


def profile(values):
    counts = {"numeric": 0, "text": 0, "missing": 0, "nonfinite_numeric": 0}
    for value in values:
        kind, parsed = parse_cell(value)
        counts[kind] += 1
        if kind == "numeric" and not math.isfinite(parsed):
            counts["nonfinite_numeric"] += 1
    present = [key for key in ("numeric", "text") if counts[key]]
    counts["type"] = present[0] if len(present) == 1 else "mixed" if present else "missing_only"
    return counts


def compare_column(name, old, new, rtol, atol, examples=3):
    result = {"column": name, "archived_profile": profile(old), "fresh_profile": profile(new),
              "compared_rows": min(len(old), len(new)), "numeric_pairs": 0,
              "numeric_exact_differences": 0, "numeric_tolerance_failures": 0,
              "missing_pairs": 0, "missing_mismatches": 0,
              "text_pairs": 0, "text_mismatches": 0, "type_mismatches": 0,
              "max_absolute_difference": 0., "max_relative_difference_to_archived": 0.,
              "archived_zero_to_nonzero_count": 0, "examples": []}
    maxabs = maxrel = 0.
    for index, (a, b) in enumerate(zip(old, new)):
        ka, va = parse_cell(a)
        kb, vb = parse_cell(b)
        reason = None
        delta = rel = None
        if ka == kb == "missing":
            result["missing_pairs"] += 1
        elif "missing" in (ka, kb):
            result["missing_mismatches"] += 1
            reason = "missingness_changed"
        elif ka == kb == "numeric":
            result["numeric_pairs"] += 1
            if va != vb:
                result["numeric_exact_differences"] += 1
                delta = abs(vb-va)
                # Infinity is explicitly serialized; no epsilon hides zero-reference changes.
                rel = delta/abs(va) if va != 0 and math.isfinite(va) else math.inf
                maxabs, maxrel = max(maxabs, delta), max(maxrel, rel)
                if va == 0 and vb != 0:
                    result["archived_zero_to_nonzero_count"] += 1
                if not (math.isfinite(va) and math.isfinite(vb) and delta <= atol+rtol*abs(va)):
                    result["numeric_tolerance_failures"] += 1
                    reason = "numeric_outside_tolerance"
        elif ka == kb == "text":
            result["text_pairs"] += 1
            if va != vb:
                result["text_mismatches"] += 1
                reason = "text_changed"
        else:
            result["type_mismatches"] += 1
            reason = "numeric_text_type_changed"
        if reason and len(result["examples"]) < examples:
            result["examples"].append({"row_index_zero_based": index, "archived": a, "fresh": b,
                                       "reason": reason,
                                       "absolute_difference": finite_json(delta) if delta is not None else None,
                                       "relative_difference_to_archived": finite_json(rel) if rel is not None else None})
    result["max_absolute_difference"] = finite_json(maxabs)
    result["max_relative_difference_to_archived"] = finite_json(maxrel)
    failures = sum(result[k] for k in ("numeric_tolerance_failures", "missing_mismatches", "text_mismatches", "type_mismatches"))
    result["status"] = "differ" if failures else "tolerant" if result["numeric_exact_differences"] else "identical"
    return result


def compare_table(relative, archived_root, fresh_root, rtol, atol):
    old, new = archived_root/relative, fresh_root/relative
    record = {"source": relative, "archived_sha256": digest(old), "fresh_sha256": digest(new)}
    if not old.is_file() or not new.is_file():
        return {**record, "status": "missing", "missing_sides": [s for s, p in (("archived", old), ("fresh", new)) if not p.is_file()]}
    try:
        ac, ar, ah = read_table(old)
        fc, fr, fh = read_table(new)
        schema = ac == fc and ah == fh
        record.update(archived_rows=len(ar), fresh_rows=len(fr), row_count_match=len(ar) == len(fr),
                      archived_columns=ac, fresh_columns=fc, schema_match=schema,
                      archived_headerless=ah, fresh_headerless=fh)
        pairs = [(name, ac.index(name), fc.index(name)) for name in ac if name in fc]
        columns = [compare_column(name, [r[ia] for r in ar], [r[ib] for r in fr], rtol, atol)
                   for name, ia, ib in pairs]
        record["columns"] = columns
        record["differing_columns"] = [c["column"] for c in columns if c["status"] == "differ"]
        record["tolerant_columns"] = [c["column"] for c in columns if c["status"] == "tolerant"]
        record["numeric_pairs"] = sum(c["numeric_pairs"] for c in columns)
        record["numeric_tolerance_failures"] = sum(c["numeric_tolerance_failures"] for c in columns)
        record["numeric_exact_differences"] = sum(c["numeric_exact_differences"] for c in columns)
        record["missing_mismatches"] = sum(c["missing_mismatches"] for c in columns)
        record["text_mismatches"] = sum(c["text_mismatches"] for c in columns)
        record["type_mismatches"] = sum(c["type_mismatches"] for c in columns)
        if not schema or len(ar) != len(fr) or record["differing_columns"]:
            record["status"] = "differ"
        elif record["archived_sha256"] == record["fresh_sha256"]:
            record["status"] = "identical"
        else:
            record["status"] = "tolerant"
    except (OSError, ValueError, csv.Error) as exc:
        record.update(status="audit_fail", error=f"{type(exc).__name__}: {exc}")
    return record


def execution_audit(fresh_root, execution=None):
    path = Path(execution).resolve() if execution is not None else fresh_root.parent/"data_omnibus"/"execution.json"
    # Default execution.json lives in a run's data_omnibus subdirectory; the
    # optional explicit record defines a self-contained publication directory.
    log_root = path.parent if execution is not None else fresh_root.parent
    result = {"path": str(path), "sha256": digest(path), "log_root": str(log_root), "tasks": []}
    if not path.is_file():
        result.update(status="missing", missing_tasks=list(range(1, 11)))
        return result
    records = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(records, list):
        raise ValueError("execution.json must contain a list")
    seen = {}
    for record in records:
        task = record.get("task")
        if isinstance(task, bool) or not isinstance(task, int) or task not in range(1, 11) or task in seen:
            raise ValueError(f"Invalid or duplicate execution task: {task!r}")
        seen[task] = record
    missing = []
    for task in range(1, 11):
        if task not in seen:
            missing.append(task)
            result["tasks"].append({"task": task, "status": "missing"})
        else:
            record = dict(seen[task])
            log_value = record.get("log")
            log = (log_root/log_value).resolve() if isinstance(log_value, str) else None
            if log is not None and not log.is_relative_to(log_root.resolve()):
                raise ValueError(f"Execution log escapes run tree: {log_value}")
            record["log_exists"] = log is not None and log.is_file()
            record["log_sha256"] = digest(log) if log is not None else None
            result["tasks"].append(record)
    result["missing_tasks"] = missing
    result["status"] = "missing" if missing else "all_passed" if all(r.get("status") == "passed" and r.get("returncode") == 0 and r["log_exists"] for r in result["tasks"]) else "not_all_passed_or_logs_missing"
    return result


def audit(archived_root, fresh_root, provenance, rtol=1e-7, atol=1e-12, execution=None):
    result = {"schema": "ai4pharm-omnibus-recompute-audit-1", "status": "audit_fail",
              "archived_root": str(archived_root), "fresh_root": str(fresh_root),
              "figure_provenance": str(provenance), "figure_provenance_sha256": digest(provenance),
              "rtol": rtol, "atol": atol,
              "comparison_rule": "abs(fresh - archived) <= atol + rtol * abs(archived)",
              "row_alignment": "position; no sorting or join; reordered rows count as differences",
              "missing_tokens_case_insensitive": sorted(MISSING),
              "scope": "CSV files directly cited by ten figure records; no images, clock fields, or nonnumeric JSON metadata equality", "tasks": []}
    try:
        records = json.loads(provenance.read_text(encoding="utf-8-sig"))
        if not isinstance(records, list):
            raise ValueError("Figure provenance must contain a list")
        by_task = {}
        for record in records:
            task = record.get("task")
            if isinstance(task, bool) or not isinstance(task, int) or task not in range(1, 11) or task in by_task:
                raise ValueError(f"Invalid or duplicate provenance task: {task!r}")
            by_task[task] = record
        for task in range(1, 11):
            record = by_task.get(task)
            entry = {"task": task, "files": []}
            if record is None:
                entry.update(status="missing", reason="Task absent from figure provenance")
            else:
                sources = record.get("sources")
                if not isinstance(sources, list) or not all(isinstance(s, str) for s in sources):
                    raise ValueError(f"Invalid sources list for Task {task}")
                for relative in dict.fromkeys(sources):
                    if Path(relative).suffix.lower() != ".csv":
                        continue
                    for root in (archived_root, fresh_root):
                        full = (root/relative).resolve()
                        if not full.is_relative_to(root/"projects"):
                            raise ValueError(f"CSV source escapes projects tree: {relative}")
                    entry["files"].append(compare_table(relative, archived_root, fresh_root, rtol, atol))
                if not entry["files"]:
                    entry.update(status="missing", reason="No selected CSV sources for task")
                else:
                    entry["status"] = worst_status([f["status"] for f in entry["files"]])
            result["tasks"].append(entry)
        result["execution"] = execution_audit(fresh_root, execution)
        files = [f for t in result["tasks"] for f in t["files"]]
        result["file_status_counts"] = {s: sum(f["status"] == s for f in files) for s in STATUSES}
        result["selected_csv_files"] = len(files)
        result["numeric_pairs"] = sum(f.get("numeric_pairs", 0) for f in files)
        result["numeric_tolerance_failures"] = sum(f.get("numeric_tolerance_failures", 0) for f in files)
        result["table_comparison_status"] = worst_status([t["status"] for t in result["tasks"]])
        execution_status = result["execution"]["status"]
        result["status"] = worst_status([result["table_comparison_status"],
                                         "identical" if execution_status == "all_passed" else
                                         "missing" if execution_status == "missing" else "differ"])
        result["interpretation"] = "Reproducibility audit of selected tables only; equality does not validate the scientific models or their assumptions."
    except (OSError, ValueError, TypeError, KeyError, csv.Error) as exc:
        result.update(status="audit_fail", error=f"{type(exc).__name__}: {exc}")
    return result


def worst_status(statuses):
    return next((s for s in reversed(STATUSES) if s in statuses), "missing")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fresh-root", type=Path, required=True, help="Canonical fresh tree containing projects/")
    parser.add_argument("--out", type=Path, required=True, help="NEW audit JSON path; existing files are refused")
    parser.add_argument("--archive-root", type=Path, default=ROOT)
    parser.add_argument("--provenance", type=Path, help="Defaults to ARCHIVE_ROOT/data_omnibus/figure_provenance.json")
    parser.add_argument("--execution", type=Path, help="Execution JSON; explicit record resolves relative logs from its parent. Default: FRESH_ROOT.parent/data_omnibus/execution.json, logs from FRESH_ROOT.parent")
    parser.add_argument("--rtol", type=float, default=1e-7)
    parser.add_argument("--atol", type=float, default=1e-12)
    args = parser.parse_args(argv)
    if not all(math.isfinite(v) and v >= 0 for v in (args.rtol, args.atol)):
        parser.error("Tolerances must be finite and nonnegative")
    archived, fresh, out = args.archive_root.resolve(), args.fresh_root.resolve(), args.out.resolve()
    provenance = args.provenance.resolve() if args.provenance else archived/"data_omnibus/figure_provenance.json"
    if out.suffix.lower() != ".json":
        parser.error("--out must be a .json path")
    if out.exists() or any(out.is_relative_to(root/"projects") for root in (archived, fresh)):
        parser.error("--out must be new and outside both input projects trees")
    result = audit(archived, fresh, provenance, args.rtol, args.atol, args.execution)
    result["audit_script_sha256"] = digest(Path(__file__))
    out.parent.mkdir(parents=True, exist_ok=True)
    # Exclusive creation also prevents a race from overwriting another artifact.
    with out.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(result, handle, ensure_ascii=False, allow_nan=False, indent=2)
        handle.write("\n")
    print(json.dumps({"status": result["status"], "selected_csv_files": result.get("selected_csv_files"),
                      "file_status_counts": result.get("file_status_counts"), "out": str(out)}, ensure_ascii=False))
    return 2 if result["status"] == "audit_fail" else 0


if __name__ == "__main__":
    raise SystemExit(main())
