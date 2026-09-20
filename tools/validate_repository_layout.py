"""Read-only validation of this repository's seven-project migration.

Only the explicitly requested JSON report is written. No scientific computation,
network access, Git mutation, or change to archived evidence is performed.
"""
from __future__ import annotations

import argparse
import hashlib
import html
import json
from pathlib import Path, PurePosixPath
import re
import subprocess
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
BASELINE = "3785df88d244bf276de7f9b6e4d1d74995b530b9"
PATH_MAP = ROOT / "docs/reorganization/path_map.json"
PROJECTS = ("task01_lead_developability", "task02_tpd", "task03_covalent_kinetics",
            "task04_pbpk", "task05_cyp_ddi", "task06_asd", "task07_qsp")
# Manifest within project, artifact map, artifact base, driver.
MANIFESTS = (
    ("results_task1/run_manifest.json", "artifacts", ".", "run_task1_mpo_admet_developability.py"),
    ("manifest_task2.json", "sha256", ".", "run_task2_tpd_ternary_cooperativity.py"),
    ("data_task3/run_manifest.json", "sha256", ".", "run_task3_covalent_kinetics_residence_time.py"),
    ("manifest.json", "files_sha256", ".", "run_task4_pbpk_pharmacokinetics_dose_prediction.py"),
    ("results/manifest.json", "files", "results", "run_task5_cyp_ddi_mechanism_based_inhibition.py"),
    ("manifest.json", "files", ".", "run_task6_asd_formulation_supersaturation_kinetics.py"),
    ("manifest.json", "files_sha256", ".", "run_task7_qsp_tumor_immune_pkpd_synergy.py"),
)
SCIENCE_SUFFIXES = {".csv", ".png", ".svg", ".sdf", ".xyz"}
SCIENCE_JSON_NAMES = {
    "reference_panel.json", "hibit_hill_fit.json", "lfer_regression.json",
    "parameters_template.json", "developability_results.json", "compound_example.json",
    "compound_panel_example.json", "parameters_used.json", "parameters.json",
    "results_summary.json", "summary.json", "survival_statistics.json", "sources.json",
    "verification.json", "release_verification.json",
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path):
    return path.relative_to(ROOT).as_posix()


def bounded_path(base, name):
    candidate = (base / name).resolve()
    if not candidate.is_relative_to(ROOT):
        raise ValueError(f"Path escapes this repository: {name}")
    return candidate


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT, stderr=subprocess.PIPE)


def inventory():
    actual = Path(git("rev-parse", "--show-toplevel").decode().strip()).resolve()
    if actual != ROOT:
        raise ValueError("Validator must be in this repository's tools directory")
    tracked = set(git("ls-files", "-z").decode("utf-8").split("\0")) - {""}
    listed = set(git("ls-files", "--cached", "--others", "--exclude-standard", "-z").decode("utf-8").split("\0")) - {""}
    return tracked, {p for p in listed if bounded_path(ROOT, p).is_file()}


def load_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def error(errors, check, **details):
    errors.append({"check": check, **details})


def migration_checks(tracked, included, errors):
    mapping = load_json(PATH_MAP)
    if mapping.get("baseline_commit") != BASELINE:
        error(errors, "baseline_commit", expected=BASELINE, actual=mapping.get("baseline_commit"))
    records = mapping.get("records", [])
    if not isinstance(records, list) or not records:
        raise ValueError("path_map.json must contain a nonempty records list")
    old_paths, new_paths, checked, missing, duplicates, unchanged = [], [], [], [], [], []
    for record in records:
        old, new, previous = (record[k] for k in ("old_path", "new_path", "previous_sha256"))
        if not re.fullmatch(r"[0-9a-f]{64}", previous):
            error(errors, "path_map_sha256_format", old_path=old)
        old_paths.append(old)
        new_paths.append(new)
        target = bounded_path(ROOT, new)
        if not target.is_file():
            missing.append(new)
            error(errors, "migration_target_missing", old_path=old, new_path=new)
            continue
        if new not in included:
            error(errors, "migration_target_not_in_git_inventory", new_path=new,
                  note="Target is ignored and untracked; ordinary staging would omit it")
        checked.append(new)
        if old != new and old in tracked and bounded_path(ROOT, old).is_file():
            duplicates.append(old)
            error(errors, "tracked_old_path_duplicate", old_path=old, new_path=new)
        suffix = PurePosixPath(new).suffix.lower()
        scientific = suffix in SCIENCE_SUFFIXES or (suffix == ".json" and PurePosixPath(new).name in SCIENCE_JSON_NAMES)
        archive = new.startswith("docs/archive/")
        prompt = PurePosixPath(new).name == "TASK.md"
        if scientific or archive or prompt:
            actual = sha(target)
            reason = "scientific_artifact" if scientific else "frozen_archive" if archive else "original_prompt"
            unchanged.append({"path": new, "category": reason, "sha256": actual, "byte_identical": actual == previous})
            if actual != previous:
                error(errors, "migrated_bytes_changed", path=new, category=reason, expected=previous, actual=actual)
    if len(set(old_paths)) != len(old_paths):
        error(errors, "duplicate_old_paths_in_map")
    if len(set(new_paths)) != len(new_paths):
        error(errors, "multiple_old_paths_map_to_one_target")
    try:
        baseline_files = set(git("ls-tree", "-r", "--name-only", "-z", BASELINE).decode("utf-8").split("\0")) - {""}
        omitted, extra = sorted(baseline_files - set(old_paths)), sorted(set(old_paths) - baseline_files)
        if omitted or extra:
            error(errors, "baseline_mapping_coverage", omitted_baseline_paths=omitted, unexpected_old_paths=extra)
    except subprocess.CalledProcessError:
        error(errors, "baseline_tree_unavailable", commit=BASELINE,
              action="Supply local baseline history; CI checkout should use fetch-depth: 0. Validator does not fetch.")
        baseline_files = set()
    return {"baseline_commit": BASELINE, "baseline_files": len(baseline_files), "map_records": len(records),
            "targets_checked": len(checked), "missing_targets": missing, "tracked_old_duplicates": duplicates,
            "preserved_artifact_checks": unchanged}


def unfenced_markdown(text):
    lines, fence = [], None
    for line in text.splitlines():
        match = re.match(r"^ {0,3}(`{3,}|~{3,})", line)
        if match:
            token = match.group(1)
            if fence is None:
                fence = token
            elif token[0] == fence[0] and len(token) >= len(fence):
                fence = None
            continue
        if fence is None:
            lines.append(line)
    return "\n".join(lines)


def destination(raw):
    raw = raw.strip()
    if raw.startswith("<"):
        return raw[1:raw.find(">")].strip() if ">" in raw else raw[1:]
    return re.split(r'\s+[\"\']', raw, maxsplit=1)[0].strip()


def markdown_destinations(text):
    """Inline/image links with balanced parentheses, references, HTML href/src."""
    text = re.sub(r"(`+).*?\1", "", unfenced_markdown(text))
    found = []
    for match in re.finditer(r"!?\[[^\]\n]*\]\(", text):
        start, depth, escaped, angle = match.end(), 1, False, False
        for pos in range(start, len(text)):
            char = text[pos]
            if escaped:
                escaped = False
            elif char == "\\":
                escaped = True
            elif char == "<":
                angle = True
            elif char == ">":
                angle = False
            elif not angle and char == "(":
                depth += 1
            elif not angle and char == ")":
                depth -= 1
                if depth == 0:
                    found.append(destination(text[start:pos]))
                    break
    found.extend(destination(m.group(1)) for m in re.finditer(r"^ {0,3}\[[^\]\n]+\]:\s*(.+)$", text, re.M))
    found.extend(m.group(2) for m in re.finditer(r"\b(?:href|src)\s*=\s*([\"'])(.*?)\1", text, re.I))
    return found


def markdown_anchors(text):
    text = unfenced_markdown(text)
    anchors = set(re.findall(r"\b(?:id|name)\s*=\s*[\"']([^\"']+)[\"']", text, flags=re.I))
    counts, headings = {}, []
    lines = text.splitlines()
    for index, line in enumerate(lines):
        heading = re.match(r"^ {0,3}#{1,6}\s+(.+?)\s*#*\s*$", line)
        if heading:
            headings.append(heading.group(1))
        elif index + 1 < len(lines) and line.strip() and re.fullmatch(r" {0,3}(?:=+|-+)\s*", lines[index + 1]):
            headings.append(line.strip())
    for title in headings:
        title = html.unescape(re.sub(r"<[^>]*>", "", title)).lower()
        title = re.sub(r"\[([^\]]+)\]\([^)]+\)", r"\1", title)
        slug = re.sub(r"[^\w\- ]", "", title).replace(" ", "-")
        duplicate = counts.get(slug, 0)
        counts[slug] = duplicate + 1
        anchors.add(slug if duplicate == 0 else f"{slug}-{duplicate}")
    return anchors


def check_markdown_links(included, errors):
    documents, checked, excluded, anchors, cli_scripts = [], [], [], {}, []
    for name in sorted(included):
        if not name.lower().endswith(".md"):
            continue
        if name.startswith("docs/archive/") or PurePosixPath(name).name == "TASK.md":
            excluded.append(name)
            continue
        source = bounded_path(ROOT, name)
        documents.append(name)
        content = source.read_text(encoding="utf-8-sig")
        # Commands are documentation too: inspect script arguments even in fences.
        # Output/config paths can intentionally be new and are not checked here.
        for match in re.finditer(r'''\b(?:python(?:\d+(?:\.\d+)*)?|py)\s+([^\s`"'<>]+\.py)\b''', content):
            script = match.group(1)
            try:
                root_target = bounded_path(ROOT, script)
                local_target = bounded_path(source.parent, script)
                target = root_target if root_target.is_file() else local_target
                cli_scripts.append({"document": name, "script": script, "resolved": relative(target)})
                if not target.is_file():
                    error(errors, "documented_python_script_missing", document=name, script=script)
            except ValueError:
                error(errors, "documented_python_script_outside_repository", document=name, script=script)
        for raw in markdown_destinations(content):
            raw = html.unescape(re.sub(r"\\([() ])", r"\1", raw))
            parsed = urlsplit(raw)
            if not raw or parsed.scheme or parsed.netloc:
                continue
            try:
                if not parsed.path:
                    target = source
                elif parsed.path.startswith("/"):
                    target = bounded_path(ROOT, unquote(parsed.path).lstrip("/"))
                else:
                    target = bounded_path(source.parent, unquote(parsed.path))
            except ValueError:
                error(errors, "local_markdown_link_outside_repository", document=name, target=raw)
                continue
            record = {"document": name, "target": raw, "resolved": relative(target)}
            checked.append(record)
            if not target.exists():
                error(errors, "local_markdown_target_missing", **record)
            elif target.is_file() and target.suffix.lower() == ".md" and parsed.fragment:
                if target not in anchors:
                    anchors[target] = markdown_anchors(target.read_text(encoding="utf-8-sig"))
                if unquote(parsed.fragment) not in anchors[target]:
                    error(errors, "local_markdown_anchor_missing", **record)
    return {"documents_checked": documents, "local_links_checked": len(checked),
            "links": checked, "python_cli_script_references": cli_scripts,
            "excluded_frozen_documents": excluded}


def check_project_manifests(errors, selected=PROJECTS):
    checked = []
    for project, (manifest_name, key, artifact_base, driver) in zip(PROJECTS, MANIFESTS):
        if project not in selected:
            continue
        folder = ROOT / "projects" / project
        path = folder / manifest_name
        record = {"project": project, "manifest": relative(path), "artifact_hashes_checked": 0,
                  "current_source_hashes_checked": 0}
        checked.append(record)
        if not path.is_file():
            error(errors, "manifest_missing", path=relative(path))
            continue
        manifest = load_json(path)
        hashes = manifest.get(key)
        if not isinstance(hashes, dict) or not hashes:
            error(errors, "manifest_artifact_map_missing", path=relative(path), key=key)
            continue
        for name, expected in hashes.items():
            try:
                target = bounded_path(folder / artifact_base, name)
                if not target.is_relative_to(folder):
                    raise ValueError("Artifact path outside its project")
            except ValueError as exc:
                error(errors, "manifest_artifact_path_invalid", manifest=relative(path), target=name, detail=str(exc))
                continue
            if not target.is_file():
                error(errors, "manifest_artifact_missing", manifest=relative(path), target=name)
                continue
            actual = sha(target)
            record["artifact_hashes_checked"] += 1
            if actual != expected:
                error(errors, "manifest_artifact_hash_mismatch", path=relative(target), expected=expected, actual=actual)
        migration = manifest.get("layout_migration", {})
        sources = migration.get("current_sources_sha256")
        if not isinstance(sources, dict) or not sources:
            error(errors, "manifest_current_sources_missing", manifest=relative(path),
                  note="Original script_sha256 remains execution provenance; layout_migration records current sources")
            continue
        if migration.get("baseline_commit") != BASELINE:
            error(errors, "manifest_migration_baseline", manifest=relative(path))
        required = {relative(folder / driver)}
        if project == "task05_cyp_ddi":
            required.add("projects/task04_pbpk/run_task4_pbpk_pharmacokinetics_dose_prediction.py")
        if not required.issubset(sources):
            error(errors, "manifest_required_sources_missing", manifest=relative(path), missing=sorted(required - set(sources)))
        for name, expected in sources.items():
            try:
                target = bounded_path(ROOT, name)
            except ValueError as exc:
                error(errors, "manifest_source_path_invalid", manifest=relative(path), target=name, detail=str(exc))
                continue
            if not target.is_file():
                error(errors, "manifest_current_source_missing", manifest=relative(path), target=name)
                continue
            record["current_source_hashes_checked"] += 1
            actual = sha(target)
            if actual != expected:
                error(errors, "manifest_current_source_hash_mismatch", path=name, expected=expected, actual=actual)
    return checked


def inspect():
    errors = []
    result = {"validator": "validate_repository_layout", "scope": "Current local repository only; no scientific validity claim"}
    try:
        tracked, included = inventory()
        result["migration"] = migration_checks(tracked, included, errors)
        readmes = [f"projects/{project}/README.md" for project in PROJECTS]
        result["project_readmes"] = readmes
        for name in readmes:
            if name not in included:
                error(errors, "project_readme_missing_or_ignored", path=name)
        result["markdown"] = check_markdown_links(included, errors)
        result["manifests"] = check_project_manifests(errors)
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as exc:
        error(errors, "validator_cannot_complete", detail=str(exc))
    result.update(status="passed" if not errors else "failed", errors=errors)
    return result


def write_report(output, result):
    path = bounded_path(ROOT, output)
    if path.is_relative_to(ROOT / "docs/archive"):
        raise ValueError("Archived evidence is frozen and cannot be a report destination")
    if path.exists():
        if not path.is_file():
            raise ValueError("--out must name a FILE, not a directory")
        try:
            previous = load_json(path)
        except (ValueError, UnicodeError):
            raise ValueError("Refusing to replace an existing non-validation file") from None
        if not isinstance(previous, dict) or previous.get("validator") != result["validator"]:
            raise ValueError("Refusing to replace another artifact; choose a new filename")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True, help="Repository-relative JSON report FILE")
    args = parser.parse_args()
    result = inspect()
    try:
        output = write_report(args.out, result)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(json.dumps({"status": result["status"], "report": relative(output), "errors": result["errors"]}, ensure_ascii=False, indent=2))
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
