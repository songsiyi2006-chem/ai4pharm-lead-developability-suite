"""Verify the integrated atlas against its data, code and artifact manifests."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from omnibus.catalog import FIGURES, PROJECTS, DRIVERS


def load(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_map(base, mapping, label, errors):
    for name, expected in mapping.items():
        path = (base/name).resolve()
        if not path.is_relative_to(base.resolve()) or not path.is_file():
            errors.append(f"{label}: missing or escaping path {name}")
        elif digest(path) != expected:
            errors.append(f"{label}: checksum mismatch {name}")


def inspect(bundle=ROOT, source=ROOT, reports=ROOT):
    errors = []
    try:
        manifest = load(bundle/"manifest_omnibus.json")
        if manifest.get("schema") != "ai4pharm-omnibus-1":
            errors.append("Unsupported manifest schema")
        for label, base, key in (("artifact", bundle, "files_sha256"), ("code", ROOT, "source_code_sha256"),
                                 ("driver", ROOT, "scientific_driver_sha256")):
            mapping = manifest.get(key, {})
            if not mapping:
                errors.append(f"Empty {key}")
            verify_map(base, mapping, label, errors)
        expected_drivers = {f"projects/{p}/{d}" for p, d in zip(PROJECTS, DRIVERS)}
        if set(manifest.get("scientific_driver_sha256", {})) != expected_drivers:
            errors.append("Scientific driver provenance must cover exactly ten modules")
        source_hashes = load(bundle/"data_omnibus/source_hashes.json")
        verify_map(source, source_hashes, "input", errors)
        records = load(bundle/"data_omnibus/figure_provenance.json")
        if sorted(r["task"] for r in records) != list(range(1, 11)):
            errors.append("Ten distinct task figures required")
        inputs = set()
        image_records = []
        for r in records:
            i = r["task"]-1
            if not 0 <= i < 10 or r["file"] != FIGURES[i]:
                errors.append(f"Wrong figure task mapping: {r['file']}")
                continue
            name = "figures_omnibus/"+r["file"]
            path = bundle/name
            if manifest["files_sha256"].get(name) != r["sha256"]:
                errors.append(f"Figure provenance/manifest mismatch: {name}")
            with Image.open(path) as im:
                dpi = im.info.get("dpi", ())
                if len(dpi) != 2 or any(abs(v-300) > .02 for v in dpi):
                    errors.append(f"Not 300 DPI: {name}")
                if min(im.size) < 1500:
                    errors.append(f"Insufficient multipanel image dimensions: {name}")
                image_records.append({"task": i+1, "pixels": list(im.size), "dpi": list(dpi)})
                im.verify()
            if len(r.get("panels", [])) < 3 or not r.get("caveats") or not r.get("sources"):
                errors.append(f"Incomplete panels, caveats or sources: Task {i+1}")
            inputs.update(r["sources"])
            for derived in r.get("derived_files", [])+r.get("derived_data", []):
                if derived not in manifest["files_sha256"]:
                    errors.append(f"Unmanifested derived table: {derived}")
        if inputs != set(source_hashes):
            errors.append("Source hash coverage must equal all figure inputs")
        required = {"index.html", "data_omnibus/figure_provenance.json", "data_omnibus/source_hashes.json",
                    "data_omnibus/run_summary.json", *("figures_omnibus/"+x for x in FIGURES)}
        if required-set(manifest["files_sha256"]):
            errors.append("Manifest omits required atlas artifacts")
        for lang in ("ZH", "EN"):
            name = f"AI4PHARM_DECADE_TREATISE_{lang}.md"
            report = (reports/name).read_text(encoding="utf-8-sig")
            for filename in FIGURES:
                if "figures_omnibus/"+filename not in report:
                    errors.append(f"Missing figure link in {name}: {filename}")
            if len(report) < 15000:
                errors.append(f"Master report unexpectedly short: {name}")
        result = {"validator": "validate_omnibus", "status": "passed" if not errors else "failed",
                  "figure_count": len(records), "input_files": len(inputs), "figure_metadata": image_records,
                  "artifact_hashes_checked": len(manifest["files_sha256"]), "errors": errors,
                  "scope": "Artifact/source integrity, figure metadata and report coverage; not biological validation"}
    except (OSError, ValueError, KeyError, TypeError) as exc:
        result = {"validator": "validate_omnibus", "status": "failed", "errors": errors+[str(exc)]}
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bundle", type=Path, default=ROOT)
    p.add_argument("--source-root", type=Path, default=ROOT)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    destination = a.out.resolve()
    if destination.exists():
        if destination.suffix != ".json" or load(destination).get("validator") != "validate_omnibus":
            p.error("Refusing to overwrite a non-validator artifact")
    protected = [ROOT/"projects", ROOT/"docs/archive", a.bundle.resolve()/"data_omnibus", a.bundle.resolve()/"figures_omnibus"]
    if destination == a.bundle.resolve()/"manifest_omnibus.json" or any(destination.is_relative_to(x.resolve()) for x in protected):
        p.error("Validation output must not overwrite scientific inputs or bundle artifacts")
    result = inspect(a.bundle.resolve(), a.source_root.resolve())
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, indent=2, ensure_ascii=False)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps({k: v for k, v in result.items() if k != "figure_metadata"}, ensure_ascii=False))
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
