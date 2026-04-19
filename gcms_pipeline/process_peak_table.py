#!/usr/bin/env python3
"""Normalize GC-MS peak tables for the cuticular compound review interface."""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path


FIELD_ALIASES = {
    "rt": ["rt", "retention_time", "retention time", "time"],
    "area": ["area", "intensity", "height", "abundance"],
    "name": ["name", "compound", "candidate"],
    "formula": ["formula"],
    "class": ["class", "compound_class", "type"],
    "match": ["match", "score", "similarity"],
    "spectrum": ["spectrum", "ions", "mass_spectrum"],
    "notes": ["notes", "comment"],
}


def sniff_dialect(path: Path) -> csv.Dialect:
    sample = path.read_text(encoding="utf-8-sig")[:4096]
    try:
        return csv.Sniffer().sniff(sample, delimiters=",\t")
    except csv.Error:
        return csv.excel_tab if "\t" in sample.splitlines()[0] else csv.excel


def normalized_key(name: str) -> str:
    return name.strip().lower()


def field_lookup(headers: list[str]) -> dict[str, str | None]:
    normalized = {normalized_key(header): header for header in headers}
    lookup: dict[str, str | None] = {}
    for canonical, aliases in FIELD_ALIASES.items():
        lookup[canonical] = next((normalized[a] for a in aliases if a in normalized), None)
    return lookup


def parse_spectrum(text: str) -> list[tuple[float, float]]:
    ions: list[tuple[float, float]] = []
    for mz, intensity in re.findall(r"(\d+(?:\.\d+)?)\s*[:=]\s*(\d+(?:\.\d+)?)", text or ""):
        ions.append((float(mz), float(intensity)))
    return ions


def has_ion(ions: list[tuple[float, float]], target: float, threshold: float = 25) -> bool:
    return any(abs(mz - target) <= 1 and intensity >= threshold for mz, intensity in ions)


def infer_class(name: str, supplied_class: str, spectrum: str) -> str:
    supplied = supplied_class.strip().lower()
    if supplied and supplied not in {"unknown", "unidentified"}:
        return supplied_class.strip()

    ions = parse_spectrum(spectrum)
    lowered_name = name.lower()

    if has_ion(ions, 74) and has_ion(ions, 43):
        return "ester"
    if has_ion(ions, 55) and has_ion(ions, 69):
        return "alkene"
    if has_ion(ions, 57) and has_ion(ions, 71) and has_ion(ions, 85):
        if "methyl" in lowered_name or "branched" in lowered_name:
            return "methyl-branched alkane"
        return "alkane"
    if has_ion(ions, 43) and has_ion(ions, 58):
        return "other"
    return supplied_class.strip() or "unknown"


def confidence_from_score(score: str, compound_class: str) -> str:
    try:
        numeric = float(score)
    except ValueError:
        numeric = 0

    if compound_class == "unknown":
        return "unreviewed"
    if numeric >= 85:
        return "probable compound"
    if numeric >= 70:
        return "tentative compound"
    return "tentative class"


def read_rows(path: Path) -> list[dict[str, str]]:
    dialect = sniff_dialect(path)
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, dialect=dialect)
        if not reader.fieldnames:
            return []
        lookup = field_lookup(reader.fieldnames)
        output = []
        for index, row in enumerate(reader, start=1):
            get = lambda field: row.get(lookup[field] or "", "").strip()
            spectrum = get("spectrum")
            name = get("name") or "Unknown"
            compound_class = infer_class(name, get("class"), spectrum)
            match = get("match")
            notes = get("notes")
            if compound_class != (get("class") or "unknown"):
                note = f"Class hint from fragment pattern: {compound_class}."
                notes = f"{notes} {note}".strip()
            output.append(
                {
                    "id": f"P{index:02d}",
                    "rt": get("rt") or str(index),
                    "area": get("area") or "0",
                    "name": name,
                    "formula": get("formula"),
                    "class": compound_class,
                    "match": match,
                    "spectrum": spectrum,
                    "confidence": confidence_from_score(match, compound_class),
                    "notes": notes,
                }
            )
        return output


def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    fields = ["id", "rt", "area", "name", "formula", "class", "match", "spectrum", "confidence", "notes"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print("Usage: process_peak_table.py input.csv output.csv", file=sys.stderr)
        return 2

    input_path = Path(argv[1])
    output_path = Path(argv[2])
    rows = read_rows(input_path)
    write_rows(output_path, rows)
    print(f"Wrote {len(rows)} peaks to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))

