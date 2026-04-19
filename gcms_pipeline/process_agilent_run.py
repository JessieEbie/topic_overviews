#!/usr/bin/env python3
"""Build a browser-workstation bundle from an Agilent .D folder or zipped .D run."""

from __future__ import annotations

import json
import math
import re
import struct
import sys
import zipfile
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Scan:
    time_min: float
    tic: float
    ions: list[tuple[float, float]]


def read_agilent_string(data: bytes, offset: int, gap: int = 1) -> str:
    length = data[offset] * gap
    chars: list[str] = []
    for i in range(0, length, gap):
        code = data[offset + 1 + i]
        if code:
            chars.append(chr(code))
    return "".join(chars).strip()


def parse_agilent_ms_buffer(data: bytes) -> tuple[list[Scan], dict[str, str]]:
    if struct.unpack_from(">I", data, 0)[0] != 0x01320000:
        raise ValueError("This DATA.MS file does not match the expected ChemStation header.")
    type_name = read_agilent_string(data, 0x4, 1)
    num_times = struct.unpack_from(">I", data, 0x116)[0] if type_name == "MSD Spectral File" else struct.unpack_from("<H", data, 0x142)[0]
    offset = struct.unpack_from(">H", data, 0x10A)[0] * 2 - 2
    scans: list[Scan] = []
    for _ in range(num_times):
        if offset + 18 >= len(data):
            break
        offset += 2
        time = struct.unpack_from(">I", data, offset)[0]
        offset += 4
        offset += 6
        pair_count = struct.unpack_from(">H", data, offset)[0]
        offset += 2
        offset += 4
        ions: list[tuple[float, float]] = []
        tic = 0.0
        for _ in range(pair_count):
            if offset + 4 > len(data):
                break
            mz = struct.unpack_from(">H", data, offset)[0] / 20
            offset += 2
            enc = struct.unpack_from(">H", data, offset)[0]
            offset += 2
            intensity = (8 ** (enc >> 14)) * (enc & 0x3FFF)
            ions.append((mz, float(intensity)))
            tic += intensity
        offset += 10
        scans.append(Scan(time_min=time / 60000, tic=tic, ions=ions))
    metadata = {
        "type": type_name,
        "date": read_agilent_string(data, 0xB2, 1),
        "method": read_agilent_string(data, 0xE4, 1),
    }
    return scans, metadata


def parse_ini_text(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    section = ""
    for line in text.splitlines():
        trimmed = line.strip()
        if not trimmed or trimmed.startswith(";"):
            continue
        section_match = re.match(r"^\[(.+)\]$", trimmed)
        if section_match:
            section = section_match.group(1)
            continue
        if "=" not in trimmed:
            continue
        key, value = trimmed.split("=", 1)
        result[f"{section}.{key.strip()}" if section else key.strip()] = value.strip()
    return result


def parse_runstart_text(text: str) -> dict[str, str]:
    meta: dict[str, str] = {}
    for line in text.splitlines():
        match = re.match(r"^\s*([^=]+?)\s*=\s*(.+)\s*$", line)
        if not match:
            continue
        key, value = match.group(1).strip(), match.group(2).strip()
        if "Methfile" in key:
            meta["methodFile"] = value
        elif "Sample Name" in key:
            meta["sampleName"] = value
        elif "Miscellaneous" in key:
            meta["miscInfo"] = value
        elif "Datafile" in key:
            meta["dataFile"] = value
    return meta


def parse_acqmeth_text(text: str) -> dict[str, str]:
    meta: dict[str, str] = {}
    method_path = re.search(r"^[ \t]*([A-Z]:\\[^\r\n]+\.M)\s*$", text, re.MULTILINE)
    if method_path:
        meta["methodPath"] = method_path.group(1).strip()
    run_time = re.search(r"Run time:\s+([0-9.]+)\s+min", text, re.IGNORECASE)
    if run_time:
        meta["runTimeMin"] = run_time.group(1)
    column = re.search(r"DB-1ms[^\r\n]*", text, re.IGNORECASE)
    if column:
        meta["column"] = column.group(0).strip()
    return meta


def parse_associated_metadata(file_text_map: dict[str, str], run_name: str, ms_metadata: dict[str, str]) -> dict[str, str]:
    result = {
        "runName": run_name,
        "acquiredAt": ms_metadata.get("date", ""),
        "methodName": ms_metadata.get("method", ""),
        "sampleName": "",
        "miscInfo": "",
        "methodFile": "",
        "scanRange": "",
        "sourceTemp": "",
        "quadTemp": "",
        "pressure": "",
        "runTimeMin": "",
    }
    for path, text in file_text_map.items():
        lowered = path.lower()
        if lowered.endswith("runstart.txt"):
            parsed = parse_runstart_text(text)
            result["methodFile"] = parsed.get("methodFile", result["methodFile"])
            result["sampleName"] = parsed.get("sampleName", result["sampleName"])
            result["miscInfo"] = parsed.get("miscInfo", result["miscInfo"])
        elif lowered.endswith("acqmeth.txt"):
            parsed = parse_acqmeth_text(text)
            result["methodFile"] = parsed.get("methodPath", result["methodFile"])
            result["runTimeMin"] = parsed.get("runTimeMin", result["runTimeMin"])
            result["column"] = parsed.get("column", result.get("column", ""))
        elif lowered.endswith("pre_post.ini"):
            parsed = parse_ini_text(text)
            result["miscInfo"] = parsed.get("POSTRUN.Miscinfo", result["miscInfo"])
            result["acquiredAt"] = parsed.get("POSTRUN.Date", result["acquiredAt"])
            low = parsed.get("Scan Parameters.lowmass1")
            high = parsed.get("Scan Parameters.highmass1")
            if low and high:
                result["scanRange"] = f"m/z {low}-{high}"
            if parsed.get("POSTRUN.SourceTemp"):
                result["sourceTemp"] = f"{parsed['POSTRUN.SourceTemp']} C"
            if parsed.get("POSTRUN.QuadTemp"):
                result["quadTemp"] = f"{parsed['POSTRUN.QuadTemp']} C"
            result["pressure"] = parsed.get("POSTRUN.Pressure", result["pressure"])
    if not result["methodName"] and result["methodFile"]:
        result["methodName"] = re.split(r"[\\/]", result["methodFile"])[-1].removesuffix(".M")
    return result


def smooth_values(values: list[float], radius: int = 7) -> list[float]:
    smoothed: list[float] = []
    for i in range(len(values)):
        lo = max(0, i - radius)
        hi = min(len(values) - 1, i + radius)
        chunk = values[lo : hi + 1]
        smoothed.append(sum(chunk) / len(chunk))
    return smoothed


def spectrum_from_scan_window(scans: list[Scan], left: int, right: int, center_index: int) -> str:
    lo = max(left, center_index - 2)
    hi = min(right, center_index + 2)
    bins: dict[int, float] = {}
    for i in range(lo, hi + 1):
        distance = abs(center_index - i)
        weight = 1.0 if distance == 0 else 0.7 if distance == 1 else 0.45
        for mz, intensity in scans[i].ions:
            key = round(mz)
            bins[key] = bins.get(key, 0.0) + intensity * weight
    ranked = sorted(bins.items(), key=lambda item: item[1], reverse=True)
    top_by_intensity = ranked[:22]
    high_mz_tail = sorted(bins.items(), key=lambda item: item[0], reverse=True)[:8]
    keep: dict[int, float] = {}
    for mz, intensity in [*top_by_intensity, *high_mz_tail]:
        keep[int(mz)] = max(intensity, keep.get(int(mz), 0.0))
    ions = list(keep.items())
    if not ions:
        return ""
    max_intensity = max(intensity for _, intensity in ions)
    normalized = [(mz, max(1, round(intensity / max_intensity * 100))) for mz, intensity in ions]
    normalized.sort(key=lambda item: item[0])
    return " ".join(f"{mz}:{intensity}" for mz, intensity in normalized)


def parse_spectrum(text: str) -> list[tuple[float, float]]:
    ions: list[tuple[float, float]] = []
    for mz, intensity in re.findall(r"(\d+(?:\.\d+)?)\s*[:=]\s*(\d+(?:\.\d+)?)", text or ""):
        ions.append((float(mz), float(intensity)))
    return ions


def has_ion(ions: list[tuple[float, float]], target: float, threshold: float = 25) -> bool:
    return any(abs(mz - target) <= 1 and intensity >= threshold for mz, intensity in ions)


def infer_class(spectrum: str) -> str:
    ions = parse_spectrum(spectrum)
    if has_ion(ions, 74) and has_ion(ions, 43):
        return "ester"
    if has_ion(ions, 55) and has_ion(ions, 69):
        return "alkene"
    if has_ion(ions, 57) and has_ion(ions, 71) and has_ion(ions, 85):
        return "alkane"
    if has_ion(ions, 43) and has_ion(ions, 58):
        return "other"
    return "unknown"


def detect_raw_peaks(scans: list[Scan], run_meta: dict[str, str]) -> list[dict[str, object]]:
    if not scans:
        return []
    tic_values = [scan.tic for scan in scans]
    smooth = smooth_values(tic_values, 7)
    max_value = max(smooth)
    candidates: list[dict[str, float | int]] = []
    for i in range(4, len(smooth) - 4):
        current = smooth[i]
        if current > smooth[i - 1] and current >= smooth[i + 1] and current > max_value * 0.045:
            left, right = i, i
            while left > 2 and smooth[left - 1] <= smooth[left]:
                left -= 1
            while right < len(smooth) - 3 and smooth[right + 1] <= smooth[right]:
                right += 1
            left_floor = min(smooth[max(0, left - 10) : left + 1])
            right_floor = min(smooth[right : min(len(smooth), right + 11)])
            prominence = current - max(left_floor, right_floor)
            if prominence < max_value * 0.012:
                continue
            if candidates and i - int(candidates[-1]["index"]) < 18:
                if current > float(candidates[-1]["height"]):
                    candidates[-1] = {"index": i, "height": current, "left": left, "right": right, "prominence": prominence}
            else:
                candidates.append({"index": i, "height": current, "left": left, "right": right, "prominence": prominence})

    candidates = sorted(candidates, key=lambda item: (-float(item["prominence"]), -float(item["height"])))[:36]
    candidates.sort(key=lambda item: int(item["index"]))
    peaks: list[dict[str, object]] = []
    for idx, candidate in enumerate(candidates, start=1):
        left = int(candidate["left"])
        right = int(candidate["right"])
        center = int(candidate["index"])
        baseline = max(smooth[left], smooth[right])
        area = sum(max(0.0, tic_values[i] - baseline) for i in range(left, right + 1))
        weighted = [max(0.0, scans[i].tic - baseline) for i in range(left, right + 1)]
        weight_total = sum(weighted)
        rt = (
            sum(scans[left + offset].time_min * weight for offset, weight in enumerate(weighted)) / weight_total
            if weight_total
            else scans[center].time_min
        )
        spectrum = spectrum_from_scan_window(scans, left, right, center)
        peak_class = infer_class(spectrum)
        peaks.append(
            {
                "id": f"R{idx:02d}",
                "rt": round(rt, 4),
                "area": round(area),
                "name": "Unknown raw peak",
                "formula": "",
                "class": peak_class,
                "match": 0,
                "spectrum": spectrum,
                "annotation": {
                    "name": "Unknown raw peak",
                    "class": peak_class,
                    "confidence": "tentative class",
                    "notes": f"Imported directly from Agilent DATA.MS{f' using {run_meta.get('methodName')}' if run_meta.get('methodName') else ''}. Spectrum is averaged around the peak apex; peak picking remains approximate and should be reviewed manually.",
                },
            }
        )
    return peaks


def simplify_trace(scans: list[Scan], max_points: int = 2200) -> list[dict[str, float]]:
    if not scans:
        return []
    if len(scans) <= max_points:
        return [{"rt": scan.time_min, "value": scan.tic} for scan in scans]
    step = math.ceil(len(scans) / max_points)
    trace: list[dict[str, float]] = []
    for i in range(0, len(scans), step):
        chunk = scans[i : i + step]
        avg_rt = sum(scan.time_min for scan in chunk) / len(chunk)
        max_val = max(scan.tic for scan in chunk)
        trace.append({"rt": avg_rt, "value": max_val})
    return trace


def collect_from_directory(path: Path) -> tuple[str, bytes, dict[str, str]]:
    if path.is_file() and path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path) as archive:
            names = archive.namelist()
            data_name = next((name for name in names if name.lower().endswith("data.ms") and not name.startswith("__MACOSX/")), None)
            if not data_name:
                raise FileNotFoundError("No DATA.MS was found inside the ZIP file.")
            folder = "/".join(data_name.split("/")[:-1])
            texts: dict[str, str] = {}
            for name in names:
                if not name.startswith(folder) or name.startswith("__MACOSX/") or not re.search(r"\.(txt|ini)$", name, re.IGNORECASE):
                    continue
                try:
                    texts[name] = archive.read(name).decode("utf-8")
                except UnicodeDecodeError:
                    texts[name] = archive.read(name).decode("latin1", "replace")
            return Path(folder).name or path.stem, archive.read(data_name), texts

    if path.is_dir() and path.suffix.lower() == ".d":
        data_path = path / "DATA.MS"
        if not data_path.exists():
            raise FileNotFoundError("No DATA.MS was found in that .D folder.")
        texts: dict[str, str] = {}
        for sidecar in path.iterdir():
            if sidecar.suffix.lower() in {".txt", ".ini"}:
                try:
                    texts[str(sidecar)] = sidecar.read_text(encoding="utf-8")
                except UnicodeDecodeError:
                    texts[str(sidecar)] = sidecar.read_text(encoding="latin1", errors="replace")
        return path.name, data_path.read_bytes(), texts

    raise FileNotFoundError("Input must be an Agilent .D folder or a zipped .D folder.")


def build_bundle(input_path: Path) -> dict[str, object]:
    run_name, ms_data, text_map = collect_from_directory(input_path)
    scans, ms_meta = parse_agilent_ms_buffer(ms_data)
    run_meta = parse_associated_metadata(text_map, run_name, ms_meta)
    peaks = detect_raw_peaks(scans, run_meta)
    return {
        "version": 1,
        "bundleType": "gcms-analysis-bundle",
        "dataSourceLabel": f"Processed raw import from {run_name}",
        "dataSourceKind": "raw",
        "selectedId": peaks[0]["id"] if peaks else None,
        "peaks": peaks,
        "compareRuns": [],
        "currentTrace": simplify_trace(scans),
        "currentRunMeta": run_meta,
    }


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print("Usage: process_agilent_run.py input.D|input.zip output.json", file=sys.stderr)
        return 2
    input_path = Path(argv[1])
    output_path = Path(argv[2])
    bundle = build_bundle(input_path)
    output_path.write_text(json.dumps(bundle, indent=2), encoding="utf-8")
    print(f"Wrote workstation bundle with {len(bundle['peaks'])} peaks to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
