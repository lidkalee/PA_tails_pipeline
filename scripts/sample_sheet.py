"""Load one or more paired FASTQ samples from Snakemake configuration."""

import csv
from pathlib import Path

from scripts.fastq_pairs import SAFE_SAMPLE, discover_pairs


def _validate(sample, r1, r2, location):
    if not SAFE_SAMPLE.fullmatch(sample):
        raise ValueError(
            f"{location}: sample name {sample!r} must contain only ASCII letters, "
            "digits, '_' or '-' and start with a letter or digit"
        )
    if not r1 or not r2:
        raise ValueError(f"{location}: both R1 and R2 paths are required")
    if r1 == r2:
        raise ValueError(f"{location}: R1 and R2 must be different files")
    return {"R1": r1, "R2": r2}


def load_samples(config):
    """Return ordered sample-to-FASTQ mapping; sheet takes precedence over legacy keys."""
    if config.get("samples"):
        path = Path(config["samples"])
        with path.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None or not {"sample", "R1", "R2"}.issubset(reader.fieldnames):
                raise ValueError(f"{path}: expected tab-separated columns sample, R1, R2")
            samples = {}
            for number, row in enumerate(reader, start=2):
                location = f"{path}:{number}"
                sample = (row.get("sample") or "").strip()
                r1 = (row.get("R1") or "").strip()
                r2 = (row.get("R2") or "").strip()
                if sample in samples:
                    raise ValueError(f"{location}: duplicate sample {sample!r}")
                samples[sample] = _validate(sample, r1, r2, location)
        if not samples:
            raise ValueError(f"{path}: sample sheet is empty")
        return samples
    if config.get("fastq_dir"):
        if config.get("r1_dir") or config.get("r2_dir"):
            raise ValueError("Use fastq_dir or both r1_dir and r2_dir, not both modes")
        pairs = discover_pairs(config["fastq_dir"], config["fastq_dir"])
        return {sample: _validate(sample, str(r1), str(r2), "fastq_dir")
                for sample, (r1, r2) in pairs.items()}
    if config.get("r1_dir") or config.get("r2_dir"):
        if not config.get("r1_dir") or not config.get("r2_dir"):
            raise ValueError("Both r1_dir and r2_dir are required")
        pairs = discover_pairs(config["r1_dir"], config["r2_dir"])
        return {sample: _validate(sample, str(r1), str(r2), "read directories")
                for sample, (r1, r2) in pairs.items()}
    sample = str(config["sample"])
    return {sample: _validate(sample, str(config["R1"]), str(config["R2"]), "config")}
