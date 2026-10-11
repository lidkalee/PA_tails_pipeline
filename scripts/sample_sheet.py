"""Load one or more paired FASTQ samples from Snakemake configuration."""

import csv
import re
from pathlib import Path

SAFE_SAMPLE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]*$")


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
    sample = str(config["sample"])
    return {sample: _validate(sample, str(config["R1"]), str(config["R2"]), "config")}
