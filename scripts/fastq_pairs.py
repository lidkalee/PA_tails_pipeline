"""Identify exact R1/R2 FASTQ pairs in one or two directories."""

import re
from pathlib import Path

SAFE_SAMPLE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]*$")

DEFAULT_R1 = r"^(?P<sample>.+)[_.-]R1(?P<suffix>[_-]\d+)?\.(?:fastq|fq)(?:\.gz)?$"
DEFAULT_R2 = r"^(?P<sample>.+)[_.-]R2(?P<suffix>[_-]\d+)?\.(?:fastq|fq)(?:\.gz)?$"
FASTQ = re.compile(r"\.(?:fastq|fq)(?:\.gz)?$", re.IGNORECASE)


def _scan(folder, expression, mate, ignore_expression=None):
    folder = Path(folder)
    if not folder.is_dir():
        raise ValueError(f"{mate} directory does not exist: {folder}")
    pattern = re.compile(expression, re.IGNORECASE)
    ignore_pattern = re.compile(ignore_expression, re.IGNORECASE) if ignore_expression else None
    if "sample" not in pattern.groupindex:
        raise ValueError(f"{mate} regex must have a named (?P<sample>...) group")
    found = {}
    unmatched = []
    for path in sorted(folder.iterdir()):
        if not path.is_file() or not FASTQ.search(path.name):
            continue
        match = pattern.fullmatch(path.name)
        if not match:
            if ignore_pattern and ignore_pattern.fullmatch(path.name):
                continue
            unmatched.append(path.name)
            continue
        sample = match.group("sample")
        if not SAFE_SAMPLE.fullmatch(sample):
            raise ValueError(f"Unsafe sample ID {sample!r} from {path.name}")
        if sample in found:
            raise ValueError(f"Multiple {mate} FASTQs for {sample!r}: {found[sample]} and {path}")
        found[sample] = (path.resolve(), match.groupdict().get("suffix"))
    if unmatched:
        raise ValueError(
            f"{mate} FASTQs did not match the regex: {', '.join(unmatched[:5])}. "
            "Use --r1-regex/--r2-regex or write the sample sheet manually."
        )
    if not found:
        raise ValueError(f"No {mate} FASTQs found in {folder}")
    return found


def discover_pairs(r1_dir, r2_dir, r1_regex=DEFAULT_R1, r2_regex=DEFAULT_R2):
    same_folder = Path(r1_dir).resolve() == Path(r2_dir).resolve()
    r1 = _scan(r1_dir, r1_regex, "R1", r2_regex if same_folder else None)
    r2 = _scan(r2_dir, r2_regex, "R2", r1_regex if same_folder else None)
    if r1.keys() != r2.keys():
        only_r1 = ", ".join(sorted(r1.keys() - r2.keys())) or "none"
        only_r2 = ", ".join(sorted(r2.keys() - r1.keys())) or "none"
        raise ValueError(f"Unpaired samples (R1 only: {only_r1}; R2 only: {only_r2})")
    for sample in r1:
        if r1[sample][0] == r2[sample][0]:
            raise ValueError(f"{sample!r}: R1 and R2 resolved to the same file")
        if r1[sample][1] != r2[sample][1]:
            raise ValueError(f"{sample!r}: R1 and R2 numeric read suffixes differ")
    return {sample: (r1[sample][0], r2[sample][0]) for sample in sorted(r1)}


