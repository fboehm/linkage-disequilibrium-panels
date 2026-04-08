#!/usr/bin/env python3
"""Count sample sizes by ancestry group in the HAPNEST public dataset.

Usage:
    python3 scripts/count_hapnest_ancestries.py

Reads resources/hapnest_public/synthetic_v1.sample and reports the number
of individuals per ancestry label, plus the total.
"""

from collections import Counter
from pathlib import Path

sample_file = Path("resources/hapnest_public/synthetic_v1.sample")

if not sample_file.exists():
    raise FileNotFoundError(f"Not found: {sample_file}")

labels = [line.strip() for line in sample_file.read_text().splitlines() if line.strip()]
counts = Counter(labels)

print(f"{'Ancestry':<12} {'N':>8}")
print("-" * 22)
for ancestry, n in sorted(counts.items()):
    print(f"{ancestry:<12} {n:>8,}")
print("-" * 22)
print(f"{'Total':<12} {sum(counts.values()):>8,}")
