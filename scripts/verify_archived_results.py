#!/usr/bin/env python3
"""Verify key values in the archived study outputs without requiring the raw MaleCNS connectome."""
from pathlib import Path
import csv
import statistics

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"

def read_rows(name):
    with (RESULTS / name).open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))

direct = read_rows("direct_targets_all.csv")
shared = read_rows("direct_convergence_shared.csv")
matched = read_rows("matched_dan_null_10000.csv")
perm = read_rows("permutation_null_10000.csv")
lesions = read_rows("virtual_lesion_summary.csv")

assert len(direct) == 34172, len(direct)
assert len(shared) == 2636, len(shared)
assert len(matched) == 10000, len(matched)
assert len(perm) == 10000, len(perm)
assert len(lesions) == 14, len(lesions)

m_shared = statistics.fmean(float(r["shared_targets"]) for r in matched)
m_all4 = statistics.fmean(float(r["all4_shared"]) for r in matched)
m_balance = statistics.fmean(float(r["mean_balance"]) for r in matched)
m_syn = statistics.fmean(float(r["shared_total_synapses"]) for r in matched)

checks = [
    ("matched shared-target mean", m_shared, 1260.2658, 1e-9),
    ("matched all-four mean", m_all4, 2.7816, 1e-9),
    ("matched balance mean", m_balance, 0.43577963030535577, 1e-12),
    ("matched synaptic-mass mean", m_syn, 40805.5276, 1e-9),
]
for label, got, expected, tol in checks:
    if abs(got - expected) > tol:
        raise AssertionError(f"{label}: got {got}, expected {expected}")

print("Archived-output verification passed.")
print(f"Direct downstream targets: {len(direct):,}")
print(f"Shared direct targets: {len(shared):,}")
print(f"Matched-DAN permutations: {len(matched):,}")
print(f"Within-seed label permutations: {len(perm):,}")
print(f"Virtual-lesion candidate types: {len(lesions)}")
