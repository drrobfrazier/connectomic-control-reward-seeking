#!/usr/bin/env python3
"""
Master pipeline for:
Connectomic Control of Maladaptive Reward Seeking:
Whole-CNS Identification of Circuit Vulnerabilities in Drosophila

Runs the remaining computational study in one pass:
1) degree/strength-matched DAN null model (fixed reward network; matched 4-DAN punishment null)
2) whole-CNS pathway-to-descending-neuron structural reach
3) candidate interface ranking
4) FlyBrain in-silico whole-type lesion screen across parameter settings
5) robustness aggregation + paper-ready CSVs + summary/manifest

Expected project layout:
<repository root>/
FlyBrain checkout supplied via FLYBRAIN_DIR or ~/flybrain/
"""
from __future__ import annotations

import csv
import math
import os
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.ipc as ipc

# ------------------------------- paths/config
PROJECT = Path(__file__).resolve().parents[1]
FLYBRAIN = Path(os.environ.get("FLYBRAIN_DIR", str(Path.home() / "flybrain"))).expanduser().resolve()
RESULTS = PROJECT / "results"
LOGS = PROJECT / "logs"
TABLES = PROJECT / "tables"
SEEDS_FILE = PROJECT / "data" / "seeds" / "seed_registry_v1.csv"
WEIGHTS_FILE = FLYBRAIN / "connectome-weights-male-cns-v1.0-minconf-0.5.feather"
NEURONS_FILE = FLYBRAIN / "neurons.csv"
GRAPH_FILE = FLYBRAIN / "malecns_graph.npz"

for p in (RESULTS, LOGS, TABLES):
    p.mkdir(parents=True, exist_ok=True)

RNG_SEED = 20261008
N_MATCHED_PERM = 10_000
MATCH_POOL = 30
MATCH_DRAW_TOP = 12
TOP_STRUCTURAL = 100
TOP_LESION_TYPES = 12

# Simulation settings. These are intentionally modest enough for a local machine
# while still providing a parameter-robust screen.
SIM_HZ = (30.0, 60.0)
SIM_WEIGHT_SCALE = (0.35, 0.50)
SIM_MS = 300.0
SIM_STIM_MS = 200.0
SIM_TRIALS = 3
SIM_DT = 0.1

# ------------------------------- utilities

def log(msg=""):
    print(msg, flush=True)


def empirical_upper(null, obs):
    null = np.asarray(null, float)
    return (1 + np.sum(null >= obs)) / (len(null) + 1)


def empirical_lower(null, obs):
    null = np.asarray(null, float)
    return (1 + np.sum(null <= obs)) / (len(null) + 1)


def percentile_rank(x):
    x = np.asarray(x, float)
    n = len(x)
    if n <= 1:
        return np.ones(n)
    order = np.argsort(x, kind="stable")
    ranks = np.empty(n, float)
    ranks[order] = np.arange(n, dtype=float) / (n - 1)
    return ranks


def write_csv(path: Path, rows, fieldnames=None):
    rows = list(rows)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        if fieldnames:
            w.writeheader()
            w.writerows(rows)


def exact_type_indices(C, typ):
    return np.array([i for i, t in enumerate(C.type) if t == typ], dtype=np.int32)


def safe_float(x, default=0.0):
    try:
        return float(x)
    except Exception:
        return default


# ------------------------------- validate inputs
required = [SEEDS_FILE, WEIGHTS_FILE, NEURONS_FILE, GRAPH_FILE]
missing = [str(p) for p in required if not p.exists()]
if missing:
    raise SystemExit("Missing required files:\n" + "\n".join(missing))

log("=" * 78)
log("MASTER STUDY PIPELINE")
log("Connectomic Control of Maladaptive Reward Seeking")
log("=" * 78)
log(f"Project: {PROJECT}")
log(f"FlyBrain: {FLYBRAIN}")
log()

# ------------------------------- load annotations
log("[1/7] Loading annotations and frozen seed registry...")
ann = {}
ann_rows = []
with NEURONS_FILE.open(newline="", encoding="utf-8") as f:
    for r in csv.DictReader(f):
        b = int(r["body"])
        rr = dict(r)
        rr["body"] = b
        ann[b] = rr
        ann_rows.append(rr)

seed_rows = []
with SEEDS_FILE.open(newline="", encoding="utf-8") as f:
    seed_rows = list(csv.DictReader(f))

reward_bodies = np.array([int(r["body_id"]) for r in seed_rows if r["study_role"] == "reward"], dtype=np.int64)
punish_bodies = np.array([int(r["body_id"]) for r in seed_rows if r["study_role"] == "punishment"], dtype=np.int64)
all_seed_bodies = set(reward_bodies.tolist()) | set(punish_bodies.tolist())

if len(reward_bodies) != 115 or len(punish_bodies) != 4:
    raise SystemExit(f"Seed registry mismatch: reward={len(reward_bodies)}, punishment={len(punish_bodies)}")

# Broad DAN pool from the same FlyBrain annotation table.
dan_bodies = sorted({
    int(r["body"]) for r in ann_rows
    if str(r.get("class", "")).strip().lower() == "dan"
    and str(r.get("nt", "")).strip().lower() == "dopamine"
})
dan_bodies = np.array(dan_bodies, dtype=np.int64)
log(f"  Reward seeds: {len(reward_bodies)}")
log(f"  Punishment seeds: {len(punish_bodies)}")
log(f"  Annotated dopamine DAN pool: {len(dan_bodies)}")

# ------------------------------- scan raw graph for all DAN outgoing edges
log("\n[2/7] Extracting outgoing edges for all annotated DANs from raw MaleCNS...")
source = pa.memory_map(str(WEIGHTS_FILE), "r")
reader = ipc.open_file(source)
schema = reader.schema
pre_col_i = schema.get_field_index("body_pre")
post_col_i = schema.get_field_index("body_post")
weight_col_i = schema.get_field_index("weight")
if min(pre_col_i, post_col_i, weight_col_i) < 0:
    raise SystemExit(f"Required raw-connectome columns missing. Schema:\n{schema}")

dan_arrow = pa.array(dan_bodies.tolist(), type=schema.field(pre_col_i).type)
pre_chunks, post_chunks, wt_chunks = [], [], []
for bi in range(reader.num_record_batches):
    batch = reader.get_batch(bi)
    m = pc.is_in(batch.column(pre_col_i), value_set=dan_arrow)
    if not pc.any(m).as_py():
        continue
    pre_chunks.extend(pc.filter(batch.column(pre_col_i), m).to_pylist())
    post_chunks.extend(pc.filter(batch.column(post_col_i), m).to_pylist())
    wt_chunks.extend(pc.filter(batch.column(weight_col_i), m).to_pylist())

edge_pre = np.asarray(pre_chunks, dtype=np.int64)
edge_post = np.asarray(post_chunks, dtype=np.int64)
edge_wt = np.asarray(wt_chunks, dtype=np.float64)
del pre_chunks, post_chunks, wt_chunks
log(f"  DAN-originating raw edges: {len(edge_pre):,}")

# Compress DAN and target IDs.
dan_sorted = np.sort(dan_bodies)
pre_di = np.searchsorted(dan_sorted, edge_pre)
if not np.all(dan_sorted[pre_di] == edge_pre):
    raise SystemExit("DAN indexing failure")

target_bodies = np.unique(edge_post)
post_ti = np.searchsorted(target_bodies, edge_post)

# Per-DAN sparse target arrays, aggregated just in case multiple rows exist.
order = np.argsort(pre_di, kind="stable")
pre_di_s = pre_di[order]
post_ti_s = post_ti[order]
wt_s = edge_wt[order]
bounds = np.searchsorted(pre_di_s, np.arange(len(dan_sorted) + 1))

dan_targets = []
dan_weights = []
dan_degree = np.zeros(len(dan_sorted), dtype=np.int32)
dan_strength = np.zeros(len(dan_sorted), dtype=np.float64)
for di in range(len(dan_sorted)):
    a, b = bounds[di], bounds[di + 1]
    t = post_ti_s[a:b]
    w = wt_s[a:b]
    if len(t):
        u, inv = np.unique(t, return_inverse=True)
        sw = np.bincount(inv, weights=w, minlength=len(u)).astype(np.float64)
    else:
        u = np.zeros(0, dtype=np.int64)
        sw = np.zeros(0, dtype=np.float64)
    dan_targets.append(u)
    dan_weights.append(sw)
    dan_degree[di] = len(u)
    dan_strength[di] = sw.sum()

body_to_di = {int(b): i for i, b in enumerate(dan_sorted)}
reward_di = np.array([body_to_di[int(b)] for b in reward_bodies], dtype=np.int32)
punish_di = np.array([body_to_di[int(b)] for b in punish_bodies], dtype=np.int32)

# Fixed reward arrays on the DAN-target universe.
def aggregate_group(dis):
    ts, ws = [], []
    for di in dis:
        ts.append(dan_targets[int(di)])
        ws.append(dan_weights[int(di)])
    if not ts:
        return np.zeros(0, np.int64), np.zeros(0, float), np.zeros(0, np.int32)
    cat_t = np.concatenate(ts) if ts else np.zeros(0, np.int64)
    cat_w = np.concatenate(ws) if ws else np.zeros(0, float)
    if not len(cat_t):
        return np.zeros(0, np.int64), np.zeros(0, float), np.zeros(0, np.int32)
    u, inv, counts = np.unique(cat_t, return_inverse=True, return_counts=True)
    sw = np.bincount(inv, weights=cat_w, minlength=len(u)).astype(float)
    return u, sw, counts.astype(np.int32)

r_t, r_w_sparse, r_src_sparse = aggregate_group(reward_di)
reward_w = np.zeros(len(target_bodies), dtype=float)
reward_src = np.zeros(len(target_bodies), dtype=np.int32)
reward_w[r_t] = r_w_sparse
reward_src[r_t] = r_src_sparse
reward_present = reward_src > 0


def punishment_metrics(dis):
    p_t, p_w_sparse, p_src_sparse = aggregate_group(dis)
    shared_mask_sparse = reward_present[p_t]
    st = p_t[shared_mask_sparse]
    pw = p_w_sparse[shared_mask_sparse]
    ps = p_src_sparse[shared_mask_sparse]
    if not len(st):
        return dict(shared_targets=0, all4_shared=0, mean_balance=np.nan, shared_total_synapses=0.0)
    rw = reward_w[st]
    mean_r = rw / len(reward_di)
    mean_p = pw / len(dis)
    den = mean_r + mean_p
    dom = np.divide(mean_r - mean_p, den, out=np.zeros_like(den), where=den > 0)
    bal = 1.0 - np.abs(dom)
    return dict(
        shared_targets=int(len(st)),
        all4_shared=int(np.sum(ps == len(dis))),
        mean_balance=float(np.mean(bal)),
        shared_total_synapses=float(np.sum(rw + pw)),
    )

obs_matched = punishment_metrics(punish_di)
log("  Observed (recomputed from raw DAN subset):")
for k, v in obs_matched.items():
    log(f"    {k}: {v}")

# ------------------------------- matched 4-DAN null model
log("\n[3/7] Running degree/strength-matched DAN null model (10,000 draws)...")
rng = np.random.default_rng(RNG_SEED)
nonseed_di = np.array([i for i, b in enumerate(dan_sorted) if int(b) not in all_seed_bodies], dtype=np.int32)
if len(nonseed_di) < 4:
    raise SystemExit("Not enough non-seed DANs for matched null")

# Match on log1p outgoing unique-target degree and total outgoing synaptic strength.
X = np.column_stack([np.log1p(dan_degree), np.log1p(dan_strength)])
mu = X[nonseed_di].mean(axis=0)
sd = X[nonseed_di].std(axis=0, ddof=1)
sd[sd == 0] = 1.0
Z = (X - mu) / sd
match_orders = []
for odi in punish_di:
    dist = np.sqrt(np.sum((Z[nonseed_di] - Z[int(odi)]) ** 2, axis=1))
    ord_local = nonseed_di[np.argsort(dist, kind="stable")]
    match_orders.append(ord_local[: min(MATCH_POOL, len(ord_local))])

null_shared = np.zeros(N_MATCHED_PERM, float)
null_all4 = np.zeros(N_MATCHED_PERM, float)
null_balance = np.zeros(N_MATCHED_PERM, float)
null_syn = np.zeros(N_MATCHED_PERM, float)
selected_quartets = []

for b in range(N_MATCHED_PERM):
    used = set()
    chosen = []
    # Randomize slot order to avoid systematic collisions in overlapping match pools.
    for slot in rng.permutation(len(match_orders)):
        pool = [int(x) for x in match_orders[int(slot)] if int(x) not in used]
        if not pool:
            pool = [int(x) for x in nonseed_di if int(x) not in used]
        top = pool[: min(MATCH_DRAW_TOP, len(pool))]
        pick = int(rng.choice(top))
        used.add(pick)
        chosen.append((int(slot), pick))
    chosen = [p for _, p in sorted(chosen)]
    m = punishment_metrics(np.array(chosen, dtype=np.int32))
    null_shared[b] = m["shared_targets"]
    null_all4[b] = m["all4_shared"]
    null_balance[b] = m["mean_balance"]
    null_syn[b] = m["shared_total_synapses"]
    if b < 100:
        selected_quartets.append([int(dan_sorted[x]) for x in chosen])
    if (b + 1) % 1000 == 0:
        log(f"  matched null: {b + 1:,}/{N_MATCHED_PERM:,}")

matched_summary = {
    "shared_targets_observed": obs_matched["shared_targets"],
    "shared_targets_null_mean": float(np.nanmean(null_shared)),
    "shared_targets_p_upper": empirical_upper(null_shared, obs_matched["shared_targets"]),
    "all4_shared_observed": obs_matched["all4_shared"],
    "all4_shared_null_mean": float(np.nanmean(null_all4)),
    "all4_shared_p_upper": empirical_upper(null_all4, obs_matched["all4_shared"]),
    "mean_balance_observed": obs_matched["mean_balance"],
    "mean_balance_null_mean": float(np.nanmean(null_balance)),
    "mean_balance_p_lower": empirical_lower(null_balance, obs_matched["mean_balance"]),
    "shared_total_synapses_observed": obs_matched["shared_total_synapses"],
    "shared_total_synapses_null_mean": float(np.nanmean(null_syn)),
    "shared_total_synapses_p_upper": empirical_upper(null_syn, obs_matched["shared_total_synapses"]),
}

write_csv(
    RESULTS / "matched_dan_null_10000.csv",
    ({
        "permutation": i + 1,
        "shared_targets": null_shared[i],
        "all4_shared": null_all4[i],
        "mean_balance": null_balance[i],
        "shared_total_synapses": null_syn[i],
    } for i in range(N_MATCHED_PERM)),
    ["permutation", "shared_targets", "all4_shared", "mean_balance", "shared_total_synapses"],
)

with (RESULTS / "matched_dan_null_summary.txt").open("w", encoding="utf-8") as f:
    f.write("DEGREE/STRENGTH-MATCHED DAN NULL MODEL\n")
    f.write("======================================\n")
    f.write("Reward population held fixed at the 115 preregistered PAM reward DANs.\n")
    f.write("Each observed punishment DAN was replaced by a non-seed dopamine DAN matched\n")
    f.write("on log outgoing target degree and log outgoing synaptic strength.\n\n")
    for k, v in matched_summary.items():
        f.write(f"{k}: {v}\n")

log("  Matched-null summary:")
for k, v in matched_summary.items():
    log(f"    {k}: {v}")

# ------------------------------- load FlyBrain engine and compute pathway reach
log("\n[4/7] Computing whole-CNS structural reach to descending neurons...")
sys.path.insert(0, str(FLYBRAIN))
from engine import Connectome, resolve_params, simulate  # noqa: E402

C = Connectome(str(FLYBRAIN))
N = C.N
pre = C.pre
post = C.indices
w = C.nsyn.astype(np.float64)
row_sum = np.bincount(pre, weights=w, minlength=N)

dn_idx = np.flatnonzero(np.array([str(s).startswith("descending") for s in C.superclass], dtype=bool)).astype(np.int32)
dn0 = np.zeros(N, dtype=np.float64)
dn0[dn_idx] = 1.0

def backward_step(prev):
    num = np.bincount(pre, weights=w * prev[post], minlength=N)
    return np.divide(num, row_sum, out=np.zeros(N, dtype=float), where=row_sum > 0)

reach1 = backward_step(dn0)
reach2 = backward_step(reach1)
reach3 = backward_step(reach2)
dn_reach = 0.50 * reach1 + 0.30 * reach2 + 0.20 * reach3

# ------------------------------- build candidate ranking from observed direct convergence
log("\n[5/7] Ranking structural reward/punishment interface candidates...")
p_t, p_w_sparse, p_src_sparse = aggregate_group(punish_di)
punish_w = np.zeros(len(target_bodies), dtype=float)
punish_src = np.zeros(len(target_bodies), dtype=np.int32)
punish_w[p_t] = p_w_sparse
punish_src[p_t] = p_src_sparse
shared_ti = np.flatnonzero((reward_src > 0) & (punish_src > 0))

candidate_rows = []
for ti in shared_ti:
    body = int(target_bodies[ti])
    a = ann.get(body, {})
    rw = float(reward_w[ti]); pw = float(punish_w[ti])
    rs = int(reward_src[ti]); ps = int(punish_src[ti])
    mr = rw / len(reward_di); mp = pw / len(punish_di)
    den = mr + mp
    dom = (mr - mp) / den if den > 0 else 0.0
    rfrac = rs / len(reward_di); pfrac = ps / len(punish_di)
    cross_source = math.sqrt(max(rfrac * pfrac, 0.0))
    cross_weight = math.sqrt(max(mr * mp, 0.0))
    ci = C.body_to_idx.get(body)
    dr = float(dn_reach[ci]) if ci is not None else 0.0
    candidate_rows.append({
        "body_post": body,
        "type": a.get("type", ""),
        "instance": a.get("instance", ""),
        "superclass": a.get("superclass", ""),
        "class": a.get("class", ""),
        "nt": a.get("nt", ""),
        "reward_synapses": int(round(rw)),
        "punishment_synapses": int(round(pw)),
        "reward_sources": rs,
        "punishment_sources": ps,
        "reward_source_fraction": rfrac,
        "punishment_source_fraction": pfrac,
        "mean_reward_synapses_per_seed": mr,
        "mean_punishment_synapses_per_seed": mp,
        "dominance": dom,
        "balance": 1.0 - abs(dom),
        "cross_source": cross_source,
        "cross_weight": cross_weight,
        "dn_reach_1hop": float(reach1[ci]) if ci is not None else 0.0,
        "dn_reach_2hop": float(reach2[ci]) if ci is not None else 0.0,
        "dn_reach_3hop": float(reach3[ci]) if ci is not None else 0.0,
        "dn_reach_weighted": dr,
    })

# Composite rank uses only structural quantities, not simulation outcomes.
pr_cross_source = percentile_rank([r["cross_source"] for r in candidate_rows])
pr_cross_weight = percentile_rank([math.log1p(r["cross_weight"]) for r in candidate_rows])
pr_dn = percentile_rank([r["dn_reach_weighted"] for r in candidate_rows])
for i, r in enumerate(candidate_rows):
    r["interface_score"] = 0.35 * pr_cross_source[i] + 0.35 * pr_cross_weight[i] + 0.30 * pr_dn[i]
    r["polarity"] = "reward-dominant" if r["dominance"] > 0.10 else ("punishment-dominant" if r["dominance"] < -0.10 else "balanced")

candidate_rows.sort(key=lambda r: r["interface_score"], reverse=True)
write_csv(RESULTS / "structural_candidate_ranking.csv", candidate_rows)
write_csv(TABLES / "top100_structural_candidates.csv", candidate_rows[:TOP_STRUCTURAL])

# Aggregate by cell type.
type_agg = defaultdict(list)
for r in candidate_rows:
    if r["type"]:
        type_agg[r["type"]].append(r)

type_rows = []
for typ, rows in type_agg.items():
    scores = np.array([r["interface_score"] for r in rows], float)
    doms = np.array([r["dominance"] for r in rows], float)
    inds = exact_type_indices(C, typ)
    type_rows.append({
        "type": typ,
        "shared_target_neurons": len(rows),
        "whole_type_size": len(inds),
        "max_interface_score": float(scores.max()),
        "mean_interface_score": float(scores.mean()),
        "median_dominance": float(np.median(doms)),
        "max_dn_reach": float(max(r["dn_reach_weighted"] for r in rows)),
    })
type_rows.sort(key=lambda r: (r["max_interface_score"], r["mean_interface_score"]), reverse=True)
write_csv(RESULTS / "structural_candidate_types.csv", type_rows)

log("  Top structural candidates:")
for r in candidate_rows[:15]:
    log(f"    body={r['body_post']:<8} type={r['type']:<18} score={r['interface_score']:.3f} dom={r['dominance']:+.3f} DNreach={r['dn_reach_weighted']:.4f}")

# ------------------------------- choose candidate types for virtual-lesion screen
primary_types = []
for tr in type_rows:
    typ = tr["type"]
    if not typ or typ.startswith("KC") or typ in {"APL", "DPM"}:
        continue
    if 1 <= tr["whole_type_size"] <= 50:
        primary_types.append(typ)
    if len(primary_types) >= TOP_LESION_TYPES:
        break

# Force biologically salient types already observed, if they exist, then add broad-hub controls.
for typ in ("MBON05", "MBON11", "SMP108"):
    if exact_type_indices(C, typ).size and typ not in primary_types:
        primary_types.append(typ)
controls = [typ for typ in ("APL", "DPM") if exact_type_indices(C, typ).size]
lesion_types = primary_types + [t for t in controls if t not in primary_types]

write_csv(
    RESULTS / "virtual_lesion_candidate_types.csv",
    ({"type": t, "n_neurons": int(exact_type_indices(C, t).size), "role": "control" if t in controls else "candidate"} for t in lesion_types),
    ["type", "n_neurons", "role"],
)
log(f"\n[6/7] Virtual lesion screen: {len(lesion_types)} cell types")
log("  " + ", ".join(lesion_types))

# Readout definition grounded in canonical MBON-valence literature:
# glutamatergic MBONs -> avoidance-associated; cholinergic/GABAergic MBONs -> approach-associated.
mbon_idx = np.array([i for i, t in enumerate(C.type) if str(t).startswith("MBON")], dtype=np.int32)
approach_idx = np.array([i for i in mbon_idx if str(C.nt[i]).lower() in {"acetylcholine", "gaba"}], dtype=np.int32)
avoid_idx = np.array([i for i in mbon_idx if str(C.nt[i]).lower() == "glutamate"], dtype=np.int32)
if len(approach_idx) == 0 or len(avoid_idx) == 0:
    raise SystemExit(f"MBON readout definition failed: approach={len(approach_idx)}, avoidance={len(avoid_idx)}")

reward_idx = np.array([C.body_to_idx[int(b)] for b in reward_bodies if int(b) in C.body_to_idx], dtype=np.int32)
punish_idx = np.array([C.body_to_idx[int(b)] for b in punish_bodies if int(b) in C.body_to_idx], dtype=np.int32)
if len(reward_idx) != 115 or len(punish_idx) != 4:
    raise SystemExit(f"Simulation seed mapping failed: reward={len(reward_idx)}, punishment={len(punish_idx)}")

approach_mask = np.zeros(N, bool); approach_mask[approach_idx] = True
avoid_mask = np.zeros(N, bool); avoid_mask[avoid_idx] = True


def run_sim(stims, q, silence=None):
    runs = simulate(C, stims, q, silence=silence)
    a_spikes = 0
    v_spikes = 0
    total = 0
    active_union = np.zeros(N, bool)
    for _, idx in runs:
        total += len(idx)
        if len(idx):
            a_spikes += int(approach_mask[idx].sum())
            v_spikes += int(avoid_mask[idx].sum())
            active_union[idx] = True
    sec = q["trials"] * (q["t_ms"] / 1000.0)
    arate = a_spikes / max(len(approach_idx), 1) / sec
    vrate = v_spikes / max(len(avoid_idx), 1) / sec
    rdi = (arate - vrate) / (arate + vrate + 1e-12)
    return {
        "approach_rate_hz": arate,
        "avoidance_rate_hz": vrate,
        "rdi": rdi,
        "total_spikes": total,
        "active_neurons": int(active_union.sum()),
    }

sim_rows = []
for wscale in SIM_WEIGHT_SCALE:
    for hz in SIM_HZ:
        q = resolve_params({
            "model": "shiu",
            "t_ms": SIM_MS,
            "trials": SIM_TRIALS,
            "seed": RNG_SEED,
            "w_scale": wscale,
            "dt": SIM_DT,
        })
        conditions = {
            "reward": [dict(idx=reward_idx, hz=hz, t_on=0.0, t_off=SIM_STIM_MS)],
            "punishment": [dict(idx=punish_idx, hz=hz, t_on=0.0, t_off=SIM_STIM_MS)],
            "conflict": [
                dict(idx=reward_idx, hz=hz, t_on=0.0, t_off=SIM_STIM_MS),
                dict(idx=punish_idx, hz=hz, t_on=0.0, t_off=SIM_STIM_MS),
            ],
        }
        log(f"  intact simulations: weight_scale={wscale}, hz={hz}")
        intact = {name: run_sim(stims, q, silence=None) for name, stims in conditions.items()}
        for name, out in intact.items():
            sim_rows.append({
                "weight_scale": wscale, "hz": hz, "condition": name,
                "lesion_type": "INTACT", "lesion_n": 0,
                **out,
                "delta_rdi_vs_intact": 0.0,
            })

        for li, typ in enumerate(lesion_types, 1):
            sil = exact_type_indices(C, typ)
            log(f"    lesion {li:02d}/{len(lesion_types):02d}: {typ} (n={len(sil)})")
            # Reward-only tests preservation; conflict tests selective control.
            for name in ("reward", "conflict"):
                out = run_sim(conditions[name], q, silence=sil)
                sim_rows.append({
                    "weight_scale": wscale, "hz": hz, "condition": name,
                    "lesion_type": typ, "lesion_n": len(sil),
                    **out,
                    "delta_rdi_vs_intact": out["rdi"] - intact[name]["rdi"],
                })

write_csv(RESULTS / "virtual_lesion_results.csv", sim_rows)

# Aggregate simulation robustness by lesion type.
lesion_summary = []
for typ in lesion_types:
    cr = [r for r in sim_rows if r["lesion_type"] == typ and r["condition"] == "conflict"]
    rr = [r for r in sim_rows if r["lesion_type"] == typ and r["condition"] == "reward"]
    if not cr or not rr:
        continue
    cdelta = np.array([safe_float(r["delta_rdi_vs_intact"]) for r in cr])
    rdelta = np.array([safe_float(r["delta_rdi_vs_intact"]) for r in rr])
    med_c = float(np.median(cdelta))
    med_r = float(np.median(rdelta))
    sign = 1 if med_c > 0 else (-1 if med_c < 0 else 0)
    consistency = float(np.mean(np.sign(cdelta) == sign)) if sign else 0.0
    specificity = float(np.median(np.abs(cdelta)) - np.median(np.abs(rdelta)))
    ratio = float(np.median(np.abs(cdelta)) / (np.median(np.abs(rdelta)) + 0.02))
    robust_score = float(np.median(np.abs(cdelta)) * consistency / (np.median(np.abs(rdelta)) + 0.02))
    lesion_summary.append({
        "type": typ,
        "n_neurons": int(exact_type_indices(C, typ).size),
        "median_conflict_delta_rdi": med_c,
        "median_reward_delta_rdi": med_r,
        "median_abs_conflict_delta": float(np.median(np.abs(cdelta))),
        "median_abs_reward_delta": float(np.median(np.abs(rdelta))),
        "conflict_effect_sign_consistency": consistency,
        "conflict_specificity_absdiff": specificity,
        "conflict_to_reward_effect_ratio": ratio,
        "robust_lesion_score": robust_score,
    })
lesion_summary.sort(key=lambda r: r["robust_lesion_score"], reverse=True)
write_csv(RESULTS / "virtual_lesion_summary.csv", lesion_summary)

# ------------------------------- final manifest + summary
log("\n[7/7] Writing study summary and reproducibility manifest...")
try:
    git_hash = subprocess.check_output(["git", "-C", str(FLYBRAIN), "rev-parse", "HEAD"], text=True).strip()
except Exception:
    git_hash = "unavailable"

manifest = {
    "pipeline_random_seed": RNG_SEED,
    "flybrain_git_commit": git_hash,
    "python": sys.version.replace("\n", " "),
    "numpy": np.__version__,
    "pyarrow": pa.__version__,
    "matched_null_permutations": N_MATCHED_PERM,
    "simulation_model": "shiu",
    "simulation_hz": ",".join(map(str, SIM_HZ)),
    "simulation_weight_scale": ",".join(map(str, SIM_WEIGHT_SCALE)),
    "simulation_ms": SIM_MS,
    "simulation_stim_ms": SIM_STIM_MS,
    "simulation_trials": SIM_TRIALS,
    "simulation_dt_ms": SIM_DT,
    "approach_readout_neurons": len(approach_idx),
    "avoidance_readout_neurons": len(avoid_idx),
}
with (RESULTS / "reproducibility_manifest.txt").open("w", encoding="utf-8") as f:
    for k, v in manifest.items():
        f.write(f"{k}: {v}\n")

summary_path = RESULTS / "STUDY_COMPLETE_SUMMARY.txt"
with summary_path.open("w", encoding="utf-8") as f:
    f.write("CONNECTOMIC CONTROL OF MALADAPTIVE REWARD SEEKING\n")
    f.write("Whole-CNS Identification of Circuit Vulnerabilities in Drosophila\n")
    f.write("=" * 74 + "\n\n")
    f.write("MATCHED DAN NULL MODEL\n")
    f.write("----------------------\n")
    for k, v in matched_summary.items():
        f.write(f"{k}: {v}\n")
    f.write("\nTOP STRUCTURAL INTERFACE CANDIDATES\n")
    f.write("-----------------------------------\n")
    for r in candidate_rows[:20]:
        f.write(
            f"body={r['body_post']} type={r['type']} score={r['interface_score']:.4f} "
            f"dominance={r['dominance']:+.4f} polarity={r['polarity']} "
            f"DNreach={r['dn_reach_weighted']:.6f}\n"
        )
    f.write("\nVIRTUAL LESION ROBUSTNESS RANKING\n")
    f.write("---------------------------------\n")
    for r in lesion_summary[:20]:
        f.write(
            f"type={r['type']} n={r['n_neurons']} robust_score={r['robust_lesion_score']:.4f} "
            f"median_conflict_delta={r['median_conflict_delta_rdi']:+.4f} "
            f"median_reward_delta={r['median_reward_delta_rdi']:+.4f} "
            f"sign_consistency={r['conflict_effect_sign_consistency']:.2f}\n"
        )
    f.write("\nINTERPRETATION NOTE\n")
    f.write("-------------------\n")
    f.write(
        "The MaleCNS wiring is empirical, but FlyBrain neural dynamics are a simplified LIF model. "
        "Simulation results therefore test propagation and perturbation hypotheses, not drug pharmacology, "
        "dopamine receptor biology, or a literal addicted-fly state.\n"
    )

log("\n" + "=" * 78)
log("PIPELINE COMPLETE")
log("=" * 78)
log(f"Main summary: {summary_path}")
log(f"Matched DAN null: {RESULTS / 'matched_dan_null_summary.txt'}")
log(f"Structural ranking: {RESULTS / 'structural_candidate_ranking.csv'}")
log(f"Lesion screen: {RESULTS / 'virtual_lesion_results.csv'}")
log(f"Lesion summary: {RESULTS / 'virtual_lesion_summary.csv'}")
log(f"Manifest: {RESULTS / 'reproducibility_manifest.txt'}")
