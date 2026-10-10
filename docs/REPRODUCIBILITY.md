# Reproducibility

## Recorded environment

The archived reproducibility manifest records:

- Pipeline random seed: `20261008`
- FlyBrain Git commit: `c582cfdb90c49cfa7da71d0f5ec958dc17730223`
- Python: `3.14.4`
- NumPy: `2.5.1`
- PyArrow: `25.0.0`
- Matched-DAN permutations: `10,000`
- Simulation model: `shiu`
- Stimulation rates: `30` and `60 Hz`
- Synaptic-weight scales: `0.35` and `0.50`
- Simulation duration: `300 ms`
- Stimulation duration: `200 ms`
- Trials per setting: `3`
- Simulation time step: `0.1 ms`
- Approach readout neurons: `71`
- Avoidance readout neurons: `26`

## External data and software

The raw MaleCNS connectome and FlyBrain installation are intentionally not redistributed here. The master pipeline expects a FlyBrain checkout containing the MaleCNS v1.0 files used by the study. Set `FLYBRAIN_DIR` to that checkout. If the variable is not set, the script defaults to `~/flybrain`.

Expected files under the FlyBrain directory include:

- `connectome-weights-male-cns-v1.0-minconf-0.5.feather`
- `neurons.csv`
- `malecns_graph.npz`
- FlyBrain's `engine.py`

## Python reproduction

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
export FLYBRAIN_DIR=/path/to/flybrain
python scripts/finish_study_master.py
```

On PowerShell:

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
$env:FLYBRAIN_DIR="C:\path\to\flybrain"
python scripts\finish_study_master.py
```

The master pipeline reproduces the primary degree/strength-matched DAN null, structural reach/ranking, and virtual-lesion screen from the frozen seed registry and raw MaleCNS/FlyBrain inputs.

The repository also archives direct-convergence outputs and the secondary within-seed label-permutation outputs used in the manuscript.

## Verify the archived package

This verification does not require the raw MaleCNS data:

```bash
python scripts/verify_archived_results.py
```

It checks key row counts and matched-null means against the archived study results.

## Figure regeneration

Install the packages listed in `R_REQUIREMENTS.txt`, then run:

```bash
Rscript scripts/make_figures.R
```

Figures are written to `results/figures/`.
