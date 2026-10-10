# Connectomic Control of Maladaptive Reward Seeking

## Whole-CNS Identification of Circuit Vulnerabilities in *Drosophila*

Reproducibility materials and archived computational outputs accompanying the study by **Robert M. Frazier, PhD** (ORCID: **0009-0003-0861-4939**).

### Study overview

This study maps the downstream interface between prespecified reward-related and punishment-related dopaminergic populations across the complete adult male *Drosophila* central nervous system. The analysis begins with a frozen registry of **119 dopamine neurons**: 115 reward-related PAM neurons and 4 punishment-related PPL1 neurons. The seed populations generated **85,857 outgoing edges** reaching **34,172 unique downstream targets**, including **2,636 direct shared targets** receiving input from both systems.

The primary degree/strength-matched dopamine-neuron null model showed that the shared interface was substantially larger and more synaptically enriched than expected from generic DAN projection breadth and strength. The observed shared-target count was **2,636** versus a matched-null mean of **1,260.2658** (`p = 0.00009999`), and **1,181** shared targets were contacted by all four punishment neurons versus a matched-null mean of **2.7816** (`p = 0.00009999`). Mean source-normalized balance was **0.185645** versus **0.435780** under the matched null (`p = 0.009499`, lower-tail), indicating stronger polarization of the shared interface.

A secondary within-seed label permutation provides a complementary sensitivity analysis. Under that null, the total shared-target count was not individually extreme (`p = 0.13349`), whereas common targeting by all four punishment neurons and the reduced balance statistic remained highly unusual. These distinct nulls address different inferential questions and are retained separately in the repository.

Structural ranking and connectome-constrained virtual lesions were then used to prioritize experimentally testable candidate control points. MBON11 emerged as a particularly selective candidate in the virtual-lesion screen, while APL produced the largest aggregate robustness score but a less selective effect profile.

### Repository contents

```text
.
├── README.md
├── CITATION.cff
├── LICENSE
├── requirements.txt
├── R_REQUIREMENTS.txt
├── data/
│   └── seeds/
│       └── seed_registry_v1.csv
├── scripts/
│   ├── finish_study_master.py
│   ├── make_figures.R
│   └── verify_archived_results.py
├── results/
│   ├── direct_targets_all.csv
│   ├── direct_convergence_shared.csv
│   ├── direct_convergence_dominance.csv
│   ├── matched_dan_null_10000.csv
│   ├── permutation_null_10000.csv
│   ├── structural_candidate_ranking.csv
│   ├── structural_candidate_types.csv
│   ├── virtual_lesion_candidate_types.csv
│   ├── virtual_lesion_results.csv
│   ├── virtual_lesion_summary.csv
│   ├── summaries/
│   └── figures/
├── logs/
└── docs/
    ├── DATA_DICTIONARY.md
    └── REPRODUCIBILITY.md
```

### Reproducibility

The recorded analysis environment and exact FlyBrain Git commit are documented in [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md). The raw MaleCNS/FlyBrain source files are not redistributed in this repository. The master pipeline uses the frozen seed registry in `data/seeds/` and expects the source MaleCNS/FlyBrain resources to be supplied separately.

To verify the archived result package without downloading the raw connectome:

```bash
python scripts/verify_archived_results.py
```

### Statistical controls

Two null strategies are retained because they answer different questions:

1. **Degree/strength-matched DAN null:** holds the 115 reward DANs fixed and replaces each punishment DAN with a non-seed dopamine neuron matched on outgoing target degree and outgoing synaptic strength.
2. **Within-seed label permutation:** retains the same 119 DANs and all anatomical connections while randomly assigning 4 neurons as punishment and the remaining 115 as reward.

Empirical probabilities use the plus-one correction, `p = (k + 1) / (B + 1)`, with `B = 10,000` permutations.

### Interpretation boundary

The MaleCNS wiring used for structural inference is empirical. FlyBrain neural dynamics are a simplified leaky integrate-and-fire model. Virtual-lesion results are therefore interpreted as propagation and perturbation hypotheses, not as direct models of drug pharmacology, dopamine receptor biology, or a literal addicted-organism state.

### Citation

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). The publication DOI can be added there after acceptance/publication.

### Author

**Robert M. Frazier, PhD**  
This is Living Ministries, Research Divison
ORCID: 0009-0003-0861-4939
