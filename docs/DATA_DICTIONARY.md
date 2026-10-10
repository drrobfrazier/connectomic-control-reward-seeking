# Data dictionary

This file describes the archived tabular materials in the repository. Column names are preserved from the study outputs.

## `data/seeds/seed_registry_v1.csv`

Rows: **119**. Frozen study seed registry containing 115 reward DANs and 4 punishment DANs.

**Columns:** `body_id`, `flybrain_idx`, `type`, `instance`, `class`, `side`, `neurotransmitter`, `flybrain_sign`, `study_role`, `mb_compartment`

## Result tables

### `direct_convergence_dominance.csv`

Rows: **2,636**. Shared targets with source-normalized reward/punishment input, dominance, and balance metrics.

**Columns:** `body_post`, `type`, `instance`, `superclass`, `class`, `side`, `nt`, `reward_synapses`, `punishment_synapses`, `reward_sources`, `punishment_sources`, `mean_reward_synapses_per_seed`, `mean_punishment_synapses_per_seed`, `dominance_synaptic`, `balance_synaptic`, `reward_source_fraction`, `punishment_source_fraction`, `dominance_sources`, `balance_sources`, `kc_or_mbon`, `apl_or_dpm`

### `direct_convergence_shared.csv`

Rows: **2,636**. Direct downstream targets receiving input from both reward and punishment seed populations.

**Columns:** `body_post`, `target_type`, `target_instance`, `target_superclass`, `target_class`, `target_side`, `target_nt`, `identified`, `is_seed_neuron`, `reward_synapses`, `punishment_synapses`, `total_seed_synapses`, `reward_seed_neurons`, `punishment_seed_neurons`, `reward_edges`, `punishment_edges`, `reward_types`, `punishment_types`, `reward_fraction`, `minimum_convergent_synapses`, `direct_convergence`

### `direct_targets_all.csv`

Rows: **34,172**. All unique downstream targets reached by the prespecified reward and/or punishment seed populations.

**Columns:** `body_post`, `target_type`, `target_instance`, `target_superclass`, `target_class`, `target_side`, `target_nt`, `identified`, `is_seed_neuron`, `reward_synapses`, `punishment_synapses`, `total_seed_synapses`, `reward_seed_neurons`, `punishment_seed_neurons`, `reward_edges`, `punishment_edges`, `reward_types`, `punishment_types`, `reward_fraction`, `minimum_convergent_synapses`, `direct_convergence`

### `matched_dan_null_10000.csv`

Rows: **10,000**. 10,000 degree/strength-matched dopamine-neuron null draws used for the primary inference.

**Columns:** `permutation`, `shared_targets`, `all4_shared`, `mean_balance`, `shared_total_synapses`

### `permutation_null_10000.csv`

Rows: **10,000**. 10,000 within-seed label permutations retaining the same 119 DANs and randomly assigning 4 punishment labels.

**Columns:** `permutation`, `shared_targets`, `shared_total_synapses`, `all4_punishment_shared`, `kc_mbon_shared`, `apl_dpm_shared`, `mean_balance_shared`

### `structural_candidate_ranking.csv`

Rows: **2,636**. Cell-level structural interface ranking for the 2,636 shared downstream targets.

**Columns:** `body_post`, `type`, `instance`, `superclass`, `class`, `nt`, `reward_synapses`, `punishment_synapses`, `reward_sources`, `punishment_sources`, `reward_source_fraction`, `punishment_source_fraction`, `mean_reward_synapses_per_seed`, `mean_punishment_synapses_per_seed`, `dominance`, `balance`, `cross_source`, `cross_weight`, `dn_reach_1hop`, `dn_reach_2hop`, `dn_reach_3hop`, `dn_reach_weighted`, `interface_score`, `polarity`

### `structural_candidate_types.csv`

Rows: **240**. Cell-type aggregation of structural candidate metrics.

**Columns:** `type`, `shared_target_neurons`, `whole_type_size`, `max_interface_score`, `mean_interface_score`, `median_dominance`, `max_dn_reach`

### `virtual_lesion_candidate_types.csv`

Rows: **14**. Cell types entered into the virtual-lesion screen.

**Columns:** `type`, `n_neurons`, `role`

### `virtual_lesion_results.csv`

Rows: **124**. Condition- and parameter-specific intact and lesion simulation outputs.

**Columns:** `weight_scale`, `hz`, `condition`, `lesion_type`, `lesion_n`, `approach_rate_hz`, `avoidance_rate_hz`, `rdi`, `total_spikes`, `active_neurons`, `delta_rdi_vs_intact`

### `virtual_lesion_summary.csv`

Rows: **14**. Aggregated robustness metrics for the 14 lesion-tested cell types.

**Columns:** `type`, `n_neurons`, `median_conflict_delta_rdi`, `median_reward_delta_rdi`, `median_abs_conflict_delta`, `median_abs_reward_delta`, `conflict_effect_sign_consistency`, `conflict_specificity_absdiff`, `conflict_to_reward_effect_ratio`, `robust_lesion_score`
