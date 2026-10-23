# Method

- **Problem.** MaxCut on two graph families, 3-regular and two-community (planted partition, dense inside,
  sparse across), n in {8, 10, 12}, 4 random graphs (seeds 0-3) per family and size.
- **Algorithm.** Standard QAOA, p = 2. Angles are optimised once per graph in the noiseless simulator (COBYLA,
  2 restarts) and reused for every placement, so placements differ only through noise.
- **Device model.** Two chips of n/2 qubits, all-to-all inside a chip, one coupler gate per cross-chip edge
  (routing overhead is not modelled). Each cost-layer ZZ-phase gate is followed by two-qubit depolarising
  noise: p_intra = 0.005 inside a chip, p_inter = r * p_intra across the coupler, r in {1, 5, 20}. Mixer layers
  are noise-free. Noise is sampled as 200 Pauli trajectories per configuration (same random seed for every
  placement, so r = 1 gives identical numbers for all placements).
- **Placements.** (a) *random* balanced assignment (chip-oblivious, standard baseline); (b) *KL*
  Kernighan-Lin bisection (chip-aware, polynomial time); (c) *exact* minimum bisection by enumeration (upper reference).
  Reference lines: noiseless QAOA, and the r = 1 case, which is equivalent to a uniform-noise monolithic device.
- **Classical baseline.** Exact MaxCut value (enumeration) is the AR denominator.
- **Metric.** AR = <C>/C_max under noise; paired gains between placements; mean ± std over graphs.
- **Ablation / control.** r = 1 removes the fidelity asymmetry (placement should not matter); KL vs exact
  separates heuristic quality from the idea itself.
- **Noise sweep.** r in {1, 5, 20}. **Scaling.** n = 8, 10, 12.
