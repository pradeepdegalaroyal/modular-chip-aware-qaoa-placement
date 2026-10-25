# Research gap

**State of the art.** Du et al. (2027, doi:10.1016/j.future.2026.108830) point out that compilers built for
monolithic chips ignore the sharp fidelity asymmetry between low-error intra-chip gates and high-error
chip-to-chip coupler gates, and propose hardware-aware compilation for modular systems. Cai et al. (2023,
doi:10.1103/RevModPhys.95.045005) review error mitigation and stress that mitigation cost grows quickly with
circuit noise, so reducing the number of *bad* gates is preferable to mitigating them afterwards. Wang et al.
(2027, doi:10.15302/frontphys.2027.023301) call for a systematic reassessment of variational algorithms as
architectures move toward early fault tolerance, and Wołk et al. (2027, doi:10.1016/j.ins.2026.124187) work on
Quantum Approximate Optimization Algorithm (QAOA) design for constrained problems. What the compilation work does not quantify is how much *application-level*
QAOA quality a chip-aware placement recovers, and how that depends on the coupler/intra-chip error ratio and on
the problem structure.

**Gap.** For QAOA (a variational algorithm whose entangling-gate pattern is the problem graph itself), the
vertex-to-chip assignment fixes how many cost-layer gates cross a coupler. We test whether a cheap
chip-aware assignment (graph bisection minimising coupler gates) recovers a useful fraction of the approximation
ratio lost to coupler noise, relative to the chip-oblivious placement that a monolithic-style compiler gives.

**Hypothesis (falsifiable).** At p_inter/p_intra = 20 (p_intra = 0.005), Kernighan-Lin bisection raises the noisy
MaxCut approximation ratio (AR) of p=2 QAOA over random balanced placement by at least 0.02 (paired, mean over
graphs), for n = 8-12 qubits on two chips; at p_inter/p_intra = 1 the gain is exactly zero (control). Success
metric: paired AR gain, mean ± std over 4 graphs per family and size. Also reported: gap between Kernighan-Lin (KL) and the exact
minimum bisection.

**Scale.** Everything is statevector / Pauli-trajectory simulation with n <= 12 qubits on two chips of n/2
qubits. The study does not claim anything about real hardware or about quantum advantage.
