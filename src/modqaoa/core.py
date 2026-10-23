"""Core simulator: MaxCut QAOA with chip-dependent two-qubit Pauli noise.

Each cost-layer ZZ-phase gate acts on one graph edge. If the two endpoints sit on the
same chip the gate is followed by a two-qubit depolarising error of strength ``p_intra``;
if they sit on different chips (a coupler gate) the strength is ``p_inter``.
Noise is sampled as Pauli trajectories on a batched statevector. Mixer layers are
noise-free (identical for every placement, so they do not affect comparisons).
"""
from __future__ import annotations

from itertools import combinations

import networkx as nx
import numpy as np
from scipy.optimize import minimize

Partition = np.ndarray  # chip label (0/1) for every vertex


def make_graph(kind: str, n: int, seed: int) -> nx.Graph:
    """Random test graph on n vertices: 3-regular or two-community planted partition."""
    if kind == "regular3":
        return nx.random_regular_graph(3, n, seed=seed)
    if kind == "community":
        h = n // 2
        for k in range(100):
            g = nx.Graph(nx.stochastic_block_model([h, h], [[0.75, 0.12], [0.12, 0.75]],
                                                   seed=seed + 1000 * k))
            if nx.number_of_isolates(g) == 0:
                return g
        raise RuntimeError("could not draw a graph without isolated vertices")
    raise ValueError(kind)


def maxcut_diag(g: nx.Graph) -> np.ndarray:
    """Cut value of every bitstring (length 2**n)."""
    n = g.number_of_nodes()
    idx = np.arange(2**n)
    c = np.zeros(2**n)
    for i, j in g.edges:
        c += ((idx >> i) & 1) ^ ((idx >> j) & 1)
    return c


def exact_maxcut(g: nx.Graph) -> float:
    """Exact MaxCut value by enumeration."""
    return float(maxcut_diag(g).max())


def _mixer(psi: np.ndarray, beta: float, n: int) -> np.ndarray:
    """Apply RX(2*beta) on every qubit of a (T, 2**n) batch."""
    t = psi.shape[0]
    c, s = np.cos(beta), -1j * np.sin(beta)
    for q in range(n):
        v = psi.reshape(t, 2 ** (n - q - 1), 2, 2**q)
        a, b = v[:, :, 0, :], v[:, :, 1, :]
        psi = np.stack([c * a + s * b, s * a + c * b], axis=2).reshape(t, 2**n)
    return psi


def _parity(x: np.ndarray) -> np.ndarray:
    x = x.copy()
    for s in (16, 8, 4, 2, 1):
        x ^= x >> s
    return x & 1


def _depolarise2(psi: np.ndarray, i: int, j: int, rate: float, idx: np.ndarray,
                 rng: np.random.Generator) -> np.ndarray:
    """Two-qubit depolarising error: with prob. ``rate`` one of 15 non-identity Paulis."""
    if rate <= 0:
        return psi
    hit = np.flatnonzero(rng.random(psi.shape[0]) < rate)
    if hit.size == 0:
        return psi
    kinds = rng.integers(1, 16, size=hit.size)
    out = psi.copy()
    for k in np.unique(kinds):
        sel = hit[kinds == k]
        a, b = int(k) & 3, int(k) >> 2  # 2-bit codes per qubit: bit0 = X part, bit1 = Z part
        xm = (a & 1) << i | (b & 1) << j
        zm = ((a >> 1) & 1) << i | ((b >> 1) & 1) << j
        sign = 1 - 2 * _parity(idx & zm)
        out[sel] = (psi[sel] * sign)[:, idx ^ xm]
    return out


def qaoa_state(g: nx.Graph, params: np.ndarray, p_intra: float = 0.0, p_inter: float = 0.0,
               part: Partition | None = None, traj: int = 1,
               rng: np.random.Generator | None = None) -> np.ndarray:
    """Return a (traj, 2**n) batch of final states (noiseless if both rates are 0).

    ``params`` = [gamma_1..gamma_p, beta_1..beta_p].
    """
    n = g.number_of_nodes()
    p = len(params) // 2
    idx = np.arange(2**n)
    psi = np.full((traj, 2**n), 2 ** (-n / 2), dtype=complex)
    noisy = p_intra > 0 or p_inter > 0
    rng = rng or np.random.default_rng(0)
    for layer in range(p):
        for i, j in g.edges:
            cut = ((idx >> i) & 1) ^ ((idx >> j) & 1)
            psi = psi * np.exp(-1j * params[layer] * cut)
            if noisy:
                assert part is not None
                psi = _depolarise2(psi, i, j, p_inter if part[i] != part[j] else p_intra, idx, rng)
        psi = _mixer(psi, params[p + layer], n)
    return psi


def noiseless_expectation(g: nx.Graph, params: np.ndarray, diag: np.ndarray) -> float:
    """Ideal <C> of the QAOA state."""
    psi = qaoa_state(g, params)[0]
    return float(np.abs(psi) ** 2 @ diag)


def optimise_angles(g: nx.Graph, p: int, diag: np.ndarray, seed: int, restarts: int = 4) -> np.ndarray:
    """Noiseless angle optimisation (COBYLA with random restarts)."""
    rng = np.random.default_rng(seed)
    best, bv = None, -1.0
    for _ in range(restarts):
        x0 = rng.uniform(0, np.pi / 2, 2 * p)
        r = minimize(lambda x: -noiseless_expectation(g, x, diag), x0, method="COBYLA",
                     options={"maxiter": 120})
        if -r.fun > bv:
            best, bv = r.x, -r.fun
    return best


def noisy_expectation(g: nx.Graph, params: np.ndarray, diag: np.ndarray, part: Partition,
                      p_intra: float, p_inter: float, traj: int, seed: int) -> float:
    """Trajectory-averaged <C> under chip-dependent noise."""
    rng = np.random.default_rng(seed)
    psi = qaoa_state(g, params, p_intra, p_inter, part, traj, rng)
    return float(((np.abs(psi) ** 2) @ diag).mean())


# ---- placement strategies -------------------------------------------------------------
def cut_edges(g: nx.Graph, part: Partition) -> int:
    """Number of graph edges whose endpoints sit on different chips (coupler gates)."""
    return sum(int(part[i] != part[j]) for i, j in g.edges)


def random_partition(n: int, rng: np.random.Generator) -> Partition:
    """Random balanced assignment of vertices to two chips (chip-oblivious placement)."""
    part = np.zeros(n, dtype=int)
    part[rng.permutation(n)[: n // 2]] = 1
    return part


def kl_partition(g: nx.Graph, seed: int) -> Partition:
    """Kernighan-Lin balanced bisection (chip-aware heuristic)."""
    a, _ = nx.algorithms.community.kernighan_lin_bisection(g, seed=seed)
    part = np.zeros(g.number_of_nodes(), dtype=int)
    part[list(a)] = 1
    return part


def exact_partition(g: nx.Graph) -> Partition:
    """Exact minimum balanced bisection by enumeration (n <= 14)."""
    n = g.number_of_nodes()
    best, bc = None, 10**9
    for comb in combinations(range(1, n), n // 2 - 1):
        part = np.zeros(n, dtype=int)
        part[[0, *comb]] = 1
        c = cut_edges(g, part)
        if c < bc:
            best, bc = part, c
    return best
