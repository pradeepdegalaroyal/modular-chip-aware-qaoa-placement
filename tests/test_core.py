import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from modqaoa import (cut_edges, exact_maxcut, exact_partition, kl_partition, make_graph,  # noqa: E402
                     maxcut_diag, noisy_expectation, qaoa_state, random_partition)


def test_norm_and_noiseless_limit():
    g = make_graph("regular3", 8, 0)
    th = np.array([0.4, 0.3, 0.2, 0.1])
    part = random_partition(8, np.random.default_rng(0))
    psi = qaoa_state(g, th, 0.05, 0.1, part, 8, np.random.default_rng(1))
    assert np.allclose(np.linalg.norm(psi, axis=1), 1)
    a = noisy_expectation(g, th, maxcut_diag(g), part, 0, 0, 1, 0)
    assert 0 < a <= exact_maxcut(g)


def test_noise_degrades_good_circuit():
    g = make_graph("regular3", 8, 1)
    d = maxcut_diag(g)
    th = np.array([0.6, 0.4, 0.3, 0.2])
    part = random_partition(8, np.random.default_rng(0))
    clean = noisy_expectation(g, th, d, part, 0, 0, 1, 0)
    noisy = noisy_expectation(g, th, d, part, 0.5, 0.5, 300, 0)
    assert noisy < clean + 1e-9


def test_partitions():
    g = make_graph("community", 8, 0)
    ex = exact_partition(g)
    assert ex.sum() == g.number_of_nodes() // 2
    assert cut_edges(g, ex) <= cut_edges(g, kl_partition(g, 0))
