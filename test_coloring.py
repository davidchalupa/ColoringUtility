import os
import pytest
import networkx as nx
import coloring_utility
from pathlib import Path

from tools.graph_loader import load_from_col_file


script_dir = Path(os.path.dirname(__file__))


def test_davis_southern_women():
    G = nx.davis_southern_women_graph()

    expected_lower_bound = 2
    expected_colors = 2

    try:
        colors, lower_bound = coloring_utility.process(G, time_limit=60)
        num_colors = max(colors)
    except Exception as e:
        pytest.fail(f"An error occurred in coloring_utility: {e}")

    assert num_colors == expected_colors
    assert lower_bound == expected_lower_bound


def test_karate_club():
    G = nx.karate_club_graph()

    expected_lower_bound = 5
    expected_colors = 5

    try:
        colors, lower_bound = coloring_utility.process(G, time_limit=60)
        num_colors = max(colors)
    except Exception as e:
        pytest.fail(f"An error occurred in coloring_utility: {e}")

    assert num_colors == expected_colors
    assert lower_bound == expected_lower_bound


def test_les_miserables_graph():
    G = nx.les_miserables_graph()

    expected_lower_bound = 10
    expected_colors = 10

    try:
        colors, lower_bound = coloring_utility.process(G, time_limit=60)
        num_colors = max(colors)
    except Exception as e:
        pytest.fail(f"An error occurred in coloring_utility: {e}")

    assert num_colors == expected_colors
    assert lower_bound == expected_lower_bound


def test_smallest_hard_to_color_brelaz():
    """
    Smallest 3-colorable graph", for which Brelaz's heuristic gives 4 colors.
    Janczewski, R., Kubale, M., Manuszewski, K., & Piwakowski, K. (2001). The smallest hard-to-color graph for algorithm DSATUR. Discrete Mathematics, 236(1-3), 151-165.
    """
    G = load_from_col_file(script_dir / "data" / "badbre.col")

    expected_lower_bound = 3
    expected_colors = 3

    try:
        colors, lower_bound = coloring_utility.process(G, time_limit=60)
        num_colors = max(colors)
    except Exception as e:
        pytest.fail(f"An error occurred in coloring_utility: {e}")

    assert num_colors == expected_colors
    assert lower_bound == expected_lower_bound


def test_mycielski_graph_4():
    G = nx.mycielski_graph(4)

    expected_lower_bound = 4
    expected_colors = 4

    try:
        colors, lower_bound = coloring_utility.process(G, time_limit=60)
        num_colors = max(colors)
    except Exception as e:
        pytest.fail(f"An error occurred in coloring_utility: {e}")

    assert num_colors == expected_colors
    assert lower_bound == expected_lower_bound


def test_mycielski_graph_5():
    G = nx.mycielski_graph(5)

    expected_lower_bound = 5
    expected_colors = 5

    try:
        colors, lower_bound = coloring_utility.process(G, time_limit=60)
        num_colors = max(colors)
    except Exception as e:
        pytest.fail(f"An error occurred in coloring_utility: {e}")

    assert num_colors == expected_colors
    assert lower_bound == expected_lower_bound


def test_mycielski_graph_6():
    G = nx.mycielski_graph(6)

    expected_lower_bound = 5
    expected_colors = 6

    try:
        colors, lower_bound = coloring_utility.process(G, time_limit=60)
        num_colors = max(colors)
    except Exception as e:
        pytest.fail(f"An error occurred in coloring_utility: {e}")

    assert num_colors == expected_colors
    assert lower_bound >= expected_lower_bound


def test_leighton_graph_450_5a():
    G = load_from_col_file(script_dir / "data" / "dimacs" / "le450_5a.col")

    expected_lower_bound = 5
    expected_colors = 5

    try:
        colors, lower_bound = coloring_utility.process(G, time_limit=60)
        num_colors = max(colors)
    except Exception as e:
        pytest.fail(f"An error occurred in coloring_utility: {e}")

    assert num_colors == expected_colors
    assert lower_bound == expected_lower_bound


def test_leighton_graph_450_15c():
    G = load_from_col_file(script_dir / "data" / "dimacs" / "le450_15c.col")

    expected_lower_bound = 15
    expected_colors = 15

    try:
        colors, lower_bound = coloring_utility.process(G, time_limit=60)
        num_colors = max(colors)
    except Exception as e:
        pytest.fail(f"An error occurred in coloring_utility: {e}")

    assert num_colors == expected_colors
    assert lower_bound == expected_lower_bound


def test_leighton_graph_450_15d():
    G = load_from_col_file(script_dir / "data" / "dimacs" / "le450_15d.col")

    expected_lower_bound = 15
    expected_colors = 15

    try:
        colors, lower_bound = coloring_utility.process(G, time_limit=120)
        num_colors = max(colors)
    except Exception as e:
        pytest.fail(f"An error occurred in coloring_utility: {e}")

    assert num_colors == expected_colors
    assert lower_bound == expected_lower_bound


def test_random_graph_dsjc_500p1():
    G = load_from_col_file(script_dir / "data" / "dimacs" / "dsjc500.1.col")

    expected_lower_bound = 6
    expected_colors = 12

    try:
        colors, lower_bound = coloring_utility.process(G, time_limit=60)
        num_colors = max(colors)
    except Exception as e:
        pytest.fail(f"An error occurred in coloring_utility: {e}")

    assert num_colors == expected_colors
    assert lower_bound >= expected_lower_bound


def test_barabasi_albert_100_4_seed_142():
    n = 100
    w = 4
    G = nx.barabasi_albert_graph(n=n, m=w, seed=142)

    expected_lower_bound = 5
    expected_colors = 5

    try:
        colors, lower_bound = coloring_utility.process(G, time_limit=60)
        num_colors = max(colors)
    except Exception as e:
        pytest.fail(f"An error occurred in coloring_utility: {e}")

    assert num_colors == expected_colors
    assert lower_bound == expected_lower_bound


def test_barabasi_albert_100_3_seed_142():
    n = 100
    w = 3
    G = nx.barabasi_albert_graph(n=n, m=w, seed=142)

    expected_lower_bound = 4
    expected_colors = 4

    try:
        colors, lower_bound = coloring_utility.process(G, time_limit=60)
        num_colors = max(colors)
    except Exception as e:
        pytest.fail(f"An error occurred in coloring_utility: {e}")

    assert num_colors == expected_colors
    assert lower_bound == expected_lower_bound


def test_barabasi_albert_100_5_seed_442():
    n = 100
    w = 5
    G = nx.barabasi_albert_graph(n=n, m=w, seed=442)

    expected_lower_bound = 6
    expected_colors = 6

    try:
        colors, lower_bound = coloring_utility.process(G, time_limit=60)
        num_colors = max(colors)
    except Exception as e:
        pytest.fail(f"An error occurred in coloring_utility: {e}")

    assert num_colors == expected_colors
    assert lower_bound == expected_lower_bound


def test_erdos_renyi_100_0p1_seed_42():
    G = nx.erdos_renyi_graph(100, 0.1, seed=42)

    expected_lower_bound = 4
    expected_colors = 5

    try:
        colors, lower_bound = coloring_utility.process(G, time_limit=60)
        num_colors = max(colors)
    except Exception as e:
        pytest.fail(f"An error occurred in coloring_utility: {e}")

    assert num_colors == expected_colors
    assert lower_bound >= expected_lower_bound


def test_erdos_renyi_200_0p1_seed_42():
    G = nx.erdos_renyi_graph(200, 0.1, seed=42)

    expected_lower_bound = 5
    expected_colors = 7

    try:
        colors, lower_bound = coloring_utility.process(G, time_limit=60)
        num_colors = max(colors)
    except Exception as e:
        pytest.fail(f"An error occurred in coloring_utility: {e}")

    assert num_colors == expected_colors
    assert lower_bound >= expected_lower_bound

def test_barabasi_albert_100000_4_seed_42():
    n = 100000
    w = 4
    G = nx.barabasi_albert_graph(n=n, m=w, seed=42)

    expected_lower_bound = 5
    expected_colors = 5

    try:
        colors, lower_bound = coloring_utility.process(G, time_limit=60)
        num_colors = max(colors)
    except Exception as e:
        pytest.fail(f"An error occurred in coloring_utility: {e}")

    assert num_colors == expected_colors
    assert lower_bound == expected_lower_bound


def test_barabasi_albert_100000_4_seed_142():
    n = 100000
    w = 4
    G = nx.barabasi_albert_graph(n=n, m=w, seed=142)

    expected_lower_bound = 4
    expected_colors = 5

    try:
        colors, lower_bound = coloring_utility.process(G, time_limit=60)
        num_colors = max(colors)
    except Exception as e:
        pytest.fail(f"An error occurred in coloring_utility: {e}")

    assert num_colors == expected_colors
    assert lower_bound >= expected_lower_bound
