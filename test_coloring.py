import pytest
import networkx as nx
import coloring_utility


def test_barabasi_albert_100_4_seed_142():
    # ToDo: this seems to fail
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


#def test_barabasi_albert_100000_4_seed_142():
#    n = 100000
#    w = 4
#    G = nx.barabasi_albert_graph(n=n, m=w, seed=142)

#    expected_lower_bound = 6
#    expected_colors = 6

#    try:
#        colors, lower_bound = coloring_utility.process(G, time_limit=60)
#        num_colors = max(colors)
#    except Exception as e:
#        pytest.fail(f"An error occurred in coloring_utility: {e}")

#    assert num_colors == expected_colors
#    assert lower_bound == expected_lower_bound


def test_erdos_renyi_100_0p1_seed_42():
    G = nx.erdos_renyi_graph(100, 0.1, seed=42)

    expected_lower_bound = 5
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
