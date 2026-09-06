import os
import pytest
import networkx as nx
import coloring_utility
from pathlib import Path

from tools.graph_generator import generate_random_planar_graph


num_small_planar_graphs = 5
script_dir = Path(os.path.dirname(__file__))


@pytest.mark.parametrize("test_case_id", range(1, num_small_planar_graphs + 1))
def test_coloring_on_random_planar_graph_10(test_case_id):
    G, _ = generate_random_planar_graph(10)

    expected_max_colors = 4

    try:
        colors, lower_bound = coloring_utility.process(G, time_limit=60)
        num_colors = max(colors)
    except Exception as e:
        pytest.fail(f"An error occurred in coloring_utility: {e}")

    assert num_colors <= expected_max_colors
    assert lower_bound <= num_colors


@pytest.mark.parametrize("test_case_id", range(1, num_small_planar_graphs + 1))
def test_coloring_on_random_planar_graph_100(test_case_id):
    G, _ = generate_random_planar_graph(100)

    expected_max_colors = 4

    try:
        colors, lower_bound = coloring_utility.process(G, time_limit=60)
        num_colors = max(colors)
    except Exception as e:
        pytest.fail(f"An error occurred in coloring_utility: {e}")

    assert num_colors <= expected_max_colors
    assert lower_bound <= num_colors



def test_coloring_on_random_planar_graph_500():
    G, _ = generate_random_planar_graph(500)

    expected_max_colors = 4

    try:
        colors, lower_bound = coloring_utility.process(G, time_limit=60)
        num_colors = max(colors)
    except Exception as e:
        pytest.fail(f"An error occurred in coloring_utility: {e}")

    assert num_colors <= expected_max_colors
    assert lower_bound <= num_colors



def test_coloring_on_random_planar_graph_1000():
    G, _ = generate_random_planar_graph(1000)

    expected_max_colors = 4

    try:
        colors, lower_bound = coloring_utility.process(G, time_limit=60)
        num_colors = max(colors)
    except Exception as e:
        pytest.fail(f"An error occurred in coloring_utility: {e}")

    assert num_colors <= expected_max_colors
    assert lower_bound <= num_colors



def test_coloring_on_random_planar_graph_2000():
    G, _ = generate_random_planar_graph(2000)

    expected_max_colors = 4

    try:
        colors, lower_bound = coloring_utility.process(G, time_limit=60)
        num_colors = max(colors)
    except Exception as e:
        pytest.fail(f"An error occurred in coloring_utility: {e}")

    assert num_colors <= expected_max_colors
    assert lower_bound <= num_colors



def test_coloring_on_random_planar_graph_10000():
    G, _ = generate_random_planar_graph(10000)

    expected_max_colors = 4

    try:
        colors, lower_bound = coloring_utility.process(G, time_limit=60)
        num_colors = max(colors)
    except Exception as e:
        pytest.fail(f"An error occurred in coloring_utility: {e}")

    assert num_colors <= expected_max_colors
    assert lower_bound <= num_colors

