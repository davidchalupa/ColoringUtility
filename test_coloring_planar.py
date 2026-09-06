import os
import pytest
import networkx as nx
import coloring_utility
from pathlib import Path

from tools.graph_generator import generate_random_planar_graph


script_dir = Path(os.path.dirname(__file__))


def test_coloring_on_random_planar_graph_10():
    G, _ = generate_random_planar_graph(10)

    expected_max_lower_bound = 4
    expected_colors = 4

    try:
        colors, lower_bound = coloring_utility.process(G, time_limit=60)
        num_colors = max(colors)
    except Exception as e:
        pytest.fail(f"An error occurred in coloring_utility: {e}")

    assert num_colors == expected_colors
    assert lower_bound <= expected_max_lower_bound


def test_coloring_on_random_planar_graph_100():
    G, _ = generate_random_planar_graph(100)

    expected_max_lower_bound = 4
    expected_colors = 4

    try:
        colors, lower_bound = coloring_utility.process(G, time_limit=60)
        num_colors = max(colors)
    except Exception as e:
        pytest.fail(f"An error occurred in coloring_utility: {e}")

    assert num_colors == expected_colors
    assert lower_bound <= expected_max_lower_bound


def test_coloring_on_random_planar_graph_500():
    G, _ = generate_random_planar_graph(500)

    expected_max_lower_bound = 4
    expected_colors = 4

    try:
        colors, lower_bound = coloring_utility.process(G, time_limit=60)
        num_colors = max(colors)
    except Exception as e:
        pytest.fail(f"An error occurred in coloring_utility: {e}")

    assert num_colors == expected_colors
    assert lower_bound <= expected_max_lower_bound


def test_coloring_on_random_planar_graph_1000():
    G, _ = generate_random_planar_graph(1000)

    expected_max_lower_bound = 4
    expected_colors = 4

    try:
        colors, lower_bound = coloring_utility.process(G, time_limit=60)
        num_colors = max(colors)
    except Exception as e:
        pytest.fail(f"An error occurred in coloring_utility: {e}")

    assert num_colors == expected_colors
    assert lower_bound <= expected_max_lower_bound


def test_coloring_on_random_planar_graph_2000():
    G, _ = generate_random_planar_graph(2000)

    expected_max_lower_bound = 4
    expected_colors = 4

    try:
        colors, lower_bound = coloring_utility.process(G, time_limit=60)
        num_colors = max(colors)
    except Exception as e:
        pytest.fail(f"An error occurred in coloring_utility: {e}")

    assert num_colors == expected_colors
    assert lower_bound <= expected_max_lower_bound


def test_coloring_on_random_planar_graph_10000():
    G, _ = generate_random_planar_graph(10000)

    expected_max_lower_bound = 4
    expected_colors = 4

    try:
        colors, lower_bound = coloring_utility.process(G, time_limit=60)
        num_colors = max(colors)
    except Exception as e:
        pytest.fail(f"An error occurred in coloring_utility: {e}")

    assert num_colors == expected_colors
    assert lower_bound <= expected_max_lower_bound
