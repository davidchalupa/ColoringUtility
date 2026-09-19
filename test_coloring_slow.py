import os
import pytest
import networkx as nx
import coloring_utility
from pathlib import Path

from tools.graph_loader import load_from_col_file


script_dir = Path(os.path.dirname(__file__))


def test_random_graph_dsjc_500p1():
    G = load_from_col_file(script_dir / "data" / "dimacs" / "dsjc500.1.col")

    expected_lower_bound = 6
    expected_colors = 12

    try:
        colors, lower_bound = coloring_utility.process(G, time_limit=180)
        num_colors = max(colors)
    except Exception as e:
        pytest.fail(f"An error occurred in coloring_utility: {e}")

    assert num_colors == expected_colors
    assert lower_bound >= expected_lower_bound
