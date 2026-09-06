import networkx as nx


def load_from_col_file(file_path):
    G = nx.Graph()
    try:
        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip()
                if not line or line.startswith(('c', 'p')):
                    continue
                if line.startswith('e'):
                    parts = line.split()
                    u, v = int(parts[1]), int(parts[2])
                    if u != v:
                        G.add_edge(u, v)
    except Exception as e:
        print(f"File Load Error: {e}")
    return G
