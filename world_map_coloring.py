import urllib.request
import zipfile
import os
import networkx as nx
import geopandas as gpd
import matplotlib.pyplot as plt

# 1. Download and load the Natural Earth dataset manually
url = "https://naturalearth.s3.amazonaws.com/110m_cultural/ne_110m_admin_0_countries.zip"
zip_path = "ne_110m_admin_0_countries.zip"
extract_dir = "ne_110m_data"

# Download and extract the shapefiles if they don't exist yet
if not os.path.exists(extract_dir):
    print("Downloading Natural Earth data...")
    urllib.request.urlretrieve(url, zip_path)
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(extract_dir)
    print("Download complete.")

# Load the extracted shapefile
shp_file = os.path.join(extract_dir, "ne_110m_admin_0_countries.shp")
world = gpd.read_file(shp_file)

# Clean the data: Remove Antarctica and empty geometries
# Note: The raw Natural Earth dataset uses uppercase column names
world = world[(world.POP_EST > 0) & (world.NAME != "Antarctica")]

# 2. Reconstruct the graph
G = nx.Graph()

# Add all countries as nodes
for index, row in world.iterrows():
    G.add_node(row['NAME'])

# Add edges based on shared borders (spatial intersection)
for index, row in world.iterrows():
    # Find all countries that intersect with the current country, excluding itself
    neighbors = world[world.geometry.intersects(row.geometry) & (world.NAME != row['NAME'])]

    for neighbor_name in neighbors['NAME']:
        G.add_edge(row['NAME'], neighbor_name)

print(f"Graph reconstructed with {G.number_of_nodes()} countries and {G.number_of_edges()} shared borders.")

# # 3. Apply Graph Coloring
# color_map = nx.coloring.greedy_color(G, strategy='largest_first')

import coloring_utility

try:
    colors, lower_bound = coloring_utility.process(G, time_limit=120)
    num_colors = max(colors)
except Exception as e:
    raise "An error occurred in coloring_utility: {e}"
color_map = dict(zip(G.nodes(), colors))

num_colors_used = len(set(color_map.values()))
print(f"Graph colored using {num_colors_used} distinct colors.")

# 4. Visualize the result on a map
world['color_id'] = world['NAME'].map(color_map)

fig, ax = plt.subplots(1, 1, figsize=(16, 8))
world.plot(
    column='color_id',
    ax=ax,
    cmap='Set3',
    edgecolor='black',
    linewidth=0.5
)

ax.set_title("World Map Colored via NetworkX Graph Coloring", fontsize=16)
ax.set_axis_off()
plt.tight_layout()
plt.show()
