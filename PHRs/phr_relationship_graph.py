import sys
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import math
import seaborn as sns
import io
import pandas as pd 
from matplotlib import cm

NODE_SIZE=750

def compute_circular_positions(nodes, scale=1, wiggle_factor=0.025):
	angle = 2 * math.pi / len(nodes)  # Angle between each node
	pos = {}
	for i, node in enumerate(nodes):
		x = math.cos(i * angle)
		y = math.sin(i * angle)
		pos[node] = (x, y)
		radius = 1 + wiggle_factor if i % 2 == 0 else 1 - wiggle_factor
		pos[node] = (radius * x * scale, radius * y * scale)
	return pos

def get_order(species):
	if species == 'mGorGor1':
		nor_plus_acros = ['21', '22']
		nor_minus_acros = ['2a', '9', '13', '15', '18', 'Y']
		ancient_acros = ['2b', '14']
		custom_order=['1', '2a', '2b', '3', '4', '5x17', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17x5', '18', '19', '20', '21', '22', 'X', 'Y']
	elif species == 'HG002' or species == 'chm13':
		nor_plus_acros = ['13', '14', '15', '21', '22']
		nor_minus_acros = ['Y']
		ancient_acros = ['2a', '2b', '9', '18']
		custom_order=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
	elif species == 'mPanPan1' or species == 'mPanTro3':
		nor_plus_acros = ['13', '14', '18', '21', '22']
		nor_minus_acros = ['Y']
		ancient_acros = ['2a', '2b', '9', '15']
		custom_order=['1', '2a', '2b', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
	elif species == 'mSymSyn1':
		nor_plus_acros = ['21', 'Y']
		nor_minus_acros = []
		ancient_acros = []
		custom_order=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', 'X', 'Y']
	else:
		nor_plus_acros = ['2a', '2b', '9', '13', '14', '15', '18', '21', '22', 'Y']
		nor_minus_acros = []
		ancient_acros = []
		custom_order=['1', '2a', '2b', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
	return nor_plus_acros, nor_minus_acros, ancient_acros, custom_order

def rescale_edge_positions(pos, start, end, radius):
	"""
	Rescale edge positions to connect to node boundaries instead of centers.
	"""
	x0, y0 = pos[start]
	x1, y1 = pos[end]
	dx, dy = x1 - x0, y1 - y0
	dist = np.sqrt(dx**2 + dy**2)
	# Scale dx, dy to stop at the boundary of the node
	dx, dy = dx / dist * radius, dy / dist * radius
	return (x0 + dx, y0 + dy), (x1 - dx, y1 - dy)

def node_radius_data_units(node_size, ax, fig):
	"""
	Convert a node's size (in points^2) to data-unit radius for clipping edges.
	"""
	# node_size is area in points^2 → radius in points
	radius_pts = np.sqrt(node_size / np.pi)
	
	# Convert points to pixels
	pts_to_px = 1 / 72 * fig.dpi
	radius_px = radius_pts * pts_to_px
	
	# Convert pixels → data coordinates
	xlim = ax.get_xlim()
	ylim = ax.get_ylim()
	xmin, xmax = xlim
	ymin, ymax = ylim
	
	dx_data = (xmax - xmin) / ax.bbox.width * radius_px
	dy_data = (ymax - ymin) / ax.bbox.height * radius_px
	
	# Use the average scaling (circle)
	return (dx_data + dy_data) / 2

relationship_tsv = sys.argv[1]
species = sys.argv[2]

dtype_dict = {
	'species': 'str',
	'chrom1': 'str',
	'chrom2': 'str',
	'bp': 'float',
	'pid': 'float'
}

df = pd.read_csv(relationship_tsv, sep='\t', names=['species', 'chrom1', 'chrom2', 'bp', 'pid'])
df['chrom1'] = df['chrom1'].apply(lambda x: x.replace('chr', ''))
df['chrom2'] = df['chrom2'].apply(lambda x: x.replace('chr', ''))

mean_bp = np.mean(df['bp'])
mean_pid = np.mean(df['pid'])
df['normalized_bp'] = df['bp'] / 250000 # mean_bp
df['normalized_pid'] = df['pid'] - 99
cmap = cm.get_cmap('YlOrRd')
df['edge_color'] = df['normalized_pid'].apply(lambda pid: cmap(pid))
df['edge_weight'] = df['normalized_bp']
df_subset = df[df['species'] == species]
df_subset.sort_values('normalized_pid', inplace=True, ascending=False)

print(df_subset)

G = nx.Graph()

for _, row in df_subset.iterrows():
	chrom1 = row['chrom1']
	chrom2 = row['chrom2']
	normalized_bp = row['normalized_bp']
	normalized_pid = row['normalized_pid']
	weight = row['edge_weight']
	color = row['edge_color']
	G.add_edge(chrom1, chrom2, norm_bp=normalized_bp, norm_pid=normalized_pid, weight=weight, color=color)

nor_plus_acros, nor_minus_acros, ancient_acros, custom_order = get_order(species)

for node in custom_order:
	if node not in G.nodes:
		G.add_node(node)

pos = compute_circular_positions(custom_order, 0.5)
edge_colors = [G[u][v]['color'] for u, v in G.edges()]
edge_weights = [G[u][v]['weight'] for u, v in G.edges()]

node_colors = []
for node in G.nodes():
	if node in nor_plus_acros:
		node_colors.append('#003863')
	elif node in nor_minus_acros:
		node_colors.append('#007BD1')
	elif node in ancient_acros:
		node_colors.append('#64D6FF')
	else:
		node_colors.append('#C8C8C8')

node_labels = {node: node for node in G.nodes()}
font_colors = {node: '#D0D0D0' if node in nor_plus_acros else 'black' for node in G.nodes()}

plt.figure(figsize=(10, 10))
nx.draw_networkx_nodes(G, pos, node_color=node_colors, edgecolors=node_colors, node_size=NODE_SIZE)

edges_sorted = sorted(
	G.edges(data=True),
	key=lambda x: x[2]["norm_pid"],
	reverse=False  # lowest first → highest drawn last = on top
)

ax = plt.gca()
fig = plt.gcf()

node_size = NODE_SIZE
node_r = node_radius_data_units(node_size, ax, fig)



for u, v, data in edges_sorted:
	x0, y0 = pos[u]
	x1, y1 = pos[v]

	dx = x1 - x0
	dy = y1 - y0
	dist = np.sqrt(dx**2 + dy**2)

	# shrink the line at both ends
	shrink = node_r

	x0_new = x0 + dx/dist * shrink
	y0_new = y0 + dy/dist * shrink
	x1_new = x1 - dx/dist * shrink
	y1_new = y1 - dy/dist * shrink

	line = mlines.Line2D(
		[x0_new, x1_new],
		[y0_new, y1_new],
		linewidth=data["weight"],
		color=data["color"],
		antialiased=False,
		zorder=1
	)
	ax.add_line(line)

# nx.draw_networkx_edges(G, pos, width=edge_weights, edge_color=edge_colors)

for node, label in node_labels.items():
	nx.draw_networkx_labels(
		G, pos,
		labels={node: label},  # Only the current node
		font_size=16,
		font_color=font_colors[node],
		font_family='Arial'
	)

plt.title(species + ' Relationship Graph')
plt.axis('off')
plt.show()
plt.savefig(sys.argv[3], format='svg')
plt.clf()
