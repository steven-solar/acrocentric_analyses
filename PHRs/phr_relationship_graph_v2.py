import sys
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import math
import seaborn as sns
import io
import pandas as pd 
from matplotlib import cm

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

def get_chrom(chrom, species):
	if species == 'chm13':
		return chrom.replace('chr', '')
	elif species in ['HG002', 'mSymSyn1']:
		return chrom.split('_')[0].replace('chr', '')
	else:
		return chrom.split('_')[2].replace('hsa', '')
		
def rename_chroms(row):
	species = row['species']
	row['chrom1'] = get_chrom(row['chrom1'], species)
	row['chrom2'] = get_chrom(row['chrom2'], species)	
	return row

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

def node_size_to_radius(node_size):
    """Convert node_size (area in points²) to node_radius (radius in data units)."""
    return np.sqrt(node_size / np.pi)

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
df = df.apply(rename_chroms, axis=1)
df = df[df['chrom1'] != df['chrom2']]
df = df.drop_duplicates(subset=['species', 'chrom1', 'chrom2'], keep='first')

mean_bp = np.mean(df['bp'])
print(mean_bp)
mean_pid = np.mean(df['pid'])
df['normalized_bp'] = df['bp'] / 250000 # mean_bp
df['normalized_pid'] = df['pid'] - 99
# norm = plt.Normalize(vmin=df['normalized_pid'].min(), vmax=df['normalized_pid'].max())
cmap = cm.get_cmap('YlOrRd')
# df['edge_color'] = df['normalized_pid'].apply(lambda pid: cmap(norm(pid)))
df['edge_color'] = df['normalized_pid'].apply(lambda pid: cmap(pid))
df['edge_weight'] = df['normalized_bp']
df_subset = df[(df['species'] == species) & (df['chrom1'] != df['chrom2'])]
G = nx.Graph()
for _, row in df_subset.iterrows():
	chrom1 = row['chrom1']
	chrom2 = row['chrom2']
	normalized_bp = df_subset['normalized_bp']
	normalized_pid = df_subset['normalized_pid']
	weight = row['edge_weight']
	color = row['edge_color']
	G.add_edge(chrom1, chrom2, norm_bp=normalized_bp, norm_pid=normalized_pid, weight=weight, color=color)

nor_plus_acros, nor_minus_acros, ancient_acros, custom_order = get_order(species)

for node in custom_order:
	if node not in G.nodes:
		G.add_node(node)

if species == 'HG002':
	G.add_edge('X', 'Y', norm_bp=3290000/250000, norm_pid=0.6219, weight=3290000/250000, color=cmap(0.6219))
	G.add_edge('Y', 'X', norm_bp=3290000/250000, norm_pid=0.6219, weight=3290000/250000, color=cmap(0.6219))


pos = compute_circular_positions(custom_order, 0.5)
edge_colors = [G[u][v]['color'] for u, v in G.edges()]
edge_weights = [G[u][v]['weight'] for u, v in G.edges()]

node_colors = []
# for node in G.nodes():
# 	if node in nor_plus_acros:
# 		node_colors.append('#289ADA') #274C77
# 	elif node in nor_minus_acros:
# 		node_colors.append('#009e73') #6096BA
# 	elif node in ancient_acros:
# 		node_colors.append('#cc79a7') #A3CEF1
# 	else:
# 		node_colors.append('#C0C0C0') #E7ECEF

# for node in G.nodes():
# 	if node in nor_plus_acros:
# 		node_colors.append('#054A91')
# 	elif node in nor_minus_acros:
# 		node_colors.append('#3E7CB1')
# 	elif node in ancient_acros:
# 		node_colors.append('#81A4CD')
# 	else:
# 		node_colors.append('#DBE4EE')

# for node in G.nodes():
# 	if node in nor_plus_acros:
# 		node_colors.append('#ed2224')
# 	elif node in nor_minus_acros:
# 		node_colors.append('#059e74')
# 	elif node in ancient_acros:
# 		node_colors.append('#5bb4e5')
# 	else:
# 		node_colors.append('#C0C0C0')

# for node in G.nodes():
# 	if node in nor_plus_acros:
# 		node_colors.append('#117733')
# 	elif node in nor_minus_acros:
# 		node_colors.append('#0084c6')
# 	elif node in ancient_acros:
# 		node_colors.append('#DDCC77')
# 	else:
# 		node_colors.append('#dcdcdc')

# for node in G.nodes():
# 	if node in nor_plus_acros:
# 		node_colors.append('#012a4a')
# 	elif node in nor_minus_acros:
# 		node_colors.append('#014f86')
# 	elif node in ancient_acros:
# 		node_colors.append('#468faf')
# 	else:
# 		node_colors.append('#a9d6e5')

for node in G.nodes():
	if node in nor_plus_acros:
		node_colors.append('#003863')
	elif node in nor_minus_acros:
		node_colors.append('#007BD1')
	elif node in ancient_acros:
		node_colors.append('#64D6FF')
	else:
		node_colors.append('#C8C8C8')

# node_label_colors = {}
# for node in G.nodes():
#     if node in nor_plus_acros:
#         node_label_colors[node] = '#FFFFFF'
#     else:
#         node_label_colors[node] = '#000000'

node_labels = {node: node for node in G.nodes()}
font_colors = {node: '#D0D0D0' if node in nor_plus_acros else 'black' for node in G.nodes()}

plt.figure(figsize=(10, 10))
nx.draw_networkx_nodes(G, pos, node_color=node_colors, edgecolors=node_colors, node_size=750)
print(G.edges())

# node_radius = node_size_to_radius(750)
# i=0
# for start, end in G.edges():
#     p_start, p_end = rescale_edge_positions(pos, start, end, radius=0.08)
#     plt.plot([p_start[0], p_end[0]], [p_start[1], p_end[1]], linewidth=edge_weights[i], color=edge_colors[i])
#     i+=1

nx.draw_networkx_edges(G, pos, width=edge_weights, edge_color=edge_colors)
# nx.draw_networkx_labels(G, pos, font_size=12, font_color='black')
# nx.draw_networkx_labels(
#     G, pos, 
#     labels={node: node for node in G.nodes()}, 
#     font_size=12, 
#     font_color=[node_label_colors[node] for node in G.nodes()]
# )
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

# df_subset.to_csv(sys.stdout, sep='\t', index=False)
