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
	angle = 2 * math.pi / len(nodes)
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

df['normalized_bp'] = df['bp'] / 250000 # mean_bp
df['normalized_pid'] = df['pid'] - 99
cmap = cm.get_cmap('YlOrRd')
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
nx.draw_networkx_nodes(G, pos, node_color=node_colors, edgecolors=node_colors, node_size=750)
print(G.edges())
nx.draw_networkx_edges(G, pos, width=edge_weights, edge_color=edge_colors)
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
plt.savefig(sys.argv[3])
plt.clf()
