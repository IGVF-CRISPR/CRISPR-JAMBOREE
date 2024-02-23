import os
import pandas as pd
import numpy as np
import muon as mu
from mudata import MuData
import networkx as nx
import matplotlib.pyplot as plt
from typing import List, Tuple, Dict, Any, Optional


def plot_network(
    mdata: MuData,
    central_node: Optional[str] = None,
    source_column: str = "source",
    target_column: str = "target",
    weight_column: Optional[str] = None,
    min_weight: Optional[float] = None,
    node_size_column: Optional[str] = None,
    results_key: Optional[str] = "test_results",
):
    # Get and filter the results dataframe
    results_df = mdata.uns[results_key]
    if min_weight is not None:
        results_df = results_df[results_df[weight_column].abs() >= min_weight]
    if central_node is not None:
        results_df = results_df[results_df[source_column] == central_node]

    # Create the network 
    G = nx.DiGraph()
    for i, row in results_df.iterrows():
        G.add_edge(row[source_column], row[target_column], weight=row[weight_column])
    pos = nx.circular_layout(G)

    # Draw nodes based on node size column if provided
    if node_size_column is not None:
        node_size = results_df.set_index(target_column)[node_size_column].to_dict()
        sizes = [node_size.get(n, 10) for n in G.nodes]
        sizes = [s * 100 for s in sizes]
        nx.draw(G, pos, node_color="skyblue", node_size=sizes, edge_cmap=plt.cm.Blues, arrowsize=20)
    else:
        nx.draw(G, pos, with_labels=True, node_color="skyblue", node_size=150, edge_cmap=plt.cm.Blues, arrowsize=20)

    # Offset labels
    if central_node is not None:
        label_pos = {k: (v[0]+0.2, v[1] + 0.05) for k, v in pos.items() if k != central_node}  # Adjust 0.1 as needed for the offset
        label_pos[central_node] = pos[central_node]
        nx.draw_networkx_labels(G, label_pos, font_size=10)
    else:
        label_pos = {k: (v[0]+0.2, v[1] + 0.05) for k, v in pos.items()}
        nx.draw_networkx_labels(G, pos, font_size=8)
    
    # Draw edges
    edge_labels = {(u, v): f"{d['weight']:.2f}" for u, v, d in G.edges(data=True)}
    nx.draw_networkx_edges(G, pos, edge_color=[d["weight"] for u, v, d in G.edges(data=True)], edge_cmap=plt.cm.coolwarm, arrowsize=5)
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)
    
    # Plt
    plt.tight_layout()
    plt.show()
    