


Networkx MultiGraph API: (Adjacency dict is the main data structure)
G.adj = G.edge  # Adjacency dict with {node1: {node2: {edge_key: {edge attr dict}}}}
G.edge[node1][node2] = G[node1][node2] = {edge_key: {edge attr dict}}
G.node[node1] = {dict with node attr}
G.edges(nbunch, keys=True) => [list with (source, target, key) tuples], uses G.adj
Note: We really do not expect there to be more than one edge from source to target,
in fact that would be an error, so it might be better to use a normal DiGraph rather than a MultiDiGraph?
