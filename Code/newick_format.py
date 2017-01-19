def newick_format(tree_in):
    tree = tree_in.copy()

    if len(tree.nodes) == 1:
        return "{};".format(tree.nodes[0])
    elif len(tree.nodes) == 2:
        return "({},{});".format(*tree.nodes)
    elif len(tree.nodes) > 2:
        # reduce one of the nodes in the tree
        for candidate_node in tree.nodes:
            # ignore leaves
            if candidate_node.name is not None:
                continue

            adjacent_edges = [edge for edge in tree.edges if candidate_node in edge.nodes]
            adjacent_nodes = [node for edge in adjacent_edges for node in edge.nodes if node in edge.nodes and node is not candidate_node]
            adjacent_leaves = [node for node in adjacent_nodes if node.name is not None]

            # find a node with two leaves
            if len(adjacent_leaves) == 2 or len(adjacent_leaves) == 3:
                leaf1, leaf2 = adjacent_leaves[0: 2]
                edges_to_cut = [edge for edge in adjacent_edges if leaf1 in edge.nodes or leaf2 in edge.nodes]
                candidate_node.name = "({},{})".format(leaf1, leaf2)

                # remove leaves
                tree.nodes.remove(leaf1)
                tree.nodes.remove(leaf2)
                for edge in edges_to_cut: tree.edges.remove(edge)

                # cycle this one again
                return newick_format(tree)