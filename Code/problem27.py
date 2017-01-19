#!/usr/bin/env python2
# Enumerating Unrooted Binary Trees
from newick_format import *

class Node():
    def __init__(self, name):
        self.name = name

    def __str__(self):
        if self.name is not None:
            return self.name
        else:
            return "{}".format(id(self))

class Edge():
    def __init__(self, node1, node2):
        self.nodes = [node1, node2]

    def __str__(self):
        return "{}{}".format(*self.nodes)

class Tree():
    def __init__(self, nodes=[], edges=[]):
        self.nodes = nodes
        self.edges = edges

    def __str__(self):
        return "tree{} edges:{}".format(id(self), [str(x) for x in self.edges])

    def copy(self):
        node_conversion = {node: Node(node.name) for node in self.nodes}
        new_nodes = list(node_conversion.values())
        new_edges = [Edge(node_conversion[edge.nodes[0]], node_conversion[edge.nodes[1]]) for edge in self.edges]

        new_tree = Tree(new_nodes, new_edges)
        return new_tree


def all_trees(leaves):
    assert(len(leaves) > 1)
    
    if len(leaves) == 2:
        n1, n2 = leaves
        t = Tree()
        t.nodes = [Node(n1), Node(n2)]
        t.edges = [Edge(t.nodes[0], t.nodes[1])]
        return [t]
    elif len(leaves) > 2:
        # get the smaller tree first
        ot = all_trees(leaves[:-1])
        nln = leaves[-1]
        nt = []

        # find the ways to add the new leaf
        for old_tree in ot:
            for i in range(len(old_tree.edges)):
                new_tree = old_tree.copy()
                edge_to_split = new_tree.edges[i]
                old_node1, old_node2 = edge_to_split.nodes

                # get rid of the old edge
                new_tree.edges.remove(edge_to_split)

                # add a new internal node
                internal = Node(None)
                new_tree.nodes.append(internal)

                # add the new leaf
                new_leaf = Node(nln)
                new_tree.nodes.append(new_leaf)

                # make the three new edges
                new_tree.edges.append(Edge(old_node1, internal))
                new_tree.edges.append(Edge(old_node2, internal))
                new_tree.edges.append(Edge(new_leaf, internal))

                # put this new tree in the list
                nt.append(new_tree) 

        return nt



if __name__ == '__main__':
    data = open('/Users/krishna/Downloads/rosalind_eubt.txt').read().split()
    trees = all_trees(data)

    print '\n'.join([newick_format(tree) for tree in trees])
