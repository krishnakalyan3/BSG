#!/usr/bin/env python2
# Distances in Trees
import newick as nk

def newick_parse(tree, start, end):
    tree = nk.parse_newick(tree)
    tree, labels = tree
    kx = labels[start]
    ky = labels[end]
    kx_depths = dict(kx.ancestors())
    ky_depths = dict(ky.ancestors())
    common = set([x[0] for x in kx.ancestors()]) & set([y[0] for y in ky.ancestors()])
    best = sorted(common, key=lambda k: kx_depths[k] + ky_depths[k])[0]
    return kx_depths[best] + ky_depths[best]

def parse_data(data):
	results = []
	for d in data:
		data_split = d.split('\n')
		s,e = data_split[1].split(' ')
		result = newick_parse(data_split[0],s, e)
		results.append(str(result))
	return results

def problem(tree, ks):
    tree, labels = tree
    kx = labels[ks[0]]
    ky = labels[ks[1]]
    kx_depths = dict(kx.ancestors())
    ky_depths = dict(ky.ancestors())
    common = set([x[0] for x in kx.ancestors()]) & set([y[0] for y in ky.ancestors()])
    best = sorted(common, key=lambda k: kx_depths[k] + ky_depths[k])[0]
    return kx_depths[best] + ky_depths[best]

if __name__ == '__main__':
	dataset = open("/Users/krishna/Downloads/rosalind_nwck.txt").read()
	#dataset = open("/Users/krishna/MIRI/BSG/data/25").read()
	data_tree = dataset.split('\n\n')
	tree_rep = parse_data(data_tree)
	print ' '.join(tree_rep)
	