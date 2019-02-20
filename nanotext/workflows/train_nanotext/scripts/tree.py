from argparse import Namespace

import ete3
import hdbscan
import numpy as np
from sklearn.preprocessing import normalize
import umap

from nanotext.io import load_embedding, eprint


args = Namespace(model='nanotext_r89.model')


config_umap_dim_reduction = {
    'metric': 'cosine',
    'n_components': 10,
    # 'n_neighbors': 10,
    # 'min_dist': 0.05,
    # 'spread': 5,
    'random_state': 42,
    }

config_hdbscan = {
    'min_cluster_size': 7,
    'min_samples': 1,
    'cluster_selection_method': 'eom',
    }


model = load_embedding(args.model)
names = model.docvecs.index2entity


# https://github.com/scikit-learn-contrib/hdbscan
# because t-SNE for clustering is no good
# stats.stackexchange.com/questions/308132
m = []
for i in names:  # names .. just a list of accession IDs
    m.append(model.docvecs[i])
m = np.array(m, dtype='float64')


eprint('Projecting points (dimension reduction) ...')
# Reduce dimensions before clustering
reducer = umap.UMAP(**config_umap_dim_reduction)
m_redux = reducer.fit_transform(m)
m_norm = normalize(m_redux, norm='l2', axis=1)
clusterer = hdbscan.HDBSCAN(**config_hdbscan)  
cluster_labels = clusterer.fit_predict(m_norm)


# We now want to extract the clustering hierarchy from the clusterer object
# (i.e. the tree).
# condensed tree to newick
# stackoverflow.com/questions/46444454
g = clusterer.condensed_tree_.to_networkx()

# Find root. The root of a tree has an indegree of 0.
# stackoverflow.com/questions/4122390
root = [n for n, d in g.in_degree() if d==0][0]  # there's only one
subtrees = {node: ete3.Tree(name=node) for node in g.nodes()}

# Combine nodes in tree and add edge weight information.
# add_child(dist=...)
# g[2511][45]  # {'weight': 0.4547185572001402}
# to newick tree
# stackoverflow.com/questions/51273890
[*map(
    lambda edge:subtrees[edge[0]].add_child(
        subtrees[edge[1]], dist=g[edge[0]][edge[1]]['weight']), g.edges())]
tree = subtrees[root]

with open('tree.newick', 'w+') as out:
    out.write(tree.write(format=5))

# midpoint root?
# http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#getting-midpoint-outgroup
R = tree.get_midpoint_outgroup()
# and set it as tree outgroup
tree.set_outgroup(R)


with open('tree_midpoint.newick', 'w+') as out:
    out.write(tree.write(format=5))
# we can export in multiple formats
# the example in ggtree ("Tree Manipulation") is:
# (((((((A:4,B:4):6,C:5):8,D:6):3,E:21):10,((F:4,G:12):14,H:8):13):13,((I:5,J:2):30,(K:11,L:11):2):17):4,M:56);
# this corresponds to format=5
# see "Reading and Writing Newick Trees"
# etetoolkit.org/docs/latest/tutorial/tutorial_trees.html


# clusterer.single_linkage_tree_.to_pandas()
# TODO: validity?




'''r
library(ggtree)

library(ggtree)
tree <- read.tree('tree.newick')
p <- ggtree(tree, size=0.2, color='grey90') %<+% anno + 
    geom_tippoint(aes(color=as.factor(foo)), size=0.7) +
    theme(legend.position='right') +
    scale_x_reverse() + 
    coord_flip()

ggsave('tree.pdf', p, width=25, height=12, units='cm')

'''

