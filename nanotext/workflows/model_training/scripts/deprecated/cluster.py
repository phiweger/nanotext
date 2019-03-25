from argparse import Namespace
import csv
import warnings
import sys

import click
import numpy as np
# with warnings.catch_warnings():
#     '''
#     DeprecationWarning: the imp module is deprecated in favour of importlib;
#     see the module's documentation for alternative uses
#     import imp

#     stackoverflow.com/questions/879173
#     '''
#     warnings.filterwarnings('ignore', category=DeprecationWarning)
import hdbscan
from sklearn.preprocessing import normalize
from sklearn.manifold import TSNE
import umap

from nanotext.io import load_embedding, eprint
from nanotext.utils import strip_name


# 'Clostridia', 'Gammaproteobacteria', 'Chlamydiia', 'Oxyphotobacteria'
@click.command()
@click.option('--model', '-m',
    required=True, type=click.Path(), help='nanotext model')
@click.option('--rank',
    default='class', help='At which rank to cluster')
@click.option('--projection_method',
    default='umap',
    help='Which method to use to project high dim vectors into 2 dim')
@click.option('--name',
    required=True, help='Specific name of the rank (e.g. Clostridia)')
@click.option('--taxonomy',
    help='Path to taxonomy summary from GTDB',
    required=True, type=click.Path())
@click.option('--outfile', '-o',
    required=True, type=click.Path(), help='File to write results to')
@click.option('--soft',
    default=False, help='Apply HDBSCAN soft clustering option')
@click.option('--ecotypes',
    help='Use this file of ecotypes for validation',
    default=None, type=click.Path())
def cluster_subset(model, rank, name, taxonomy, outfile, soft, ecotypes, projection_method):
    '''
    TODO: https://github.com/lmcinnes/umap/issues/90

    Iterate over a taxonomic rank such as class and cluster using HDBSCAN. One 
    reason we believe we can do this is that at higher ranks there are clear
    boundaries between organisms. The main motivation behind it is that
    HDBSCAN clusters are rather coarse when the whole dataset is clustered at 
    once, and "soft-clustering" does not scale.


    python ~/Dropbox/repos_git/nanotext/nanotext/workflows/train_nanotext/scripts/cluster.py -m nanotext_r89.model --name Clostridia --taxonomy /Users/phi/data_local/databases/gtdb/bac_taxonomy_r83.tsv -o clusters.tsv
    '''
    # args = Namespace(
    #     model='nanotext_r89.model',
    #     rank='class',
    #     name='Oxyphotobacteria',
    #     taxonomy='/Users/phi/data_local/databases/gtdb/bac_taxonomy_r83.tsv',
    #     outfile='clusters.tsv',
    #     soft=False,
    #     )

    args = Namespace(
        model=model,
        rank=rank,
        name=name,
        taxonomy=taxonomy,
        outfile=outfile,
        soft=soft,
        ecotypes=ecotypes
        )

    '''
    umap.UMAP(
    ['n_neighbors=15', 'n_components=2', "metric='euclidean'", 'n_epochs=None', 'learning_rate=1.0', "init='spectral'", 'min_dist=0.1', 'spread=1.0', 'set_op_mix_ratio=1.0', 'local_connectivity=1.0', 'repulsion_strength=1.0', 'negative_sample_rate=5', 'transform_queue_size=4.0', 'a=None', 'b=None', 'random_state=None', 'metric_kwds=None', 'angular_rp_forest=False', 'target_n_neighbors=-1', "target_metric='categorical'", 'target_metric_kwds=None', 'target_weight=0.5', 'transform_seed=42', 'verbose=False'],)
    '''
    config_umap_visualisation = {
        'metric': 'cosine',
        'n_components': 2,
        # 'n_neighbors': 5,
        # min_dist=0.05,
        # 'spread': 5,
        'random_state': 42,
        }

    config_umap_dim_reduction = {
        'metric': 'cosine',
        'n_components': 10,
        # 'n_neighbors': 10,
        # 'min_dist': 0.05,
        # 'spread': 5,
        'random_state': 42,
        }
    
    # min_cluster_size 3 leaf works great
    config_hdbscan = {
        'min_cluster_size': 5,
        # 'min_samples': 1,
        'cluster_selection_method': 'eom',
        }


    # Filter taxonomy file for a given rank
    names = []
    with open(args.taxonomy, 'r') as file:
        for line in file:
            if f'{args.rank[0]}__{args.name}' in line:  # d__ for domain etc.
                names.append(line.strip().split('\t')[0])
    eprint(f'There are {len(names)} data points for {args.rank} {args.name}.')
    
    
    # Extract only those vectors
    model = load_embedding(args.model)
    m, found = [], []
    for i in names:
        try:
            m.append(model.docvecs[strip_name(i)])
            found.append(strip_name(i))
        except KeyError:
            continue

    
    if args.ecotypes:
        ecotypes = {}
        with open(args.ecotypes) as csvfile:
            _ = next(csvfile)
            reader = csv.reader(csvfile, delimiter='\t')
            for row in reader:
                genome, e, curator = row[0], row[-2], row[-1]
                if '_' in e:
                    subtype = e
                    e = e.split('_')[0]
                else:
                    subtype = 'NA'
                ecotypes[genome] = (e, subtype, curator)

        # Extend sample list
        # TODO: this will be unnecessary once we have the r89 tax
        for i in ecotypes.keys():
            try:
                m.append(model.docvecs[i])
                found.append(i)
            except KeyError:
                continue


    m = np.array(m, dtype='float64')
    ratio = int(round(m.shape[0]/len(names), 2)*100)
    eprint(f'Of those, {m.shape[0]} ({ratio}%) are present in the model.')        

    pm = projection_method.upper()
    eprint(f'Projecting points (visualisation) using {pm} ...')
    if projection_method == 'tsne':
        projection = TSNE(n_components=2, random_state=42).fit_transform(m)
    elif projection_method == 'umap':
        reducer = umap.UMAP(**config_umap_visualisation)
        projection = reducer.fit_transform(m)
        # projection[-10:] == reducer.embedding_[-10:]
    else:
        eprint('No valid projection method. Abort!')
        sys.exit(-1)

    
    eprint('Projecting points (dimension reduction) ...')
    # Reduce dimensions before clustering
    reducer = umap.UMAP(**config_umap_dim_reduction)
    m_redux = reducer.fit_transform(m)


    eprint('Clustering ...')
    m_norm = normalize(m_redux, norm='l2', axis=1)
    # prediction_data=True, eom/ leaf
    clusterer = hdbscan.HDBSCAN(**config_hdbscan)
    cluster_labels = clusterer.fit_predict(projection)
    # Or soft clustering: init clusterer w/ prediction_data=True
    # soft_clusters = hdbscan.all_points_membership_vectors(clusterer)
    # cluster_labels = [np.argmax(x) for x in soft_clusters]
    
    
    if args.ecotypes:
        with open(args.outfile, 'w+') as out:
            out.write('c1\tc2\tname\tclade\tsubclade\tcollection\tcluster\n')
        
            for i, j, k in zip(projection, found, cluster_labels):
                c1, c2 = i
                clade, subclade, curator = ecotypes.get(j, 2*['NA']+['GTDB'])
                out.write(
                    f'{c1}\t{c2}\t{j}\t{clade}\t{subclade}\t{curator}\t{k}\n')
    else:
        with open(args.outfile, 'w+') as out:
            out.write('c1\tc2\tname\tcluster\n')
        
            for i, j, k in zip(projection, found, cluster_labels):
                c1, c2 = i
                out.write(f'{c1}\t{c2}\t{j}\t{k}\n')

    eprint('Done.')


if __name__ == '__main__':
    cluster_subset()


'''
library(ggplot2)
library(readr)


anno <- read_tsv('clusters.tsv')
anno$cluster <- as.factor(anno$cluster)
size <- 0.2


p <- ggplot() + 
    geom_point(
        data=subset(anno, cluster==-1),  # -1
        aes(x=c1, y=c2), size=size, color='grey90') +
    geom_point(
        data=subset(anno, cluster!=-1), 
        aes(x=c1, y=c2, color=cluster), size=size) + 
    # geom_point(
    #     data=subset(anno, cluster==COI),  # -1
    #     aes(x=c1, y=c2), size=2*size, color='red') +
    # guides(color=FALSE) +
    theme_classic() +
    # scale_x_continuous(limits=c(10, 25)) +
    # scale_y_continuous(limits=c(-5, 2)) +
    # scale_color_brewer(palette='Set2') +
    theme(axis.line=element_blank())


# ecotypes
anno <- read_tsv('clusters.tsv')
size <- 1
q <- ggplot() + 
    geom_point(
        data=subset(anno, is.na(anno$clade)),  # -1
        aes(x=c1, y=c2), size=size, color='grey90') +
    geom_point(
        data=subset(anno, !is.na(anno$clade)), 
        aes(x=c1, y=c2, color=clade), size=size) + 
    # geom_point(
    #     data=subset(anno, cluster==COI),  # -1
    #     aes(x=c1, y=c2), size=2*size, color='red') +
    # guides(color=FALSE) +
    theme_classic() +
    scale_x_continuous(limits=c(10, 25)) +
    scale_y_continuous(limits=c(-5, 2)) +
    scale_color_brewer(palette='Set2') +
    theme(axis.line=element_blank())
'''