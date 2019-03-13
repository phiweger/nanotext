from argparse import Namespace
from collections import Counter, defaultdict

import click

from nanotext.io import load_taxonomy_gtdb


@click.command()
@click.option('--model_name', '--model_name',
    required=True, help='nanotext model iteration')
@click.option('--clusters',
    type=click.Path(), required=True,
    help='Cluster file')
@click.option('--taxonomy',
    help='Path to taxonomy summary from GTDB',
    required=True, type=click.Path())
@click.option('--outfile', '-o',
    required=True, type=click.Path(), help='File to write results to')
def homogeneity(outfile, model_name, taxonomy, clusters):
    '''
    Take the clusters from cluster.py and check how homogenous their taxa are
    at all ranks.
    '''
    # args = Namespace(
    #     outfile='homogeneity.tsv',
    #     model_name='19',
    #     taxonomy='/Users/phi/data_local/databases/gtdb/bac_taxonomy_r86.tsv',
    #     clusters='clusters.tsv')

    args = Namespace(
        outfile=outfile,
        model_name=model_name,
        taxonomy=taxonomy,
        clusters=clusters)

    db = load_taxonomy_gtdb(args.taxonomy)
    ranks = [
        'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    
    
    clusters = defaultdict(list)
    with open(args.clusters, 'r') as file:
        _ = next(file)  # header
        for line in file:
            c1, c2, name, cluster = line.strip().split('\t')
            clusters[cluster].append(db[name])
    
    # model cluster rank ratio
    with open(args.outfile, 'w+') as out:
        # out.write('cluster\trank\tratio\n')
    
        for k, v in clusters.items():
            if k != '-1':  # the "no cluster" cluster in HDBSCAN
                d = defaultdict(list)  # a dict for each cluster
                for i in v:
                    for rank, taxon in zip(ranks, i):
                        d[rank].append(taxon)
        
                for rank, taxa in d.items():
                    # at least one taxon is specified, i.e. not empty ''
                    if any(taxa):
                        cnt = Counter(taxa)
                        del cnt['']  # remove unknown taxa?
                        maxn = max(cnt.values())
                        # hits = [k for k, v in cnt.items() if v == maxn]
                        # pick = random.choice(hits)
                        # majority[rank] = pick
                        ratio = round(maxn/len(taxa), 2)
                        out.write(f'{k}\t{rank}\t{ratio}\t{args.model_name}\n')
                    else:
                        ratio = 0
                        out.write(f'{k}\t{rank}\t{ratio}\t{args.model_name}\n')


if __name__ == '__main__':
    homogeneity()