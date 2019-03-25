from collections import defaultdict
from itertools import product
import json

import click
import numpy as np
from sklearn.metrics import homogeneity_completeness_v_measure as hcv

from nanotext.evaluate import odd_one
from nanotext.io import eprint, load_embedding
from nanotext.utils import cosine


@click.command()
@click.option(
    '--model', '-m', 
    help='The nanotext model to evaluate', 
    type=click.Path(), required=True)
@click.option(
    '--corpus', '-i', 
    help='A subset of the corpus set aside for testing',
    type=click.Path(), required=True)
@click.option('--clusters',
    default=None, help='Use this file of ecotypes for validation',
    type=click.Path(), required=True)
@click.option(
    '--outfile', '-o', 
    help='Name of the resulting evaluation report', 
    type=click.Path(), required=True)
def evaluate(model, corpus, outfile, clusters):
    '''
    A test battery:
    
    - SOMO
    - king queen
    - various clusterings w/ associated tables
        - e coli
        - closridia
        - chlamydia
        - prochlorococcus
    
    closest genomes prochlorococcus or Tara
    '''
    results = {}


    # SOMO task -- tests word embedding quality
    eprint('SOMO task ...')
    SOMO = odd_one(corpus, model, n_not_odd=5)
    results['SOMO'] = SOMO


    # Ecotype task -- tests document embedding quality
    eprint('Ecotype task ...')
    truth, cluster_labels = [], []
    d = defaultdict(list)

    with open(clusters, 'r') as file:
        _ = next(file)  # header
        for line in file:
            row = line.strip().split('\t')
            clade = row[3] 
            cluster = row[6]
            if (clade in ['HL', 'LL']) and (cluster != -1):
                d[clade].append(row[2])  # genome UID
                truth.append(clade)
                cluster_labels.append(cluster)
    h, c, v = [round(i, 4) for i in hcv(truth, cluster_labels)]
    results.update(
        dict(zip('homogeneity completeness vscore'.split(), [h, c, v])))


    # Distance
    eprint('Median cosine distance btw/ points of different ecotypes ...')
    model = load_embedding(model)
    l = []
    for i, j in product(d['HL'], d['LL']):
        try:
            a = model.docvecs[i]
            b = model.docvecs[j]
            l.append(cosine(a, b))
        except KeyError:
            continue

    results['ecotypes_distance'] = round(
        float(1-np.median(l)), 4)


    with open(outfile, 'w+') as out:
        json.dump(results, out, indent=4)


if __name__ == '__main__':
    evaluate()
