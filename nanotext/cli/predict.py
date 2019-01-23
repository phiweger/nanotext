import click


def cosine_similarity(a, b):
    '''
    https://en.wikipedia.org/wiki/Cosine_similarity
    stackoverflow.com/questions/18424228
    '''
    from numpy import dot
    from numpy.linalg import norm

    return dot(a, b)/(norm(a)*norm(b))


def cosim2(query, db, topn=10):
    '''
    cosim2(vv_media['1'], vv_media, 10)
    '''
    import operator 

    result = {}
    for k, v in db.items():
        result[k] = cosine_similarity(v, query)


    # sort in descending order
    sorted_result = sorted(
        result.items(), key=operator.itemgetter(1), reverse=True)
    return sorted_result[:topn]


@click.command()
@click.option(
    '--genome',
    help='Protein domain annotation',
    type=click.Path())
@click.option(
    '--embedding',
    help='Genome embedding model')
@click.option(
    '--db',
    help='Media vectors (i.e. the sum across ingredient word vectors).')
@click.option(
    '--model', 
    help='Predictive model', 
    type=click.Path())
@click.option(
    '--topn', 
    help='Top n hits to return.',
    default=10)
@click.option(
    '--out',
    help='Output path. If not specified, write to stdout.',
    default='-')
def predict(genome, embedding, db, model, out, topn):
    '''
    From a <genome> w/ annotated protein domains predict a phenotype. Requires
    the learned <model> (genotype-phenotype mapping) as well as a genome
    <embedding>. Return the closest <topn> vectors from a database <db>.

    Usage:

    \b
    nanotext predict \\
        --out - \\
        --model data/media_prediction.h5 \\
        --db data/embedding.media.json \\
        --embedding data/embedding.genomes.model \\
        --genome data/TARA_ION_MAG_00012.domtbl.tsv \\
        --topn 3
    '''
    import json
    import os

    import numpy as np
    import tensorflow as tf
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  # turn off debugging info
    from keras.models import load_model

    from nanotext.io import load_embedding, smart_open, eprint
    from nanotext.utils import infer_genome_vector


    eprint('Loading embedding model for genomes ...')
    e = load_embedding(embedding)
    eprint('Inferring genome vector ...')
    v = infer_genome_vector(genome, e)


    eprint('Loading media vector database ...')
    with open(db, 'r') as file:
        vv = json.load(file)  # vv .. vectors

    eprint('Loading predictive model ...')
    nn = load_model(model)

    y_hat = nn.predict(np.array([v]))[0]  # [0] .. only single genome for now
    sim = cosim2(y_hat, vv, topn)
    with smart_open(out) as fh:
        # fh.write('\nmedium\tcosine\n')
        for name, cos in sim:
            fh.write(f'{name}\t{round(cos, 4)}\n')
    eprint('Done.')


