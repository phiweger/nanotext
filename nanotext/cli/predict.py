import click

@click.command()
@click.option(
    '--genome',
    help='...',
    type=click.Path())
@click.option(
    '--embedding',
    help='...')
@click.option(
    '--embedding_target',
    help='...')
@click.option(
    '--model', 
    help='Predictive model.', 
    type=click.Path())
@click.option(
    '--topn', 
    help='...',
    default=10)
@click.option(
    '--out',
    help='Output path. If not specifies, write to stdout.')
def predict(genome, embedding, model, out, topn, embedding_target):
    '''
    Usage:

    nanotext predict --genome x --model medium.model --out -
    nanotext predict --genome x --model water.model --out -
    '''
    import os

    import tensorflow as tf
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  # turn off debugging info
    from keras.models import load_model

    from nanotext.io import load_embedding, smart_open, eprint
    from nanotext.utils import infer_genome_vector

        
    eprint('Loading embedding model ...')
    e = load_embedding(embedding)

    eprint('Loading predictive model ...')
    nn = load_model(model)

    if embedding_target:
        # et = load_embedding(embedding_target)
        
        from gensim.models.doc2vec import Doc2Vec
        print(embedding_target)
        et = Doc2Vec.load(embedding_target)
        
        print(et.docvecs.most_similar('1'))
        v = infer_genome_vector(genome, e)
        pred = get_target([v], nn, et, topn=topn)
        print(pred)
        
        # for k, v in pred:
        #     print(k)
        #     print(et.docvecs.most_similar(k))


def get_target(vv, nn, model, topn=10):
    '''
    vv .. a list of vectors
    '''
    import numpy as np

    y_hat = nn.predict(np.array(vv))
    return model.docvecs.most_similar(y_hat, topn=topn)



