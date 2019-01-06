import click

@click.command()
@click.option(
    '--model',
    help='Params for training the embedding.')
@click.option(
    '--topn', 
    help='Corpus of protein domains to train on.', 
    type=click.Path())
@click.option(
    '--out',
    help='Model name (output).')
def predict(genome, topn, out):
    '''
    Usage:

    nanotext predict --genome x --model medium.model --out -
    nanotext predict --genome x --model water.model --out -
    '''
    print('pred')
    # predict medium