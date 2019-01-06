import click

# TODO: nanotext index
# https://github.com/facebookresearch/faiss
# https://github.com/facebookresearch/faiss/wiki/Index-IO,-index-factory,-cloning-and-hyper-parameter-tuning#io-and-deep-copying-indexes


@click.command()
@click.option(
    '--genome',
    help='Params for training the embedding.',
    type=click.Path())
@click.option(
    '--topn', 
    help='Corpus of protein domains to train on.',
    default=10)
@click.option(
    '--embedding',
    help='Model name (output).',
    required=True, type=click.Path())
@click.option(
    '--out',
    help='Output path. If not specifies, write to stdout.',
    default='-')
def search(genome, topn, embedding, out):
    '''
    Usage:

    \b
    nanotext search \\
        --embedding embedding.genomes.model --topn 3 --out - \\
        --genome .../tara/annotation/TARA_ION_MAG_00012/orfs.domtbl.tsv
    # Loading model ...
    # Inferring genome vector ...
    # GCA_000634215.1 0.93
    # GCF_000759935.1 0.93
    # GCF_000759855.1 0.93
    # Done.
    '''
    from nanotext.io import load_embedding, smart_open, eprint
    from nanotext.utils import infer_genome_vector

    eprint('Loading model ...')
    model = load_embedding(embedding)
    eprint('Inferring genome vector ...')
    v = infer_genome_vector(genome, model)
    sim = model.docvecs.most_similar([v], topn=topn)

    with smart_open(out) as fh:
        for name, cos in sim:
            fh.write(f'{name}\t{round(cos, 2)}\n')
    eprint('Done.')


@click.command()
@click.option(
    '--genome',
    help='Params for training the embedding.')
@click.option(
    '--topn', 
    help='Corpus of protein domains to train on.', 
    type=click.Path())
@click.option(
    '--out',
    help='Model name (output).')
def compare(genome, topn, out):
    print('foo')
    # cosine distance btw/ 2 genomes


@click.command()
@click.option(
    '--genome',
    help='Params for training the embedding.')
@click.option(
    '--topn', 
    help='Corpus of protein domains to train on.', 
    type=click.Path())
@click.option(
    '--out',
    help='Model name (output).')
def taxonomy(genome, topn, out):
    print('cello')
    # get similar genomes and their taxon -- try assign taxon


@click.command()
@click.option(
    '--genome',
    help='Params for training the embedding.')
@click.option(
    '--topn', 
    help='Corpus of protein domains to train on.', 
    type=click.Path())
@click.option(
    '--out',
    help='Model name (output).')
def lookup(genome, topn, out):
    print('hello')
    # given 16S return genome vector
