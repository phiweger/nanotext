import click

# TODO: nanotext index
# https://github.com/facebookresearch/faiss
# https://github.com/facebookresearch/faiss/wiki/Index-IO,-index-factory,-cloning-and-hyper-parameter-tuning#io-and-deep-copying-indexes


@click.command()
@click.option(
    '--genome',
    help='Protein domain annotation',
    type=click.Path())
@click.option(
    '--topn', 
    help='Top n hits to return',
    default=10)
@click.option(
    '--embedding',
    help='Genome embedding model',
    required=True, type=click.Path())
@click.option(
    '--out',
    help='Output path. If not specified, write to stdout.',
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
    # GCA_000634215.1 0.9344
    # GCF_000759935.1 0.9282
    # GCF_000759855.1 0.9276
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
            fh.write(f'{name}\t{round(cos, 4)}\n')
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
    '--topn', 
    help='Number of top similar vectors', default=10)
@click.option(
    '--outfile', '-o',
    help='Where to write results', required=True, type=click.Path())
@click.option(
    '--embedding',
    help='Genome embedding model', required=True, type=click.Path())
@click.option(
    '--taxonomy',
    help='File path to GTDB taxonomy', required=True, type=click.Path())
@click.option(
    '--query',
    help='File path genome annotation', required=True, type=click.Path())
@click.option(
    '--fmt',
    help='Query fmt (pfamscan or hmmer)', default='hmmer')
@click.option(
    '--steps',
    help='How many epochs for vector inference', default=200)
def taxonomy(query, taxonomy, embedding, topn, outfile, fmt, steps):
    '''
    Given a query vector, get the <n> closest vectors and their taxonomy and
    then report their <raw> taxonomy or use <majority vote> to identify the
    most likely (?) one.

    Usage:

    \b
    nanotext taxonomy \\
        --embedding nanotext_r89.model --taxonomy bac_taxonomy_r86.tsv \\
        --query JFOD01_pfam.tsv --fmt pfamscan --topn 10 -o results.json

    '''
    from collections import Counter, defaultdict
    import json
    import random

    from nanotext.io import load_taxonomy_gtdb, load_embedding, eprint
    from nanotext.utils import infer_genome_vector


    ranks = [
        'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    results = {}
    results['topn'] = topn
    results['raw'] = {}
    notfound = []
    
    db = load_taxonomy_gtdb(taxonomy)
    model = load_embedding(embedding)

    v = infer_genome_vector(query, model, fmt=fmt, steps=steps)
    sim = model.docvecs.most_similar([v], topn=topn)
    
    vote = defaultdict(list)
    names = {}
    for name, cos_sim in sim:
        names[name] = round(cos_sim, 4)
        try:
            d = {k: v for k, v in zip(ranks, db[name])}
            for k, v in d.items():
                vote[k].append(v)
        except KeyError:
            # eprint(f'{name} has no taxonomy record')
            notfound.append(name)
            continue

    results['notfound'] = notfound
    results['similarity'] = names

    majority, majority_ratio = {}, {}

    for rank, taxa in vote.items():
        results['raw'][rank] = taxa
        cnt = Counter(taxa)
        maxn = max(cnt.values())
        hits = [k for k, v in cnt.items() if v == maxn]
        pick = random.choice(hits)
        majority[rank] = pick
        majority_ratio[rank] = round(cnt[pick]/len(taxa), 2)

    results['majority'] = majority
    results['ratio'] = majority_ratio
    
    with open(outfile, 'w+') as out:
        json.dump(dict(sorted(results.items())), out, indent=4)


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
