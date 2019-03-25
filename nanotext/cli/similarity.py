import click

# TODO: nanotext index
# https://github.com/facebookresearch/faiss
# https://github.com/facebookresearch/faiss/wiki/Index-IO,-index-factory,-cloning-and-hyper-parameter-tuning#io-and-deep-copying-indexes


@click.command()
@click.option(
    '--annotation',
    help='Protein domain annotation from pfam_scan.pl',
    type=click.Path())
@click.option(
    '--topn', 
    help='Top n hits to return',
    default=10)
@click.option(
    '--models',
    help='Genome embedding models',
    required=True, type=click.Path())
@click.option(
    '--mode',
    help='Model w/ focus on core/ accessory/ an ensemble of domains',
    default='core')
@click.option(
    '--taxonomy',
    help='Lookup GTDB taxonomy for hits',
    default=False)
@click.option(
    '--out',
    help='Output path (tsv format). If not specified, write to stdout.',
    default='-')
def search(annotation, topn, models, mode, taxonomy, out):
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
    from nanotext.classes import GenomeModel
    from nanotext.io import load_embedding, smart_open, eprint
    from nanotext.utils import infer_genome_vector, get_taxa_from_names

    # fp, mode, 

    eprint('Loading model ...')
    model = GenomeModel(models, mode=mode, norm='l2')
    v = model.infer(annotation, fmt='pfamscan', steps=1000)
    sim = model.search(v, topn)
    with smart_open(out) as fh:
        for name, cos in sim:
            fh.write(f'{name}\t{round(float(cos), 4)}\n')

    if taxonomy:
        names, _ = zip(*sim)
        df = get_taxa_from_names(taxonomy, names)
        df.to_csv('taxonomy.tsv', sep='\t', index=None)

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

    '''
    TODO: new fmt

    name, cos, ranks

    2nd output

    majority vote across columns
    '''
    from collections import Counter, defaultdict
    import json
    import pdb
    import random

    import numpy as np
    from sklearn.manifold import TSNE
    import umap

    from nanotext.io import load_taxonomy_gtdb, load_embedding, eprint
    from nanotext.utils import infer_genome_vector, strip_name


    config_umap_visualisation = {
        'metric': 'cosine',
        'n_components': 2,
        # 'repulsion_strength': 5,
        # 'n_neighbors': 5,
        # min_dist=0.05,
        # 'spread': 5,
        'random_state': 42,
        }


    ranks = [
        'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    notfound = []
    
    db = load_taxonomy_gtdb(taxonomy)
    model = load_embedding(embedding)

    v_query = infer_genome_vector(query, model, fmt=fmt, steps=steps)
    sim = model.docvecs.most_similar([v_query], topn=topn)


    taxcollector = {i: [] for i in ranks}
    distance = {}

    for name, cos_sim in sim:
        distance[name] = round(cos_sim, 4)
        try:
            for k, v in zip(ranks, db[name]):
                taxcollector[k].append(v)

        except KeyError:
            eprint(f'{name} has no taxonomy record')
            continue


    # What is the last uniform rank?
    cache = ()
    for i in ranks:
        if (len(set(taxcollector[i])) == 1) and (taxcollector[i][0] != ''):
            cache = (i, taxcollector[i][0])
            continue
        else:
            pass
    

    # p__Firmicutes_A
    # Collect the UIDs for this rank.
    eprint(f'Will collect all vectors for {cache[0]} {cache[1]} ...')
    names = []
    with open(taxonomy, 'r') as file:
        for line in file:
            # if 'c__Clostridia' in line:
            # if 'f__Pseudomonadaceae' in line:
            # if ('p__Firmicutes_A' in line) or ('p__Firmicutes_B' in line): 
            if f'{cache[0][0]}__{cache[1]}' in line:  # d__ for domain etc.
                names.append(line.strip().split('\t')[0])


    # Collect the associated document vector for each UID.
    m, found = [], []
    for name in names:
        try:
            m.append(model.docvecs[strip_name(name)])
            found.append(strip_name(name))
        except KeyError:
            continue


    # Project into 2D.
    eprint(f'Projecting with UMAP ...')
    m = np.array(m, dtype='float64')
    reducer = umap.UMAP(**config_umap_visualisation)
    # projection = reducer.fit_transform(m)
    eprint(f'Projecting with TSNE ...')
    projection = TSNE(n_components=2, random_state=42).fit_transform(m)

    results = defaultdict(list)
    # results['query'].extend(reducer.transform([v_query])[0])
    # results['query'].extend(7*['query'])

    for i, j in zip(found, projection):
        results[i].extend(j)
        results[i].extend(db[i])


    # Add distance info.
    for k, v in results.items():
        results[k].append(distance.get(k, 'NA'))

    # majority, majority_ratio = {}, {}

    # for rank, taxa in vote.items():
    #     results['raw'][rank] = taxa
    #     cnt = Counter(taxa)
    #     maxn = max(cnt.values())
    #     hits = [k for k, v in cnt.items() if v == maxn]
    #     pick = random.choice(hits)
    #     majority[rank] = pick
    #     majority_ratio[rank] = round(cnt[pick]/len(taxa), 2)

    # results['majority'] = majority
    # results['ratio'] = majority_ratio
    
    with open(outfile, 'w+') as out:
        out.write('\t'.join(
            'name c1 c2 domain phylum class order family genus species cos'.split())+'\n')
        # json.dump(dict(sorted(results.items())), out, indent=4)
        for k, v in results.items():
            line = '\t'.join([str(i) for i in [k]+v])+'\n'
            out.write(line)


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
