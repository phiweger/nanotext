import click

@click.command()
@click.option(
    '--config',
    help='Params for training the embedding.')
@click.option(
    '--corpus', 
    help='Corpus of protein domains to train on.', 
    type=click.Path())
@click.option(
    '--out',
    help='Model name (output).')
def train(corpus, config, out):
    '''
    Usage:

    \b
    nanotext train --corpus corpus.txt --out embedding.genomes.model

    '''
    import random
    
    from gensim.models.doc2vec import TaggedDocument
    from gensim.models import Doc2Vec
    from gensim.models.callbacks import CallbackAny2Vec
    from tqdm import tqdm


    tagged_data = []
    cnt = 0
    print('Reading corpus ...')
    with open(corpus, 'r') as file:
        for line in tqdm(file):
            genome, contig, domains = line.strip().split('\t')
            domains = domains.split(',')
    
            tagged_data.append(TaggedDocument(words=domains, tags=[genome]))
            cnt += len(domains)
    
    # TODO: add genome count
    print(f'Corpus has {len(tagged_data)} contigs.')
    # Corpus has 2,782,580 contigs.
    print(f'Corpus has {cnt} domains.')
    # Corpus has 144,988,433 domains.


    epochs = 10  # will be time <epochs>, by default 5 per iteration
    random.seed(42)
    _ = random.shuffle(tagged_data)
    
    
    model = Doc2Vec(
        vector_size=100,
        hs=0, negative=5,  # neg10
        # from word2vec, defaults, e.g. Doc2Vec().negative is 5
        # 1 for archaea
        min_count=3,
        sample=0.001,  # default 0.001
        workers=8,
        window=10,  #w5
        # dm=1,
        dm=0, dbow_words=1,
        epochs=epochs,
        alpha=0.025,
        min_alpha=0.0001,
        seed=42,
        )

    model.build_vocab(tagged_data)

    print('Training starts ...')
    _ = model.train(
        tagged_data, total_examples=model.corpus_count, epochs=model.epochs)

    model.save(out)


    