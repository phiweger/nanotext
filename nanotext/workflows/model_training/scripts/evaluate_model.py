import click

from nanotext.evaluate import odd_one, ecotype_task
from nanotext.io import load_ecotypes, load_embedding, eprint


@click.command()
@click.option(
    '--model', '-m', 
    help='The nanotext model to evaluate', 
    type=click.Path(), required=True)
@click.option(
    '--corpus', '-i', 
    help='A subset of the corpus set aside for testing',
    type=click.Path(), required=True)
@click.option(
    '--ecotypes', '-e', 
    help='Path to the ecotype labels', 
    type=click.Path(), required=True)
@click.option(
    '--outfile', '-o', 
    help='Name of the resulting evaluation report', 
    type=click.Path(), required=True)
def evaluate(model, corpus, ecotypes, outfile):
    '''
    A test battery:
    
    - SOMO task -- find the domain in a sequence that does not fit
    - Ecotype task -- separate niche-specific subpopulations of misc species

    Returns a list of accuracy - (sub)task pairs.
    '''
    with open(outfile, 'w+') as out:

        # (1)
        eprint('SOMO task ...')
        SOMO = odd_one(corpus, model, n_not_odd=5)
        out.write(f'{str(SOMO)}\tSOMO\n')

        # (2)
        eprint('Ecotype task ...')
        m = load_embedding(model)
        for task in ['vibrio', 'prochlorococcus', 'pseudomonas']:
            rank, eco = load_ecotypes(f'{ecotypes}/{task}.tsv')
            d = {k: v2 for k, (v1, v2) in eco.items()}
            for k, v in ecotype_task(d, m).items():
                if k != 'NA':
                    out.write(f'{v}\t{k}\n')

        eprint('Evaluation done.')


if __name__ == '__main__':
    evaluate()
