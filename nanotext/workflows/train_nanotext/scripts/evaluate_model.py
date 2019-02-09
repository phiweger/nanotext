import click

from nanotext.evaluate import odd_one


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
    '--outfile', '-o', 
    help='Name of the resulting evaluation report', 
    type=click.Path(), required=True)
def evaluate(model, corpus, outfile):
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
    SOMO = odd_one(corpus, model, n_not_odd=5)
    with open(outfile, 'w+') as out:
        out.write(str(SOMO)+'\n')


if __name__ == '__main__':
    evaluate()
