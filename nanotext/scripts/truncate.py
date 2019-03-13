import click


@click.command()
@click.option(
    '--infile', help='Fasta file to be truncated',
    type=click.Path(), required=True)
@click.option(
    '--suffix', help='Suffix of input file, e.g. ".fna"',
    type=click.Path(), default='.fna')
@click.option(
    '--outdir', help='Where to store truncated files',
    type=click.Path(), required=True)
@click.option(
    '--stepsize', help='Truncate from <low> to <high> in this stepsize',
    default=0.05)
@click.option(
    '--low', help='Minimum of stepsize range',
    default=0)
@click.option(
    '--high', help='Maximum of stepsize range',
    default=1)
@click.option(
    '--mode', help='remove_sequence or remove_contigs?',
    default='remove_sequence')
def truncate(infile, suffix, outdir, stepsize, low, high, mode):
    '''Truncate a fasta file

    Usage:

    \b
    mkdir test
    python truncate.py \\
        --infile GCA_900130355.1_10625_5_68_genomic.fna --suffix .fna \\
        --outdir test --stepsize 0.1
    '''
    import os

    from Bio import SeqIO
    import numpy as np
    
    from nanotext.utils import truncate


    fa = SeqIO.index(infile, 'fasta')
    # trunc = {}

    outfile_prefix = os.path.basename(infile).replace(suffix, '')
    
    # Loop over fasta and randomly discard contigs.
    if mode == 'remove_contigs':
        print('Not implemented yet, nothing will happen.')

    # Loop over fasta entries and successively remove part of the sequences.
    elif mode == 'remove_sequence':

        for by in np.arange(low, high, stepsize):
            
            fp_out = f'{outdir}/{outfile_prefix}_truncatedby{int(by*100)}.fna'
            with open(fp_out, 'w+') as out:
                for i in fa:
                    seq = fa[i].seq.__str__()  # truncate() fn takes a list
                    for j, fragment in enumerate(truncate([seq], by=by)):
                        name = f'{i}_truncatedby{int(by*100)}_fragment{j}'
                        out.write(f'>{name}\n{fragment}\n')
                        # trunc[name] = fragment

    # eprint('Done.')


if __name__ == '__main__':
    truncate()