from glob import glob
import os

import click
from pybedtools import BedTool
from tqdm import tqdm

from nanotext.io import load_domains
from nanotext.utils import create_domain_sequence


@click.command()
@click.option(
    '--indir', help='Folder w/ input (annotation) files',
    type=click.Path(), required=True)
@click.option(
    '--outfile', help='Where to store corpus',
    type=click.Path(), required=True)
@click.option(
    '--minlen', help='Minimum number of domains per contig',
    default=0)
def ingest(indir, outfile, minlen):
    '''
    Ingest a collection of protein domain annotations generated using
    `pfam_scan.pl` into a corpus for later use with Doc2Vec.
    '''
    files = glob(f'{indir}/*')  # 145817

    # the next step will take roughly 10 hours
    with open(outfile, 'w+') as out:
        for file in tqdm(files):
        
            result = load_domains(file, fmt='pfamscan')
            dom = BedTool(list(result.values()))
        
            # make sure Pfam ID is truncated: PF00815.1 -> PF00815
            seq = create_domain_sequence(
                dom, keep_unknown=True, fmt_fn=lambda x: x.split('.')[0])
        
            text = list(seq.values())
            genome = os.path.basename(file).strip('_pfam.tsv')
            # e.g. ...
            # UBA9934_pfam.tsv
            # GB_GCA_001790445.1_pfam.tsv
            # RS_GCF_000012865.1_pfam.tsv

            # adjust name if pfam table name from GTDB
            if (not 'UBA' in genome) and \
               (any(x in genome for x in ['RS_', 'GB_'])):
                genome = '_'.join(genome.split('_')[1:])

    
            for k, v in seq.items():
                if len(v) > minlen:
                    out.write(f'{genome}\t{k}\t{",".join(v)}\n')


if __name__ == '__main__':
    ingest()