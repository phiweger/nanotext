import argparse
from collections import defaultdict
from glob import glob
import json
import random

from tqdm import tqdm

from nanotext.io import tokenize_ensembl


def get_name_from_fp(fp):
    import os

    prefix = os.path.basename(fp).replace('.json', '')
    # acinetobacter_baumannii_aye
    organism = '_'.join(prefix.split('_')[:2])
    return organism


def main():
    '''
    Collect Ensembl annotations, balance them and return a training corpus as
    well as a list of files not included in the corpus for testing.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--annotations', 
        help='Directory to Ensembl annotations -- will recursively collect all .json files.')
    parser.add_argument('--cap', type=int, 
        help='Maximum number of samples for each species (to balance data set).')
    parser.add_argument('--outdir', 
        help='Output directory.')
    # flag (switch), stackoverflow.com/questions/8259001
    parser.add_argument('--keep_unknown', action='store_true', 
        help='Include "unknown" label for ORFs w/o domains.')
    args = parser.parse_args()


    # stackoverflow.com/questions/2186525
    files = glob(f'{args.annotations}/**/*.json', recursive=True)
    seen_organisms = defaultdict(int)


    with open(f'{args.outdir}/corpus.ensemble.train.txt', 'w+') as train, \
         open(f'{args.outdir}/corpus.ensemble.test.txt', 'w+') as test: 
        
        for fp in tqdm(files):
            organism = get_name_from_fp(fp)
            # print(organism)
        
            if seen_organisms[organism] < 2*args.cap:  
            # <cap> for training and testing
                
                with open(fp, 'r') as file:
                    anno = json.load(file)
                
                genome, domains = tokenize_ensembl(
                    anno, keep_unknown=args.keep_unknown)
                out = random.choice([train, test])
    
                for k, v in domains.items():  
                    # k is contig, v is list of domains
                    out.write(f'{genome}\t{k}\t{",".join(v)}\n')
                    out.flush()  # stackoverflow.com/questions/3167494
    
                seen_organisms[organism] += 1


if __name__ == '__main__':
    main()