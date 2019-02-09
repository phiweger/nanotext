import argparse
from collections import defaultdict
from glob import glob
import json

import numpy as np
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
    parser.add_argument('--test_ratio', type=float, default=0.1, 
        help='Size of test set in the range [0, 1].')
    # flag (switch), stackoverflow.com/questions/8259001
    parser.add_argument('--keep_unknown', action='store_true', 
        help='Include "unknown" label for ORFs w/o domains.')
    args = parser.parse_args()


    print('Recursively parsing annotation directory ...')
    # stackoverflow.com/questions/2186525
    files = glob(f'{args.annotations}/**/*.json', recursive=True)
    # There are 43914 files.
    _ = np.random.shuffle(files)  # should there be any meaning in the ordering
    seen_organisms = defaultdict(list)


    for fp in files:
        organism = get_name_from_fp(fp)
        seen_organisms[organism].append(fp)


    print(f'Balancing samples to a maximum of n={args.cap} per species ...')
    balanced_sample = []
    for v in seen_organisms.values():
        try:
            balanced_sample.extend(
                np.random.choice(v, size=args.cap, replace=False))
        except ValueError:  
            # Cannot take a larger sample than population when 'replace=False'
            balanced_sample.extend(v)  # take all
    print(f'Corpus (train, test) holds {len(balanced_sample)} samples.')


    print('Parsing annotations ...')
    with open(f'{args.outdir}/corpus.ensemble.train.txt', 'w+') as train, \
         open(f'{args.outdir}/corpus.ensemble.test.txt', 'w+') as test: 
        
        for fp in tqdm(balanced_sample):
            organism = get_name_from_fp(fp)

            with open(fp, 'r') as file:
                anno = json.load(file)
            
            genome, domains = tokenize_ensembl(
                anno, keep_unknown=args.keep_unknown)

            out = np.random.choice(
                [train, test], p=[1-args.test_ratio, args.test_ratio])
            # Counter([np.random.choice([1, 2]) for _ in range(1000)])
            # Counter({2: 494, 1: 506})
    
            for k, v in domains.items():
                # k is contig, v is list of domains
                out.write(f'{genome}\t{k}\t{",".join(v)}\n')
                out.flush()  # stackoverflow.com/questions/3167494


if __name__ == '__main__':
    main()