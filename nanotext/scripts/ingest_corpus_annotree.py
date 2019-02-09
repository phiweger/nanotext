from collections import defaultdict
from glob import glob
import os
import random

from gensim.models.doc2vec import TaggedDocument
from gensim.models import Doc2Vec
from pybedtools import BedTool
from tqdm import tqdm

from nanotext.io import load_domains
from nanotext.utils import remove_overlap, deduplicate, create_domain_sequence


# fp = '/Users/phi/data_local/databases/annotree/pfam_archaea/RS_GCF_001282785.1_pfam.tsv'


# with open('corpus_annotree_bacteria.txt', 'w+') as out:
#     for fp in tqdm(glob('pfam_bacteria/*')):
#         text = tokenize_pfam_scan(fp)
#         for line in text:
#             out.write(f'{name},{line}\n')


# fp = '/Users/phi/data_local/databases/annotree/pfam_archaea/*'
fp = '/Users/phi/data_local/databases/annotree/pfam_bacteria/*'
tagged_data = []
corpus = {}

with open('/Users/phi/tmp/corpus.annotree.train.txt', 'w+') as out:
    for file in tqdm(glob(fp)):
    
        # TODO: skip overlap removal and deduplication for now, we need strand
        # info
        # for this
        # Rich Hickey would be proud ...
        # dom = deduplicate(remove_overlap(load_domains(file, fmt='pfamscan')))
    
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
        if not 'UBA' in genome:
            genome = '_'.join(genome.split('_')[1:])

        for k, v in seq.items():
            out.write(f'{genome}\t{k}\t{",".join(v)}\n')

    # corpus[name] = text  # This is not necessary, here for exploring tokens later.



