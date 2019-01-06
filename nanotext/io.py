import contextlib

@contextlib.contextmanager
def smart_open(filename=None):
    '''
    stackoverflow.com/questions/17602878
    '''
    import sys
    if filename and filename != '-':
        fh = open(filename, 'w')
    else:
        fh = sys.stdout

    try:
        yield fh
    finally:
        if fh is not sys.stdout:
            fh.close()


def eprint(*args, **kwargs):
    '''
    stackoverflow.com/questions/5574702
    '''
    import sys
    print(*args, file=sys.stderr, **kwargs)


def load_domains(fp, fmt='pfamscan'):
    '''
    fmt .. format can be 'pfamscan' or 'hmmer'

    TODO: filter, recommendation bacteria:
    E-value 1e-18, cov 0.35
    for bacteria, use E-value < 1e-18 and coverage > 0.35


    > Our suggestion is that for plants, use E-value < 1e-23 and coverage > 0.2; for bacteria, use E-value < 1e-18 and coverage > 0.35; and for fungi, use E-value < 1e-17 and coverage > 0.45.

    http://csbl.bmb.uga.edu/dbCAN/download/readme.txt

    data based on CAZyDB (v6) released on 07/20/2017
    '''
    from nanotext.utils import to_interval

    if fmt == 'pfamscan':
        num_lines_in_header = 29
        header = 'seq_id alignment_start alignment_end envelope_start envelope_end hmm_acc hmm_name type hmm_start hmm_end hmm_length bit_score E-value significance clan'.split()

        intervals = {}
        with open(fp, 'r') as file:
        
            for i in range(num_lines_in_header):
                _ = next(file)  # skip header
        
            for line in file:
                d = {k:v for k, v in zip(header, line.strip().split())}
                # Interval(chrom, start, end, name=".", score=".", strand=".",
                # otherfields=None)
                uid =  d['seq_id']
                start = int(d['envelope_start'])
                end = int(d['envelope_end'])
                intervals[(uid, start, end)] = to_interval(
                    uid, start, end, d['hmm_acc'], d['E-value'])
        return intervals

    elif fmt == 'hmmer':
        try:
            intervals = {}
            with open(fp, 'r') as file:
                try:
                    header = next(file).strip().split('\t')
                except StopIteration:
                    print('Empty file. Abort!')
                    return None

                header[4] = 'accession_query'  
                # originally named "accession" -- would create 2 columns w/ 
                # same name, bad
                for line in file:
                    d = {k:v for k, v in zip(header, line.split('\t'))}
                    uid = d['query name']
                    start, end = int(d['from']), int(d['to'])
                    intervals[(uid, start, end)] = to_interval(
                        uid, start, end, d['accession'], d['E-value'])
            return intervals
        except KeyError:
            print('Did you reformat the hmmscan result using HmmPy.py?')
            print('Yes, I thought so. Abort!')
            return None

    else:
        print(f'Only formats "pfamscan" and "hmmer" are supported. Abort!')
        return None


def load_orfs(fp, fmt='gff'):
    '''
    orfs = load_orfs('orfs.gff')
    '''
    from pybedtools import Interval

    if fmt == 'gff':
        gff_header = 'contig x1 coding start end x2 strand x3 attrs'.split()
        orfs = {}

        with open(fp, 'r') as file:
            for line in file:
                if not line[0] == '#':  # skip header
                    d = {k:v for k,v in zip(gff_header, line.strip().split('\t'))}
                    _, ix = d['attrs'].split(';')[0].replace('ID=', '').split('_')

                    contig = d['contig']
                    orfs[(d['contig'], int(ix))] = Interval(
                        d['contig'], int(d['start']), int(d['end']), 
                        name=ix, strand=d['strand'])
    elif fmt == 'fasta':
        pass

    else:
        print(f'Only formats "gff" and "fasta" are supported. Abort!')
        return None
    
    return orfs


def tokenize_ensembl(anno, keep_unknown=True):
    '''
    Take an Ensembl annotation loaded from json and extract the domains, like
    words in a document. Collapse the double strand into one string of domains.

    Omit all domains that are not on the chromosome, e.g. plasmids such
    as p3ABAYE. The GTDB omits plasmids, bc/ they behave differently from the
    chromosome and follow their own "life history" (think the selfish gene).

    Return the assembly ID and a dict of contigs (keys) and associated
    domains (values), the latter linearly ordered. If <keep_unknown>, ORFs
    without Pfam annotations will receive the "unknown" label.
    '''
    from collections import defaultdict

    genome = anno['assembly']['accession']

    d = {}

    for g in anno['genes']:  # g .. gene
        start = int(g['start'])
        end = int(g['end'])
        strand = g['strand']
        
        # TODO: How to filter plasmids? <region> too heterogenous, not 
        # informative.
        # region = g['seq_region_name']  

        try:
            contig = g['ENA_FEATURE_GENE'][0].split(':')[0]
        except KeyError:
            contig = 'NA'
        # print(g['transcripts']['exons'])

        uid = (contig, start, end)  
        # <strand> not needed for now, bc/ we won't disambiguate/ deduplicate
        # the domains for now (no E-value for domain plausibility in Ensembl
        # annotation).
        
        try:
            domains = g['Pfam']
        except KeyError:
            if keep_unknown:
                domains = ['unknown']
            else:
                pass

        if (contig != 'NA'):  # and (region == 'Chromosome') -- see TODO above
            d[uid] = domains

    # By sorting starts and ends, we recreate the contiguity of ORFs.
    # We can then linearly write the domains, whose context is preserved.
    result = defaultdict(list)
    
    for k, v in sorted(d.items()):
        # [(('CU459141.1', 170, 1567), ['PF11638', 'PF00308', 'PF08299']), ...
        result[k[0]].extend(v)

    return genome, result


def load_embedding(fp, algorithm='doc2vec'):
    from gensim.models import Doc2Vec

    if algorithm == 'doc2vec':
        model = Doc2Vec.load(fp)
        return model
    else:
        print('Not implemented yet, sorry.')
        return None


def save_embedding(prefix, model):
    '''
    Save gensim model in GloVe format.

    https://radimrehurek.com/gensim/scripts/glove2word2vec.html
    https://nlp.stanford.edu/projects/glove/

    word1 0.123 0.134 0.532 0.152
    word2 0.934 0.412 0.532 0.159
    word3 0.334 0.241 0.324 0.188
    ...
    word9 0.334 0.241 0.324 0.188
    '''
    # word vectors
    with open(f'{prefix}.wv.txt', 'w+') as out:
        for k in model.wv.vocab.keys():
            v = ' '.join([str(i) for i in model.wv[k]])
            out.write(f'{k} {v}\n')
    # document vectors
    with open(f'{prefix}.dv.txt', 'w+') as out:
        for k in model.docvecs.doctags.keys():
            v = ' '.join([str(i) for i in model.docvecs[k]])
            out.write(f'{k} {v}\n')