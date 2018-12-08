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
                d = {k:v for k, v in zip(header, line.split())}
                # Interval(chrom, start, end, name=".", score=".", strand=".",
                # otherfields=None)
                uid =  d['seq_id']
                start = int(d['envelope_start'])
                end = int(d['envelope_end'])
        
                intervals[(uid, start, end)] = to_interval(
                    uid, start, end, d['hmm_acc'], d['E-value'])

        return intervals

    elif fmt == 'hmmer':
        print('Not yet implemented, sorry. Abort!')
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
        region = g['seq_region_name']  

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
            
            # unknown label in list bc/ domains is a list which we'll flatten
            # RESULT: does not make a difference to the odd one out metric
        if (contig != 'NA') and (region == 'Chromosome'):
            d[uid] = domains

    # By sorting starts and ends, we recreate the contiguity of ORFs.
    # We can then linearly write the domains, whose context is preserved.
    result = defaultdict(list)
    
    for k, v in sorted(d.items()):
        # [(('CU459141.1', 170, 1567), ['PF11638', 'PF00308', 'PF08299']), ...
        result[k[0]].extend(v)

    return genome, result