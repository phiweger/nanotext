def strip_name(name):
    '''
    Takes an ID from GTDB and retrieves the original (NCBI) ID -- basically
    sometimes removes a prefix.

    name = 'RS_GCF_000372645.1'
    strip_name(i['uid'])
    # 'GCA_002387705.1'
    '''
    if 'UBA' in name:
        return name
    elif 'U_' in name:
        return name
    else:
        return '_'.join(name.split('_')[1:])


def cosine(a, b):
    '''
    cosine similarity, range(-1, 1)
    '''
    from numpy import dot
    from numpy.linalg import norm

    return dot(a, b)/(norm(a)*norm(b))


def get_interval_uid(interval):
    contig, start, end = interval.fields[:3]
    return (contig, int(start), int(end))


def to_interval(uid, start, end, name='', score='.', strand='.'):
    '''
    Convenience fn to contruct interval w/ all-positional args.

    Main advantage: If working w/ different domain formats (HMMER domtbl
    vs. pfam_scan.pl output), this is the only fn we need to change.
    '''
    from pybedtools import Interval

    return Interval(
        uid, int(start), int(end), name=name, score=score, strand=strand)


def is_nested(u, v):
    '''This fn is not symmetric: (u in v) != (v in u). Fn asks former.'''
    if (u.start > v.start) and (u.end < v.end):
        return True
    else:
        return False


def deduplicate(intervals):
    '''
    Given a BedTool() object of intervals, remove successions of intervals
    w/ the same name field (3) -- e.g. consecutive domains of the same name.
    For those, return only the first occurrence, i.e.

    A-B-C-C-C-D becomes A-B-C-D

    Example data:
    
    NZ_CVUA01000001.1_2711  1037    1107    PF01131.15  6.3e-14 .
    NZ_CVUA01000001.1_2711  1142    1180    PF01396.14  1.5e-14 .
    NZ_CVUA01000001.1_2711  1227    1270    PF01396.14  1.4e-08 .
    NZ_CVUA01000001.1_2711  1313    1372    PF14520.1   5.8e-10 .
    
    There are weird edge cases:
    
    NZ_CVUA01000001.1_953   80      382     PF03104.14  1e-33   .
    NZ_CVUA01000001.1_953   891     993     PF00136.16  4.5e-05 .
    NZ_CVUA01000001.1_953   991     1101    PF13403.1   1.9e-08 .
    NZ_CVUA01000001.1_953   1402    1718    PF00136.16  3.7e-33 .
    NZ_CVUA01000001.1_957   5       224     PF13476.1   4.5e-25 .
    '''
    from pybedtools import BedTool

    result = []
    cache = 0
    for i in intervals:
        if not cache:
            cache = i
            result.append(i)
        else:
            if \
                (i.fields[3] == cache.fields[3]) and \
                (i.fields[0] == cache.fields[0]):
                # first condition: same Pfam ID
                # second condition: same contig
                pass  # omit domain, i.e. deduplicate
            else:
                cache = i
                result.append(i)

    return BedTool(result)


def remove_overlap(intervals, deduplicate=True):
    '''
    Take a dict of the form:

    {(contig, start, end): Interval(...)}

    and remove those features that overlap based on their plausibility --
    namely their E-value in the score column. Smaller E-values are more 
    plausible.
    '''
    from pybedtools import BedTool
    from nanotext.utils import to_interval, is_nested, get_interval_uid

    result = intervals.copy()  # do not modify input in place
    dom = BedTool(list(result.values()))

    for i in dom.intersect(dom, wa=True, wb=True):  # contig order is preserved
        u = to_interval(*i.fields[:6])
        v = to_interval(*i.fields[6:])
        # features can (1) overlap partly, (2) be nested or (3) not overlap
        try:
            # (1) overlap but not w/ self (54/2964 of test cases) -- 
            # drop smaller E
            if not u == v:
                evalues = [float(j.fields[4]) for j in [u, v]]
                # get the interval w/ the higher E-value and delete
                ix = evalues.index(max(evalues))
                k = get_interval_uid([u, v][ix])
                _ = result.pop(k, None)
    
            # (2) all features overlap w/ themselves, trivial -- skip
            else:
                pass
        
        # (3) nested domains (78/2964 test cases) -- get the larger interval
        except NotImplementedError:
            # "Features are nested -- comparison undefined"
            # e.g. helix-loop-helix in larger domain
            if is_nested(u, v):  # asks: u in v?
                _ = result.pop(get_interval_uid(u), None)
            else:
                _ = result.pop(get_interval_uid(v), None)

    return BedTool(list(result.values()))


def split_orf_uid(s):
    '''
    in:  NZ_CVUA01000001.1_993
    out: (NZ_CVUA01000001.1, 993)

    If there is no suffix of ORF uid, return 1
    '''
    *contig, orf = s.split('_')
    try:
        return '_'.join(contig), int(orf)
    except ValueError:  # int() of string
        return s, 1



def create_domain_sequence(domains, keep_unknown=True, fmt_fn=lambda x: x):
    '''
    We can pass a <fmt_fn> to reformat domain names (e.g. turn PF07005.14
    into PF07005). By default, the identity function is used, which just
    returns the input as is.

    To e.g. reformat PF00815.1 -> PF00815: fmt_fn=lambda x: x.split('.')[0]

    If <keep_unknown>, all ORFs w/o domain calls are turned into "unknown"
    domains, 1 per ORF.
    '''
    from collections import defaultdict
    from nanotext.utils import split_orf_uid

    d = defaultdict(list)
    
    for i in domains:
        d[split_orf_uid(i.fields[0])].append(fmt_fn(i.fields[3]))
    
    result = defaultdict(list)
    cache_orf = 0  # prodigal ORFs are indexed starting at 1
    cache_contig = ''

    for k, v in sorted(d.items()):  # only sorts keys
        contig, orf = k
        
        # this block makes the routine contig aware
        if cache_contig != contig:
            cache_contig = contig
            cache_orf = 0  # reset the ORF counter when on each new contig
    
        if keep_unknown and (cache_orf+1 != orf):
            result[contig].extend(['unknown']*(orf-cache_orf-1))
    
        result[contig].extend(v)
        cache_orf = orf

    return result


def infer_genome_vector(fp, model, steps=200, fmt='hmmer', truncate_by=0):
    '''
    From a genome annotation either from Pfam or HMMER (formatted w/ HMMPy.py)
    infer a genome vector.

    Note that the genomes contigs are concatenated, which at the boundaries of
    the contigs creates aritificial context, i.e. domains that do not normally
    co-occur. However, in our experiments this seemed to not affect the
    accuracy of the resulting vectors much, and for reasonably fragmented
    genome assemblies this should not affect results much.

    We infer the vector w/ 200 steps (iterations) which during testing provided
    a mean cosine distance of the estimates < 0.01 -- more steps will reduce this variance and increase time needed to compute the vector.
    '''
    from pybedtools import BedTool

    from nanotext.io import load_embedding, load_domains
    from nanotext.utils import create_domain_sequence

    result = load_domains(fp, fmt=fmt)
    dom = BedTool(list(result.values()))
    seq = create_domain_sequence(
        dom, keep_unknown=True, fmt_fn=lambda x: x.split('.')[0])
    # fmt_fn here splits version number from Pfam domain PF00001.1 -> PF00001  
    
    # concatenate protein domain sequences from contigs
    if truncate_by:
        flat = [item for sublist in truncate(seq.values(), truncate_by) for item in sublist]
    else:
        flat = [item for sublist in seq.values() for item in sublist]
    
    # 200 epochs inference gives a varience < 0.01 cosine distance on
    # when repeatedly inferring vectors (from our experiments)
    return model.infer_vector(flat, steps=steps)  
    
    
def truncate(sequences, by=0.5):
    '''Given several sequences, truncate them <by> a given fraction.'''
    import random

    for seq in sequences:
        cut = int(by*len(seq))
        start = random.choice(range(0, len(seq)-cut))
        seq1 = seq[:start]
        seq2 = seq[start+cut:]
        for j in [seq1, seq2]:
            if j:
                yield j


def subset_taxonomy(query, taxa):
    '''
    Given a query in the GTDB format of <first letter rank>__<name> (e.g.
    "f__Pseudomonadaceae"), return all IDs from the taxonomy database. 
    '''
    taxon, name = query.split('__')
    groups = 'd p c o f g s'.split()
    result = []
    
    for k, v in taxa.items():
        d = dict(zip(groups, v))
        if d[taxon] == name:
            result.append(k)
    return result


def subset_model_by_rank(model, db, taxon):
    '''
    From a taxonomy database, select those records in the model that match.
    '''
    entries = model.docvecs.index2entity
    ranks = 'd p c o f g s'.split()
    # domain phylum class order family genus species
    rank, name = taxon.split('__')
    ix = [i for i, j in enumerate(ranks) if j == rank][0]
    
    cnt, found = 0, []
    for entry in entries:
        record = db.get(entry, None)
        try:
            if record[ix] == name:
                found.append(entry)
        except TypeError:
            cnt += 1
            continue

    return found


def get_vectors(l, model, normalized=False):
    '''Get array of document vectors given an ID list and a model'''
    import numpy as np
    from sklearn.preprocessing import normalize

    from nanotext.io import eprint

    m, found = [], []
    cnt = 0

    for i in l:
        try:
            m.append(model.docvecs[i])
            found.append(i)
        except KeyError:
            cnt += 1

    m = np.array(m, dtype='float64')
    if normalized:
        m = normalize(m, norm='l2', axis=1)
    if cnt:
        eprint(f'{cnt} entries missing ({round(cnt/len(l), 4)} %)')
    return found, m



def flatten(l):
    '''
    Given a nested iterable return unnested items in list.
    '''
    # import itertools
    # stackoverflow.com/questions/952914
    # return [item for sublist in l for item in sublist]
    # return list(itertools.chain.from_iterable(l))
    # Problem: both split strings
    # stackoverflow.com/questions/5286541
    for x in l:
        if hasattr(x, '__iter__') and not isinstance(x, str):
            for y in flatten(x):
                yield y
        else:
            yield x


def get_names_from_taxon(db, taxon):
    '''Given a taxon in GTDB format (e.g. c__Clostridia) return model UIDs'''
    from nanotext.utils import strip_name
    from nanotext.io import dbopen

    groups = {
        'c': 'class',
        's': 'species',
        'p': 'phylum',
        'f': 'family',
        'd': 'domain',
        'o': 'order',
        'g': 'genus',}
    rank = groups[taxon.split('__')[0]]  # c__Clostridia

    with dbopen(db) as cursor:
        tablename = 'metadata'
        t = (taxon,)
        cursor.execute(
            f'SELECT accession FROM {tablename} WHERE gtdb_{rank}=?', t)
        l = [strip_name(i[0]) for i in cursor.fetchall()]  
        # ('GB_GCA_002409805.1',),
    return l


def index_model(names, models, norm='l2'):
    '''
    To normalize or not to normalize:
    
    - stats.stackexchange.com/questions/177905
    - stackoverflow.com/questions/36034454

    Usage:

    fp = f'{base}/models/{n}/nanotext_r89.model'
    model3 = load_embedding(fp)
    m3 = subtract_mean(model3)
    found, m, index = index_model(names, [m1, m2, m3], norm=norm)
    '''
    import faiss
    import numpy as np
    from sklearn.preprocessing import normalize
    
    from nanotext.io import eprint
    
    m = []
    found, notfound = [], 0

    # first take mean of vectors ...
    for i in names:
        model_vv = []
        try:
            for model in models:
                model_vv.append(model[i])
        except KeyError:
            notfound += 1
            continue
        
        sum_ = np.sum(model_vv, axis=0)/len(model_vv)
        found.append(i)
        m.append(sum_)
        # if only one model is present, this will return the original vector
    
    db = np.array(m, dtype='float32')
    dim = db.shape[1]  # dimensions
    
    # ... then normalize
    if not norm:
        index = faiss.IndexFlatL2(dim)
    elif norm == 'l2':
        index = faiss.IndexFlatIP(dim)
        db = normalize(db, norm=norm, axis=1)
        # the inner product IP of two unit length vectors = cosine similarity
    else:
        raise ValueError('This norm is not supported, abort!')

    index.add(db)
    if notfound > 0:
        fraction = round(notfound/len(names), 4)
        eprint(f'{notfound} entries ({fraction}) not found.')
    return found, db, index


def subtract_mean(model, names=None, dtype='float32'):
    '''Subtract mean vector from model and return vector collection (matrix)

    If no <names> (IDs) are provided, then all vectors in the model will be
    averaged and subtracted. Otherwise this operation is only performed on the
    subset defined by <names>.

    Suggestion from "All-but-the-top" paper (https://arxiv.org/abs/1702.01417).

    Usage:

    m_ = subtract_mean(model)  # m_ .. m minus
    '''
    import numpy as np
    from nanotext.io import eprint
    
    vv = []
    if not names:
        names = model.docvecs.index2entity
    
    notfound = 0
    for i in names:
        try:
            vv.append(model.docvecs[i])
        except KeyError:
            notfound += 1
            continue
    
    # eprint(f'{notfound} names {round(notfound/len(names), 4)} not found ...')

    m = np.array(vv, dtype=dtype)
    mu = np.mean(m, axis=0)
    return dict(zip(names, m-mu))


def query_model(v, index, names, topn=1, norm='l2'):
    '''k-nearest neighbor search'''
    import numpy as np
    from sklearn.preprocessing import normalize

    v = np.array([v], dtype='float32')
    if not norm:
        pass
    elif norm == 'l2':
        v = normalize(v, norm=norm, axis=1)
    else:
        raise ValueError('This norm is not supported, abort!')

    D, I = index.search(v, topn)
    
    result = {}
    for d, i in zip(D[0], I[0]):
        result[names[i]] = d

    return result


def get_taxa_from_names(db, names):
    '''Given a list of IDs return GTDB taxonomy

    The order of <names> is preserved, e.g. if they are sorted by distance.
    '''
    import pandas as pd

    from nanotext.utils import strip_name
    from nanotext.io import dbopen, eprint

    
    with dbopen(db) as cursor:
        # cannot use placeholders here
        # stackoverflow.com/questions/31277027
        
        # needs to be a tuple otherwise OperationalError: no such table ...
        n = tuple(names)
        # hack that covers case of only one query
        # problem is that tuple([1]) -> tuple(1,)
        # -- trailing comma causes error
        if len(n) == 1:
            n = n*2
        
        statement = f'SELECT accession_redux, gtdb_taxonomy FROM metadata WHERE accession_redux IN {n}'
        # statement = f'SELECT accession_redux, gtdb_taxonomy FROM metadata WHERE accession_redux = {n}'
        
        cursor.execute(statement)
        l = cursor.fetchall()
    
    # order of names is preserved, e.g. when ordered by distance
    taxa = {}
    for name, taxon in l:
        # 'd__Bacteria;p__Cyanobacteriota;c__Cyanobacteriia;[...]'
        taxa[name] = [name] + [j.split('__')[1] for j in taxon.split(';')]

    found, notfound = [], []
    for i in names:
        taxon = taxa.get(i, None)
        if taxon:
            found.append(taxon)
        else:
            notfound.append(i)

    # df = pd.DataFrame.from_records([taxa[i] for i in names])  
    # preserves order
    if notfound:
        eprint(f'Did not find {len(notfound)} out of {len(names)} queries:')
        eprint(notfound)
    df = pd.DataFrame.from_records(found)
    if df.empty:  # no query found -- len(pd.DataFrame.from_records([]))
        return df
    else:
        df.columns = 'name domain phylum class order family genus species'.split()
        return df

