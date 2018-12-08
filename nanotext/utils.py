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
    '''
    *contig, orf = s.split('_')
    return '_'.join(contig), int(orf)


def create_domain_sequence(domains, keep_unknown=True, fmt_fn=lambda x: x):
    '''
    We can pass a <fmt_fn> to reformat domain names (e.g. turn PF07005.14
    into PF07005). By default, the identity function is used, which just
    returns the input as is.

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





