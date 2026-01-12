import sys
import time
from concurrent import futures
from functools import partial
import pandas as pd
import re
import numpy as np


from ..utils.CHARMio import parse_pairs, parse_gtf, write_pairs

def cli(args)->int:
    filename, gtf_file, out_name, thread = \
        args.filename[0], args.gtf_filename, args.out_name, args.num_thread
    # parsing ref gtf and pairs file
    pairs = parse_pairs(filename)
    # build in-mem exon index
    gtf = parse_gtf(gtf_file)
    try:
        typegenome = re.search(r'rabbit', gtf_file).group(0)
    except:
        typegenome = "else"
    if typegenome == "rabbit":
        ref = build_in_memory_index(get_exon_rab(gtf))
    else:
        ref = build_in_memory_index(get_exon(gtf))
    #do search
    cleaned = clean_splicing(pairs, ref, thread)
    write_pairs(cleaned, out_name)
def get_exon(gtf:pd.DataFrame) -> pd.DataFrame:
    # extract exon-gene_name from gtf table
    relevant = gtf.query('feature == "exon" & source == "HAVANA"') #using HAVANA only
    gene_id = relevant["group"].str.extract('gene_id "(ENSG[0-9]{11}.[0-9])";') #extract gene name from group
    gene_id.columns = ["gene_id"] # extract returns dataframe rather than series
    # don't mind strand
    return pd.concat([relevant.drop(["group","feature","source","score","strand","frame"],axis=1),gene_id],axis=1)
def build_in_memory_index(exons:pd.DataFrame) -> dict:
    # split by chr and using IntervalIndex to enable searching
    ref_dict = {key : value for key, value in exons.groupby("seqname")}
    # build index by chromosome
    for chromosome in ref_dict:
        # using start, end attrs as index
        by_chr_table = ref_dict[chromosome]
        bed_tuple = by_chr_table.set_index(['start','end']).index 
        bed_interval_index = pd.IntervalIndex.from_tuples(bed_tuple)
        by_chr_table.index = bed_interval_index
        ref_dict[chromosome] = by_chr_table.drop(["start","end","seqname"],axis=1)
    sys.stderr.write("CHARMtools::clean_splicing: index done.\n")
    return ref_dict


def _build_interval_cache(ref_dict: dict) -> dict:
    """
    Prepare start positions and prefix max ends per chromosome for fast interval lookup.
    """
    cache = {}
    for chromosome, table in ref_dict.items():
        if table.empty:
            continue
        starts = table.index.left.to_numpy()
        ends = table.index.right.to_numpy()
        order = np.argsort(starts, kind="mergesort")
        starts = starts[order]
        ends = ends[order]
        prefix_end = np.maximum.accumulate(ends)
        cache[chromosome] = (starts, prefix_end)
    return cache


def _search_chromosome_cached(item: tuple, cache: dict) -> pd.DataFrame:
    """
    Find contacts whose two legs fall into the same exon on a chromosome.
    Interval membership matches the original IntervalIndex.contains logic
    (open on the left, closed on the right).
    """
    chromosome, contacts = item
    if contacts.empty:
        return contacts
    interval_cache = cache.get(chromosome)
    if interval_cache is None:
        return contacts.iloc[0:0]
    starts, prefix_end = interval_cache
    pos1 = contacts["pos1"].to_numpy()
    pos2 = contacts["pos2"].to_numpy()
    pos_min = np.minimum(pos1, pos2)
    pos_max = np.maximum(pos1, pos2)
    idx = np.searchsorted(starts, pos_min, side="left") - 1
    valid = idx >= 0
    if not np.any(valid):
        return contacts.iloc[0:0]
    safe_idx = np.where(valid, idx, 0)
    start_at_idx = starts[safe_idx]
    end_at_idx = prefix_end[safe_idx]
    hits = valid & (pos_min > start_at_idx) & (pos_max <= end_at_idx)
    return contacts.loc[hits]


def clean_splicing(pairs:pd.DataFrame, ref:dict, thread:int)->pd.DataFrame:
    '''
    clean contacts from splicing
    '''
    t0 = time.time()
    #intra = pairs.query('chrom1 == chrom2') # only search for intra
    intra = pairs.loc[pairs['chrom1'] == pairs['chrom2']]
    interval_cache = _build_interval_cache(ref)
    input_data = [(chromosome, per_chr_contacts) for chromosome, per_chr_contacts in intra.groupby("chrom1", sort=False)]
    sys.stderr.write("CHARMtools::clean_splicing: input parsed, search in %d thread\n" % thread)
    working_func = partial(_search_chromosome_cached, cache=interval_cache)
    if thread > 1:
        with futures.ThreadPoolExecutor(thread) as pool:
            result_groups = list(pool.map(working_func, input_data))
    else:
        result_groups = [working_func(item) for item in input_data]
    result = pd.concat(result_groups, axis=0) if len(result_groups) > 0 else intra.iloc[0:0]
    cleaned = pairs.drop(result.index, axis=0) # clean contacts
    cleaned.attrs.update(pairs.attrs)
    print("clean_splicing: %d contacts removed in %s\n" %(len(result), pairs.attrs["name"]) )
    sys.stderr.write("clean_splicing: finished in %.2fs\n" % (time.time() - t0))
    return cleaned

def get_exon_rab(gtf:pd.DataFrame) -> pd.DataFrame:
    # extract exon-gene_name from gtf table
    relevant = gtf.query('feature == "exon"') #using HAVANA only
    gene_id = relevant["group"].str.extract('gene_id "([A-Za-z0-9_-]+)";') #extract gene name from group
    gene_id.columns = ["gene_id"] # extract returns dataframe rather than series
    # don't mind strand
    return pd.concat([relevant.drop(["group","feature","source","score","strand","frame"],axis=1),gene_id],axis=1)

if __name__ == "__main__":
    # infile = "/shareb/zliu/project/202212/hires_pipe/result/cleaned_pairs/c12/OrgfE951001.pairs.gz"
    # gtf_file = "/share/Data/public/ref_genome/mouse_ref/M23/raw_data/annotation.gtf"
    # outfile = "/shareb/zliu/project/202212/hires_pipe/result/cleaned_pairs/c123/OrgfE951001.pairs.gz"
    infile = sys.argv[1]
    gtf_file = sys.argv[2]
    outfile = sys.argv[3]
    typegenome = sys.argv[4]
    thread = 6

    pairs = parse_pairs(infile)
    # build in-mem exon index
    gtf = parse_gtf(gtf_file)
    if typegenome == "rab":
        ref = build_in_memory_index(get_exon_rab(gtf))
    else:
        ref = build_in_memory_index(get_exon(gtf))
    # do search
    cleaned = clean_splicing(pairs, ref, thread)
    write_pairs(cleaned, outfile)
