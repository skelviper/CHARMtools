import numpy as np
import pandas as pd
import gzip
import os
from pkgutil import get_data
from io import StringIO
# import ref
from .. import ref
from functools import partial

import re
import h5py
import scipy

# Functions for CHARM/HiRES pipeline 
def divide_name(filename):
    #home-made os.path.splitext, for it can't handle "name.a.b.c" properly
    basename = os.path.basename(filename)
    parts = basename.split(".") #split return >= 1 length list
    if len(parts) == 1:
        return parts[0], ""
    else:
        return parts[0], "."+".".join(parts[1:])

def parse_pairs(filename:str)->"Cell":
    '''
    read from 4DN's standard .pairs format
    compatible with all hickit originated pairs-like format 
    '''
    #comment lines are stored in dataframe.attrs["comment"]
    name_array = "readID chrom1 pos1 chrom2 pos2 strand1 strand2 phase0 phase1 phase_prob00 phase_prob01 phase_prob10 phase_prob11".split()
    #read comment line
    with gzip.open(filename,"rt") as f:
        comments = []
        for line in f.readlines():
            if line[0] != "#":
                break
            comments.append(line)
    #read table format data
    pairs = pd.read_table(filename, header=None, comment="#",low_memory=False)
    pairs.attrs["comments"] = comments
    pairs.attrs["name"], _ = divide_name(filename) # get real sample name
    #assign column names
    pairs.columns = name_array[0:pairs.shape[1]]
    #sys.stderr.write("pairs_parser: %s parsed \n" % filename)
    return pairs

def check_index_binsize(df: pd.DataFrame) -> int:
    """
    Calculate the resolution of a DataFrame with a MultiIndex (chr, pos).
    Resolution is defined as the smallest non-zero difference between the 'pos' values for each 'chr'.
    Args:
        df: A pandas DataFrame with a MultiIndex where the first level is 'chr' and the second level is 'pos'.
    Returns:
        The resolution as an integer.
    """
    # Initialize a large number to find the minimum resolution
    min_resolution = float('inf')

    # Group by the chromosome level of the index
    for chromosome, group in df.groupby(level=0):
        # Calculate differences between adjacent positions
        pos_diffs = group.index.get_level_values(1).to_series().diff().dropna()

        # Exclude zero differences and find the minimum
        non_zero_diffs = pos_diffs[pos_diffs != 0]
        if not non_zero_diffs.empty:
            chrom_min_res = non_zero_diffs.min()
            # Update overall minimum resolution if this chromosome's min is smaller
            min_resolution = min(min_resolution, chrom_min_res)

    # If min_resolution is still infinity, there were no non-zero resolutions found
    if min_resolution == float('inf'):
        raise ValueError("No non-zero position differences found in the DataFrame.")

    return int(min_resolution)

def write_pairs(pairs:pd.DataFrame, out_name:str):
    '''
    write dataframe to tab delimited zipped file
    reserve comment lines, no dataframe index
    headers store in last comment line
    need to sort with upper triangle label
    '''
    #sys.stderr.write("write to %s\n" % out_name)
    with gzip.open(out_name,"wt") as f:
        pairs.attrs["comments"].pop()
        pairs.attrs["comments"].append("#columns:" + "\t".join(pairs.columns) + "\n")
        f.write("".join(pairs.attrs["comments"]))
        pairs.to_csv(f, sep="\t", header=False, index=False, mode="a")

def parse_3dg(filename:str,s2m=False,sorting=False)->pd.DataFrame:
    # read in hickit 3dg file(or the .xyz file)
    # norm chr name alias

    ## get alias file in package
    ## reference is a "data module" with its own __init__.py
    dat = get_data(ref.__name__, "chrom_alias.csv")
    dat_f = StringIO(dat.decode())
    norm_chr = fill_func_ref(
                    converter_template,
                    dat_f,
                    "alias")
    ## read comments
    with open(filename) as f:
        comments = []
        for line in f.readlines():
            if line[0] != "#":
                break
            comments.append(line)
    ## read real positions
    s = pd.read_table(filename, 
                      comment="#",header=None,
                     index_col=[0,1],
                     converters={0:norm_chr})
    s.columns = "x y z".split()
    s.index.names = ["chr","pos"]
    ## assume last comment is backbone_unit
    if len(comments) > 0:
        s.attrs["comments"] = comments
        backbone_unit = float(comments[-1].split(":")[1].strip())
        s.attrs["backbone_unit"] = backbone_unit    
    if sorting:
        s.sort_index(inplace=True)
    if s2m:
        binsize = check_index_binsize(s)
        s.index = pd.MultiIndex.from_arrays(
            [s.index.get_level_values(0), s.index.get_level_values(1) + binsize//2],
            names=["chr","pos"]
        )
    return s

def m2s_index(_3dg:pd.DataFrame)->pd.DataFrame:
    '''
    Convert from middle-as-pos to start-as-pos
    Input:
        a structure dataframe with 2-level multiindex
    Output:
        a structure dataframe with 2-level multiindex
    '''
    binsize = check_index_binsize(_3dg)
    _3dg.index = pd.MultiIndex.from_arrays(
        [_3dg.index.get_level_values(0), _3dg.index.get_level_values(1) - binsize//2],
        names=["chr","pos"]
    )
    if _3dg.index.get_level_values(1).min() < 0:
        raise ValueError("Negative position found in the dataframe, perhaps the dataframe has been converted to start-as-pos already.")
    return _3dg
def write_3dg(_3dg:pd.DataFrame, outname:str, m2s=False):
    if m2s:
        _3dg = m2s_index(_3dg)
    _3dg.to_csv(outname, sep="\t",header=None)
    return 0

#Functions for analysis
def getMatrixFromMCOOLs(filepath:str, genome_coord1:str,genome_coord2=None, resolution=40000, balance=False)->np.ndarray:
    """
    intput: mcool filepath ,
            genome_coord(e.g. 'chr1:35,900,000-40,900,000'), 
            resolution(should be included in mcoolfile)
    output: numpy 2d array
    """
    import cooler
    cool = filepath+"::/resolutions/"+str(resolution)

    c = cooler.Cooler(cool)
    if(genome_coord2 == None):
        genome_coord2 = genome_coord1
    matrix = c.matrix(balance=balance).fetch(genome_coord1,genome_coord2).astype("double")

    return matrix

def getMatrixFromCooler(filepath:str, genome_coord1:str,genome_coord2=None, resolution=40000, balance=False)->np.ndarray:
    """
    intput: cooler or mcool filepath, file type determined by last extension name.
            genome_coord(e.g. 'chr1:35,900,000-40,900,000'), 
            resolution(should be included in mcoolfile)
    output: numpy 2d array
    """
    import cooler
    
    if filepath.split(sep='.')[-1] == "mcool":
        return getMatrixFromMCOOLs(filepath, genome_coord1,genome_coord2, resolution, balance)
    
    c = cooler.Cooler(filepath)
    if(genome_coord2 == None):
        genome_coord2 = genome_coord1
    matrix = c.matrix(balance=balance).fetch(genome_coord1,genome_coord2).astype("double")

    return matrix

def cooltoolsGetObsExp(filepath:str, genome_coord1:str,genome_coord2=None, resolution=40000, balance=False)->np.ndarray:
    """
    input: cooler or mcool path
    output: observe / expected matrix as np.ndarry 
    """
    import cooltools
    pass
    # if filepath.split(sep='.')[-1] == "mcool":

def parse_gtf(filename:str) -> pd.DataFrame:
    # read gtf, get exons
    gencode = pd.read_table(filename, comment="#", header=None)
    gencode.columns="seqname source feature start end score strand frame group".split()
    return gencode

## norm name
def converter_template(c_in:str,ref_dict:pd.DataFrame):
    # a reat_table converter function
    #print(ref_dict)
    return ref_dict[c_in]
def fill_func_ref(template_func:callable, ref_file:str, index_col:str)->callable:
    # read in ref_file for template_fucn, generate new func
    # hope will boost new func's speed
    
    ## read in ref_file, get ref_dict in memory
    ref_df = pd.read_csv(ref_file, index_col=index_col)
    ref_dict = ref_df.iloc[:,0].to_dict()
    working_func = partial(template_func, ref_dict=ref_dict)
    return working_func

def pairs_describe(pairs):
    """%
    Generate basic statistics of pairs including contacts number, inter_pairs_ratio and non_contact_ratio.
    """
    contacts_number = pairs.shape[0]
    inter_pairs_ratio = pairs[pairs['chrom1'] != pairs['chrom2']].shape[0] / contacts_number
    non_contact_ratio = pairs[(pairs['chrom1'] == pairs['chrom2']) & (pairs["pos2"] - pairs["pos1"] < 1000)].shape[0] / contacts_number
    return [pairs.attrs["name"],contacts_number, inter_pairs_ratio, non_contact_ratio]


# Functions for simpleDiff-python
def load_sparse_form_h5(filename):
    with h5py.File(filename, 'r') as f:
        data = f['data'][:]
        indices = f['indices'][:]
        indptr = f['indptr'][:]
        shape = f['shape'][:]
        csr = scipy.sparse.csr_matrix((data, indices, indptr), shape=shape)
        return csr

def read_subregion_from_h5(filename, row_start, row_end, col_start, col_end):
    csr = load_sparse_form_h5(filename)
    subregion = csr[row_start:row_end, col_start:col_end].todense()
    subregion = subregion + subregion.T - np.diag(np.diag(subregion))
    return subregion

def read_mat_h5(filepath,genome_coord,resolution = 10000):
    split = re.split("[:-]",genome_coord.replace(",",""))
    if len(split) == 3:
        chrom,start,end = split
        start = int(start)
        end = int(end)
        mat = read_subregion_from_h5(filepath,start//resolution,end//resolution,start//resolution,end//resolution)
        mat = np.log2(mat+1)
    else:
        chrom = split[0] 
        csr = load_sparse_form_h5(filepath)
        mat = csr.todense()
        mat = np.log2(mat+1)
    return mat