import numpy as np
import pandas as pd
import gzip
import os
from pkgutil import get_data
from io import StringIO
import ref
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
    name_array = "readID chr1 pos1 chr2 pos2 strand1 strand2 phase0 phase1 phase_prob00 phase_prob01 phase_prob10 phase_prob11".split()
    #read comment line
    with gzip.open(filename,"rt") as f:
        comments = []
        for line in f.readlines():
            if line[0] != "#":
                break
            comments.append(line)
    #read table format data
    pairs = pd.read_table(filename, header=None, comment="#")
    pairs.attrs["comments"] = comments
    pairs.attrs["name"], _ = divide_name(filename) # get real sample name
    #assign column names
    pairs.columns = name_array[0:pairs.shape[1]]
    #sys.stderr.write("pairs_parser: %s parsed \n" % filename)
    return pairs

def check_index_binsize(s: pd.DataFrame) -> int:
    """
    Check if the index of s, considering its first level grouping, has sorted and consistent binsizes.
    
    The function assumes that the second level of the MultiIndex is numerical (e.g., genomic coordinates or timestamps).
    It calculates the binsize for each group defined by the first level and returns a series with the binsizes.
    If any inconsistency in binsize within a group is found, raises ValueError.
    Input:
        s: A pandas DataFrame with a MultiIndex where the first level represents chromosome names and the second level
              represents positions or other values that should have a consistent difference (binsize).
    Output:
        binsize
    """

    # Ensure the index is sorted
    s = s.sort_index()
    # Calculate binsizes per group based on the first level
    # last 2 element is not used, because the last one is NaN and in some situations the second last one is partial binsize
    # negative period is used to ensure first element is not NaN
    result_dfs = []
    for name, group in s.groupby(level=0):
        new_df = pd.Series(
            -group.index.get_level_values(1).diff(-1),
            index = group.index
            ).rename("binsizes").iloc[:-2]
        result_dfs.append(new_df)
    binsizes = pd.concat(result_dfs, axis=0).dropna().astype(int)
    if binsizes.empty:
        print("Warning: No binsize found.")
        # just use the first binsize
        binsize = -s.index.get_level_values(1).diff()[0]
    elif len(binsizes.unique()) > 1:
        print("Warning: Inconsistent binsizes found in the input file %s" % binsizes.unique())
        binsize = binsizes.dropna().unique()[0]
    else:
        binsize = binsizes.dropna().unique()[0]
    return binsize

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

def parse_3dg(file:str, sorting=False, s2m=False)->pd.DataFrame:
    """
    Read in hickit 3dg file(or the .xyz file)
    Norm chr name alias
    Read into dataframe.attrs if has comments, treat last comment line as backbone_unit
    Input:
        filename: file path
        sorting: whether to sort chromosome and positions
        s2m: whether to use mid point of bin as position
    Output:
        3 col dataframe with 2-level multindex: (chrom, pos) x, y, z
    """

    ## read comments
    comments = read_comments(file, next_line=False)
    ## read real positions
    s = pd.read_table(file, 
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

def parquet2pairs(parquet_path,pairs_path):
    cell = pd.read_parquet(parquet_path)
    temp = cell#[["read_name","align1_fragment_start","align1_fragment_end","align1_chrom","align1_haplotype","align2_fragment_start","align2_fragment_end","align2_chrom","align2_haplotype"]]
    temp = temp.assign(pos1 = (temp.align1_fragment_start + temp.align1_fragment_end) // 2) 
    temp = temp.assign(pos2 = (temp.align2_fragment_start + temp.align2_fragment_end) // 2)
    temp = temp[["read_name","align1_chrom","pos1","align2_chrom","pos2","align1_strand","align2_strand","align1_haplotype","align2_haplotype"]]
    temp["align1_strand"] = temp["align1_strand"].replace(True, '+')
    temp["align1_strand"] = temp["align1_strand"].replace(False, '-')
    temp["align2_strand"] = temp["align2_strand"].replace(True, '+')
    temp["align2_strand"] = temp["align2_strand"].replace(False, '-')
    temp['align1_haplotype'] = temp['align1_haplotype'].replace(-1, '.')
    temp['align2_haplotype'] = temp['align2_haplotype'].replace(-1, '.')
    temp['align1_haplotype'] = temp['align1_haplotype'].replace(1, '0')
    temp['align2_haplotype'] = temp['align2_haplotype'].replace(2, '1')
    temp.columns = ["readID","chr1","pos1","chr2","pos2","strand1","strand2","phase0","phase1"]

    # only keep autosomes and chrX and chrY
    temp = temp[temp.chr1.str.contains("chr[0-9]+|chrX|chrY")]
    temp = temp[temp.chr2.str.contains("chr[0-9]+|chrX|chrY")]
    unique_chroms = set(list(temp['chr1'].unique()) + list(temp['chr2'].unique())) 
    unique_chroms = sorted(unique_chroms, key=lambda x: (int(x[3:]) if x[3:].isdigit() else float('inf'), x))
    temp["chr1"] = pd.Categorical(temp["chr1"], categories=unique_chroms)
    temp["chr2"] = pd.Categorical(temp["chr2"], categories=unique_chroms)
    temp = temp.sort_values(["chr1","chr2","pos1","pos2"])

        #write tsv with header
        # for GRCh38
    header = """## pairs format v1.0  
#sorted: chr1-chr2-pos1-pos2 
#shape: upper triangle
#chromosome: chr1 248956422
#chromosome: chr2 242193529
#chromosome: chr3 198295559
#chromosome: chr4 190214555
#chromosome: chr5 181538259
#chromosome: chr6 170805979
#chromosome: chr7 159345973
#chromosome: chr8 145138636
#chromosome: chr9 138394717
#chromosome: chr10 133797422
#chromosome: chr11 135086622
#chromosome: chr12 133275309
#chromosome: chr13 114364328
#chromosome: chr14 107043718
#chromosome: chr15 101991189
#chromosome: chr16 90338345
#chromosome: chr17 83257441
#chromosome: chr18 80373285
#chromosome: chr19 58617616
#chromosome: chr20 64444167
#chromosome: chr21 46709983
#chromosome: chr22 50818468
#chromosome: chrX 156040895
#chromosome: chrY 57227415
"""

    temp.attrs['comments'] = [header,""]


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

    #temp.to_csv("test.pairs.gz",sep="\t",index=False,header=None)
    write_pairs(temp, pairs_path)
    return None

def pairs_describe(pairs):
    """%
    Generate basic statistics of pairs including contacts number, inter_pairs_ratio and non_contact_ratio.
    """
    contacts_number = pairs.shape[0]
    inter_pairs_ratio = pairs[pairs['chr1'] != pairs['chr2']].shape[0] / contacts_number
    non_contact_ratio = pairs[(pairs['chr1'] == pairs['chr2']) & (pairs["pos2"] - pairs["pos1"] < 1000)].shape[0] / contacts_number
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