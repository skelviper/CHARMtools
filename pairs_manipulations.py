# funcitons for manipulate 4DN pairs format files
from . import CHARMio
import glob
import pandas as pd

from joblib import Parallel, delayed
import multiprocessing
import numpy as np
import matplotlib.pyplot as plt
import cooltools.lib.plotting

# IO
# CHARMio.parse_pairs
# CHARMio.write_pairs

def sort_pairs_in_parallel(combined_pairs,num_cores = None,keylist = ["chr1","pos1","chr2","pos2"]):
    if(num_cores is None):
        num_cores = multiprocessing.cpu_count() - 5
    num_partitions = num_cores * 2
    partitions = np.array_split(combined_pairs, num_partitions)
    results = Parallel(n_jobs=num_cores)(delayed(_sort_partition)(partition,keylist=keylist) for partition in partitions)
    sorted_pairs = pd.concat(results)
    sorted_pairs = sorted_pairs.sort_values(keylist)
    sorted_pairs.reset_index(drop=True, inplace=True)
    return sorted_pairs

def _sort_partition(partition,keylist = ["chr1","pos1","chr2","pos2"]):
    return partition.sort_values(keylist)

def parse_and_combine_pairs(file_paths, num_cores = None) -> pd.DataFrame:
    if isinstance(file_paths, str):
            file_paths = glob.glob(file_paths)
    if(num_cores is None):
        num_cores = multiprocessing.cpu_count() - 5
    results = Parallel(n_jobs=num_cores)(delayed(CHARMio.parse_pairs)(i) for i in file_paths)
    combined_pairs = pd.concat(results)
    combined_pairs = sort_pairs_in_parallel(combined_pairs,num_cores = num_cores,keylist = ["chr1","pos1","chr2","pos2"])
    combined_pairs.reset_index(drop=True,inplace=True)
    return combined_pairs

def pairs_to_bedgraph(pairs, binsize=500, num_cores=None,normalize_method = "cpm"):
    left = pairs[["chr1","pos1"]]
    left.columns = ["chr","pos"]
    right = pairs[["chr2","pos2"]]
    right.columns = ["chr","pos"]
    pairs = pd.concat([left,right])
    bined_pairs = pairs.assign(bin = pairs.pos // binsize * binsize).groupby(["chr","bin"]).aggregate("count").reset_index()
    bined_pairs.columns = ["chr","start","count"]
    bined_pairs = bined_pairs.assign(end = bined_pairs.start + binsize)
    bined_pairs = bined_pairs[["chr","start","end","count"]]
    if(normalize_method == "cpm"):
        bined_pairs["count"] = bined_pairs["count"] / (bined_pairs["count"].sum() / 1e6)
    return bined_pairs

def _generate_genomic_intervals(chrom_sizes:str,resolutions=500): 
    #return a bed-like dataframe
    chrom_sizes = pd.read_csv(chrom_sizes,sep="\t",header=None)
    chrom_sizes.columns = ["chr","size"]
    #keep writing
    intervals = []
    for index, row in chrom_sizes.iterrows():
        chr_name = row["chr"]
        chr_size = row["size"]
        for i in range(0,chr_size,resolutions):
            intervals.append([chr_name,i,i+resolutions])
    intervals_df = pd.DataFrame(intervals,columns=["chr","start","end"])
    return intervals_df

def _bedgraph_to_dense(bedgraph,chrom_sizes:str,resolution = 500):
    genome_intervals = _generate_genomic_intervals(chrom_sizes,resolution)
    dense = pd.merge(genome_intervals,bedgraph,how="left",on=["chr","start","end"])
    dense = dense.fillna(0)
    dense = dense[["chr","start","end","count"]]
    return dense

# def _pileup_partition(bed,bedgraph,flank,resolution):
#     bed.columns = ["chr","start","end"]
#     bed = bed.assign(midpoint = (bed["start"] + bed["end"]) // 2 // resolution * resolution)
#     bed = bed.assign(start = bed["midpoint"] - flank)
#     bed = bed.assign(end = bed["midpoint"] + flank + resolution)
#     bed = bed[["chr","start","end"]]

#     bedgraph_dict = bedgraph.groupby("chr")
#     result = []
#     for index, row in bed.iterrows():
#         chr_name = row["chr"]
#         start = row["start"]
#         end = row["end"]
#         if chr_name in bedgraph_dict.groups:
#             chr_group = bedgraph_dict.get_group(chr_name)
#             chr_group = chr_group[(chr_group["start"] >= start) & (chr_group["end"] <= end)]
#             #print(chr_group["count"].values)
#             result.append(chr_group["count"].values)
#         else:
#             result.append(np.zeros((int((end-start)/resolution),)))
#     return np.array(result)
        
# def pileup_bedgraph_on_bed(bed, bedgraph, flank=20000, resolution = 500,num_cores = None):
#     bed = bed.iloc[:,0:3]
#     bed.columns = ["chr","start","end"]
#     bed = bed.assign(midpoint = (bed["start"] + bed["end"]) // 2 // resolution * resolution)
#     bed = bed.assign(start = bed["midpoint"] - flank)
#     bed = bed.assign(end = bed["midpoint"] + flank + resolution)
#     bed = bed[["chr","start","end"]]

#     if(num_cores is None):
#         num_cores = multiprocessing.cpu_count() - 5
#     num_partitions = num_cores * 2
#     partitions = np.array_split(bed, num_partitions)
#     results = Parallel(n_jobs=num_cores)(delayed(_pileup_partition)(partition,bedgraph,flank,resolution) for partition in partitions)
#     results = np.array(results)
#     results = results.reshape(-1,results.shape[2])
#     return results

def pileup_bedgraph_on_bed(bed, bedgraph, flank=20000, resolution = 500):
    bed = bed.iloc[:,0:3]
    bed.columns = ["chr","start","end"]
    bed = bed.assign(midpoint = (bed["start"] + bed["end"]) // 2 // resolution * resolution)
    bed = bed.assign(start = bed["midpoint"] - flank)
    bed = bed.assign(end = bed["midpoint"] + flank + resolution)
    bed = bed[["chr","start","end"]]

    bedgraph_dict = bedgraph.groupby("chr")
    result = []
    for index, row in bed.iterrows():
        chr_name = row["chr"]
        start = row["start"]
        end = row["end"]
        if chr_name in bedgraph_dict.groups:
            chr_group = bedgraph_dict.get_group(chr_name)
            chr_group = chr_group[(chr_group["start"] >= start) & (chr_group["end"] <= end)]
            #print(chr_group["count"].values)
            result.append(chr_group["count"].values)
        else:
            result.append(np.zeros((int((end-start)/resolution),)))
    return np.array(result)


def plot_enrichment(matrix):
    fig, axs = plt.subplots(2, 1, figsize=(4, 4), dpi=150, gridspec_kw={'height_ratios': [1, 2]},sharex=True,constrained_layout=True)
    plt.sca(axs[0])
    axs[0]= plt.plot(np.mean(matrix,axis=0))
    plt.sca(axs[1])
    matrix = matrix[np.argsort(np.sum(matrix,axis=1))[::-1],:]
    ax = axs[1].imshow(matrix, cmap="fall", aspect="auto",vmin = 0,vmax= 0.5)

    plt.colorbar(ax)
    plt.yticks([])
    plt.xticks([0,40,80],["-20kb","center","+20kb"])

    plt.show()
    return(fig)

def pairs_enrichment_pipeline(pairspath:str,regions:str,flank:int,resolution:int,num_cores:int,chrom_sizes:str):
    pairs = parse_and_combine_pairs(pairspath,num_cores=num_cores)
    bed = pd.read_csv(regions,sep="\t",header=None).sample(10000)
    bedgraph = pairs_to_bedgraph(pairs,binsize=resolution,num_cores=num_cores,normalize_method="cpm")
    bedgraph = _bedgraph_to_dense(bedgraph,chrom_sizes=chrom_sizes,resolution=resolution)
    pileup = pileup_bedgraph_on_bed(bed,bedgraph,flank=flank,resolution=resolution)#num_cores=num_cores)
    fig = plot_enrichment(pileup)
    return fig
