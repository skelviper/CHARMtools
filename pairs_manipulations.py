# funcitons for manipulate 4DN pairs format files
from . import CHARMio
import glob
import pandas as pd

from joblib import Parallel, delayed
import multiprocessing
import numpy as np
import matplotlib.pyplot as plt
import tqdm
import cooltools.lib.plotting
import gzip

def _bedgraph_to_dense(bedgraph,chrom_sizes:str,resolution = 500):
    genome_intervals = _generate_genomic_intervals(chrom_sizes,resolution)
    dense = pd.merge(genome_intervals,bedgraph,how="left",on=["chr","start","end"])
    dense = dense.fillna(0)
    dense = dense[["chr","start","end","count"]]
    return dense

def bed_to_bedgraph(bed,chrom_sizes,binsize = 100):
    bed.columns = ["chr","start","end"]
    bed=bed.assign(pos = lambda x: (x["start"] + x["end"]) // 2)
    bed = bed[["chr","pos"]]
    bined_bed = bed.assign(bin = bed.pos // binsize * binsize).groupby(["chr","bin"]).aggregate("count").reset_index()
    bined_bed.columns = ["chr","start","count"]
    bined_bed = bined_bed.assign(end = bined_bed.start + binsize)
    bined_bed = bined_bed[["chr","start","end","count"]]
    dense = _bedgraph_to_dense(bined_bed,chrom_sizes=chrom_sizes,resolution=binsize)
    return dense
    
def pairs_to_bedgraph(pairspath,chrom_sizes,binsize=500, num_cores=1,chunksize = None,normalize_method = "cpm"):
    combined_bedgraph = pd.DataFrame()
    if chunksize is None:
        with gzip.open(pairspath, 'rt') as f:
            #print("here"
            results = _process_chunk(pd.read_csv(f, sep='\t',comment="#",header=None),binsize,chrom_sizes)
            #print(results)
            combined_bedgraph = results
    else:
        with gzip.open(pairspath, 'rt') as f:
            results = Parallel(n_jobs=num_cores)(delayed(_process_chunk)(chunk,binsize,chrom_sizes) for chunk in tqdm.tqdm(pd.read_csv(f, sep='\t', chunksize=chunksize)))
            combined_bedgraph = pd.concat(results).groupby(["chr","start","end"]).agg({"count":"sum"}).reset_index()

    if(normalize_method == "cpm"):
        combined_bedgraph["count"] = combined_bedgraph["count"] / (combined_bedgraph["count"].sum() / 1e6)
    return combined_bedgraph

def multiple_pairs_to_bedgraph(pairs_paths,chrom_sizes,binsize=500, num_cores=10,normalize_method = "cpm"):
    combined_bedgraph = pd.DataFrame()
    #pairs_paths = glob.glob(pairs_paths)
    results = Parallel(n_jobs=num_cores)(delayed(pairs_to_bedgraph)(i,chrom_sizes,binsize=binsize,chunksize = None,num_cores=10,normalize_method=None) for i in pairs_paths)
    #combined_bedgraph = pd.concat(results).groupby(["chr","start","end"]).agg({"count":"sum"}).reset_index()
    summed_count = np.sum([table["count"].values for table in results],axis=0)
    combined_bedgraph = results[0][["chr","start","end"]]
    combined_bedgraph = combined_bedgraph.assign(count = summed_count)
    if(normalize_method == "cpm"):
        combined_bedgraph["count"] = combined_bedgraph["count"] / (combined_bedgraph["count"].sum() / 1e6)
    return combined_bedgraph

def _process_chunk(chunk,binsize,chrom_sizes):
    chunk = chunk.iloc[:,0:5]
    chunk.columns = ["readid","chr1","pos1","chr2","pos2"]
    left = chunk[["chr1","pos1"]]
    left.columns = ["chr","pos"]
    right = chunk[["chr2","pos2"]]
    right.columns = ["chr","pos"]
    pairs = pd.concat([left,right])
    bined_pairs = pairs.assign(bin = pairs.pos // binsize * binsize).groupby(["chr","bin"]).aggregate("count").reset_index()
    bined_pairs.columns = ["chr","start","count"]
    bined_pairs = bined_pairs.assign(end = bined_pairs.start + binsize)
    bined_pairs = bined_pairs[["chr","start","end","count"]]
    dense = _bedgraph_to_dense(bined_pairs,chrom_sizes=chrom_sizes,resolution=binsize)
    return dense
    
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

# def pairs_to_bedgraph(pairs, binsize=500, num_cores=None,normalize_method = "cpm"):
#     left = pairs[["chr1","pos1"]]
#     left.columns = ["chr","pos"]
#     right = pairs[["chr2","pos2"]]
#     right.columns = ["chr","pos"]
#     pairs = pd.concat([left,right])
#     bined_pairs = pairs.assign(bin = pairs.pos // binsize * binsize).groupby(["chr","bin"]).aggregate("count").reset_index()
#     bined_pairs.columns = ["chr","start","count"]
#     bined_pairs = bined_pairs.assign(end = bined_pairs.start + binsize)
#     bined_pairs = bined_pairs[["chr","start","end","count"]]
#     if(normalize_method == "cpm"):
#         bined_pairs["count"] = bined_pairs["count"] / (bined_pairs["count"].sum() / 1e6)
#     return bined_pairs

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
            if len(chr_group["count"].values) != 2 * flank // resolution + 1:
                continue
            result.append(chr_group["count"].values)
        else:
            result.append(np.zeros((int((end-start)/resolution),)))
    return np.array(result)


def plot_enrichment(matrix,flank):
    fig, axs = plt.subplots(2, 1, figsize=(4, 4), dpi=150, gridspec_kw={'height_ratios': [1, 2]},sharex=True,constrained_layout=True)
    plt.sca(axs[0])
    axs[0]= plt.plot(np.mean(matrix,axis=0))
    plt.sca(axs[1])
    matrix = matrix[np.argsort(np.sum(matrix,axis=1))[::-1],:]
    #ax = axs[1].imshow(matrix, cmap="fall", aspect="auto",vmin = 0,vmax= 0.5)
    # 95% data as vmax and vmin
    vmax = np.percentile(matrix, 95)
    vmin = np.percentile(matrix, 5)
    ax = axs[1].imshow(matrix, cmap="fall", aspect="auto",vmin = vmin,vmax= vmax)

    plt.colorbar(ax)
    plt.yticks([])
    plt.xticks([0,(matrix.shape[1]-1)/2,matrix.shape[1]-1],["-"+str(flank/1000)+"k","center","+"+str(flank/1000)+"k"])

    plt.show()
    return(fig)

def pairs_enrichment_pipeline(pairspath:str,regions:str,flank:int,resolution:int,num_cores:int,chrom_sizes:str):
    import time 
    start = time.time()
    bed = pd.read_csv(regions,sep="\t",header=None)
    bedgraph = multiple_pairs_to_bedgraph(pairspath,chrom_sizes = chrom_sizes,binsize=resolution,num_cores=num_cores,normalize_method="cpm")
    #print("bedgraph generated at {} seconds".format(time.time() - start))
    pileup = pileup_bedgraph_on_bed(bed,bedgraph,flank=flank,resolution=resolution)#num_cores=num_cores)
    #return(pileup)
    # print(pileup)
    #print("pileup generated at {} seconds".format(time.time() - start))
    fig = plot_enrichment(pileup,flank=flank)
    #print("plot generated at {} seconds".format(time.time() - start))
    print("Enrichment at peak: " + str(np.mean(pileup[:,flank//resolution]) * 2 / np.mean((pileup[:,0] + pileup[:,-1]))))
    return bedgraph,fig

def remove_regions_from_pairs(pairs_path_in,pairs_path_out,regions):
    pairs = CHARMio.parse_pairs(pairs_path_in)
    pairs["readID"] = "read" + pairs.index.astype(str)
    bedpe = pairs[["readID","chr1","pos1","chr2","pos2"]]
    bedpe = bedpe.assign(start1=bedpe.pos1, end1=bedpe.pos1+1, start2=bedpe.pos2, end2=bedpe.pos2+1)[["chr1","start1","end1","chr2","start2","end2","readID"]]
    bedpe.columns = ["chrom1","start1","end1","chrom2","start2","end2","readID"]

    bedpe_bedtool = BedTool.from_dataframe(bedpe)
    bed_bedtool = BedTool.from_dataframe(regions)

    bedpe_keep = BedTool.pair_to_bed(bedpe_bedtool, bed_bedtool, type = "neither").to_dataframe()[["thickStart"]]
    bedpe_keep.columns = ["readID"]

    pairs_write = pd.merge(pairs,bedpe_keep,how="inner",on="readID")
    pairs_write.attrs["comments"] = pairs.attrs["comments"]
    print("Percent of pairs removed: %.4f" % (1-len(pairs_write)/len(pairs)))
    CHARMio.write_pairs(pairs_write,pairs_path_out)


import multiprocessing
from tqdm import tqdm

def process_bed_chunk(bed_chunk, bedgraph, flank, resolution):
    bedgraph_dict = bedgraph.groupby("chr")
    result = []

    for index, row in bed_chunk.iterrows():
        chr_name = row["chr"]
        start = row["start"]
        end = row["end"]

        if chr_name in bedgraph_dict.groups:
            chr_group = bedgraph_dict.get_group(chr_name)
            chr_group = chr_group[(chr_group["start"] >= start) & (chr_group["end"] <= end)]
            if len(chr_group["count"].values) != 2 * flank // resolution + 1:
                result.append(None)
            else:
                result.append(chr_group["count"].values)
        else:
            result.append(np.zeros((int((end-start)/resolution),)))

    return result

def pileup_bedgraph_on_bed_parallel(bed, bedgraph, flank=20000, resolution=500, chunk_size=1000):
    bed = bed.iloc[:, 0:3]
    bed.columns = ["chr", "start", "end"]
    bed = bed.assign(midpoint=(bed["start"] + bed["end"]) // 2 // resolution * resolution)
    bed = bed.assign(start=bed["midpoint"] - flank)
    bed = bed.assign(end=bed["midpoint"] + flank + resolution)
    bed = bed[["chr", "start", "end"]]

    result = []
    chunks = [bed_chunk for _, bed_chunk in bed.groupby(np.arange(len(bed)) // chunk_size)]

    with multiprocessing.Pool(20) as pool, tqdm(total=len(chunks)) as pbar:
        results = []

        def update(*a):
            pbar.update()

        for chunk in chunks:
            results.append(pool.apply_async(process_bed_chunk, (chunk, bedgraph, flank, resolution), callback=update))

        for res in results:
            result.extend(res.get())

    return np.array(result)

