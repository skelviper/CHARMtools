import math 
import gzip
import pandas as pd
import numpy as np
from sklearn.cluster import DBSCAN

# adapt from Picard's EstimateLibraryComplexity.java
def f(X, C, N):
    """
    Function representing the Lander-Waterman equation:
    C/X = 1 - exp(-N/X)
    Rearranged to:
    f(X) = C/X - (1 - exp(-N/X))
    """
    try:
        return C / X - (1 - math.exp(-N / X))
    except ZeroDivisionError:
        return float('inf')

def estimate_library_size(read_pairs, unique_read_pairs):
    """
    Estimates the size of a library based on the number of paired end molecules observed
    and the number of unique pairs observed.

    :param read_pairs: total number of read pairs (N)
    :param unique_read_pairs: number of distinct fragments observed in read pairs (C)
    :return: estimated number of distinct molecules in the library (X) or None if invalid input
    """
    read_pair_duplicates = read_pairs - unique_read_pairs
    try:
        m = 1.0
        M = 100.0

        # Check initial condition for the bounds
        if unique_read_pairs >= read_pairs or f(m * unique_read_pairs, unique_read_pairs, read_pairs) < 0:
            raise ValueError("Invalid values for pairs and unique pairs: {}, {}".format(read_pairs, unique_read_pairs))

        # Find value of M, large enough to act as other side for bisection method
        while f(M * unique_read_pairs, unique_read_pairs, read_pairs) > 0:
            M *= 10.0

        # Use bisection method (no more than 40 times) to find solution
        for _ in range(40):
            r = (m + M) / 2.0
            u = f(r * unique_read_pairs, unique_read_pairs, read_pairs)
            if u == 0:
                break
            elif u > 0:
                m = r
            else:
                M = r

        return int(unique_read_pairs * (m + M) / 2.0)
    except:
        return None
    
def seg2pairs(file_path,mapq_threshold):
    # read
    with gzip.open(file_path, 'rt') as file:
        rows = []
        for line in file:

            fields = line.strip().split('\t')
            if len(fields) < 2:
                continue  
            
            # second column
            second_col = fields[1].split('!')
            chr1, start1, end1 = second_col[:3]
            
            # last column
            last_col = fields[-1].split('!')
            chr2, start2, end2 = last_col[:3]
            
            # mapq
            mapq1 = int(second_col[5])
            mapq2 = int(last_col[5])
            
            rows.append({
                'chrom1': chr1, 'start1': start1, 'end1': end1,
                'chrom2': chr2, 'start2': start2, 'end2': end2,
                'mapq1': mapq1, 'mapq2': mapq2
            })
            
    columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'mapq1', 'mapq2']
    df = pd.DataFrame(rows, columns=columns)
    df=df.query('mapq1 > @mapq_threshold and mapq2 > @mapq_threshold').copy()

    df["start1"] = df["start1"].astype(int)
    df["end1"] = df["end1"].astype(int)
    df["start2"] = df["start2"].astype(int)
    df["end2"] = df["end2"].astype(int)

    df["pos1"] = (df["start1"] + df["end1"]) // 2
    df["pos2"] = (df["start2"] + df["end2"]) // 2

    df=df.sort_values(by=['chrom1', 'pos1', 'chrom2', 'pos2'])[['chrom1', 'pos1', 'chrom2', 'pos2']]
    df.reset_index(drop=True, inplace=True)

    return df

def DBSCAN_wrapper(column):
    try:
        return DBSCAN(eps=500, min_samples=1).fit(column.values.reshape(-1, 1)).labels_
    except:
        pass

def dedup_pairs(rawpairs:pd.DataFrame) -> pd.DataFrame:
    """
    dedup pairs by DBSCAN clustering
    """
    rawpairs = rawpairs.copy()
    raw_pairs_num = rawpairs.shape[0]
    # try chrom1 or chr1 (keeping backward compatibility)
    if "chrom1" in rawpairs.columns:
        rawpairs["cluster1"] = rawpairs.groupby("chrom1")["pos1"].transform(DBSCAN_wrapper)
        rawpairs["cluster2"] = rawpairs.groupby("chrom2")["pos2"].transform(DBSCAN_wrapper)
        rawpairs = rawpairs.groupby(["chrom1","cluster1","chrom2","cluster2"]).head(n=1)
    elif "chr1" in rawpairs.columns:
        rawpairs["cluster1"] = rawpairs.groupby("chr1")["pos1"].transform(DBSCAN_wrapper)
        rawpairs["cluster2"] = rawpairs.groupby("chr2")["pos2"].transform(DBSCAN_wrapper)
        rawpairs = rawpairs.groupby(["chr1","cluster1","chr2","cluster2"]).head(n=1)
    rawpairs = rawpairs.drop(["cluster1","cluster2"], axis=1)
    dedup_pairs_num = rawpairs.shape[0]
    rate = 100*(raw_pairs_num-dedup_pairs_num)/raw_pairs_num

    return rawpairs, rate

def pairs2coverage(pairs, resolution=100000):
    # calculate type of each pair, if the two legs are on different chromosomes or abs(pos2-pos1) > 1000 append both legs,
    # else append the middle point of the two legs
    cross_chrom_or_distant = (pairs['chrom1'] != pairs['chrom2']) | (np.abs(pairs['pos1'] - pairs['pos2']) > 1000)
    
    chrom1_distant = pairs.loc[cross_chrom_or_distant, 'chrom1']
    pos1_distant = pairs.loc[cross_chrom_or_distant, 'pos1']
    chrom2_distant = pairs.loc[cross_chrom_or_distant, 'chrom2']
    pos2_distant = pairs.loc[cross_chrom_or_distant, 'pos2']
    
    chrom_middle = pairs.loc[~cross_chrom_or_distant, 'chrom1']
    pos_middle = ((pairs.loc[~cross_chrom_or_distant, 'pos1'] + pairs.loc[~cross_chrom_or_distant, 'pos2']) // 2)

    # merge craete DataFrame and calc coverage at given resolution
    #chroms = pd.concat([chrom1_distant, chrom2_distant, chrom_middle])
    #positions = pd.concat([pos1_distant, pos2_distant, pos_middle])
    chroms = pd.concat([chrom1_distant, chrom_middle])
    positions = pd.concat([pos1_distant, pos_middle])
    covs = pd.DataFrame({'chrom': chroms, 'pos': positions})
    covs['pos'] = (covs['pos'] // resolution) * resolution
    covs = covs.groupby(['chrom', 'pos']).size().reset_index(name='count')

    return covs
