# python version of simple diff
# zliu 2024-03-10
# use schiclsuter for imputation 
import sys

import numpy as np
import pandas as pd
import scipy
import tqdm

import seaborn as sns
import matplotlib.pyplot as plt
from numba import jit

import concurrent.futures
import gc
import math
from . import CHARMio


# for reference: numba version takes 2.5s, numpy version takes 7.9s
@jit(nopython=True)
def calculate_band_z_scores(matrix):
    """
    Band-ly calculate z-scores of a matrix, used as input for SimpleDiff,
    ignoring NaN values in the calculation.
    
    input: n*n numpy 2d array
    output: band-zscore-normalized n*n numpy 2d array
    """
    n = matrix.shape[0]
    z_score_matrix = np.zeros((n, n), dtype=float)
    
    for k in range(-n + 1, n):
        band = []
        valid_count = 0
        sum_of_band = 0.0
        for i in range(max(0, -k), min(n, n - k)):
            j = i + k
            if not math.isnan(matrix[i, j]):
                band.append(matrix[i, j])
                sum_of_band += matrix[i, j]
                valid_count += 1

        band_length = len(band)
        if valid_count > 1:
            mean = sum_of_band / valid_count
            
            sum_of_squares = 0.0
            for x in band:
                sum_of_squares += (x - mean) ** 2
            
            std = math.sqrt(sum_of_squares / valid_count)
            
            if std != 0:
                for i in range(band_length):
                    row_idx = max(0, -k) + i
                    col_idx = row_idx + k
                    z_score_matrix[row_idx, col_idx] = (band[i] - mean) / std
                
    return z_score_matrix

def calculate_band_z_scores_numpy(matrix):
    """
    Band-ly calculate z-scores of a matrix, used as input for SimpleDiff,
    ignoring NaN values in the calculation.
    
    input: n*n numpy 2d array
    output: band-zscore-normalized n*n numpy 2d array, NaN preserved
    """
    n = matrix.shape[0]
    z_score_matrix = np.zeros((n,n), dtype=float)
    nan_mask = np.isnan(matrix)  # Create a mask of where NaN values are located

    for k in range(-n+1, n):
        band = np.diagonal(matrix, offset=k).copy()
        if len(band) > 1:
            mean = np.nanmean(band)  
            std = np.nanstd(band)  
            if std == 0:
                np.fill_diagonal(z_score_matrix[max(-k,0):, max(k,0):], 0)
            else:
                z_scores = (band - mean) / std
                # Ignore NaN in original data by setting z-scores to NaN where original data was NaN
                z_scores[np.isnan(band)] = np.nan
                np.fill_diagonal(z_score_matrix[max(-k,0):, max(k,0):], z_scores)

    # Preserve original NaN positions by setting NaNs in z_score_matrix where they were in the input matrix
    z_score_matrix[nan_mask] = np.nan

    return z_score_matrix

def center_band_numpy(matrix):
    """
    Band-ly center band, used as input for SimpleDiff,
    ignoring NaN values in the calculation.
    
    input: n*n numpy 2d array
    output: centered n*n numpy 2d array, NaN preserved
    """
    n = matrix.shape[0]
    centered_matrix = np.zeros((n,n), dtype=float)
    nan_mask = np.isnan(matrix) 

    for k in range(-n+1, n):
        band = np.diagonal(matrix, offset=k).copy()
        if len(band) > 1:
            mean = np.nanmean(band)  
            centered = band - mean

            centered[np.isnan(band)] = np.nan
            np.fill_diagonal(centered_matrix[max(-k,0):, max(k,0):], centered)

    centered_matrix[nan_mask] = np.nan

    return centered_matrix

def getBandVec(mat:np.array,bins:int):
    return np.concatenate([np.diagonal(mat,offset = i) for i in range(bins)],axis=0)

def multiple_testing_correction(pvalues, correction_type="FDR"):
    """
    Consistent with R - print
    correct_pvalues_for_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05,
                                          0.069, 0.07, 0.071, 0.09, 0.1])
    from https://github.com/CoBiG2/cobig_misc_scripts/blob/master/FDR.py
    """
    from numpy import array, empty
    pvalues = array(pvalues)
    sample_size = pvalues.shape[0]
    qvalues = empty(sample_size)
    if correction_type == "Bonferroni":
        # Bonferroni correction
        qvalues = sample_size * pvalues
    elif correction_type == "Bonferroni-Holm":
        # Bonferroni-Holm correction
        values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
        values.sort()
        for rank, vals in enumerate(values):
            pvalue, i = vals
            qvalues[i] = (sample_size-rank) * pvalue
    elif correction_type == "FDR":
        # Benjamini-Hochberg, AKA - FDR test
        values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
        values.sort()
        values.reverse()
        new_values = []
        for i, vals in enumerate(values):
            rank = sample_size - i
            pvalue, index = vals
            new_values.append((sample_size/rank) * pvalue)
        for i in range(0, int(sample_size)-1):
            if new_values[i] < new_values[i+1]:
                new_values[i+1] = new_values[i]
        for i, vals in enumerate(values):
            pvalue, index = vals
            qvalues[index] = new_values[i]
    return qvalues

def generate_chrom_vec(path,chrom,blacklist_df,resolution=10000,bins=200):
    mat = CHARMio.read_mat_h5(path,genome_coord=chrom)
    mat = np.log2(mat + 1)
    chrom_blacklist_df = blacklist_df.query('chrom == @chrom')
    indices=[]
    for index,row in chrom_blacklist_df.iterrows():
        start_row = math.floor(row['start'] / resolution)
        end_row = math.ceil(row['end'] / resolution)
        indices += [i for i in range(start_row,end_row)]
    indices = list(set(indices))
    
    mat[indices,:]=np.nan
    mat[:,indices]=np.nan
    np.fill_diagonal(mat,np.nan)

    z_score_band_mat = center_band_numpy(mat)
    z_score_band_vec = getBandVec(z_score_band_mat, bins)
    return z_score_band_vec

def generate_chrom_mat(dir,cellnames,chrom,blacklist_df,cores=20,resolution=10000):
    results = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=cores) as executor:
        futures = [executor.submit(generate_chrom_vec,dir+"/"+cellname+".impute."+chrom+".h5",chrom,blacklist_df,resolution) for cellname in cellnames]
        progress_bar = tqdm.tqdm(concurrent.futures.as_completed(futures),total=len(futures),desc="generate mat for compare")

        for future in progress_bar:
            results.append(future.result())
    results = np.vstack(results)
    return results

def test_diff(test_matrix_celltype1,test_matrix_celltype2,method="t",n_cores=10,chunks=50):
    n_features = test_matrix_celltype1.shape[1]
    chunk_size = (n_features + chunks - 1) // chunks  
    chunks = [list(range(i, min(i+chunk_size, n_features))) for i in range(0, n_features, chunk_size)]

    results = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_cores) as executor:
        if method == "t":
            futures = [executor.submit(scipy.stats.ttest_ind,test_matrix_celltype1[:, chunk], test_matrix_celltype2[:, chunk], alternative='two-sided') for chunk in chunks]
        else:
            futures = [executor.submit(scipy.stats.mannwhitneyu,test_matrix_celltype1[:, chunk], test_matrix_celltype2[:, chunk], alternative='two-sided') for chunk in chunks]
        progress_bar = tqdm.tqdm(range(len(futures)), total=len(futures), desc="Processing chunks")
        for i, future in enumerate(futures):
            try:
                results.append(future.result())
                progress_bar.update(1)
            except Exception as e:
                print(f"An error occurred: {e}")
                results.append(None)

    pvs = np.array([i[1] for i in results]).reshape(-1)
    stats = np.array([i[0] for i in results]).reshape(-1)

    return pvs,stats


def post_process(pvs,stats,test_matrix_celltype1,test_matrix_celltype2,reconstructMatrixShapeLength,bins=200):
    mat1_vec = np.nanmean(test_matrix_celltype1,axis=0)
    mat2_vec = np.nanmean(test_matrix_celltype2,axis=0)

    diagonals = []
    diagonals_stats = []
    diagonals_mat1 = []
    diagonals_mat2 = []

    start=0

    for i in range(bins):
        diagonals.append(pvs[start:start+reconstructMatrixShapeLength-i])
        diagonals_stats.append(stats[start:start+reconstructMatrixShapeLength-i])
        diagonals_mat1.append(mat1_vec[start:start+reconstructMatrixShapeLength-i])
        diagonals_mat2.append(mat2_vec[start:start+reconstructMatrixShapeLength-i])

        start += reconstructMatrixShapeLength-i
        
    reMat_pvs = scipy.sparse.diags(np.array(diagonals, dtype=object),[i for i in range(bins)]).toarray()
    reMat_stats = scipy.sparse.diags(np.array(diagonals_stats, dtype=object),[i for i in range(bins)]).toarray()
    reMat_mat1 = scipy.sparse.diags(np.array(diagonals_mat1, dtype=object),[i for i in range(bins)]).toarray()
    reMat_mat2 = scipy.sparse.diags(np.array(diagonals_mat2, dtype=object),[i for i in range(bins)]).toarray()
    reMat_diff = reMat_mat1 - reMat_mat2

    return reMat_pvs,reMat_stats,reMat_mat1,reMat_mat2,reMat_diff

def generate_res_df(reMat_pvs,reMat_stats,reMat_diff,chrom,resolution=10000):

    where = np.where(reMat_pvs > 0)
    where_list = (np.array(where)*resolution).tolist()
    where_list.append(reMat_stats[where].tolist())
    where_list.append(reMat_pvs[where].tolist())
    where_list.append(reMat_diff[where].tolist())

    df = pd.DataFrame(where_list).T
    df.columns = ["start1","start2","stats","pv","diff"]
    df["end1"] = df["start1"] + resolution
    df["end2"] = df["start2"] + resolution

    df["start1"] = df["start1"].astype(int)
    df["end1"] = df["end1"].astype(int)
    df["start2"] = df["start2"].astype(int)
    df["end2"] = df["end2"].astype(int)

    df["chrom1"] = chrom
    df["chrom2"] = chrom
    df = df[["chrom1","start1","end1","chrom2","start2","end2","stats","pv","diff"]]

    return df