# Hi-C matrix imputation relatred code

import torch
import os
import glob
import numpy as np
import pandas as pd
import scipy as sp
import cooler
from scipy import ndimage
from scipy import signal
from scipy import spatial
import tqdm

# compare matrix
# from skimage.metrics import structural_similarity as ssim

# sturctural 
def calc_distance(tdg:pd.DataFrame,bins = 50):
    nparray= np.array(tdg.values.tolist(),dtype=np.str_)
    res = []
    for entry_index in tqdm.tqdm(range(len(nparray))):
        for i in range(bins):
            # only calculate distance intra chromosome
            if entry_index+i >= len(nparray) or nparray[entry_index][0] != nparray[entry_index+i][0]:
                break
            temp = list(nparray[entry_index][:2]) + list(nparray[entry_index + i][:2])
            distance = np.linalg.norm(nparray[entry_index][2:].astype(np.float64)-nparray[entry_index+i][2:].astype(np.float64))
            temp.append(distance)
            res.append(temp)
    return res

def distances_to_dense_mat(distances_vector, resolution = 20000):
    # list to dataframe
    df = pd.DataFrame(distances_vector, columns=['Chromosome1', 'Position1', 'Chromosome2', 'Position2', 'Distance'])
    df['Position1'] = df['Position1'].astype(int)
    df['Position2'] = df['Position2'].astype(int)

    # relocate 
    df['Position1_bin'] = df['Position1'] // resolution - min(df['Position1']) // resolution
    df['Position2_bin'] = df['Position2'] // resolution - min(df['Position1']) // resolution

    # create empty matrix
    size = max(df['Position1_bin'].max(), df['Position2_bin'].max()) + 1
    matrix = pd.DataFrame(np.zeros((size, size)), index=range(size), columns=range(size))

    # fill dense matrix with sparse distances 
    for _, row in df.iterrows():
        matrix.loc[row['Position1_bin'], row['Position2_bin']] = row['Distance']
        matrix.loc[row['Position2_bin'], row['Position1_bin']] = row['Distance'] # symmetric

    return matrix.to_numpy()

# schicluster
def solve_rwr_inverse(stoch_matrix, alpha = 0.05):
    m = stoch_matrix*(1-alpha)
    m = m.transpose()
    #y = sp.sparse.spdiags([1] * m.shape[0], 0, m.shape[0], m.shape[0], format = "csc")
    y = np.eye(m.shape[0])
    A = y - m

    s = None
    #A = A.todense()
    #y = y.todense()
    s = sp.linalg.solve(A, y)

    s *= alpha
    s += s.transpose()
    
    if y is not None:
        del y
    if A is not None:
        del A
    if m is not None:
        del m
    return s

def schicluster_imputation_for_mat(mat,alpha=0.05,kernel_size=3,sigma=2,if_convolve=True,scale_factor=1e5):
    gauss_kernel_1d = signal.gaussian(kernel_size, std=sigma)
    gauss_kernel_2d = np.outer(gauss_kernel_1d, gauss_kernel_1d)

    if if_convolve:
        # add if since snapHi-C did not convolve the matrix
        mat = ndimage.convolve(mat, gauss_kernel_2d, mode='constant', cval=0.0)

    np.fill_diagonal(mat[1:,:-1], mat[1:,:-1] + 1)
    np.fill_diagonal(mat[:-1,1:], mat[:-1,1:] + 1)
    # mat to stochastic matrix
    mat = mat / np.nansum(mat, axis = 0)

    mat = solve_rwr_inverse(mat,alpha)
    mat = mat / np.nansum(mat) * scale_factor
    return mat    

def normalize_matrix(matrix):
    """
    z-score normalization for band
    """
    from scipy.stats import zscore
    normalized_matrix = np.zeros_like(matrix)
    for i in range(-matrix.shape[0] + 1, matrix.shape[1]):
        band = matrix.diagonal(i)
        normalized_band = zscore(band)
        
        if i >= 0:
            np.fill_diagonal(normalized_matrix[i:], normalized_band)
        else:
            np.fill_diagonal(normalized_matrix[:, -i:], normalized_band)
    
    return normalized_matrix

## scHiCluster GPU
def solve_rwr_inverse(stoch_matrix, alpha=0.05, device='cuda'):
    with torch.no_grad():
        m = stoch_matrix * (1 - alpha)
        m = m.T
        y = torch.eye(m.shape[0], device=device, dtype=torch.float32)
        A = y - m

        s = torch.linalg.solve(A, y)

        s *= alpha
        s = s + s.T.clone()  

        return s

def schicluster_imputation_for_mat(mat, alpha=0.05, kernel_size=3, sigma=2, if_convolve=True, scale_factor=1e5, device='cuda'):
    gauss_kernel_1d = signal.windows.gaussian(kernel_size, std=sigma)
    gauss_kernel_2d = np.outer(gauss_kernel_1d, gauss_kernel_1d)

    if if_convolve:
        mat = ndimage.convolve(mat, gauss_kernel_2d, mode='constant', cval=0.0)

    np.fill_diagonal(mat[1:, :-1], mat[1:, :-1] + 1)
    np.fill_diagonal(mat[:-1, 1:], mat[:-1, 1:] + 1)

    mat = mat / np.nansum(mat, axis=0)

    mat = torch.tensor(mat, dtype=torch.float32, device=device)
    mat = solve_rwr_inverse(mat, alpha, device=device)
    mat = mat / torch.nansum(mat) * scale_factor
    return mat.cpu().numpy()



## CtG
# implentation of the CtG imputaion algorithm for Hi-C data and compared with others.
# ref: Computational Enhanced Hi-C data reveals the function of structural geometry in genomic regulation  
#      https://doi.org/10.1101/2022.07.12.499232
# original authos: Yueying He et al.
def ctg_impute(W, lambda_=4):
    # Calculate the degree matrix D
    # since the inverse of D is used, we need to make sure that the diagonal elements of D are not 0
    W = np.nan_to_num(W)
    diag_indices = np.diag_indices_from(W)
    W[diag_indices] = np.where(W[diag_indices] == 0, 1, W[diag_indices])

    D = np.diag(np.sum(W, axis=1))
    # Calculate the 1-step transition probability matrix P
    P = np.linalg.inv(D) @ W 

    # Diagonalize P
    eigenvalues, eigenvectors = np.linalg.eig(P)
    Lambda = np.diag(eigenvalues)
    U = eigenvectors

    # Calculate the k-step transition probability matrix P_k and the transition propensity matrix S
    #P_k = np.linalg.matrix_power(P, k)
    #S_k = np.sum([np.exp(-lambda_ * t) * np.linalg.matrix_power(P, t) for t in range(1, k+1)], axis=0)

    # Calculate the limit of S as k approaches infinity
    S = U @ (Lambda / (np.exp(lambda_) - Lambda)) @ np.linalg.inv(U)
    
    CTG_distance_matrix = spatial.distance.squareform(spatial.distance.pdist(S, metric='cityblock'))

    return CTG_distance_matrix

def ctg_impute_v2(W, iterations=20, alpha=0.3):
    W = W.astype(np.float64)
    def normalize_matrix(W):
        """ Normalize the matrix W by column. """
        n = W.shape[0]
        for col in range(n):
            col_sum = np.sum(W[:, col])
            if np.abs(col_sum) > 1e-6:
                W[:, col] /= col_sum
        return W

    def matrix_iteration(W, iterations, alpha):
        """ Perform matrix iterations. """
        n = W.shape[0]
        sum_matrix = np.zeros_like(W)
        k_matrix = np.copy(W)

        for i in range(iterations):
            factor = np.exp(-alpha * (i + 1))
            sum_matrix += factor * k_matrix
            k_matrix = np.dot(k_matrix, W)
        return sum_matrix

    def compute_L1_distance(matrix):
        """ Compute L1 distance matrix. """
        n = matrix.shape[0]
        L1_matrix = np.zeros((n, n))

        for i in range(n):
            for j in range(n):
                L1_matrix[i, j] = np.sum(np.abs(matrix[:, i] - matrix[:, j]))
        return L1_matrix
    
        #L1_matrix = np.sum(np.abs(matrix[:, np.newaxis] - matrix[np.newaxis, :]), axis=2)
        return L1_matrix

    normalized_W = normalize_matrix(W)
    iterated_matrix = matrix_iteration(normalized_W, iterations, alpha)
    L1_distance_matrix = compute_L1_distance(iterated_matrix)

    return L1_distance_matrix

## BandNorm
def readSingleCellCount(file):
    contactMatrix = pd.read_csv(file,sep="\t",names=["chr1","bin1","chr2","bin2","diag","count"])
    # treat all inter contacts as one band
    thisCellInterCount = contactMatrix.query('diag == "-1"')['count'].sum(axis=0)

    # calc intra band per chromosome
    thisCellIntraCount = contactMatrix.query('chr1 == chr2')[['chr1','diag','count']].groupby(by = ['chr1','diag']).sum()
    return [thisCellIntraCount,thisCellInterCount]

def generateBandnormFactor(filepath:list,ncpus=40):
    cellNum = len(filepath)
    totalCellInterCount = 0
    totalCellIntraCount = []

    result = []

    with Pool(ncpus) as p:
        result += p.map(readSingleCellCount,filepath)

    for i in result:
        totalCellIntraCount.append(i[0])
        totalCellInterCount+=i[1]

    totalCellIntraCount = pd.concat(totalCellIntraCount).groupby(by = ['chr1','diag']).sum()
    totalCellIntraCount['count'] = totalCellIntraCount['count'] / cellNum 
    totalCellInterCount = totalCellInterCount / cellNum
    
    return totalCellIntraCount,totalCellInterCount

def normCell(intraCount:pd.DataFrame, interCount:int, contactMatrix:pd.DataFrame) -> pd.DataFrame:
    """
    input: inter/intraCount represent average count of the experiment,
        contactMatrix is contactMatrix of the current cell.
    """
    # calc intra part
    thisCellIntra = contactMatrix.query('chr1 == chr2')[['chr1','diag','count']].groupby(by = ['chr1','diag']).sum() 
    normFactorIntra = pd.merge(thisCellIntra,intraCount,how = 'left',on = ["chr1","diag"])
    normFactorIntra = normFactorIntra.assign(normFactor = lambda line: line.count_y / line.count_x)[["normFactor"]]

    intraDF = pd.merge(contactMatrix.query('chr1 == chr2'),normFactorIntra,on=['chr1','diag'])
    intraDF['normCount'] = intraDF['count'] * intraDF['normFactor']
    intraDF = intraDF[["chr1","bin1","chr2","bin2","normCount"]]

    # calc inter part
    interDF = contactMatrix.query('chr1 != chr2')
    normFactorInter = interCount / interDF['count'].sum(axis=0)
    interDF['normCount'] = interDF['count'] * normFactorInter
    interDF = interDF[["chr1","bin1","chr2","bin2","normCount"]]
    # combine result and return
    return pd.concat([interDF,intraDF])
    
def cell2csr_mat(cellMat:pd.DataFrame,chrAdder:pd.DataFrame,resolution:int) -> pd.DataFrame:
    """
    chrAdder is a dataframe generate by :
        chrAdder = pd.read_csv("/share/home/zliu/project/gradProject/BubbleCluster/otherFiles/chr.len.hg19.tsv",
                sep="\t",names=["chr","len","adding"])[["chr","adding"]]

    output: 	bin1abs	bin2abs	normCount
                0	0	282	1.909429
                1	0	286	0.954715
    """
    chrAdderDict = chrAdder.set_index('chr').to_dict()['adding']
 
    cellMat['bin1abs'] = cellMat.apply(lambda row: math.floor((row.bin1 + chrAdderDict[row['chrom1']])/resolution), axis=1)
    cellMat['bin2abs'] = cellMat.apply(lambda row: math.floor((row.bin2 + chrAdderDict[row['chrom2']])/resolution), axis=1)

    return cellMat[["bin1abs","bin2abs","normCount"]]
