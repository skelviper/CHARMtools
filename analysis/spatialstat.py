# reimplementation of spatial autocorrelation in https://github.com/ammondongp/3D_ATAC_PALM/blob/master/SpatialAutocorrelation3D.m

import numpy as np
from scipy.spatial import distance_matrix
from scipy.spatial import ConvexHull, cKDTree

def Ripleys_K_3D(points,r):
    v = ConvexHull(points).volume
    n = len(points)
    tree = cKDTree(points)
    count = tree.query_ball_tree(tree, r)
    count = [len(i) for i in count]
    return (np.sum(count) -n) / n / n * v

def Ripleys_H_3D(points,r):
    k = Ripleys_K_3D(points,r)
    return np.cbrt(k*3/4/np.pi) - r

def spatialxcorr_3D_without_edge(x, y, dr, vol,max_dist=10):
    # pair-wise distance
    c = distance_matrix(x, y)
    
    # max width
    width = max(np.max(y, axis=0) - np.min(y, axis=0))

    Density_y = y.shape[0] / vol
    
    m, n = x.shape
    
    steps = int(np.floor(max_dist / dr))
    
    r = np.zeros(steps)
    Corr = np.zeros(steps)
    
    #for i in tqdm.tqdm(range(steps)):
    for i in range(steps):
        r[i] = dr * (i + 1) - dr / 2
        logic = (c <= r[i] + dr / 2) & (c > r[i] - dr / 2)
        Corr[i] = np.count_nonzero(logic) / (y.shape[0] - 1) / (Density_y * np.pi * (4 * (r[i]**2) * dr + 1/3 * dr**3))
    
    return Corr, r, Density_y, width


def spatial_autocorr_3d(X, Nsim=4, dr=0.5,max_dist=10):
    """
    Calculate the spatial autocorrelation curve for 3D data
    Parameters
    ----------
    X : array
        3D data
    Nsim : int
        Number of simulations for calculating 3D shape autocorrelation
    dr : float
        Bin size
    """
    # Normalize X
    X[:, 0] -= np.min(X[:, 0])
    X[:, 1] -= np.min(X[:, 1])
    X[:, 2] -= np.min(X[:, 2])
    
    len_X = len(X)
    Xs = X.copy()

    width = np.max(np.max(X, axis=0))
    hull = ConvexHull(Xs)
    vol = hull.volume
    PointsNum = len_X
    number = round(PointsNum * width**3 / vol)

    RandomAll = []
    CorrR = []
    rR = []
    lengths = np.zeros(Nsim)
    
    #for i in tqdm.tqdm(range(Nsim)):
    for i in range(Nsim):
        random = np.random.rand(number, 3) * width
        in_hull = np.all(np.dot(hull.equations[:, :3], random.T) + hull.equations[:, 3].reshape(-1, 1) <= 0, axis=0)
        random = random[in_hull, :]
        corr, r, _, _ = spatialxcorr_3D_without_edge(random, random, dr, vol,max_dist)
        CorrR.append(corr)
        rR.append(r)
        lengths[i] = len(corr)
        RandomAll.append(random)

    Corr_R = np.mean(CorrR, axis=0)
    r_R = np.mean(rR, axis=0)
    
    # Generate shape autocorrelation curve for the real data
    Corr, r, _, _ = spatialxcorr_3D_without_edge(X, X, dr, vol,max_dist)
    
    # normalize shape autocorrelation curve
    FCorr = Corr / Corr_R
    Fr = r_R

    return FCorr, Fr

def spatial_crosscorr_3d(X,Y,Nsim=5,dr=0.5,max_dist=10):
    # normalize X, Y
    xx1 = min(np.min(X[:, 0]), np.min(Y[:, 0]))
    xx2 = min(np.min(X[:, 1]), np.min(Y[:, 1]))
    xx3 = min(np.min(X[:, 2]), np.min(Y[:, 2]))

    X[:, 0] -= xx1
    Y[:, 0] -= xx1
    X[:, 1] -= xx2
    Y[:, 1] -= xx2
    X[:, 2] -= xx3
    Y[:, 2] -= xx3

    # generate shape autocorrelation curve for XY
    width = max(np.max(np.max(X, axis=0)), np.max(np.max(Y, axis=0)))
    
    hull = ConvexHull(X)
    vol = hull.volume
    PointsNum = len(X)
    number = round(PointsNum * width**3 / vol)

    RandomAll = []
    CorrR = []
    rR = []
    lengths = np.zeros(Nsim)

    for i in range(Nsim):
        random1 = np.random.rand(number, 3) * width
        random2 = np.random.rand(number, 3) * width

        in_hull1 = np.all(np.dot(hull.equations[:, :3], random1.T) + hull.equations[:, 3].reshape(-1, 1) <= 0, axis=0)
        in_hull2 = np.all(np.dot(hull.equations[:, :3], random2.T) + hull.equations[:, 3].reshape(-1, 1) <= 0, axis=0)

        random1 = random1[in_hull1, :]
        random2 = random2[in_hull2, :]

        corr, r, _, _ = spatialxcorr_3D_without_edge(random1, random2, dr, vol, max_dist)
        CorrR.append(corr)
        rR.append(r)
        lengths[i] = len(corr)

    Corr_R = np.mean(CorrR, axis=0)
    r_R = np.mean(rR, axis=0)

    # Generate shape autocorrelation curve for the real data
    Corr, r, _, _ = spatialxcorr_3D_without_edge(X, Y, dr, vol, max_dist)
    lenC = len(Corr)

    # normalize shape autocorrelation curve
    FCorr = Corr / Corr_R
    Fr = r_R

    return FCorr, Fr