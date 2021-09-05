"""
Compartment level analysis and plot funcitons
@author zliu
@data 20210902
"""
#global dependence
import numpy as np
import pandas as pd
from . import CHARMio

# Decay profile
# for p(s) curve use log_bins=True , otherwise(e.g. normalize distance for Hi-C matrix ) use log_bins=False
def psDataFromMat(matrix, indices=None, log_bins=True, base=1.1):
    """
    ***FUNCTION COPY FROM HICSTAFF***
    Compute distance law as a function of the genomic coordinate aka P(s).
    Bin length increases exponentially with distance if log_bins is True. Works
    on dense and sparse matrices. Less precise than the one from the pairs.
    Parameters
    ----------
    matrix : numpy.array or scipy.sparse.coo_matrix
        Hi-C contact map of the chromosome on which the distance law is
        calculated.
    indices : None or numpy array
        List of indices on which to compute the distance law. For example
        compartments or expressed genes.
    log_bins : bool
        Whether the distance law should be computed on exponentially larger
        bins.
    Returns
    -------
    numpy array of floats :
        The start index of each bin.
    numpy array of floats :
        The distance law computed per bin on the diagonal
    """

    n = min(matrix.shape)
    included_bins = np.zeros(n, dtype=bool)
    if indices is None:
        included_bins[:] = True
    else:
        included_bins[indices] = True
    D = np.array(
        [
            np.average(matrix.diagonal(j)[included_bins[: n - j]])
            for j in range(n)
        ]
    )
    if not log_bins:
        return np.array(range(len(D))), D
    else:
        n_bins = int(np.log(n) / np.log(base) + 1)
        logbin = np.unique(
            np.logspace(0, n_bins - 1, num=n_bins, base=base, dtype=np.int)
        )
        logbin = np.insert(logbin, 0, 0)
        logbin[-1] = min(n, logbin[-1])
        if n < logbin.shape[0]:
            print("Not enough bins. Increase logarithm base.")
            return np.array(range(len(D))), D
        logD = np.array(
            [
                np.average(D[logbin[i - 1] : logbin[i]])
                for i in range(1, len(logbin))
            ]
        )
        return logbin[:-1], logD

def getPearsonCorrMatrix(matrix:np.ndarray)->np.ndarray:
    """
    get decay profile normalized pearson correlation matrix
    """
    n=matrix.shape[0]
    dist_matrix = np.zeros((n, n))
    _, dist_vals = psDataFromMat(matrix, log_bins=False)
    for i in range(n):
        for j in range(n):
            dist_matrix[i, j] = dist_vals[abs(j - i)]
    matrix /= dist_matrix

    matrix = np.corrcoef(matrix)
    return matrix

def addVec(a,b):
    if len(a) < len(b):
        c = b.copy()
        c[:len(a)] += a
    else:
        c = a.copy()
        c[:len(b)] += b
    return c

def getPsData(mcoolPath,chromlist,resolution=10000,celltype="unknown")->pd.DataFrame:
    matlist = [CHARMio.getMatrixFromMCOOLs(mcoolPath,genome_coord=chrom,resolution=resolution) for chrom in chromlist]
    bin = psDataFromMat(matlist[0])[0]
    value = np.array([])
    for mat in matlist:
        value = addVec(value,psDataFromMat(mat)[1])
    return pd.DataFrame({"bin":bin,"aveCount":value,"celltype":celltype})

# quick visualize functions
def plotPsCurve(mcoolsPath:list,celltypeNames:list,chroms:list,resolution=100000,title="P(s) curve",plotType="interaction"):
    """
    plotPsCurve function take bin and value
    """
    import plotly.express as px
    from IPython.display import Image

    #Calculate P(s) data, get a 3 column pd.DataFrame with (bin,resolution,celltype)
    psDataAll = []
    for i in range(len(mcoolsPath)):
        psDataAll.append(getPsData(mcoolsPath[i],["chr"+str(i+1) for i in range(len(chroms))],resolution=resolution,celltype=celltypeNames[i])) 
    merged = pd.concat(psDataAll)

    data =  pd.merge(merged,merged.groupby("celltype").sum(),how="left",on="celltype").assign(prob= lambda df: df.aveCount_x/df.aveCount_y)

    fig = px.line(x=data["bin_x"]*resolution,y=data["prob"],color=data["celltype"],title=title,log_x=True,log_y=True).update_layout(template='simple_white')
    fig.update_layout(width=800,height=600)
    fig.update_layout(xaxis_title="Genomic Distance(bp)",
                    yaxis_title="Contact Probability")
    if(plotType == "interaction"):
        return fig
    else : return Image(fig.to_image(format="png", engine="kaleido"))

def plotMatrix(matrix:np.ndarray,if_log=False,title="Matrix"):
    """
    plotMatrix function for plot hic contact matrix
    """
    import plotly.express as px 
    fig = px.imshow(matrix,color_continuous_scale=px.colors.sequential.Bluered)
    fig = fig.update_layout(template='simple_white').update_layout(width=800,height=600)

    return fig