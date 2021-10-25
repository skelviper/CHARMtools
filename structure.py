from . import CHARMio
import numpy as np


def get3dProximityStackMatrix(listCoolPath:list,genome_coord1:str,genome_coord2=None,resolution=40000):
    """
    input: 
    output: percentage matrix
    """
    matList = [CHARMio.getMatrixFromCooler(cool,genome_coord1,genome_coord2,resolution) for cool in listCoolPath]
    return np.sum(np.array(matList),0)/len(matList)/2