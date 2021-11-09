from . import CHARMio
import numpy as np
from os import path


def get3dProximityStackMatrix(listCoolPath:list,genome_coord1:str,genome_coord2=None,resolution=20000):
    """
    input: 
    output: percentage matrix
    """
    matList = []
    for cool in listCoolPath:
        if (path.exists(cool) ):
            matList.append(CHARMio.getMatrixFromCooler(cool,genome_coord1,genome_coord2,resolution))
    return np.sum(np.array(matList),0)/len(matList)/2