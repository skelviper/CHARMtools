from . import CHARMio
import numpy as np
from os import path


def get3dProximityStackMatrix(listCoolPath:list,genome_coord1:str,genome_coord2=None,resolution=20000):
    """
    input: 
    output: percentage matrix
    """

    count=0
    mat = 0
    for cool in listCoolPath:
        if (path.exists(cool) ):
            count +=1
            mat = np.add(mat,CHARMio.getMatrixFromCooler(cool,genome_coord1,genome_coord2,resolution))
    return mat/count/2

#split numpy array into different length numpy arrays
def split_array(array,length):
    return [array[i:i+length] for i in range(0, len(array), length)]