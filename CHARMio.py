#global dependence
import numpy as np
import pandas as pd

#Hi-C
def getMatrixFromMCOOLs(filepath: "str", genome_coord, resolution=40000, balance=False)->np.ndarray:
    """
    intput: mcool filepath ,
            genome_coord(e.g. 'chr1:35,900,000-40,900,000'), 
            resolution(should be included in mcoolfile)
    output: numpy 2d array
    """
    import cooler
    cool = filepath+"::/resolutions/"+str(resolution)

    c = cooler.Cooler(cool)
    matrix = c.matrix(balance=balance).fetch(genome_coord).astype("double")

    return matrix

