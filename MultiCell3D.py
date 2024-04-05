import numpy as np
import pandas as pd
from CHARMtools import Cell3D
import tqdm

class MultiCell3D:
    def __init__(self, cells):
        self.cells_dict = {cell.cellname: cell for cell in cells}
        self.num_cells = len(cells)
        self.cellnames = list(self.cells_dict.keys())
        self.resolutions = list(set([cell.resolution for cell in cells]))
        self.features = list(set(sum([cell.features for cell in cells], [])))

        metadata = []
        for cell in cells:
            metadata_dict = {"cellname": cell.cellname, **cell.metadata}
            metadata.append(metadata_dict)

        self.metadata_df = pd.DataFrame(metadata)

    def __repr__(self):
        self.get_info()
        return ""

    def get_info(self):
        print("CHARMtools MultiCell3D object v0.1")
        print(f"Object contains {self.num_cells} cells starting with: {self.cellnames[:3]}")
        print(f"Resolutions: {self.resolutions}")
        print(f"Features: {self.features}")

    def get_cell(self, cellnames):
        return [self.cells_dict[cellname] for cellname in cellnames]
    
    def get_distance_matrix(self,genome_coord,cells=None):
        """
        Calculate the distance matrix between cells for a given genomic coordinate.
        if cells is None, all cells will be used.
        """
        if cells is None:
            cells = self.get_cell(self.cellnames)
        mats = []
        for cell in tqdm.tqdm(cells):
            mats.append(cell.calc_distance_matrix(genome_coord))
        return np.nanmean(mats,axis=0)

    def get_3dproximity_matrix(self,genome_coord,distance_threshold=3,cells=None):
        """
        Calculate the 3D proximity matrix between cells.
        if cells is None, all cells will be used.
        """
        pass