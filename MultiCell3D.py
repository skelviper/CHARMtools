import numpy as np
import pandas as pd
from CHARMtools import Cell3D
import tqdm
import concurrent.futures
import sys
from functools import partial
import pybedtools
#from concurrent.futures import ThreadPoolExecutor, as_completed

# functions used 
def _concat_in_chunks(data_chunk):
    return pd.concat(data_chunk, axis=1)

def _parallel_concat(data, nproc=10):
    chunk_size = len(data) // nproc
    data_chunks = [data[i:i + chunk_size] for i in range(0, len(data), chunk_size)]

    with concurrent.futures.ProcessPoolExecutor(nproc) as executor:
        concatenated_chunks = list(tqdm.tqdm(executor.map(_concat_in_chunks, data_chunks), total=len(data_chunks)))
    
    final_result = pd.concat(concatenated_chunks, axis=1)
    
    return final_result


# handle bedtools warnings 
class FilteredStderr(object):
    def __init__(self, target):
        self.target = target
        self.ignore_strings = ["has inconsistent naming convention for record:"]

    def write(self, s):
        # Check if any of the ignore strings are in the message
        if not any(ignore_string in s for ignore_string in self.ignore_strings):
            self.target.write(s)
    
    def flush(self):
        self.target.flush()

# io
def _process_cell(enrich_cellname,path, resolution, CpG_path=None, peaks_atac=None, peaks_ct=None, flank=200):
    cellname = enrich_cellname.replace("EN", "")
    cell = Cell3D.Cell3D(cellname=cellname,
                            tdg_path=path + 'hic/processed/{i}/3d_info/clean.20k.0.3dg'.format(i=cellname),
                            resolution=resolution)
    if CpG_path is not None:
        cell.add_bedGraph_data(CpG_path, column_name="CpG", resolution=20000, type="all")
    if peaks_atac is not None:
        cell.add_bed_data(path=path + "enrich/processed/atac_all/{i}.atac.frag.bed.gz".format(i=enrich_cellname),
                            column_name="ATAC", type="all", peaks=peaks_atac, flank=200)
    else:
        cell.add_bed_data(path=path + "enrich/processed/atac_all/{i}.atac.frag.bed.gz".format(i=enrich_cellname),
                            column_name="ATAC", type="all")
    if peaks_ct is not None:
        cell.add_bed_data(path=path + "enrich/processed/ct_all/{i}.ct.frag.bed.gz".format(i=enrich_cellname),
                            column_name="CT", type="all", peaks=peaks_ct, flank=200)
    else:
        cell.add_bed_data(path=path + "enrich/processed/ct_all/{i}.ct.frag.bed.gz".format(i=enrich_cellname),
                            column_name="CT", type="all")

    return cell

def load_CHARM(enrich_cellnames, path, resolution, CpG_path=None, peaks_atac=None, peaks_ct=None, flank=200,num_cores=30):
    """
    Construct a MultiCell3D object from the CHARM dataset.
    """
    original_stderr = sys.stderr
    sys.stderr = FilteredStderr(sys.stderr)

    try:
        with concurrent.futures.ProcessPoolExecutor(num_cores) as executor:
            process_cell_partial = partial(_process_cell, path=path, resolution=resolution, CpG_path=CpG_path,
                                           peaks_atac=peaks_atac, peaks_ct=peaks_ct, flank=flank)
            cells = list(tqdm.tqdm(executor.map(process_cell_partial, enrich_cellnames), total=len(enrich_cellnames)))
    finally:
        sys.stderr = original_stderr
    pybedtools.helpers.cleanup(remove_all=True)
    #return MultiCell3D(cells)
    # for dev and debugging
    return cells

# analysis
def chromatic(mat,vec):
    # resij = matij *veci * vecj
    res = np.zeros(mat.shape)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if i != j :
                res[i,j] = mat[i,j] * vec[i] * vec[j]
            else:
                res[i,j] = mat[i,j] * vec[i]
    return res

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
        self.matrices = {}

    def __repr__(self):
        self.get_info()
        return ""

    def get_info(self):
        print("CHARMtools MultiCell3D object v0.1")
        print(f"Object contains {self.num_cells} cells starting with: {self.cellnames[:3]}")
        print(f"Resolutions: {self.resolutions}")
        print(f"Features: {self.features}")
        print(f"Object contains matrices: {list(self.matrices.keys())}")

    def get_cell(self, cellnames):
        """
        cellnames can be a list of cellnames or a single cellname.
        """
        if isinstance(cellnames, str):
            return self.cells_dict[cellnames]
        elif isinstance(cellnames, list):
            return [self.cells_dict[cellname] for cellname in cellnames]
        else:
            raise ValueError("cellnames should be a list or a string.")
    
    def calc_distance_matrix(self,genome_coord,cellnames=None):
        """
        Calculate the distance matrix between cells for a given genomic coordinate.
        if cells is None, all cells will be used.
        """
        if cellnames is None:
            cellnames = self.cellnames
        mats = []
        for cell in tqdm.tqdm(self.get_cell(cellnames)):
            mats.append(cell.calc_distance_matrix(genome_coord))
        return np.nanmean(mats,axis=0)

    def calc_3dproximity_matrix(self, genome_coord, distance_threshold=3, cellnames=None):
        """
        Calculate the 3D proximity matrix between cells.
        if cells is None, all cells will be used.
        """
        if cellnames is None:
            cellnames = self.cellnames
        mats = []
        for cell in tqdm.tqdm(self.get_cell(cellnames)):
            distance_matrix = cell.calc_distance_matrix(genome_coord)
            mats.append(np.where(np.isnan(distance_matrix), np.nan, distance_matrix <= distance_threshold))
        return np.nanmean(mats, axis=0)
    

    def calc_feature_matrix(self, genome_coord,feature, cells=None):
        """
        Get the feature matrix for a given feature.
        if cells is None, all cells will be used.
        feature should be present in the cell object feature list.
        """
        if cells is None:
            cells = self.get_cell(self.cellnames)
        feature_mats = []
        feature_vecs = []
        for cell in tqdm.tqdm(cells):
            feature_mat, feature_vec = cell.calc_feature_matrix(genome_coord, feature)
            feature_mats.append(feature_mat)
            feature_vecs.append(feature_vec)

        return np.nanmean(feature_mats,axis=0), np.array(feature_vecs)
    
    def calc_feature_proximity_matrix(self, genome_coord, feature, distance_threshold=3, cells=None):
        """
        Calculate the feature proximity matrix between cells.
        if cells is None, all cells will be used.
        feature should be present in the cell object feature list.
        """
        if cells is None:
            cells = self.get_cell(self.cellnames)
        feature_mats = []
        feature_vecs = []
        for cell in tqdm.tqdm(cells):
            feature_mat, feature_vec = cell.calc_feature_proximity_matrix(genome_coord, feature,distance_threshold)
            feature_mats.append(feature_mat)
            feature_vecs.append(feature_vec)
        return np.nanmean(feature_mats, axis=0), np.array(feature_vecs)
    
    def calc_radial_position_matrix(self,key="radial_position",**kwargs):
        """
        Calculate radial position in single cell if radial_position is not present in the features for each cell.
        Combine the radial position for all cells, add to the matrix dictionary

        Params: 
        key: str, default "radial_position"
        **kwargs: additional arguments for calc_radial_position
        """
        if key not in self.features:
            for cell in self.get_cell(self.cellnames):
                cell.calc_radial_position(**kwargs)
            self.features.append(key)
        radial_positions = []
        for cell in tqdm.tqdm(self.get_cell(self.cellnames)):
            radial_positions.append(cell.get_data()[["chrom","pos",key]].set_index(["chrom","pos"]))
        #mat = pd.concat(radial_positions,axis=1)
        mat = _parallel_concat(radial_positions)
        mat.columns = self.cellnames
        self.matrices[key] = mat
        return None