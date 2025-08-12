import numpy as np
import pandas as pd
import tqdm
import concurrent.futures
import sys
from functools import partial
import pybedtools
from .Cell3D import Cell3D
from .multicell3d_utils import _parallel_concat

# Handle bedtools warnings 
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

# IO functions
def _convert_int_resolution_to_string(num: int):
    """
    Simplify a number to a shorter format using 'k' for thousands, 'm' for millions, etc.
    
    Parameters:
    num (int): The number to be simplified.

    Returns:
    str: A simplified version of the number with 'k', 'm', etc.
    """
    if num < 1000:
        return str(num)
    elif num < 1000000:
        # Thousands - round to nearest thousand
        return f"{num // 1000}k"
    elif num < 1000000000:
        # Millions - round to nearest million
        return f"{num // 1000000}m"
    else:
        # Billions - round to nearest billion
        return f"{num // 1000000000}b"
    
def _process_cell_CHARM(enrich_cellname, path, resolution, CpG_path=None, peaks_atac=None, peaks_ct=None, flank=200, rep=1):
    try:
        cellname = enrich_cellname.replace("EN", "")
        cell = Cell3D.Cell3D(cellname=cellname,
                                tdg_path=path + f'hic/processed/{cellname}/3d_info/clean.{_convert_int_resolution_to_string(resolution)}.{rep}.3dg',
                                resolution=resolution)
        if CpG_path is not None:
            cell.add_bedGraph_data(CpG_path, column_name="CpG", resolution=resolution, type="all")
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
    except Exception as e:
        print(f"Error processing cell {enrich_cellname}: {e}")
        return None

def load_CHARM(enrich_cellnames, path, resolution, CpG_path=None, peaks_atac=None, peaks_ct=None, flank=200, num_cores=30):
    """
    Construct a MultiCell3D object from the CHARM dataset.
    """
    original_stderr = sys.stderr
    sys.stderr = FilteredStderr(sys.stderr)

    try:
        with concurrent.futures.ProcessPoolExecutor(num_cores) as executor:
            process_cell_partial = partial(_process_cell_CHARM, path=path, resolution=resolution, CpG_path=CpG_path,
                                           peaks_atac=peaks_atac, peaks_ct=peaks_ct, flank=flank)
            cells = list(tqdm.tqdm(executor.map(process_cell_partial, enrich_cellnames), total=len(enrich_cellnames)))
    finally:
        sys.stderr = original_stderr
    pybedtools.helpers.cleanup(remove_all=True)
    # return MultiCell3D(cells)
    # for dev and debugging
    return cells

def _process_cell(cellname, path, resolution):
    cell = Cell3D.Cell3D(cellname=cellname, tdg_path=path, resolution=resolution)
    return cell

def load_cells(cellnames, path, resolution, num_cores=20):
    """
    cellnames: list of cell names
    path: list of path of the tdg files
    resolution: resolution of the tdg files, e.g. 20000
    num_cores: number of cores to use
    """
    with concurrent.futures.ProcessPoolExecutor(num_cores) as executor:
        cells = list(tqdm.tqdm(executor.map(_process_cell, cellnames, path, [resolution]*len(cellnames)), total=len(cellnames)))
    from .MultiCell3D import MultiCell3D
    return MultiCell3D(cells)

class MultiCell3DIO:
    """
    Input/Output functionality for MultiCell3D class.
    """
    
    def combine_feature(self, key=None, **kwargs):
        """
        Combine feature column in single-cells into matrix store in the multicell obj.
        """
        # 检查key是否在features中
        if key not in self.features:
            raise ValueError(f"Feature {key} not found in cell features")
            
        # 收集所有细胞的特征数据
        feature_data = []
        for cell in tqdm.tqdm(self.get_cell(self.cellnames)):
            feature_data.append(cell.get_data()[["chrom", "pos", key]].set_index(["chrom", "pos"]))
            
        # 使用并行方式合并数据
        mat = _parallel_concat(feature_data)
        mat.columns = self.cellnames
        self.matrices[key] = mat
        return None

    def calc_radial_position_matrix(self, key="radial_position", **kwargs):
        """
        Calculate radial position in single cell if radial_position is not present in the features for each cell.
        Combine the radial position for all cells, add to the matrix dictionary

        Params: 
        key: str, default "radial_position"
        **kwargs: additional arguments for calc_radial_position
        """

        if key not in self.features:
            for cell in self.get_cell(self.cellnames):
                cell.calc_radial_position(key, **kwargs)
            self.features.append(key)
            
        # 调用combine_feature来合并数据
        self.combine_feature(key=key)
        return None