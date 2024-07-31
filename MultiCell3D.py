import numpy as np
import pandas as pd
from CHARMtools import Cell3D
import tqdm
import concurrent.futures
import sys
from functools import partial
import pybedtools
import warnings
from scipy import stats
from statsmodels.stats.multitest import multipletests
from functools import partial

def dev_only(func):
    import functools
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        if __name__ == "__main__":
            return func(*args, **kwargs)
        else:
            print(f"Skipping {func.__name__} as it's being imported")
    
    return wrapper

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
def _convert_int_resolution_to_string(num:int):
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
    
def _process_cell_CHARM(enrich_cellname,path, resolution, CpG_path=None, peaks_atac=None, peaks_ct=None, flank=200,rep=1):
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

def load_CHARM(enrich_cellnames, path, resolution, CpG_path=None, peaks_atac=None, peaks_ct=None, flank=200,num_cores=30):
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
    #return MultiCell3D(cells)
    # for dev and debugging
    return cells

def _process_cell(cellname,path, resolution):
    cell = Cell3D.Cell3D(cellname = cellname,tdg_path = path,resolution = resolution)
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
    return MultiCell3D(cells)

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

# 这个部分需要添加两个函数，一个是线性回归，看对于ATAC/CT 两个矩阵的情况各自能解释多少结构
# 另一个是画图，比较一下ATAC mat CT mat 和 Hi-C mat 
def chromatic_diff():
    pass

# visualization
def plot_diff(diff,chrom_plot):
    """
    TODO 
    Plot the difference between two groups of cells.
    """
    import matplotlib.pyplot as plt
    df = diff.query('chrom == @chrom_plot')
    plt.figure(figsize=(6, 4))
    plt.plot(df['pos'], df['mean_group1'], label='Group 1', color='blue', alpha=1)
    plt.plot(df['pos'], df['mean_group2'], label='Group 2', color='red', alpha=1)

    significant_points = df[df['p_value_adj'] < 0.05]
    for index, row in significant_points.iterrows():
        plt.axvspan(row['pos'], row['pos'] + 1000000, color='grey', alpha=0.2)

    plt.ylim(0.7,1.3)
    plt.xlabel(chrom_plot)
    plt.ylabel('Radial position')
    plt.legend()

    plt.show()

class MultiCell3D:
    def __init__(self, cells):
        # remove "None" in cells
        cells = [cell for cell in cells if cell is not None]
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
    
    def __getitem__(self, key):
        if isinstance(key, slice):
            selected_cellnames = self.cellnames[key]
            selected_cells = [self.cells_dict[name] for name in selected_cellnames]
            return MultiCell3D(selected_cells)
        elif isinstance(key, int):
            cellname = self.cellnames[key]
            return self.cells_dict[cellname]
        else:
            raise TypeError("Invalid argument type.")

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
        # elif is list or ndarray
        elif isinstance(cellnames, np.ndarray):
            return [self.cells_dict[cellname] for cellname in cellnames]
        elif isinstance(cellnames, list):
            return [self.cells_dict[cellname] for cellname in cellnames]
        else:
            raise ValueError("cellnames should be a list or a string.")
        
    def subset(self,cellnames):
        """
        INPUT: list of cellnames 
        OUTPUT: a new MultiCell3D object 
        """
        cells = self.get_cell(cellnames)
        obj = MultiCell3D(cells)
        for key in self.matrices.keys():
            obj.matrices[key] = self.matrices[key].loc[:,cellnames]
        return obj

    def calc_distance_matrix(self, genome_coord, cellnames=None,allele=True,combine=True):
        """
        Calculate the distance matrix between cells for a given genomic coordinate.
        If cellnames is None, all cells will be used.
        """
        if cellnames is None:
            cellnames = self.cellnames
        temp_cells = self.get_cell(cellnames)
        if allele:
            with concurrent.futures.ProcessPoolExecutor(nproc) as executor:
                results = list(tqdm.tqdm(executor.map(partial(_calc_distance,genome_coord=genome_coord), temp_cells), total=len(temp_cells)))
        else:
            with concurrent.futures.ProcessPoolExecutor(nproc) as executor:
                results1 = list(tqdm.tqdm(executor.map(partial(_calc_distance,genome_coord=genome_coord.replace(":","a:")), temp_cells), total=len(temp_cells)))
                results2 = list(tqdm.tqdm(executor.map(partial(_calc_distance,genome_coord=genome_coord.replace(":","b:")), temp_cells), total=len(temp_cells)))
            results = results1 + results2
        if combine:
            return np.nanmean(results, axis=0)
        else:
            return np.array(results)
            


    def calc_3dproximity_matrix(self, genome_coord, distance_threshold=3, cellnames=None, allele=True, combine=True, nproc=20):
        if cellnames is None:
            cellnames = self.cellnames
        temp_cells = self.get_cell(cellnames)

        if allele:
            with concurrent.futures.ProcessPoolExecutor(nproc) as executor:
                results = list(tqdm.tqdm(executor.map(partial(_calc_proximity,genome_coord=genome_coord,distance_threshold=distance_threshold), temp_cells), total=len(temp_cells)))
        else:
            with concurrent.futures.ProcessPoolExecutor(nproc) as executor:
                results1 = list(tqdm.tqdm(executor.map(partial(_calc_proximity,genome_coord=genome_coord.replace(":","a:"),distance_threshold=distance_threshold), temp_cells), total=len(temp_cells)))
                results2 = list(tqdm.tqdm(executor.map(partial(_calc_proximity,genome_coord=genome_coord.replace(":","b:"),distance_threshold=distance_threshold), temp_cells), total=len(temp_cells)))
            results = results1 + results2
        if combine:
            return np.nanmean(results, axis=0)
        else:
            return np.array(results)


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
                cell.calc_radial_position(key,**kwargs)
            self.features.append(key)
        radial_positions = []
        for cell in tqdm.tqdm(self.get_cell(self.cellnames)):
            radial_positions.append(cell.get_data()[["chrom","pos",key]].set_index(["chrom","pos"]))
        #mat = pd.concat(radial_positions,axis=1)
        mat = _parallel_concat(radial_positions)
        mat.columns = self.cellnames
        self.matrices[key] = mat
        return None

    def get_feature_vec(self,genome_coord,column_name,cellnames=None,combine=True,allele=True):
        """
        Get the feature vector for a given genomic coordinate.
        """
        if cellnames is None:
            cellnames = self.cellnames
        if allele:
            results = [cell.get_feature_vec(genome_coord,column_name) for cell in tqdm.tqdm(self.get_cell(cellnames))]
        else:
            genome_coorda = genome_coord.replace(":","a:")
            genome_coordb = genome_coord.replace(":","b:")
            results = []
            for cell in tqdm.tqdm(self.get_cell(cellnames)):
                results.append(cell.get_feature_vec(genome_coorda,column_name))
                results.append(cell.get_feature_vec(genome_coordb,column_name))
        if combine:
            return np.nanmean(np.array(results), axis=0)
        else:
            return np.array(results)

    def FindMarkers(self, matrix_key,cellnames_group1, cellnames_group2, method = "manwhitneyu"):
        """
        Find markers between two groups of cells.
        """
        df = self.matrices[matrix_key]
        diff = []
        for position in tqdm.tqdm(df.groupby(["chrom","pos"])):
            index, temp_df = position
            data1 = temp_df.loc[:,cellnames_group1].values.reshape(1,-1)[0]
            data2 = temp_df.loc[:,cellnames_group2].values.reshape(1,-1)[0]
            
            if len(data1) == 0 or len(data2) == 0:
                continue
            if method == "manwhitneyu":
                stat, p = stats.mannwhitneyu(data1, data2, alternative="two-sided")
            elif method == "ttest":
                stat, p = stats.ttest_ind(data1, data2)
            else:
                raise ValueError("method should be manwhitneyu or ttest, others not implented yet.")

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                mean_data1 = np.nanmean(data1)
                mean_data2 = np.nanmean(data2)
                mean_diff = mean_data1 - mean_data2
            temp_res = pd.DataFrame({"chrom":index[0], "pos":index[1],"mean_group1":mean_data1, "mean_group2":mean_data2,
                                      "mean_diff":mean_diff, "stat":stat, "p_value":p}, index=[0])
            diff.append(temp_res)

        diff = pd.concat(diff).dropna()
        diff["p_value_adj"] = multipletests(diff["p_value"], method="fdr_bh")[1]

        return diff

    def zoomify_matrix(self,matrix_key,resolution,combine_allele=True,inplace=False,new_key=None):
        """
        Zoomify the matrix to a given resolution.
        """
        df=self.matrices[matrix_key].reset_index()
        df['pos'] = (df['pos'] // resolution) * resolution
        df.set_index(['chrom', 'pos'], inplace=True)
        df = df.groupby(level=['chrom', 'pos']).mean()
        df = df.reset_index()

        if combine_allele:
            df["type"] = df["chrom"].str[-1]
            df["chrom"] = df["chrom"].str[:-1]
            df = df.set_index(["chrom","pos","type"])
        else:
            df = df.set_index(["chrom","pos"])

        if inplace:
            self.matrices[matrix_key] = df
        else:
            if new_key is None:
                new_key = matrix_key + "_zoom" + str(resolution)
            self.matrices[new_key] = df

    def simple_diff(self,cellnames_group1,cellnames_group2,genome_coord = None):
        """
        Calculate 
        """
        pass


# helper functions
def _calc_proximity(cell,genome_coord,distance_threshold):
    mat = cell.calc_distance_matrix(genome_coord) < distance_threshold
    return mat
def _calc_distance(cell,genome_coord):
    mat = cell.calc_distance_matrix(genome_coord)
    return mat