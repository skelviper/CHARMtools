import numpy as np
import pandas as pd
import tqdm
import concurrent.futures
from functools import partial
from scipy import stats
import warnings
from .multicell3d_utils import _auto_genome_coord

class MultiCell3DCore:
    """
    Core functionality for MultiCell3D class.
    """
    
    def __init__(self, cells):
        # remove "None" in cells
        cells = [cell for cell in cells if cell is not None]
        self.cells_dict = {cell.cellname: cell for cell in cells}
        self.num_cells = len(cells)
        self.cellnames = list(self.cells_dict.keys())
        self.resolutions = list(set([cell.resolution for cell in cells]))
        self.features = list(set(sum([cell.features for cell in cells], [])))
        self.metadata = None

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
            from .MultiCell3D import MultiCell3D
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
        return ""

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

    def subset(self, cellnames):
        """
        INPUT: list of cellnames 
        OUTPUT: a new MultiCell3D object 
        """
        cells = self.get_cell(cellnames)
        from .MultiCell3D import MultiCell3D
        obj = MultiCell3D(cells)
        for key in self.matrices.keys():
            obj.matrices[key] = self.matrices[key].loc[:,cellnames]
        return obj

    def get_data(self, nproc=20, cellnames=None, **kwargs):
        """
        Get the data for a given cellname or list of cellnames.
        Params:
            cellnames: list of cellnames
            nproc: number of processes to use
            **kwargs: additional arguments for Cell3D.get_data
        Return:
            list of dataframe
        """
        if cellnames is None:
            cellnames = self.cellnames
        temp_cells = self.get_cell(cellnames)

        if nproc == 1:
            results = [_get_data(cell=cell, **kwargs) for cell in temp_cells]
        else:
            with concurrent.futures.ProcessPoolExecutor(nproc) as executor:
                results = list(executor.map(partial(_get_data, **kwargs), [cell for cell in temp_cells]))

        return results

    def get_feature_vec(self, genome_coord, column_name, cellnames=None, combine=True, allele=True, nproc=20):
        """
        Get the feature vector for a given genomic coordinate.
        """
        if cellnames is None:
            cellnames = self.cellnames
        if allele:
            with concurrent.futures.ProcessPoolExecutor(nproc) as executor:
                results = list(executor.map(partial(_get_feature_vec, genome_coord=genome_coord, column_name=column_name), self.get_cell(cellnames)))
        else:
            genome_coorda = genome_coord.replace(":", "a:")
            genome_coordb = genome_coord.replace(":", "b:")
            with concurrent.futures.ProcessPoolExecutor(nproc) as executor:
                results1 = list(executor.map(partial(_get_feature_vec, genome_coord=genome_coorda, column_name=column_name), self.get_cell(cellnames)))
                results2 = list(executor.map(partial(_get_feature_vec, genome_coord=genome_coordb, column_name=column_name), self.get_cell(cellnames)))
            results = results1 + results2
        results = np.array(results)
        results[np.isnan(results)] = 0
        if combine:
            return np.mean(results, axis=0)
        else:
            return results

    def FindMarkers(self, matrix_key, cellnames_group1, cellnames_group2, method="manwhitneyu"):
        """
        Find markers between two groups of cells.
        """
        from statsmodels.stats.multitest import multipletests
        
        df = self.matrices[matrix_key]
        diff = []
        for position in tqdm.tqdm(df.groupby(["chrom", "pos"])):
            index, temp_df = position
            data1 = temp_df.loc[:, cellnames_group1].values.reshape(1, -1)[0]
            data2 = temp_df.loc[:, cellnames_group2].values.reshape(1, -1)[0]
            
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
            temp_res = pd.DataFrame({"chrom": index[0], "pos": index[1], "mean_group1": mean_data1, "mean_group2": mean_data2,
                                      "mean_diff": mean_diff, "stat": stat, "p_value": p}, index=[0])
            diff.append(temp_res)

        diff = pd.concat(diff).dropna()
        diff["p_value_adj"] = multipletests(diff["p_value"], method="fdr_bh")[1]

        return diff

    def zoomify_matrix(self, matrix_key, resolution, combine_allele=True, inplace=False, new_key=None):
        """
        Zoomify the matrix to a given resolution.
        """
        df = self.matrices[matrix_key].reset_index()
        df['pos'] = (df['pos'] // resolution) * resolution
        df.set_index(['chrom', 'pos'], inplace=True)
        df = df.groupby(level=['chrom', 'pos']).mean()
        df = df.reset_index()

        if combine_allele:
            df["type"] = df["chrom"].str[-1]
            df["chrom"] = df["chrom"].str[:-1]
            df = df.set_index(["chrom", "pos", "type"])
        else:
            df = df.set_index(["chrom", "pos"])

        if inplace:
            self.matrices[matrix_key] = df
        else:
            if new_key is None:
                new_key = matrix_key + "_zoom" + str(resolution)
            self.matrices[new_key] = df


# Helper functions
def _get_data(cell, **kwargs):
    df = cell.get_data(**kwargs)
    df["cellname"] = cell.cellname
    return df

def _get_feature_vec(cell, genome_coord, column_name):
    return cell.get_feature_vec(genome_coord, column_name)