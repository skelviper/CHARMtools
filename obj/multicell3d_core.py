import numpy as np
import pandas as pd
import tqdm
import concurrent.futures
from functools import partial
from scipy import stats
import warnings
import anndata as ad
from typing import List, Dict, Optional, Union, Tuple
from ..utils.helper import auto_genome_coord

class MultiCell3DCore:
    """
    Core functionality for MultiCell3D class.
    """
    
    def __init__(self, cells):
        # remove "None" in cells
        cells = [cell for cell in cells if cell is not None]
        self.cells = cells
        self.cells_dict = {cell.cellname: cell for cell in cells}
        self.num_cells = len(cells)
        self.cellnames = list(self.cells_dict.keys())
        self.resolutions = list(set([cell.resolution for cell in cells]))
        self.features = list(set(sum([cell.features for cell in cells], [])))
        self.metadata = None

        # Build metadata DataFrame from cells
        self._build_metadata()
        
        # Initialize storage
        self.matrices = {}  # Legacy matrix storage for backward compatibility
        self.adata_dict = {}  # Dictionary to store different AnnData objects
        
    def _build_metadata(self):
        """
        Build metadata DataFrame from Cell3D objects.
        """
        metadata_list = []
        for cell in self.cells:
            cell_metadata = {'cellname': cell.cellname}
            
            if hasattr(cell, 'metadata') and cell.metadata is not None:
                if isinstance(cell.metadata, dict):
                    cell_metadata.update(cell.metadata)
                elif isinstance(cell.metadata, pd.Series):
                    cell_metadata.update(cell.metadata.to_dict())
                else:
                    # Try to convert to dict
                    try:
                        cell_metadata.update(dict(cell.metadata))
                    except:
                        pass
            
            metadata_list.append(cell_metadata)
        
        if metadata_list:
            self.metadata_df = pd.DataFrame(metadata_list)
            self.metadata_df.set_index('cellname', inplace=True)
        else:
            # Create empty metadata DataFrame
            self.metadata_df = pd.DataFrame(index=self.cellnames)
            self.metadata_df.index.name = 'cellname'
        
        # Ensure all cells are represented
        missing_cells = set(self.cellnames) - set(self.metadata_df.index)
        if missing_cells:
            missing_df = pd.DataFrame(index=list(missing_cells))
            self.metadata_df = pd.concat([self.metadata_df, missing_df])

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
        print("CHARMtools MultiCell3D object v2.0.0-anndata")
        print(f"Object contains {self.num_cells} cells starting with: {self.cellnames[:3]}")
        print(f"Resolutions: {self.resolutions}")
        print(f"Features: {self.features}")
        print(f"Object contains matrices: {list(self.matrices.keys())}")
        print(f"Object contains AnnData matrices: {list(self.adata_dict.keys())}")
        print(f"Metadata columns: {list(self.metadata_df.columns)}")
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
    
    def add_matrix_to_anndata(self, matrix_name: str, matrix_data: Union[np.ndarray, pd.DataFrame], 
                             obs_metadata: Optional[pd.DataFrame] = None,
                             var_metadata: Optional[pd.DataFrame] = None,
                             **kwargs):
        """
        Add a matrix to AnnData storage.
        
        Parameters:
        -----------
        matrix_name : str
            Name of the matrix (e.g., 'contact_matrix', 'distance_matrix')
        matrix_data : np.ndarray or pd.DataFrame
            The matrix data
        obs_metadata : pd.DataFrame, optional
            Metadata for observations (cells)
        var_metadata : pd.DataFrame, optional
            Metadata for variables (features/positions)
        **kwargs
            Additional arguments to store in uns
        """
        # Use provided metadata or default to class metadata
        if obs_metadata is None:
            obs_metadata = self.metadata_df.copy()
        
        # Create AnnData object
        if isinstance(matrix_data, pd.DataFrame):
            # Ensure obs and matrix rows match
            common_obs = list(set(obs_metadata.index) & set(matrix_data.index))
            if len(common_obs) != len(matrix_data.index):
                warnings.warn(f"Not all matrix rows have corresponding metadata. Using {len(common_obs)} common observations.")
            
            obs_metadata = obs_metadata.loc[common_obs]
            matrix_data = matrix_data.loc[common_obs]
            
            # Handle MultiIndex columns by joining with '-'
            if isinstance(matrix_data.columns, pd.MultiIndex):
                new_columns = ['-'.join(map(str, col)) for col in matrix_data.columns]
                matrix_data.columns = new_columns
            
            # Handle MultiIndex in obs_metadata index
            if isinstance(obs_metadata.index, pd.MultiIndex):
                obs_metadata = obs_metadata.reset_index()
                obs_metadata.index = obs_metadata.index.astype(str)
            
            adata = ad.AnnData(X=matrix_data.values, obs=obs_metadata, var=var_metadata)
            adata.var.index = matrix_data.columns
        else:
            # numpy array
            if matrix_data.shape[0] != obs_metadata.shape[0]:
                warnings.warn(f"Matrix rows ({matrix_data.shape[0]}) don't match metadata rows ({obs_metadata.shape[0]})")
                min_rows = min(matrix_data.shape[0], obs_metadata.shape[0])
                matrix_data = matrix_data[:min_rows]
                obs_metadata = obs_metadata.iloc[:min_rows]
            
            adata = ad.AnnData(X=matrix_data, obs=obs_metadata, var=var_metadata)
        
        # Add additional information to uns
        adata.uns['matrix_name'] = matrix_name
        for key, value in kwargs.items():
            adata.uns[key] = value
        
        # Store in dictionary
        self.adata_dict[matrix_name] = adata
        
        # Also store in legacy format for backward compatibility
        if isinstance(matrix_data, np.ndarray):
            self.matrices[matrix_name] = pd.DataFrame(matrix_data, 
                                                    index=obs_metadata.index,
                                                    columns=[f"feature_{i}" for i in range(matrix_data.shape[1])])
        else:
            self.matrices[matrix_name] = matrix_data
    
    def get_anndata(self, matrix_name: str) -> ad.AnnData:
        """
        Get AnnData object for a specific matrix.
        
        Parameters:
        -----------
        matrix_name : str
            Name of the matrix
            
        Returns:
        --------
        anndata.AnnData
            The AnnData object
        """
        if matrix_name not in self.adata_dict:
            raise ValueError(f"Matrix '{matrix_name}' not found in AnnData storage")
        return self.adata_dict[matrix_name]
    
    def list_anndata_matrices(self) -> List[str]:
        """
        List all available AnnData matrices.
        
        Returns:
        --------
        List[str]
            List of matrix names
        """
        return list(self.adata_dict.keys())
    
    def remove_anndata_matrix(self, matrix_name: str):
        """
        Remove an AnnData matrix.
        
        Parameters:
        -----------
        matrix_name : str
            Name of the matrix to remove
        """
        if matrix_name in self.adata_dict:
            del self.adata_dict[matrix_name]
        if matrix_name in self.matrices:
            del self.matrices[matrix_name]


# Helper functions
def _get_data(cell, **kwargs):
    df = cell.get_data(**kwargs)
    df["cellname"] = cell.cellname
    return df

def _get_feature_vec(cell, genome_coord, column_name):
    return cell.get_feature_vec(genome_coord, column_name)