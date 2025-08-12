import numpy as np
import pandas as pd
import tqdm
import concurrent.futures
from functools import partial
from scipy import stats
from .multicell3d_utils import _auto_genome_coord

class MultiCell3DAnalysis:
    """
    Analysis methods for MultiCell3D class.
    """
    
    def calc_insilico_GAM(self, genome_coord, cellnames=None, slice_width=3, num_slices=10):
        """
        Calculate the insilico GAM matrix for a given genomic coordinate.
        """
        if cellnames is None:
            cellnames = self.cellnames
        temp_cells = self.get_cell(cellnames)
        slices = []
        for cell in tqdm.tqdm(temp_cells):
            for i in range(num_slices):
                slices.append(cell.get_data_slice(genome_coord=genome_coord, if_rotate=True, if_full=True, slice_width=slice_width)["in_slice"])

        slices_temp = np.array(slices).T
        slices_temp = slices_temp.astype(int)
        mean_slices = slices_temp.mean(axis=1)
        mean_both = slices_temp @ slices_temp.T / slices_temp.shape[1]
        mean_product = np.outer(mean_slices, mean_slices)
        result = mean_both - mean_product

        return result
    
    def calc_insilico_SPRITE(self, genome_coord, cellnames=None, n_sphere=2000, sample_frac=0.5, radius=3):
        """
        Calculate the insilico SPRITE matrix for a given genomic coordinate.
        """
        if cellnames is None:
            cellnames = self.cellnames
        temp_cells = self.get_cell(cellnames)
        spheres = []
        for i in tqdm.tqdm(range(n_sphere)):
            # random select a cell
            cell = np.random.choice(temp_cells)
            spheres.append(cell.get_data_sphere(genome_coord, sample_frac=sample_frac)["in_ball"].values, radius=3)
        spheres_temp = np.array(spheres).T
        spheres_temp = spheres_temp.astype(int)
        mean_spheres = spheres_temp.mean(axis=1)
        mean_both = spheres_temp @ spheres_temp.T / spheres_temp.shape[1]
        mean_product = np.outer(mean_spheres, mean_spheres)
        result = mean_both - mean_product

        return result
 
    def calc_distance_matrix(self, genome_coord, cellnames=None, allele=True, combine=True, nproc=20):
        """
        Calculate the distance matrix between cells for a given genomic coordinate.
        If cellnames is None, all cells will be used.
        """
        if cellnames is None:
            cellnames = self.cellnames
        temp_cells = self.get_cell(cellnames)
        if allele:
            with concurrent.futures.ProcessPoolExecutor(nproc) as executor:
                results = list(executor.map(partial(_calc_distance, genome_coord=genome_coord), temp_cells))
        else:
            with concurrent.futures.ProcessPoolExecutor(nproc) as executor:
                results1 = list(executor.map(partial(_calc_distance, genome_coord=genome_coord.replace(":", "a:")), temp_cells))
                results2 = list(executor.map(partial(_calc_distance, genome_coord=genome_coord.replace(":", "b:")), temp_cells))
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
                results = list(executor.map(partial(_calc_proximity, genome_coord=genome_coord, distance_threshold=distance_threshold), temp_cells))
        else:
            with concurrent.futures.ProcessPoolExecutor(nproc) as executor:
                results1 = list(executor.map(partial(_calc_proximity, genome_coord=genome_coord.replace(":", "a:"), distance_threshold=distance_threshold), temp_cells))
                results2 = list(executor.map(partial(_calc_proximity, genome_coord=genome_coord.replace(":", "b:"), distance_threshold=distance_threshold), temp_cells))
            results = results1 + results2
        if combine:
            return np.nanmean(results, axis=0)
        else:
            return np.array(results)

    def calc_scABC_pred_gene(self, tss_genome_coord, flank=2000000, expression_key="UMIs_tss", distance_type="3d",
                             activity_keys=["atac_sum_in_radius_2", "ct_sum_in_radius_2"], cellnames=None, nproc=20, allele=True, return_data=False):
        """
        Calculate single-cell ABC prediction for gene expression.
        """
        if cellnames is None:   
            cellnames = self.cellnames
        temp_cells = self.get_cell(cellnames)

        if allele:
            with concurrent.futures.ProcessPoolExecutor(nproc) as executor:
                results = list(tqdm.tqdm(executor.map(partial(_calc_scABC_pred_gene,
                    tss_genome_coord=tss_genome_coord,
                    flank=flank,
                    expression_key=expression_key,
                    distance_type=distance_type,
                    activity_keys=activity_keys), temp_cells), 
                total=len(temp_cells)))
            expressions = np.array([result[0] for result in results])
            abc = np.array([result[1] for result in results])
        
        else:
            with concurrent.futures.ProcessPoolExecutor(nproc) as executor:
                results1 = list(tqdm.tqdm(executor.map(partial(_calc_scABC_pred_gene,
                    tss_genome_coord=tss_genome_coord.replace(":", "a:"),
                    flank=flank,
                    expression_key=expression_key,
                    distance_type=distance_type,
                    activity_keys=activity_keys), temp_cells), 
                total=len(temp_cells)))
                results2 = list(tqdm.tqdm(executor.map(partial(_calc_scABC_pred_gene,
                    tss_genome_coord=tss_genome_coord.replace(":", "b:"),
                    flank=flank,
                    expression_key=expression_key,
                    distance_type=distance_type,
                    activity_keys=activity_keys), temp_cells), 
                total=len(temp_cells)))
            results = results1 + results2
            expressions = np.array([result[0] for result in results])
            abc = np.array([result[1] for result in results])

        cors = []
        pvs = []

        for i in range(abc.shape[1]):
            cor = stats.spearmanr(expressions, abc[:, i], nan_policy='omit')
            cors.append(cor[0])
            pvs.append(cor[1])

        df = pd.DataFrame({"cor": cors, "p_value": pvs})
        chrom, start, end = _auto_genome_coord(tss_genome_coord)

        start = start - flank
        end = end + flank
        df["chrom"] = tss_genome_coord.split(":")[0]
        df["pos"] = np.arange(start, end, temp_cells[0].resolution)
        if return_data:
            return df[["chrom", "pos", "cor", "p_value"]], expressions, abc
        else:
            return df[["chrom", "pos", "cor", "p_value"]]

    def calc_feature_matrix(self, genome_coord, feature, cells=None):
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

        return np.nanmean(feature_mats, axis=0), np.array(feature_vecs)
    
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
            feature_mat, feature_vec = cell.calc_feature_proximity_matrix(genome_coord, feature, distance_threshold)
            feature_mats.append(feature_mat)
            feature_vecs.append(feature_vec)
        return np.nanmean(feature_mats, axis=0), np.array(feature_vecs)

    def simple_diff(self, cellnames_group1, cellnames_group2, genome_coord=None):
        """
        Calculate simple difference between two groups.
        """
        pass


# Helper functions for analysis
def _calc_proximity(cell, genome_coord, distance_threshold):
    mat = cell.calc_distance_matrix(genome_coord) < distance_threshold
    return mat

def _calc_distance(cell, genome_coord):
    mat = cell.calc_distance_matrix(genome_coord)
    return mat

def _calc_scABC_pred_gene(cell, tss_genome_coord, flank, expression_key, activity_keys, distance_type):
    expression, abc = cell.calc_scABC_pred_gene(tss_genome_coord, flank, expression_key, activity_keys, distance_type)
    return [expression, abc]

def chromatic(mat, vec):
    """
    Calculate chromatic interaction.
    resij = matij * veci * vecj
    """
    res = np.zeros(mat.shape)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if i != j:
                res[i, j] = mat[i, j] * vec[i] * vec[j]
            else:
                res[i, j] = mat[i, j] * vec[i]
    return res

def chromatic_diff():
    """
    TODO: Calculate chromatic difference.
    """
    pass