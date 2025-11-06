# Cell3D Data Module - Data processing and retrieval
import warnings
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
from scipy.spatial.distance import pdist, squareform
from sklearn.neighbors import NearestNeighbors
from scipy.spatial import cKDTree
from scipy.stats import rankdata
from ..utils.helper import auto_genome_coord

class Cell3DData:
    """Data processing and retrieval for Cell3D objects"""
    
    def _sparse_to_dense(self,genome_coord,sparse_df):
        chrom,start,end = auto_genome_coord(genome_coord)

        if start is None:
            start = 0
        if end is None:
            if self.chrom_length is None:
                raise ValueError("Chrom length is not available, please run add_chrom_length first")
            end = self.chrom_length.set_index("chrom").to_dict()["size"][chrom]

        positions = np.arange(start,end,self.resolution)
        dense_df = pd.DataFrame(positions,columns=['pos'])
        sparse_df['pos'] = sparse_df['pos'].astype(int)
        dense_df = pd.merge(dense_df, sparse_df, on='pos', how='left')
        dense_df['chrom'] = chrom

        columns_order = ['chrom', 'pos'] + [col for col in dense_df.columns if col not in ['chrom', 'pos']]
        dense_df = dense_df[columns_order]
        return dense_df

    def get_data(self, genome_coord=None, if_dense=False, if_sort=True):
        """Get data from the Cell3D object"""
        if self.on_disk:
            self.to_memory()
        
        if genome_coord is None:
            df = self.tdg.copy()
        else:
            df = self.tdg.query(f"chrom == '{genome_coord}'").copy()
        
        if if_sort:
            df = df.sort_values(["chrom", "pos"])

        if if_dense:
            df = self._sparse_to_dense(genome_coord, df)
        
        return df

    def get_data_slice(self, genome_coord, start_pos, end_pos):
        """Get a slice of data based on genomic coordinates"""
        if self.on_disk:
            self.to_memory()
        
        df = self.tdg.query(f"chrom == '{genome_coord}' and pos >= {start_pos} and pos <= {end_pos}").copy()
        return df.sort_values(["chrom", "pos"])

    def get_data_cluster(self, cluster_id, cluster_column="cluster"):
        """Get data for a specific cluster"""
        if self.on_disk:
            self.to_memory()
        
        if cluster_column not in self.tdg.columns:
            raise ValueError(f"Cluster column '{cluster_column}' not found in data")
        
        df = self.tdg.query(f"{cluster_column} == {cluster_id}").copy()
        return df.sort_values(["chrom", "pos"])

    def get_data_sphere(self, center_x, center_y, center_z, radius):
        """Get data points within a sphere"""
        if self.on_disk:
            self.to_memory()
        
        df = self.tdg.copy()
        distances = np.sqrt((df['x'] - center_x)**2 + (df['y'] - center_y)**2 + (df['z'] - center_z)**2)
        return df[distances <= radius].sort_values(["chrom", "pos"])

    def get_feature_vec(self, feature_name):
        """Get feature vector for a specific feature"""
        if self.on_disk:
            self.to_memory()
        
        if feature_name not in self.tdg.columns:
            raise ValueError(f"Feature '{feature_name}' not found in data")
        
        return self.tdg[feature_name].values

    def add_knn_density(self, k=10, column_name="knn_density"):
        """Add k-nearest neighbors density to the data"""
        if self.on_disk:
            self.to_memory()
        
        coords = self.tdg[['x', 'y', 'z']].values
        nbrs = NearestNeighbors(n_neighbors=k+1, algorithm='ball_tree').fit(coords)
        distances, indices = nbrs.kneighbors(coords)
        
        # Calculate density as inverse of mean distance to k nearest neighbors
        # Exclude the first neighbor (which is the point itself)
        mean_distances = np.mean(distances[:, 1:], axis=1)
        density = 1.0 / (mean_distances + 1e-10)  # Add small epsilon to avoid division by zero
        
        self.tdg[column_name] = density
        self.features.append(column_name)

    # def add_feature_in_radius(self, feature_name, radius, new_column_name=None, aggregation='sum'):
    #     """Add aggregated feature values within a radius for each point"""
    #     if self.on_disk:
    #         self.to_memory()
        
    #     if feature_name not in self.tdg.columns:
    #         raise ValueError(f"Feature '{feature_name}' not found in data")
        
    #     if new_column_name is None:
    #         new_column_name = f"{feature_name}_radius_{radius}"
        
    #     coords = self.tdg[['x', 'y', 'z']].values
    #     feature_values = self.tdg[feature_name].values
        
    #     # Calculate pairwise distances
    #     distances = squareform(pdist(coords))
        
    #     # For each point, find neighbors within radius and aggregate feature values
    #     aggregated_values = []
    #     for i in range(len(coords)):
    #         neighbors_mask = distances[i] <= radius
    #         neighbor_values = feature_values[neighbors_mask]
            
    #         if aggregation == 'sum':
    #             agg_value = np.sum(neighbor_values)
    #         elif aggregation == 'mean':
    #             agg_value = np.mean(neighbor_values)
    #         elif aggregation == 'max':
    #             agg_value = np.max(neighbor_values)
    #         elif aggregation == 'min':
    #             agg_value = np.min(neighbor_values)
    #         elif aggregation == 'count':
    #             agg_value = len(neighbor_values)
    #         else:
    #             raise ValueError(f"Unsupported aggregation method: {aggregation}")
            
    #         aggregated_values.append(agg_value)
        
    #     self.tdg[new_column_name] = aggregated_values
    #     self.features.append(new_column_name)

    def add_feature_in_radius(self, feature, radius, aggregation="sum",
                            include_self=True, rank=False, new_column_name=None,
                            nan_policy="propagate"):
        """
        KDTree 半径聚合：支持 sum/mean/max/min/count，控制是否包含自身，可选秩归一化。
        nan_policy: "propagate" | "omit"
        """

        if self.on_disk:
            self.to_memory()
        if feature not in self.tdg.columns:
            raise ValueError(f"Feature '{feature}' not found")

        if not hasattr(self, "kdtree") or self.kdtree is None:
            coords = self.tdg.filter(regex=r"^(x|y|z)$").values  # 或自定义坐标列
            self.kdtree = cKDTree(coords)

        indices_list = self.kdtree.query_ball_tree(self.kdtree, r=radius)

        if not include_self:
            indices_list = [[j for j in idxs if j != i] for i, idxs in enumerate(indices_list)]

        s = self.tdg[feature]
        vals = []
        for idxs in indices_list:
            if len(idxs) == 0:
                group = s.iloc[idxs]  # 空
            else:
                group = s.iloc[idxs]
            if nan_policy == "omit":
                group = group.dropna()

            if aggregation == "sum":
                vals.append(group.sum())
            elif aggregation == "mean":
                vals.append(group.mean())
            elif aggregation == "max":
                vals.append(group.max())
            elif aggregation == "min":
                vals.append(group.min())
            elif aggregation == "count":
                vals.append(len(group))
            else:
                raise ValueError("aggregation must be one of sum/mean/max/min/count")

        if rank:
            vals = rankdata(vals, nan_policy="omit") / float(np.isfinite(vals).sum())

        col = new_column_name or f"{feature}_{aggregation}_in_radius_{radius}"
        self.tdg[col] = vals
        self.features.append(col)



    def calc_distance_matrix(self, genome_coord=None, metric='euclidean'):
        """Calculate distance matrix between points
        
        Parameters:
            genome_coord: str, format like chrom:start-end or list/tuple of chrom,start,end.
                         If provided, calculates distance matrix for specific genomic region.
                         If None, calculates for all points using the metric parameter.
            metric: str, distance metric for scipy.spatial.distance.pdist (default: 'euclidean')
                   Only used when genome_coord is None.
        
        Returns:
            distance_matrix: numpy array of pairwise distances
        """
        if self.on_disk:
            self.to_memory()
        
        # If genome_coord is provided, filter data for that region
        if genome_coord is not None:
            coord_info = auto_genome_coord(genome_coord)
            if coord_info is None:
                raise ValueError(f"Invalid genome coordinate format: {genome_coord}")
            
            chrom, start, end = coord_info
            
            # Filter data for the specified genomic region
            data_subset = self.tdg.query("chrom == @chrom & pos >= @start & pos <= @end").copy()
            
            if len(data_subset) == 0:
                raise ValueError(f"No data found for genomic region: {genome_coord}")
            
            coords = data_subset[['x', 'y', 'z']].values
        else:
            # Use all points
            coords = self.tdg[['x', 'y', 'z']].values
        
        # Calculate distance matrix
        if metric == 'euclidean':
            distances = squareform(pdist(coords, metric='euclidean'))
        else:
            distances = squareform(pdist(coords, metric=metric))
        
        return distances

    def calc_feature_matrix(self, feature_name, metric='correlation'):
        """Calculate feature similarity/distance matrix"""
        if self.on_disk:
            self.to_memory()
        
        if feature_name not in self.tdg.columns:
            raise ValueError(f"Feature '{feature_name}' not found in data")
        
        feature_values = self.tdg[feature_name].values.reshape(-1, 1)
        
        if metric == 'correlation':
            # Calculate correlation matrix
            corr_matrix = np.corrcoef(feature_values.flatten())
            return corr_matrix
        else:
            # Use scipy's pdist for other metrics
            distances = squareform(pdist(feature_values, metric=metric))
            return distances

    def calc_feature_proximity_matrix(self, features_list, metric='euclidean'):
        """Calculate proximity matrix based on multiple features"""
        if self.on_disk:
            self.to_memory()
        
        # Check if all features exist
        missing_features = [f for f in features_list if f not in self.tdg.columns]
        if missing_features:
            raise ValueError(f"Features not found: {missing_features}")
        
        feature_matrix = self.tdg[features_list].values
        
        if metric == 'euclidean':
            distances = squareform(pdist(feature_matrix, metric='euclidean'))
        else:
            distances = squareform(pdist(feature_matrix, metric=metric))
        
        return distances