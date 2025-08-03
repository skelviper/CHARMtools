# Cell3D Data Module - Data processing and retrieval
import warnings
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
from scipy.spatial.distance import pdist, squareform
from sklearn.neighbors import NearestNeighbors

class Cell3DData:
    """Data processing and retrieval for Cell3D objects"""
    
    @staticmethod
    def _sparse_to_dense(sparse_matrix, shape):
        """Convert sparse matrix to dense format"""
        dense_matrix = np.full(shape, np.nan)
        for i, j, value in zip(sparse_matrix['i'], sparse_matrix['j'], sparse_matrix['data']):
            dense_matrix[i, j] = value
            if i != j:  # Make symmetric
                dense_matrix[j, i] = value
        return dense_matrix

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
        
        if if_dense and hasattr(self, 'hic_matrix') and self.hic_matrix is not None:
            # Convert sparse HiC matrix to dense if requested
            if genome_coord:
                chrom_indices = df.index.tolist()
                dense_hic = self._sparse_to_dense(self.hic_matrix, (len(self.tdg), len(self.tdg)))
                df['hic_dense'] = [dense_hic[i] for i in chrom_indices]
        
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

    def add_feature_in_radius(self, feature_name, radius, new_column_name=None, aggregation='sum'):
        """Add aggregated feature values within a radius for each point"""
        if self.on_disk:
            self.to_memory()
        
        if feature_name not in self.tdg.columns:
            raise ValueError(f"Feature '{feature_name}' not found in data")
        
        if new_column_name is None:
            new_column_name = f"{feature_name}_radius_{radius}"
        
        coords = self.tdg[['x', 'y', 'z']].values
        feature_values = self.tdg[feature_name].values
        
        # Calculate pairwise distances
        distances = squareform(pdist(coords))
        
        # For each point, find neighbors within radius and aggregate feature values
        aggregated_values = []
        for i in range(len(coords)):
            neighbors_mask = distances[i] <= radius
            neighbor_values = feature_values[neighbors_mask]
            
            if aggregation == 'sum':
                agg_value = np.sum(neighbor_values)
            elif aggregation == 'mean':
                agg_value = np.mean(neighbor_values)
            elif aggregation == 'max':
                agg_value = np.max(neighbor_values)
            elif aggregation == 'min':
                agg_value = np.min(neighbor_values)
            elif aggregation == 'count':
                agg_value = len(neighbor_values)
            else:
                raise ValueError(f"Unsupported aggregation method: {aggregation}")
            
            aggregated_values.append(agg_value)
        
        self.tdg[new_column_name] = aggregated_values
        self.features.append(new_column_name)

    def calc_distance_matrix(self, metric='euclidean'):
        """Calculate distance matrix between all points"""
        if self.on_disk:
            self.to_memory()
        
        coords = self.tdg[['x', 'y', 'z']].values
        
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