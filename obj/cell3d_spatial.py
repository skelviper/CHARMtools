# Cell3D Spatial Module - Spatial analysis and clustering
import warnings
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import DBSCAN, KMeans
from sklearn.neighbors import NearestNeighbors

class Cell3DSpatial:
    """
    Spatial analysis and clustering for Cell3D objects
    基于结构计算各类统计信息。
    """
    
    def calc_intermingling(self, chrom1, chrom2, radius=1.0):
        """Calculate intermingling between two chromosomes"""
        if self.on_disk:
            self.to_memory()
        
        # Get data for both chromosomes
        chrom1_data = self.tdg.query(f"chrom == '{chrom1}'")[['x', 'y', 'z']].values
        chrom2_data = self.tdg.query(f"chrom == '{chrom2}'")[['x', 'y', 'z']].values
        
        if len(chrom1_data) == 0 or len(chrom2_data) == 0:
            return 0.0
        
        # Calculate distances between all pairs
        intermingling_count = 0
        total_pairs = 0
        
        for point1 in chrom1_data:
            for point2 in chrom2_data:
                distance = np.linalg.norm(point1 - point2)
                if distance <= radius:
                    intermingling_count += 1
                total_pairs += 1
        
        return intermingling_count / total_pairs if total_pairs > 0 else 0.0

    def calc_3D_cluster(self, method='dbscan', eps=0.5, min_samples=5, n_clusters=8, **kwargs):
        """Perform 3D clustering on the data"""
        if self.on_disk:
            self.to_memory()
        
        coords = self.tdg[['x', 'y', 'z']].values
        
        if method == 'dbscan':
            clustering = DBSCAN(eps=eps, min_samples=min_samples, **kwargs)
            cluster_labels = clustering.fit_predict(coords)
        elif method == 'kmeans':
            clustering = KMeans(n_clusters=n_clusters, **kwargs)
            cluster_labels = clustering.fit_predict(coords)
        else:
            raise ValueError(f"Unsupported clustering method: {method}")
        
        self.tdg['cluster'] = cluster_labels
        self.features.append('cluster')
        
        return cluster_labels

    def point_to_3Dcluster(self, point_coords, cluster_column='cluster'):
        """Find which cluster a given point belongs to"""
        if self.on_disk:
            self.to_memory()
        
        if cluster_column not in self.tdg.columns:
            raise ValueError(f"Cluster column '{cluster_column}' not found")
        
        coords = self.tdg[['x', 'y', 'z']].values
        
        # Find the nearest point in the dataset
        distances = np.linalg.norm(coords - point_coords, axis=1)
        nearest_idx = np.argmin(distances)
        
        return self.tdg.iloc[nearest_idx][cluster_column]

    def calc_scABC_pred_gene(self, gene_coords, radius=1.0, feature_name='accessibility'):
        """Calculate single-cell ABC prediction for genes"""
        if self.on_disk:
            self.to_memory()
        
        if feature_name not in self.tdg.columns:
            raise ValueError(f"Feature '{feature_name}' not found")
        
        predictions = []
        
        for gene_coord in gene_coords:
            # Find points within radius of gene
            coords = self.tdg[['x', 'y', 'z']].values
            distances = np.linalg.norm(coords - gene_coord, axis=1)
            nearby_mask = distances <= radius
            
            # Calculate ABC score (simplified version)
            if np.any(nearby_mask):
                nearby_features = self.tdg.loc[nearby_mask, feature_name]
                abc_score = np.sum(nearby_features * (1.0 / (distances[nearby_mask] + 1e-10)))
            else:
                abc_score = 0.0
            
            predictions.append(abc_score)
        
        return np.array(predictions)

    def expandStructure(self, expansion_factor=1.5):
        """Expand the 3D structure by a given factor"""
        if self.on_disk:
            self.to_memory()
        
        # Calculate center of mass
        center_x = np.mean(self.tdg['x'])
        center_y = np.mean(self.tdg['y'])
        center_z = np.mean(self.tdg['z'])
        
        # Expand coordinates relative to center
        self.tdg['x'] = center_x + (self.tdg['x'] - center_x) * expansion_factor
        self.tdg['y'] = center_y + (self.tdg['y'] - center_y) * expansion_factor
        self.tdg['z'] = center_z + (self.tdg['z'] - center_z) * expansion_factor

    def calc_radial_position(self, center=None, returnValue = False):
        """Calculate radial position of each point from center"""
        if self.on_disk:
            self.to_memory()
        
        if self.tdg is None:
            raise ValueError("No 3D structure data available. Please ensure the tdg data is properly loaded.")
        
        if center is None:
            # Use center of mass as default center
            center = [
                np.mean(self.tdg['x']),
                np.mean(self.tdg['y']),
                np.mean(self.tdg['z'])
            ]
        
        # Calculate radial distances
        radial_distances = np.sqrt(
            (self.tdg['x'] - center[0])**2 + 
            (self.tdg['y'] - center[1])**2 + 
            (self.tdg['z'] - center[2])**2
        )
        
        self.tdg['radial_position'] = radial_distances
        self.features.append('radial_position')
        
        if returnValue:
            return radial_distances
        else:
            return None

    def feature_radial_distribution(self, feature_name, n_bins=10, center=None):
        """Calculate radial distribution of a feature"""
        if self.on_disk:
            self.to_memory()
        
        if feature_name not in self.tdg.columns:
            raise ValueError(f"Feature '{feature_name}' not found")
        
        # Calculate radial positions if not already done
        if 'radial_position' not in self.tdg.columns:
            self.calc_radial_position(center=center)
        
        # Create radial bins
        max_radius = np.max(self.tdg['radial_position'])
        bin_edges = np.linspace(0, max_radius, n_bins + 1)
        
        # Calculate feature distribution in each bin
        distribution = []
        for i in range(n_bins):
            mask = (self.tdg['radial_position'] >= bin_edges[i]) & \
                   (self.tdg['radial_position'] < bin_edges[i + 1])
            
            if np.any(mask):
                bin_mean = np.mean(self.tdg.loc[mask, feature_name])
            else:
                bin_mean = 0.0
            
            distribution.append(bin_mean)
        
        return np.array(distribution), bin_edges

    def calc_radius_gyration(self, by_chromosome=False):
        """Calculate radius of gyration"""
        if self.on_disk:
            self.to_memory()
        
        if by_chromosome:
            rg_values = {}
            for chrom in self.tdg['chrom'].unique():
                chrom_data = self.tdg.query(f"chrom == '{chrom}'")[['x', 'y', 'z']].values
                if len(chrom_data) > 0:
                    center = np.mean(chrom_data, axis=0)
                    distances_sq = np.sum((chrom_data - center)**2, axis=1)
                    rg = np.sqrt(np.mean(distances_sq))
                    rg_values[chrom] = rg
            return rg_values
        else:
            coords = self.tdg[['x', 'y', 'z']].values
            center = np.mean(coords, axis=0)
            distances_sq = np.sum((coords - center)**2, axis=1)
            rg = np.sqrt(np.mean(distances_sq))
            return rg

    def calc_feature_distances(self, feature_name, reference_points=None):
        """Calculate distances from feature-positive points to reference points"""
        if self.on_disk:
            self.to_memory()
        
        if feature_name not in self.tdg.columns:
            raise ValueError(f"Feature '{feature_name}' not found")
        
        # Get feature-positive points
        feature_positive = self.tdg[self.tdg[feature_name] > 0][['x', 'y', 'z']].values
        
        if reference_points is None:
            # Use center of mass as reference
            reference_points = [np.mean(self.tdg[['x', 'y', 'z']].values, axis=0)]
        
        distances = []
        for ref_point in reference_points:
            point_distances = np.linalg.norm(feature_positive - ref_point, axis=1)
            distances.append(point_distances)
        
        return distances

    def calc_singlecell_compartment(self, feature_name, method='kmeans', n_compartments=2):
        """Calculate single-cell compartments based on features"""
        if self.on_disk:
            self.to_memory()
        
        if feature_name not in self.tdg.columns:
            raise ValueError(f"Feature '{feature_name}' not found")
        
        feature_values = self.tdg[feature_name].values.reshape(-1, 1)
        
        if method == 'kmeans':
            clustering = KMeans(n_clusters=n_compartments, random_state=42)
            compartment_labels = clustering.fit_predict(feature_values)
        elif method == 'threshold':
            # Simple threshold-based compartmentalization
            threshold = np.median(feature_values)
            compartment_labels = (feature_values.flatten() > threshold).astype(int)
        else:
            raise ValueError(f"Unsupported compartmentalization method: {method}")
        
        self.tdg[f'{feature_name}_compartment'] = compartment_labels
        self.features.append(f'{feature_name}_compartment')
        
        return compartment_labels

    def calc_volume(self, method='convex_hull'):
        """Calculate volume of the 3D structure"""
        if self.on_disk:
            self.to_memory()
        
        coords = self.tdg[['x', 'y', 'z']].values
        
        if method == 'convex_hull':
            try:
                from scipy.spatial import ConvexHull
                hull = ConvexHull(coords)
                return hull.volume
            except ImportError:
                warnings.warn("scipy.spatial.ConvexHull not available, using bounding box method")
                method = 'bounding_box'
        
        if method == 'bounding_box':
            # Calculate bounding box volume
            min_coords = np.min(coords, axis=0)
            max_coords = np.max(coords, axis=0)
            volume = np.prod(max_coords - min_coords)
            return volume
        
        raise ValueError(f"Unsupported volume calculation method: {method}")

    def calc_structural_metrics(self):
        """Calculate various structural metrics"""
        if self.on_disk:
            self.to_memory()
        
        coords = self.tdg[['x', 'y', 'z']].values
        
        metrics = {}
        
        # Center of mass
        center_of_mass = np.mean(coords, axis=0)
        metrics['center_of_mass'] = center_of_mass
        
        # Radius of gyration
        rg = self.calc_radius_gyration()
        metrics['radius_of_gyration'] = rg
        
        # Volume
        try:
            volume = self.calc_volume()
            metrics['volume'] = volume
        except Exception as e:
            metrics['volume'] = f"Error: {str(e)}"
        
        # Coordinate ranges
        metrics['coordinate_ranges'] = {
            'x_range': np.ptp(coords[:, 0]),
            'y_range': np.ptp(coords[:, 1]),
            'z_range': np.ptp(coords[:, 2])
        }
        
        # Density (points per unit volume)
        try:
            volume_val = metrics['volume']
            if isinstance(volume_val, (int, float)):
                metrics['density'] = len(coords) / volume_val
        except:
            pass
        
        return metrics
    
    def calc_chromosome_territories_metrics(self):
        """Analyze chromosome territory organization"""
        if self.on_disk:
            self.to_memory()
        
        territory_analysis = {}
        
        for chrom in self.tdg['chrom'].unique():
            chrom_data = self.tdg.query(f"chrom == '{chrom}'")[['x', 'y', 'z']].values
            
            if len(chrom_data) > 0:
                # Territory center
                territory_center = np.mean(chrom_data, axis=0)
                
                # Territory radius of gyration
                distances_from_center = np.linalg.norm(chrom_data - territory_center, axis=1)
                territory_rg = np.sqrt(np.mean(distances_from_center**2))
                
                # Territory volume (convex hull if possible)
                try:
                    from scipy.spatial import ConvexHull
                    if len(chrom_data) >= 4:  # Minimum points for 3D convex hull
                        hull = ConvexHull(chrom_data)
                        territory_volume = hull.volume
                    else:
                        territory_volume = None
                except ImportError:
                    territory_volume = None
                
                territory_analysis[chrom] = {
                    'center': territory_center,
                    'radius_of_gyration': territory_rg,
                    'volume': territory_volume,
                    'n_points': len(chrom_data)
                }
        
        return territory_analysis