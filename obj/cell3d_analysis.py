# Cell3D Analysis Module - Advanced analysis methods
import warnings
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr
from scipy.spatial.distance import pdist, squareform

class Cell3DAnalysis:
    """Advanced analysis methods for Cell3D objects"""
    
    def calculate_RMSD(self, other_cell3d, align=True):
        """Calculate Root Mean Square Deviation between two Cell3D objects"""
        if self.on_disk:
            self.to_memory()
        if other_cell3d.on_disk:
            other_cell3d.to_memory()
        
        # Get common genomic positions
        self_coords = self.tdg[['chrom', 'pos', 'x', 'y', 'z']].copy()
        other_coords = other_cell3d.tdg[['chrom', 'pos', 'x', 'y', 'z']].copy()
        
        # Merge on common positions
        merged = pd.merge(self_coords, other_coords, on=['chrom', 'pos'], suffixes=('_self', '_other'))
        
        if len(merged) == 0:
            raise ValueError("No common genomic positions found between the two Cell3D objects")
        
        # Extract coordinates
        coords_self = merged[['x_self', 'y_self', 'z_self']].values
        coords_other = merged[['x_other', 'y_other', 'z_other']].values
        
        if align:
            # Align structures using Procrustes analysis
            coords_other_aligned = self._procrustes_align(coords_self, coords_other)
        else:
            coords_other_aligned = coords_other
        
        # Calculate RMSD
        diff = coords_self - coords_other_aligned
        rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
        
        return rmsd
    
    def _procrustes_align(self, reference, target):
        """Align target coordinates to reference using Procrustes analysis"""
        # Center both coordinate sets
        ref_centered = reference - np.mean(reference, axis=0)
        target_centered = target - np.mean(target, axis=0)
        
        # Calculate cross-covariance matrix
        H = target_centered.T @ ref_centered
        
        # Singular Value Decomposition
        U, S, Vt = np.linalg.svd(H)
        
        # Calculate rotation matrix
        R = Vt.T @ U.T
        
        # Ensure proper rotation (det(R) = 1)
        if np.linalg.det(R) < 0:
            Vt[-1, :] *= -1
            R = Vt.T @ U.T
        
        # Apply rotation and translation
        target_aligned = (R @ target_centered.T).T + np.mean(reference, axis=0)
        
        return target_aligned
    
    def calc_feature_distances_v1(self, feature_name, reference_value=None, method='euclidean'):
        """Calculate distances between feature values (legacy version)"""
        if self.on_disk:
            self.to_memory()
        
        if feature_name not in self.tdg.columns:
            raise ValueError(f"Feature '{feature_name}' not found in data")
        
        feature_values = self.tdg[feature_name].values
        
        if reference_value is None:
            reference_value = np.mean(feature_values)
        
        if method == 'euclidean':
            distances = np.abs(feature_values - reference_value)
        elif method == 'manhattan':
            distances = np.abs(feature_values - reference_value)
        elif method == 'relative':
            distances = np.abs(feature_values - reference_value) / (np.abs(reference_value) + 1e-10)
        else:
            raise ValueError(f"Unsupported distance method: {method}")
        
        return distances
    
    def _auto_genome_coord(self, genome_coord_list=None):
        """Automatically determine genome coordinates if not provided"""
        if self.on_disk:
            self.to_memory()
        
        if genome_coord_list is None:
            genome_coord_list = sorted(self.tdg['chrom'].unique())
        
        return genome_coord_list
    
    def mat_cor_with_na(self, matrix1, matrix2, method='pearson'):
        """Calculate correlation between two matrices handling NaN values"""
        # Flatten matrices and remove NaN pairs
        flat1 = matrix1.flatten()
        flat2 = matrix2.flatten()
        
        # Create mask for valid (non-NaN) values
        valid_mask = ~(np.isnan(flat1) | np.isnan(flat2))
        
        if np.sum(valid_mask) < 2:
            return np.nan
        
        valid_flat1 = flat1[valid_mask]
        valid_flat2 = flat2[valid_mask]
        
        if method == 'pearson':
            correlation, p_value = pearsonr(valid_flat1, valid_flat2)
        elif method == 'spearman':
            correlation, p_value = spearmanr(valid_flat1, valid_flat2)
        else:
            raise ValueError(f"Unsupported correlation method: {method}")
        
        return correlation
    
    def _roated_hic_mat(self, hic_matrix, rotation_angle=0):
        """Rotate HiC matrix (placeholder for matrix rotation)"""
        # This is a placeholder implementation
        # In practice, this would involve complex matrix transformations
        if rotation_angle == 0:
            return hic_matrix
        
        # For now, just return the original matrix
        warnings.warn("Matrix rotation not fully implemented")
        return hic_matrix
    
    def plotDistanceMatrix(self, genome_coord=None, method='euclidean', 
                          log_transform=False, figsize=(8, 6), save_path=None):
        """Plot distance matrix with additional options"""
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError("matplotlib is required for plotting")
        
        if self.on_disk:
            self.to_memory()
        
        # Get data
        if genome_coord is None:
            plot_data = self.tdg.copy()
        else:
            plot_data = self.tdg.query(f"chrom == '{genome_coord}'").copy()
        
        # Calculate distance matrix
        coords = plot_data[['x', 'y', 'z']].values
        distances = squareform(pdist(coords, metric=method))
        
        if log_transform:
            distances = np.log10(distances + 1e-10)
        
        # Create plot
        fig, ax = plt.subplots(figsize=figsize)
        
        im = ax.imshow(distances, cmap='viridis', aspect='auto')
        plt.colorbar(im, ax=ax, label='Distance')
        
        ax.set_title(f'Distance Matrix - {genome_coord or "All Chromosomes"}')
        ax.set_xlabel('Genomic Position Index')
        ax.set_ylabel('Genomic Position Index')
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        plt.tight_layout()
        return fig, ax
    
    def compare_structures(self, other_cell3d, metrics=['rmsd', 'correlation'], 
                         features_to_compare=None):
        """Comprehensive comparison between two Cell3D structures"""
        if self.on_disk:
            self.to_memory()
        if other_cell3d.on_disk:
            other_cell3d.to_memory()
        
        results = {}
        
        # RMSD comparison
        if 'rmsd' in metrics:
            try:
                rmsd = self.calculate_RMSD(other_cell3d)
                results['rmsd'] = rmsd
            except Exception as e:
                results['rmsd'] = f"Error: {str(e)}"
        
        # Correlation comparison
        if 'correlation' in metrics:
            # Get common positions
            self_coords = self.tdg[['chrom', 'pos', 'x', 'y', 'z']].copy()
            other_coords = other_cell3d.tdg[['chrom', 'pos', 'x', 'y', 'z']].copy()
            merged = pd.merge(self_coords, other_coords, on=['chrom', 'pos'], suffixes=('_self', '_other'))
            
            if len(merged) > 0:
                # Coordinate correlations
                x_corr = self.mat_cor_with_na(merged['x_self'].values, merged['x_other'].values)
                y_corr = self.mat_cor_with_na(merged['y_self'].values, merged['y_other'].values)
                z_corr = self.mat_cor_with_na(merged['z_self'].values, merged['z_other'].values)
                
                results['coordinate_correlations'] = {
                    'x': x_corr,
                    'y': y_corr,
                    'z': z_corr
                }
        
        # Feature comparison
        if features_to_compare:
            feature_correlations = {}
            for feature in features_to_compare:
                if feature in self.tdg.columns and feature in other_cell3d.tdg.columns:
                    # Merge on common positions
                    self_feature = self.tdg[['chrom', 'pos', feature]].copy()
                    other_feature = other_cell3d.tdg[['chrom', 'pos', feature]].copy()
                    merged_feature = pd.merge(self_feature, other_feature, on=['chrom', 'pos'], suffixes=('_self', '_other'))
                    
                    if len(merged_feature) > 0:
                        corr = self.mat_cor_with_na(merged_feature[f'{feature}_self'].values, 
                                                  merged_feature[f'{feature}_other'].values)
                        feature_correlations[feature] = corr
            
            results['feature_correlations'] = feature_correlations
        
        return results
    
    def calculate_structural_metrics(self):
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
    
    def analyze_chromosome_territories(self):
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