# Cell3D Utils Module - Utility functions and helpers
import warnings
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix

class Cell3DUtils:
    """Utility functions for Cell3D objects"""
    
    @staticmethod
    def _point_cloud_rotation(coords, rotation_matrix):
        """Apply rotation matrix to point cloud coordinates"""
        if coords.shape[1] != 3:
            raise ValueError("Coordinates must be Nx3 array")
        
        if rotation_matrix.shape != (3, 3):
            raise ValueError("Rotation matrix must be 3x3")
        
        # Apply rotation
        rotated_coords = coords @ rotation_matrix.T
        
        return rotated_coords
    
    @staticmethod
    def create_rotation_matrix(axis, angle):
        """Create rotation matrix around specified axis"""
        axis = axis.lower()
        angle_rad = np.radians(angle)
        
        if axis == 'x':
            rotation_matrix = np.array([
                [1, 0, 0],
                [0, np.cos(angle_rad), -np.sin(angle_rad)],
                [0, np.sin(angle_rad), np.cos(angle_rad)]
            ])
        elif axis == 'y':
            rotation_matrix = np.array([
                [np.cos(angle_rad), 0, np.sin(angle_rad)],
                [0, 1, 0],
                [-np.sin(angle_rad), 0, np.cos(angle_rad)]
            ])
        elif axis == 'z':
            rotation_matrix = np.array([
                [np.cos(angle_rad), -np.sin(angle_rad), 0],
                [np.sin(angle_rad), np.cos(angle_rad), 0],
                [0, 0, 1]
            ])
        else:
            raise ValueError("Axis must be 'x', 'y', or 'z'")
        
        return rotation_matrix
    
    @staticmethod
    def create_bins_genetable(genes_df, resolution=40000, genome_size_file=None):
        """Create bins gene table from genes dataframe"""
        genes = genes_df.copy()
        
        # Ensure required columns
        required_cols = ['chrom', 'start', 'end', 'gene_id']
        if not all(col in genes.columns for col in required_cols):
            raise ValueError(f"Genes dataframe must contain columns: {required_cols}")
        
        # Filter for standard chromosomes
        genes = genes[genes['chrom'].str.contains('chr')].copy()
        
        # Create bins
        bins_list = []
        
        for chrom in genes['chrom'].unique():
            chrom_genes = genes[genes['chrom'] == chrom].copy()
            
            if len(chrom_genes) == 0:
                continue
            
            # Determine chromosome length
            if genome_size_file is not None:
                try:
                    genome_sizes = pd.read_csv(genome_size_file, sep='\t', header=None, names=['chrom', 'size'])
                    chrom_size = genome_sizes[genome_sizes['chrom'] == chrom]['size'].iloc[0]
                except:
                    chrom_size = chrom_genes['end'].max()
            else:
                chrom_size = chrom_genes['end'].max()
            
            # Create bins for this chromosome
            n_bins = int(np.ceil(chrom_size / resolution))
            
            for i in range(n_bins):
                bin_start = i * resolution
                bin_end = min((i + 1) * resolution, chrom_size)
                
                # Find genes overlapping with this bin
                overlapping_genes = chrom_genes[
                    (chrom_genes['start'] < bin_end) & (chrom_genes['end'] > bin_start)
                ]
                
                # Create bin entry
                bin_entry = {
                    'chrom': chrom,
                    'start': bin_start,
                    'end': bin_end,
                    'pos': bin_start,  # Position for compatibility
                    'n_genes': len(overlapping_genes),
                    'genes': ';'.join(overlapping_genes['gene_id'].tolist()) if len(overlapping_genes) > 0 else ''
                }
                
                bins_list.append(bin_entry)
        
        bins_df = pd.DataFrame(bins_list)
        return bins_df
    
    @staticmethod
    def validate_coordinates(coords):
        """Validate 3D coordinates array"""
        if not isinstance(coords, np.ndarray):
            coords = np.array(coords)
        
        if coords.ndim != 2 or coords.shape[1] != 3:
            raise ValueError("Coordinates must be Nx3 array")
        
        # Check for NaN or infinite values
        if np.any(np.isnan(coords)):
            warnings.warn("NaN values found in coordinates")
        
        if np.any(np.isinf(coords)):
            warnings.warn("Infinite values found in coordinates")
        
        return coords
    
    @staticmethod
    def normalize_coordinates(coords, method='center'):
        """Normalize coordinates using different methods"""
        coords = Cell3DUtils.validate_coordinates(coords)
        
        if method == 'center':
            # Center coordinates around origin
            center = np.mean(coords, axis=0)
            normalized_coords = coords - center
        
        elif method == 'standardize':
            # Standardize to zero mean and unit variance
            mean = np.mean(coords, axis=0)
            std = np.std(coords, axis=0)
            normalized_coords = (coords - mean) / (std + 1e-10)
        
        elif method == 'minmax':
            # Scale to [0, 1] range
            min_vals = np.min(coords, axis=0)
            max_vals = np.max(coords, axis=0)
            normalized_coords = (coords - min_vals) / (max_vals - min_vals + 1e-10)
        
        elif method == 'unit_sphere':
            # Project onto unit sphere
            center = np.mean(coords, axis=0)
            centered_coords = coords - center
            distances = np.linalg.norm(centered_coords, axis=1, keepdims=True)
            normalized_coords = centered_coords / (distances + 1e-10)
        
        else:
            raise ValueError(f"Unsupported normalization method: {method}")
        
        return normalized_coords
    
    @staticmethod
    def calculate_pairwise_distances(coords1, coords2=None, metric='euclidean'):
        """Calculate pairwise distances between coordinate sets"""
        coords1 = Cell3DUtils.validate_coordinates(coords1)
        
        if coords2 is None:
            coords2 = coords1
        else:
            coords2 = Cell3DUtils.validate_coordinates(coords2)
        
        if metric == 'euclidean':
            # Efficient computation using broadcasting
            diff = coords1[:, np.newaxis, :] - coords2[np.newaxis, :, :]
            distances = np.sqrt(np.sum(diff**2, axis=2))
        
        elif metric == 'manhattan':
            diff = coords1[:, np.newaxis, :] - coords2[np.newaxis, :, :]
            distances = np.sum(np.abs(diff), axis=2)
        
        elif metric == 'chebyshev':
            diff = coords1[:, np.newaxis, :] - coords2[np.newaxis, :, :]
            distances = np.max(np.abs(diff), axis=2)
        
        else:
            raise ValueError(f"Unsupported distance metric: {metric}")
        
        return distances
    
    @staticmethod
    def filter_by_distance(coords, center_point, max_distance):
        """Filter coordinates by distance from center point"""
        coords = Cell3DUtils.validate_coordinates(coords)
        center_point = np.array(center_point)
        
        if center_point.shape != (3,):
            raise ValueError("Center point must be 3D coordinate")
        
        distances = np.linalg.norm(coords - center_point, axis=1)
        mask = distances <= max_distance
        
        return coords[mask], mask
    
    @staticmethod
    def create_grid_coordinates(bounds, resolution):
        """Create regular grid coordinates within bounds"""
        x_min, x_max, y_min, y_max, z_min, z_max = bounds
        
        x_coords = np.arange(x_min, x_max + resolution, resolution)
        y_coords = np.arange(y_min, y_max + resolution, resolution)
        z_coords = np.arange(z_min, z_max + resolution, resolution)
        
        # Create meshgrid
        X, Y, Z = np.meshgrid(x_coords, y_coords, z_coords, indexing='ij')
        
        # Flatten and combine
        grid_coords = np.column_stack([
            X.ravel(),
            Y.ravel(),
            Z.ravel()
        ])
        
        return grid_coords
    
    @staticmethod
    def interpolate_missing_coordinates(coords, positions, method='linear'):
        """Interpolate missing coordinates based on genomic positions"""
        coords = Cell3DUtils.validate_coordinates(coords)
        positions = np.array(positions)
        
        if len(coords) != len(positions):
            raise ValueError("Coordinates and positions must have same length")
        
        # Find missing coordinates (NaN values)
        missing_mask = np.any(np.isnan(coords), axis=1)
        
        if not np.any(missing_mask):
            return coords  # No missing coordinates
        
        valid_mask = ~missing_mask
        
        if np.sum(valid_mask) < 2:
            warnings.warn("Not enough valid coordinates for interpolation")
            return coords
        
        # Interpolate each dimension separately
        interpolated_coords = coords.copy()
        
        for dim in range(3):
            if method == 'linear':
                interpolated_coords[missing_mask, dim] = np.interp(
                    positions[missing_mask],
                    positions[valid_mask],
                    coords[valid_mask, dim]
                )
            else:
                raise ValueError(f"Unsupported interpolation method: {method}")
        
        return interpolated_coords
    
    @staticmethod
    def detect_outliers(coords, method='iqr', threshold=1.5):
        """Detect outlier coordinates"""
        coords = Cell3DUtils.validate_coordinates(coords)
        
        if method == 'iqr':
            # Interquartile range method
            outlier_mask = np.zeros(len(coords), dtype=bool)
            
            for dim in range(3):
                q1 = np.percentile(coords[:, dim], 25)
                q3 = np.percentile(coords[:, dim], 75)
                iqr = q3 - q1
                
                lower_bound = q1 - threshold * iqr
                upper_bound = q3 + threshold * iqr
                
                dim_outliers = (coords[:, dim] < lower_bound) | (coords[:, dim] > upper_bound)
                outlier_mask |= dim_outliers
        
        elif method == 'zscore':
            # Z-score method
            z_scores = np.abs((coords - np.mean(coords, axis=0)) / np.std(coords, axis=0))
            outlier_mask = np.any(z_scores > threshold, axis=1)
        
        elif method == 'distance':
            # Distance from center method
            center = np.mean(coords, axis=0)
            distances = np.linalg.norm(coords - center, axis=1)
            distance_threshold = np.mean(distances) + threshold * np.std(distances)
            outlier_mask = distances > distance_threshold
        
        else:
            raise ValueError(f"Unsupported outlier detection method: {method}")
        
        return outlier_mask
    
    @staticmethod
    def smooth_trajectory(coords, window_size=3, method='moving_average'):
        """Smooth coordinate trajectory"""
        coords = Cell3DUtils.validate_coordinates(coords)
        
        if window_size < 1 or window_size > len(coords):
            raise ValueError("Invalid window size")
        
        smoothed_coords = coords.copy()
        
        if method == 'moving_average':
            # Simple moving average
            for i in range(len(coords)):
                start_idx = max(0, i - window_size // 2)
                end_idx = min(len(coords), i + window_size // 2 + 1)
                smoothed_coords[i] = np.mean(coords[start_idx:end_idx], axis=0)
        
        elif method == 'gaussian':
            # Gaussian smoothing
            from scipy.ndimage import gaussian_filter1d
            sigma = window_size / 3.0  # Convert window size to sigma
            
            for dim in range(3):
                smoothed_coords[:, dim] = gaussian_filter1d(coords[:, dim], sigma=sigma)
        
        else:
            raise ValueError(f"Unsupported smoothing method: {method}")
        
        return smoothed_coords

    def _auto_genome_coord(genome_coord):
        """
        Automatically convert genome_coord to chrom,start,end format
        INPUT:
            genome_coord: str, list, or tuple
                Format like "chrom:start-end" or ["chrom", start, end] or "chrom"
        OUTPUT:
            chrom, start, end: tuple
        """
        import re
        # determine the genome_coord format
        if isinstance(genome_coord, str):
            if ":" in genome_coord:
                chrom, start, end = re.split(":|-", genome_coord)
                start, end = int(start), int(end)
            else:
                chrom, start, end = genome_coord, None, None
        elif isinstance(genome_coord, (list, tuple)):
            chrom, start, end = genome_coord
        else:
            raise ValueError('Genome_coord should be str or list/tuple. e.g. "chr1a:10000-20000" or ["chr1a",10000,20000] or "chr1a"')
        
        return chrom, start, end