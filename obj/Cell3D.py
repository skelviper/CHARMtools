# Cell3D - Refactored modular implementation
# Main Cell3D class that inherits from all modules

import warnings
import pandas as pd
import numpy as np
from pathlib import Path

# Import all modules
from .cell3d_core import Cell3DCore
from .cell3d_data import Cell3DData
from .cell3d_features import Cell3DFeatures
from .cell3d_spatial import Cell3DSpatial
from .cell3d_output import Cell3DOutput
from .cell3d_visualization import Cell3DVisualization
from .cell3d_analysis import Cell3DAnalysis
from .cell3d_utils import Cell3DUtils
from .cell3d_dev import Cell3DDev

class Cell3D(Cell3DCore, Cell3DData, Cell3DFeatures, Cell3DSpatial, 
             Cell3DOutput, Cell3DVisualization, Cell3DAnalysis, 
             Cell3DUtils, Cell3DDev):
    """
    Cell3D: A comprehensive class for 3D chromatin structure analysis
    
    This class provides functionality for:
    - Loading and managing 3D genomic data
    - Adding various genomic features (ChIP-seq, ATAC-seq, RNA-seq)
    - Spatial analysis and clustering
    - Visualization and plotting
    - Structure analysis and comparison
    - Output format conversion (CIF, PDB, XYZ)
    
    The class is modularly designed with functionality split across:
    - Core: Basic data management and I/O
    - Data: Data processing and retrieval
    - Features: Feature addition and management
    - Spatial: Spatial analysis and clustering
    - Output: Format conversion and export
    - Visualization: 3D plotting and visualization
    - Analysis: Advanced analysis methods
    - Utils: Utility functions
    - Dev: Development and experimental features
    """
    
    def __init__(self, tdg_path=None, cellname=None, resolution=40000, 
                 on_disk=False, **kwargs):
        """
        Initialize Cell3D object
        
        Parameters:
        -----------
        tdg_path : str, optional
            Path to 3DG file
        cellname : str, optional
            Name of the cell
        resolution : int, default 40000
            Genomic resolution in base pairs
        on_disk : bool, default False
            Whether to keep data on disk to save memory
        **kwargs : dict
            Additional parameters
        """
        # Initialize core attributes
        self.cellname = cellname
        self.resolution = resolution
        self.on_disk = on_disk
        self.features = []
        self.metadata = {}
        self.tdg = None
        self.chrom_length = None
        self.expected = None
        self.hic_matrix = None
        
        # Load data if path provided
        if tdg_path is not None:
            self._load_tdg(tdg_path)
            
            # Set cellname from filename if not provided
            if self.cellname is None:
                self.cellname = Path(tdg_path).stem
    
    def __repr__(self):
        """String representation of Cell3D object"""
        if self.tdg is not None:
            n_points = len(self.tdg)
            n_chroms = self.tdg['chrom'].nunique() if 'chrom' in self.tdg.columns else 0
            n_features = len(self.features)
            
            return (f"Cell3D(cellname='{self.cellname}', "
                   f"resolution={self.resolution}, "
                   f"n_points={n_points}, "
                   f"n_chromosomes={n_chroms}, "
                   f"n_features={n_features}, "
                   f"on_disk={self.on_disk})")
        else:
            return f"Cell3D(cellname='{self.cellname}', empty=True)"
    
    def __len__(self):
        """Return number of genomic bins"""
        if self.tdg is not None:
            return len(self.tdg)
        return 0
    
    def copy(self):
        """Create a deep copy of the Cell3D object"""
        new_cell3d = Cell3D(
            cellname=self.cellname,
            resolution=self.resolution,
            on_disk=False  # Always create copy in memory
        )
        
        if self.tdg is not None:
            if self.on_disk:
                self.to_memory()
            new_cell3d.tdg = self.tdg.copy()
        
        new_cell3d.features = self.features.copy()
        new_cell3d.metadata = self.metadata.copy()
        
        if self.chrom_length is not None:
            new_cell3d.chrom_length = self.chrom_length.copy()
        
        if self.expected is not None:
            new_cell3d.expected = self.expected.copy()
        
        if self.hic_matrix is not None:
            new_cell3d.hic_matrix = self.hic_matrix.copy()
        
        return new_cell3d
    
    @property
    def chromosomes(self):
        """Get list of chromosomes in the data"""
        if self.tdg is not None:
            if self.on_disk:
                self.to_memory()
            return sorted(self.tdg['chrom'].unique())
        return []
    
    @property
    def coordinate_ranges(self):
        """Get coordinate ranges for x, y, z dimensions"""
        if self.tdg is not None:
            if self.on_disk:
                self.to_memory()
            
            return {
                'x': (self.tdg['x'].min(), self.tdg['x'].max()),
                'y': (self.tdg['y'].min(), self.tdg['y'].max()),
                'z': (self.tdg['z'].min(), self.tdg['z'].max())
            }
        return None
    
    @property
    def center_of_mass(self):
        """Get center of mass coordinates"""
        if self.tdg is not None:
            if self.on_disk:
                self.to_memory()
            
            return {
                'x': self.tdg['x'].mean(),
                'y': self.tdg['y'].mean(),
                'z': self.tdg['z'].mean()
            }
        return None
    
    def summary(self):
        """Get summary information about the Cell3D object"""
        summary_info = {
            'cellname': self.cellname,
            'resolution': self.resolution,
            'on_disk': self.on_disk,
            'n_features': len(self.features),
            'features': self.features,
            'metadata_keys': list(self.metadata.keys())
        }
        
        if self.tdg is not None:
            if self.on_disk:
                # Get basic info without loading full data
                summary_info.update({
                    'n_points': 'on_disk',
                    'chromosomes': 'on_disk',
                    'coordinate_ranges': 'on_disk'
                })
            else:
                summary_info.update({
                    'n_points': len(self.tdg),
                    'chromosomes': self.chromosomes,
                    'coordinate_ranges': self.coordinate_ranges,
                    'center_of_mass': self.center_of_mass
                })
        else:
            summary_info.update({
                'n_points': 0,
                'chromosomes': [],
                'coordinate_ranges': None
            })
        
        return summary_info
    
    def validate(self):
        """Validate the Cell3D object for consistency"""
        issues = []
        
        # Check basic attributes
        if self.cellname is None:
            issues.append("cellname is None")
        
        if self.resolution <= 0:
            issues.append("resolution must be positive")
        
        # Check data
        if self.tdg is None:
            issues.append("No data loaded (tdg is None)")
        else:
            if self.on_disk:
                self.to_memory()
            
            # Check required columns
            required_cols = ['chrom', 'pos', 'x', 'y', 'z']
            missing_cols = [col for col in required_cols if col not in self.tdg.columns]
            if missing_cols:
                issues.append(f"Missing required columns: {missing_cols}")
            
            # Check for NaN coordinates
            if self.tdg[['x', 'y', 'z']].isnull().any().any():
                issues.append("NaN values found in coordinates")
            
            # Check feature consistency
            missing_features = [f for f in self.features if f not in self.tdg.columns]
            if missing_features:
                issues.append(f"Features listed but not found in data: {missing_features}")
        
        if issues:
            warnings.warn(f"Validation issues found: {'; '.join(issues)}")
            return False
        
        return True


# Backward compatibility: keep original function names as module-level functions
def create_bins_genetable(genes_df, resolution=40000, genome_size_file=None):
    """Create bins gene table from genes dataframe (backward compatibility)"""
    return Cell3DUtils.create_bins_genetable(genes_df, resolution, genome_size_file)


# Module metadata
__version__ = "2.0.0"
__author__ = "CHARM Tools Team"
__description__ = "Modular Cell3D implementation for 3D chromatin structure analysis"
