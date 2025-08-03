# CHARMtools obj module
# Modular Cell3D implementation

from .Cell3D import Cell3D
from .cell3d_core import Cell3DCore
from .cell3d_data import Cell3DData
from .cell3d_features import Cell3DFeatures
from .cell3d_spatial import Cell3DSpatial
from .cell3d_output import Cell3DOutput
from .cell3d_visualization import Cell3DVisualization
from .cell3d_analysis import Cell3DAnalysis
from .cell3d_utils import Cell3DUtils
from .cell3d_dev import Cell3DDev
from .MultiCell3D import MultiCell3D

# Backward compatibility
from .Cell3D import create_bins_genetable

__all__ = [
    'Cell3D',
    'Cell3DCore',
    'Cell3DData', 
    'Cell3DFeatures',
    'Cell3DSpatial',
    'Cell3DOutput',
    'Cell3DVisualization',
    'Cell3DAnalysis',
    'Cell3DUtils',
    'Cell3DDev',
    'MultiCell3D',
    'create_bins_genetable'
]

__version__ = "2.0.0"
__author__ = "CHARM Tools Team"
__description__ = "Modular Cell3D implementation for 3D chromatin structure analysis"