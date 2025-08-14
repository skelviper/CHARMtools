from .multicell3d_core import MultiCell3DCore
from .multicell3d_analysis import MultiCell3DAnalysis
from .multicell3d_io import MultiCell3DIO, load_CHARM, load_cells
from .multicell3d_utils import dev_only
from .multicell3d_visualization import MultiCell3DVisualization

class MultiCell3D(MultiCell3DCore, MultiCell3DAnalysis, MultiCell3DIO, MultiCell3DVisualization):
    """
    MultiCell3D class for handling multiple Cell3D objects.
    
    This class combines functionality from multiple modules:
    - MultiCell3DCore: Core functionality, data management, and AnnData support
    - MultiCell3DAnalysis: Analysis methods for 3D genomics data
    - MultiCell3DIO: Input/output operations and data loading
    - MultiCell3DVisualization: Visualization methods
    """
    
    def __init__(self, cells):
        """
        Initialize MultiCell3D object with integrated AnnData support.
        
        Parameters:
        -----------
        cells : list
            List of Cell3D objects
        """
        MultiCell3DCore.__init__(self, cells)

# Export main functions for backward compatibility
__all__ = ['MultiCell3D', 'load_CHARM', 'load_cells', 'dev_only']