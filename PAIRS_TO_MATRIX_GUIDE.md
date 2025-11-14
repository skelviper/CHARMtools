# Guide: Converting Pairs DataFrame to Numpy Matrix in CHARMtools

## Current State

After thorough exploration of the CHARMtools codebase, there is **NO direct function** that converts a pairs DataFrame to a numpy contact matrix.

## Available Functions

### 1. CHARMio Module (`utils/CHARMio.py`)

#### Reading/Writing Pairs
- **`parse_pairs(filename: str)`** - Reads pairs file into a pandas DataFrame
  - Returns: DataFrame with columns like `["readID", "chrom1", "pos1", "chrom2", "pos2", ...]`
  - Preserves comments as attributes
  
- **`write_pairs(pairs: pd.DataFrame, out_name: str)`** - Writes pairs DataFrame to file
  - Writes compressed pairs file in 4DN format

#### Reading Matrices from Cooler Files
- **`getMatrixFromCooler(filepath: str, genome_coord1: str, genome_coord2=None, resolution=40000, balance=False)`** 
  - Reads matrix from .cool or .mcool files (NOT from pairs DataFrame)
  - Returns: numpy.ndarray
  - Example: `matrix = getMatrixFromCooler("data.mcool", "chr1:1000000-2000000", resolution=40000)`

- **`getMatrixFromMCOOLs(filepath: str, genome_coord1: str, genome_coord2=None, resolution=40000, balance=False)`**
  - Specifically for .mcool files
  - Returns: numpy.ndarray

#### Reading Sparse Matrices from HDF5
- **`load_sparse_form_h5(filename)`** - Loads sparse matrix from HDF5 file
  - Returns: scipy.sparse.csr_matrix

- **`read_mat_h5(filepath, genome_coord, resolution=10000)`** - Reads subregion from HDF5
  - Returns: numpy.ndarray (log2 transformed)

### 2. Pairs Manipulation Module (`utils/pairs_manipulations.py`)

- **`sortPairs(pairs: pd.DataFrame)`** - Sorts pairs DataFrame
- **`parse_and_combine_pairs(file_paths, num_cores=None)`** - Combines multiple pairs files
- **`pairs_to_bedgraph(...)`** - Converts pairs to bedgraph format (1D signal)
- **`pileup_bedgraph_on_bed(...)`** - Pileup analysis (NOT a contact matrix)

## Typical Workflow

The standard bioinformatics workflow for Hi-C data analysis is:

```
1. Pairs File (.pairs.gz)
   ↓
2. Use cooler/pairtools (external tools)
   ↓
3. Cooler File (.cool or .mcool)
   ↓
4. CHARMtools: getMatrixFromCooler()
   ↓
5. Numpy Matrix
```

### Using External Tools

To convert pairs to a contact matrix, you would typically use:

```bash
# Using cooler (recommended)
cooler cload pairs \
  -c1 2 -p1 3 -c2 4 -p2 5 \
  genome.chrom.sizes:resolution \
  input.pairs.gz \
  output.cool
```

Then in Python:
```python
from utils.CHARMio import getMatrixFromCooler
matrix = getMatrixFromCooler("output.cool", "chr1:0-10000000", resolution=40000)
```

## What's Missing: Direct Pairs-to-Matrix Function

### Function Signature (if implemented)

```python
def pairs_to_matrix(
    pairs: pd.DataFrame,
    chrom: str,
    start: int,
    end: int, 
    resolution: int = 40000,
    balanced: bool = False
) -> np.ndarray:
    """
    Convert pairs DataFrame to contact matrix.
    
    Parameters:
    -----------
    pairs : pd.DataFrame
        DataFrame with columns ["chrom1", "pos1", "chrom2", "pos2"]
    chrom : str
        Chromosome name (e.g., "chr1")
    start : int
        Start position in base pairs
    end : int
        End position in base pairs
    resolution : int
        Bin size in base pairs
    balanced : bool
        Whether to apply ICE normalization
        
    Returns:
    --------
    np.ndarray
        Contact matrix as 2D numpy array
    """
    pass
```

### Implementation Location

If this function were to be implemented, it should go in:
- **Primary option**: `utils/CHARMio.py` (alongside other matrix functions)
- **Alternative**: `utils/pairs_manipulations.py` (alongside other pairs functions)

### Implementation Approach

```python
def pairs_to_matrix(pairs, chrom, start, end, resolution=40000):
    """Convert pairs DataFrame to contact matrix"""
    import numpy as np
    
    # Filter pairs for the region
    region_pairs = pairs[
        (pairs['chrom1'] == chrom) & 
        (pairs['chrom2'] == chrom) &
        (pairs['pos1'] >= start) & (pairs['pos1'] < end) &
        (pairs['pos2'] >= start) & (pairs['pos2'] < end)
    ].copy()
    
    # Bin the positions
    region_pairs['bin1'] = ((region_pairs['pos1'] - start) // resolution).astype(int)
    region_pairs['bin2'] = ((region_pairs['pos2'] - start) // resolution).astype(int)
    
    # Count contacts per bin pair
    n_bins = (end - start) // resolution
    matrix = np.zeros((n_bins, n_bins))
    
    for _, row in region_pairs.iterrows():
        b1, b2 = row['bin1'], row['bin2']
        if 0 <= b1 < n_bins and 0 <= b2 < n_bins:
            matrix[b1, b2] += 1
            if b1 != b2:  # Make symmetric for upper triangle pairs
                matrix[b2, b1] += 1
    
    return matrix
```

## Related Functions in CHARMtools

### Cell3D Object Methods
The `Cell3D` class has methods for working with 3D structures derived from Hi-C:
- `calc_distance_matrix(genome_coord)` - Calculate spatial distance matrix
- `calc_feature_matrix(feature_name, metric='correlation')` - Feature similarity matrix

### MultiCell3D Methods
- `calc_distance_matrix(genome_coord, cellnames=None)` - Multi-cell distance matrices
- `calc_feature_matrix(genome_coord, feature, cells=None)` - Feature matrices across cells

## Summary

**To answer the original question**: 

CHARMtools does NOT currently have a function that directly converts a pairs DataFrame to a numpy contact matrix. 

**Recommended approach**:
1. Use external tools (cooler) to convert pairs → cooler file
2. Use `CHARMio.getMatrixFromCooler()` to read the matrix

**Alternative**: Implement `pairs_to_matrix()` function in `utils/CHARMio.py` following the approach outlined above.
