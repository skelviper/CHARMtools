# Utilities

The `utils` module provides essential utility functions and helper tools that support the core functionality of CHARMtools. These utilities handle file I/O, data manipulation, and common operations used throughout the package.

## Overview

The utilities module includes:
- File I/O operations for various genomic data formats
- Data format conversions and manipulations
- Helper functions for common tasks
- Color generation and visualization support

## Core Modules

### 1. CHARMio.py - Input/Output Operations

**Purpose**: Comprehensive I/O functions for reading and writing genomic data formats used in 3D chromatin analysis.

#### Key Functions

##### File Parsing Functions

**`parse_pairs(filename: str) -> pd.DataFrame`**
- **Input**: Path to .pairs file (4DN standard format, gzipped or plain text)
- **Output**: pandas DataFrame with contact data
- **Features**:
  - Automatic comment line handling
  - Flexible column assignment
  - Memory-efficient reading
  - Preserves metadata in DataFrame attributes

```python
from CHARMtools.utils.CHARMio import parse_pairs

# Load pairs file
pairs_data = parse_pairs("contacts.pairs.gz")
print(f"Loaded {len(pairs_data)} contacts")
print(f"Sample name: {pairs_data.attrs['name']}")
```

**`parse_3dg(filename: str, s2m=False, sorting=False) -> pd.DataFrame`**
- **Input**: Path to .3dg file containing 3D coordinates
- **Output**: pandas DataFrame with 3D structure data
- **Parameters**:
  - `s2m`: Convert single-cell to multi-cell format
  - `sorting`: Sort coordinates by genomic position

```python
from CHARMtools.utils.CHARMio import parse_3dg

# Load 3D structure
structure = parse_3dg("cell_structure.3dg", sorting=True)
print(f"Structure contains {len(structure)} beads")
```

**`parse_gtf(filename: str) -> pd.DataFrame`**
- **Input**: Path to GTF annotation file
- **Output**: pandas DataFrame with gene annotations
- **Features**: Handles GENCODE/ENSEMBL GTF formats

##### File Writing Functions

**`write_pairs(pairs: pd.DataFrame, out_name: str)`**
- **Input**: pandas DataFrame with pairs data, output filename
- **Output**: Writes .pairs file in 4DN standard format
- **Features**: Preserves comments and metadata

**`write_3dg(_3dg: pd.DataFrame, outname: str, m2s=False)`**
- **Input**: pandas DataFrame with 3D coordinates, output filename
- **Output**: Writes .3dg file
- **Parameters**: `m2s` for multi-cell to single-cell conversion

##### Matrix Operations

**`getMatrixFromMCOOLs(filepath: str, genome_coord1: str, genome_coord2=None, resolution=40000, balance=False) -> np.ndarray`**
- **Input**: Path to .mcool file, genomic coordinates, resolution
- **Output**: Contact matrix as numpy array
- **Features**: Supports balanced and raw matrices

**`getMatrixFromCooler(filepath: str, genome_coord1: str, genome_coord2=None, resolution=40000, balance=False) -> np.ndarray`**
- **Input**: Path to .cool file, genomic coordinates, resolution
- **Output**: Contact matrix as numpy array

##### Utility Functions

**`check_index_binsize(df: pd.DataFrame) -> int`**
- **Input**: DataFrame with MultiIndex (chr, pos)
- **Output**: Detected resolution/bin size
- **Purpose**: Automatically determine data resolution

**`divide_name(filename: str) -> tuple`**
- **Input**: Filename with extensions
- **Output**: Tuple of (basename, extensions)
- **Purpose**: Robust filename parsing

### 2. helper.py - General Helper Functions

**Purpose**: Collection of utility functions for common data processing tasks.

#### Key Functions

- Data validation and type checking
- Coordinate system conversions
- Statistical helper functions
- Memory management utilities

### 3. pairs_manipulations.py - Pairs Data Processing

**Purpose**: Specialized functions for manipulating and processing contact pairs data.

#### Key Functions

- Pairs filtering and subsetting
- Distance calculations
- Contact frequency analysis
- Data quality metrics

**Example Usage**:
```python
from CHARMtools.utils.pairs_manipulations import filter_pairs_by_distance

# Filter pairs by genomic distance
filtered_pairs = filter_pairs_by_distance(pairs_data, min_dist=1000, max_dist=1000000)
```

### 4. generateColor2.py - Color Generation

**Purpose**: Generate color schemes and palettes for visualization.

#### Key Functions

**Color Palette Generation**
- Categorical color schemes
- Continuous color maps
- Custom color interpolation
- Accessibility-friendly palettes

**Example Usage**:
```python
from CHARMtools.utils.generateColor2 import generate_palette

# Generate colors for chromosomes
chrom_colors = generate_palette(n_colors=23, palette_type='categorical')
```

## Data Format Support

### Supported Input Formats

1. **.pairs** - 4DN standard contact pairs format
2. **.3dg** - 3D genome structure coordinates
3. **.gtf/.gff** - Gene annotation files
4. **.cool/.mcool** - Cooler contact matrices
5. **.h5** - HDF5 data files
6. **.parquet** - Parquet data files

### Supported Output Formats

1. **.pairs** - Contact pairs
2. **.3dg** - 3D coordinates
3. **.csv/.tsv** - Tabular data
4. **.h5** - HDF5 format
5. **.cif** - Crystallographic Information File

## Common Usage Patterns

### 1. Loading and Converting Data

```python
from CHARMtools.utils.CHARMio import parse_pairs, parse_3dg, write_pairs

# Load different data types
pairs = parse_pairs("contacts.pairs.gz")
structure = parse_3dg("structure.3dg")

# Process data
filtered_pairs = pairs[pairs['chr1'] == pairs['chr2']]  # cis contacts only

# Save processed data
write_pairs(filtered_pairs, "cis_contacts.pairs")
```

### 2. Matrix Operations

```python
from CHARMtools.utils.CHARMio import getMatrixFromCooler

# Load contact matrix
matrix = getMatrixFromCooler(
    "contacts.cool", 
    "chr1:1000000-2000000", 
    resolution=10000, 
    balance=True
)

print(f"Matrix shape: {matrix.shape}")
```

### 3. File Format Detection

```python
from CHARMtools.utils.CHARMio import divide_name, check_index_binsize

# Parse filename
basename, extension = divide_name("sample.pairs.gz")
print(f"Base: {basename}, Extension: {extension}")

# Detect resolution
resolution = check_index_binsize(data_with_coords)
print(f"Detected resolution: {resolution} bp")
```

## Performance Considerations

### Memory Management

- Use chunked reading for large files
- Leverage pandas' memory-efficient data types
- Consider using Parquet format for intermediate storage

### I/O Optimization

- Prefer compressed formats (gzip) for storage
- Use appropriate data types to reduce memory usage
- Cache frequently accessed data

### Best Practices

1. **File Handling**:
   - Always check file existence before processing
   - Use context managers for file operations
   - Handle compressed and uncompressed files uniformly

2. **Data Validation**:
   - Validate input formats before processing
   - Check for required columns in DataFrames
   - Handle missing or malformed data gracefully

3. **Error Handling**:
   - Provide informative error messages
   - Log processing steps for debugging
   - Implement fallback strategies for common issues

## Integration with CHARMtools

The utilities module is designed to work seamlessly with other CHARMtools components:

- **Cell3D Objects**: Utilities provide the I/O backend for Cell3D data loading
- **Preprocessing**: File parsing functions support preprocessing workflows
- **Analysis**: Matrix operations enable contact analysis
- **Visualization**: Color generation supports plotting functions

## Troubleshooting

### Common Issues

1. **File Format Errors**
   - Verify file format compatibility
   - Check for proper compression
   - Validate column structures

2. **Memory Issues**
   - Use chunked processing for large files
   - Consider data type optimization
   - Monitor memory usage during operations

3. **Performance Problems**
   - Profile I/O operations
   - Use appropriate file formats
   - Consider parallel processing where applicable

### Getting Help

For utility-specific issues:
- Check function docstrings for detailed parameter information
- Review example usage in other CHARMtools modules
- Consult the [Troubleshooting Guide](Troubleshooting.md)

## Next Steps

After familiarizing yourself with the utilities:
1. Explore [Data Preprocessing](Data-Preprocessing.md) workflows
2. Learn about [Cell3D Objects](Cell3D-Objects.md) for data analysis
3. Check out [Visualization](Visualization.md) tools for data exploration

---

*For more information about CHARMtools, visit the [Home](Home.md) page or check the [Installation Guide](Installation-Guide.md).*