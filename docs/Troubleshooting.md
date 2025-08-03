# Troubleshooting Guide

This guide provides solutions to common issues encountered when using CHARMtools. If you don't find a solution here, please check our [GitHub Issues](https://github.com/your-repo/CHARMtools/issues) or contact the development team.

## Installation Issues

### Python Environment Problems

#### Issue: Package conflicts during installation
```bash
ERROR: pip's dependency resolver does not currently consider all the packages that are installed
```

**Solution**:
1. Create a fresh conda environment:
```bash
conda create -n charmtools python=3.9
conda activate charmtools
```

2. Install dependencies in order:
```bash
conda install numpy pandas scipy scikit-learn matplotlib
pip install CHARMtools
```

#### Issue: ImportError for compiled extensions
```python
ImportError: cannot import name '_validate_lengths' from 'sklearn.utils.validation'
```

**Solution**:
1. Update scikit-learn:
```bash
pip install --upgrade scikit-learn
```

2. If issues persist, reinstall with conda:
```bash
conda install scikit-learn
```

### Dependency Issues

#### Issue: Missing R dependencies
```
Error: R package 'xyz' not found
```

**Solution**:
1. Install R packages manually:
```r
install.packages(c("ggplot2", "dplyr", "reshape2"))
```

2. For Bioconductor packages:
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GenomicRanges")
```

#### Issue: CUDA/GPU related errors
```
RuntimeError: CUDA out of memory
```

**Solution**:
1. Reduce batch size or data size
2. Use CPU-only mode:
```python
import os
os.environ['CUDA_VISIBLE_DEVICES'] = ''
```

## Data Loading Issues

### File Format Problems

#### Issue: Corrupted or incomplete .3dg files
```
ValueError: could not convert string to float
```

**Solution**:
1. Validate file format:
```python
from CHARMtools.utils.CHARMio import parse_3dg

try:
    data = parse_3dg("your_file.3dg")
    print("File loaded successfully")
except Exception as e:
    print(f"Error: {e}")
    # Check file manually
```

2. Check file structure:
```bash
head -n 10 your_file.3dg
# Should show: chr pos x y z format
```

#### Issue: .pairs file parsing errors
```
ParserError: Error tokenizing data
```

**Solution**:
1. Check file compression:
```bash
file your_file.pairs.gz
# Should show: gzip compressed data
```

2. Validate pairs format:
```python
from CHARMtools.utils.CHARMio import parse_pairs

# Try with error handling
try:
    pairs = parse_pairs("your_file.pairs.gz")
except Exception as e:
    print(f"Parsing error: {e}")
    # Check first few lines manually
```

### Memory Issues

#### Issue: Out of memory during data loading
```
MemoryError: Unable to allocate array
```

**Solution**:
1. Use chunked loading:
```python
from CHARMtools.obj.Cell3D import Cell3D

# Load with memory optimization
cell = Cell3D(
    tdg_path="large_file.3dg",
    on_disk=True,  # Keep data on disk
    chunk_size=100000
)
```

2. Subsample large datasets:
```python
# Load subset of data
cell = Cell3D(tdg_path="large_file.3dg")
subset_cell = cell.subset(chromosomes=['chr1', 'chr2'])
```

## Analysis Issues

### Performance Problems

#### Issue: Slow analysis on large datasets

**Solution**:
1. Use parallel processing:
```python
from CHARMtools.obj.Cell3D import Cell3D

cell = Cell3D(tdg_path="data.3dg")
# Enable multiprocessing
results = cell.spatial_analysis(n_jobs=4)
```

2. Optimize data structures:
```python
# Convert to memory-efficient format
cell.to_memory()  # If data fits in memory
# or
cell.to_disk()    # For large datasets
```

#### Issue: Analysis functions hanging or crashing

**Solution**:
1. Check data validity:
```python
# Validate Cell3D object
validation_results = cell.validate()
if not validation_results['valid']:
    print("Data validation failed:")
    print(validation_results['errors'])
```

2. Use debug mode:
```python
import logging
logging.basicConfig(level=logging.DEBUG)

# Run analysis with detailed logging
results = cell.spatial_analysis(verbose=True)
```

### Numerical Issues

#### Issue: NaN or infinite values in results
```
Warning: invalid value encountered in divide
```

**Solution**:
1. Check for missing data:
```python
# Check for NaN values
data = cell.get_data()
nan_count = data.isnull().sum().sum()
print(f"NaN values found: {nan_count}")

# Remove or impute missing values
clean_data = data.dropna()
```

2. Handle edge cases:
```python
# Use robust analysis parameters
results = cell.spatial_analysis(
    handle_missing='interpolate',
    robust=True
)
```

## Visualization Issues

### Display Problems

#### Issue: Plots not displaying in Jupyter notebooks
```
<Figure size 640x480 with 1 Axes>
# No plot shown
```

**Solution**:
1. Enable inline plotting:
```python
%matplotlib inline
import matplotlib.pyplot as plt
plt.ion()  # Turn on interactive mode
```

2. For 3D plots:
```python
%matplotlib widget
# or
%matplotlib notebook
```

#### Issue: Poor quality or pixelated plots

**Solution**:
1. Increase DPI:
```python
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
```

2. Use vector formats:
```python
cell.plot_3d_structure(save_path='plot.svg', format='svg')
```

### Interactive Visualization Issues

#### Issue: Interactive plots not responding

**Solution**:
1. Check backend:
```python
import matplotlib
print(matplotlib.get_backend())

# Switch to interactive backend
matplotlib.use('Qt5Agg')  # or 'TkAgg'
```

2. Update packages:
```bash
pip install --upgrade matplotlib plotly ipywidgets
```

## API and Integration Issues

### FastAPI Server Problems

#### Issue: Server won't start
```
ERROR: Could not import 'api:app'
```

**Solution**:
1. Check dependencies:
```bash
pip install fastapi uvicorn python-multipart
```

2. Verify file structure:
```bash
ls -la api.py
# Ensure api.py exists and is readable
```

3. Start with explicit path:
```bash
uvicorn CHARMtools.api:app --host 0.0.0.0 --port 8000
```

#### Issue: API endpoints returning errors
```
500 Internal Server Error
```

**Solution**:
1. Check server logs:
```bash
uvicorn CHARMtools.api:app --log-level debug
```

2. Test with simple requests:
```python
import requests

# Test health endpoint
response = requests.get("http://localhost:8000/health")
print(response.status_code, response.json())
```

## Platform-Specific Issues

### Windows Issues

#### Issue: Path separator problems
```
FileNotFoundError: [Errno 2] No such file or directory
```

**Solution**:
1. Use pathlib for cross-platform paths:
```python
from pathlib import Path

data_path = Path("data") / "structure.3dg"
cell = Cell3D(tdg_path=str(data_path))
```

2. Use forward slashes or raw strings:
```python
# Use forward slashes
path = "C:/data/structure.3dg"
# or raw strings
path = r"C:\data\structure.3dg"
```

### macOS Issues

#### Issue: Permission denied errors
```
PermissionError: [Errno 13] Permission denied
```

**Solution**:
1. Fix file permissions:
```bash
chmod +r your_data_file.3dg
```

2. Check directory permissions:
```bash
ls -la /path/to/data/
```

### Linux Issues

#### Issue: Missing system libraries
```
OSError: libGL.so.1: cannot open shared object file
```

**Solution**:
1. Install system dependencies:
```bash
# Ubuntu/Debian
sudo apt-get install libgl1-mesa-glx libglib2.0-0

# CentOS/RHEL
sudo yum install mesa-libGL glib2
```

## Performance Optimization

### Memory Optimization

#### Issue: High memory usage

**Solution**:
1. Use data types efficiently:
```python
# Optimize data types
data = cell.get_data()
data = data.astype({
    'chr': 'category',
    'pos': 'int32',  # Instead of int64
    'x': 'float32',  # Instead of float64
    'y': 'float32',
    'z': 'float32'
})
```

2. Process data in chunks:
```python
# Process chromosomes separately
for chrom in cell.chromosomes:
    subset = cell.subset(chromosomes=[chrom])
    results = subset.spatial_analysis()
    # Process results
    del subset  # Free memory
```

### Speed Optimization

#### Issue: Slow processing

**Solution**:
1. Use vectorized operations:
```python
import numpy as np

# Use numpy operations instead of loops
data = cell.get_data()
distances = np.sqrt(
    (data['x'].values[:, None] - data['x'].values) ** 2 +
    (data['y'].values[:, None] - data['y'].values) ** 2 +
    (data['z'].values[:, None] - data['z'].values) ** 2
)
```

2. Enable parallel processing:
```python
# Use all available cores
results = cell.spatial_analysis(n_jobs=-1)
```

## Data Quality Issues

### Validation Problems

#### Issue: Data validation failures

**Solution**:
1. Check data completeness:
```python
# Validate data structure
validation = cell.validate()
if not validation['valid']:
    for error in validation['errors']:
        print(f"Validation error: {error}")
```

2. Clean problematic data:
```python
# Remove invalid coordinates
data = cell.get_data()
clean_data = data.dropna(subset=['x', 'y', 'z'])
clean_data = clean_data[clean_data['x'].between(-1000, 1000)]
```

### Missing Data

#### Issue: Incomplete genomic coverage

**Solution**:
1. Check coverage statistics:
```python
coverage_stats = cell.get_coverage_stats()
print(f"Coverage: {coverage_stats}")
```

2. Impute missing regions:
```python
from CHARMtools.analysis.imputation import impute_missing

imputed_cell = impute_missing(
    cell,
    method='interpolation',
    min_coverage=0.5
)
```

## Getting Additional Help

### Documentation and Resources

1. **Check Documentation**: Review relevant sections:
   - [Installation Guide](Installation-Guide.md)
   - [Quick Start](Quick-Start.md)
   - [Cell3D Objects](Cell3D-Objects.md)
   - [API Documentation](API_DOCUMENTATION.md)

2. **Example Notebooks**: Look for similar use cases in example notebooks

3. **Community Support**:
   - GitHub Issues: Report bugs and request features
   - Discussion Forums: Ask questions and share experiences
   - Stack Overflow: Tag questions with `charmtools`

### Reporting Issues

When reporting issues, please include:

1. **System Information**:
```python
import sys
import platform
import CHARMtools

print(f"Python version: {sys.version}")
print(f"Platform: {platform.platform()}")
print(f"CHARMtools version: {CHARMtools.__version__}")
```

2. **Error Details**:
   - Full error traceback
   - Minimal reproducible example
   - Input data characteristics (size, format)
   - Expected vs. actual behavior

3. **Environment Details**:
   - Package versions (`pip list` or `conda list`)
   - Hardware specifications
   - Operating system version

### Contact Information

For direct support:
- **Email**: [support@charmtools.org](mailto:support@charmtools.org)
- **GitHub Issues**: [https://github.com/your-repo/CHARMtools/issues](https://github.com/your-repo/CHARMtools/issues)
- **Documentation**: [https://charmtools.readthedocs.io](https://charmtools.readthedocs.io)

---

*This troubleshooting guide is regularly updated. For the latest version, visit the [CHARMtools documentation](Home.md).*