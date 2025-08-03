# Installation Guide

This guide will walk you through the installation process for CHARMtools on different operating systems.

## 📋 System Requirements

### Minimum Requirements
- **Python**: 3.7 or higher (3.8+ recommended)
- **Operating System**: Linux, macOS, or Windows
- **Memory**: 8GB RAM minimum
- **Storage**: 1GB for installation, additional space for data
- **Internet**: Required for downloading dependencies

### Recommended Requirements
- **Python**: 3.9 or 3.10
- **Memory**: 16GB+ RAM for large datasets
- **Storage**: 10GB+ free space
- **CPU**: Multi-core processor for parallel processing

### Optional Requirements
- **R**: 4.0+ (required for some analysis modules)
- **Git**: For development and version control
- **Jupyter**: For interactive analysis

## 🐍 Python Environment Setup

### Option 1: Using Conda/Mamba (Recommended)

Conda provides the most reliable installation experience with automatic dependency management.

#### Install Conda/Mamba

If you don't have conda installed:

```bash
# Install Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# Install Mamba (faster than conda)
conda install -c conda-forge mamba
```

#### Create CHARMtools Environment

```bash
# Create a new environment with Python 3.10
mamba create -n charmtools -c conda-forge -c bioconda python=3.10

# Activate the environment
mamba activate charmtools
```

### Option 2: Using Virtual Environment

If you prefer using pip and virtual environments:

```bash
# Create virtual environment
python -m venv charmtools_env

# Activate environment
# On Linux/macOS:
source charmtools_env/bin/activate
# On Windows:
charmtools_env\Scripts\activate

# Upgrade pip
pip install --upgrade pip
```

## 📦 Installing CHARMtools

### Method 1: From Source (Recommended)

```bash
# Clone the repository
git clone https://github.com/skelviper/CHARMtools.git
cd CHARMtools

# Install dependencies using conda/mamba
mamba env update --file charmtools.yaml

# Or install using pip
pip install -r requirements.txt

# Install CHARMtools in development mode
pip install -e .
```

### Method 2: Direct Installation

If you don't need the source code:

```bash
# Install directly from GitHub
pip install git+https://github.com/skelviper/CHARMtools.git
```

## 🔧 Core Dependencies

CHARMtools requires the following core Python packages:

### Essential Dependencies
```
numpy>=1.19.0
pandas>=1.3.0
scipy>=1.7.0
scikit-learn>=1.0.0
matplotlib>=3.3.0
seaborn>=0.11.0
plotly>=5.0.0
h5py>=3.0.0
tables>=3.6.0
```

### Optional Dependencies
```
jupyter>=1.0.0
ipywidgets>=7.6.0
numba>=0.56.0  # For performance optimization
dask>=2021.0.0  # For large dataset processing
```

### R Dependencies (Optional)

For R-based analysis modules:

```r
# Install required R packages
install.packages(c("BiocManager", "devtools"))
BiocManager::install(c("GenomicRanges", "rtracklayer", "HiCcompare"))
```

## ✅ Verification

### Test Basic Installation

```python
# Test import
import charmtools
from charmtools.obj import Cell3D

# Check version
print(f"CHARMtools version: {charmtools.__version__}")

# Test basic functionality
cell = Cell3D(cellname="test", resolution=40000)
print("Installation successful!")
```

### Run Test Suite

```bash
# Run basic tests
python -m pytest tests/ -v

# Run specific module tests
python -c "from charmtools.obj import Cell3D; print('Cell3D import successful')"
python -c "from charmtools.utils import CHARMio; print('Utils import successful')"
```

### Check Command Line Interface

```bash
# Test CLI
charmtools --help
charmtools --version
```

## 🐛 Troubleshooting

### Common Issues

#### 1. Import Errors

**Problem**: `ModuleNotFoundError: No module named 'charmtools'`

**Solution**:
```bash
# Make sure you're in the correct environment
conda activate charmtools

# Reinstall in development mode
cd CHARMtools
pip install -e .
```

#### 2. Dependency Conflicts

**Problem**: Package version conflicts

**Solution**:
```bash
# Create fresh environment
mamba create -n charmtools_clean python=3.10
mamba activate charmtools_clean

# Install step by step
mamba install -c conda-forge numpy pandas scipy matplotlib
mamba install -c conda-forge scikit-learn seaborn plotly
pip install -e .
```

#### 3. Memory Issues

**Problem**: Out of memory errors

**Solution**:
- Use `on_disk=True` when creating Cell3D objects
- Process data in chunks
- Increase system swap space
- Use a machine with more RAM

#### 4. R Integration Issues

**Problem**: R modules not working

**Solution**:
```bash
# Install R and required packages
sudo apt-get install r-base r-base-dev  # Ubuntu/Debian
brew install r  # macOS

# Install rpy2 for Python-R interface
mamba install -c conda-forge rpy2
```

### Platform-Specific Issues

#### Linux
```bash
# Install system dependencies
sudo apt-get update
sudo apt-get install build-essential python3-dev
```

#### macOS
```bash
# Install Xcode command line tools
xcode-select --install

# Install Homebrew if needed
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

#### Windows
- Install Microsoft Visual C++ Build Tools
- Use Windows Subsystem for Linux (WSL) for better compatibility

## 🔄 Updating CHARMtools

### Update from Source

```bash
cd CHARMtools
git pull origin main
pip install -e . --upgrade
```

### Update Dependencies

```bash
# Update conda environment
mamba env update --file charmtools.yaml --prune

# Or update pip packages
pip install --upgrade -r requirements.txt
```

## 🗑️ Uninstallation

### Remove CHARMtools

```bash
# If installed with pip
pip uninstall charmtools

# Remove conda environment
mamba env remove -n charmtools

# Remove source directory
rm -rf CHARMtools/
```

## 🚀 Next Steps

Once installation is complete:

1. **[Quick Start](Quick-Start)** - Learn basic usage
2. **[Tutorial](Tutorial)** - Follow step-by-step examples
3. **[Cell3D Objects](Cell3D-Objects)** - Understand core data structures
4. **[Examples](Examples)** - See real-world applications

## 💡 Performance Tips

### For Large Datasets
- Use `on_disk=True` for Cell3D objects
- Install `numba` for faster computations
- Consider using `dask` for parallel processing
- Allocate sufficient swap space

### For Development
- Install in development mode: `pip install -e .`
- Use `pre-commit` hooks for code quality
- Set up IDE integration (VS Code, PyCharm)

## 📞 Getting Help

If you encounter issues during installation:

1. Check this troubleshooting section
2. Search [GitHub Issues](https://github.com/skelviper/CHARMtools/issues)
3. Create a new issue with:
   - Your operating system and version
   - Python version
   - Complete error message
   - Installation method used

---

**Installation complete?** Continue to the [Quick Start](Quick-Start) guide!