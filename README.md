# CHARMtools

**Single-cell Hi-C and Multi-omics Data Analysis and Visualization Toolkit**

[![Python](https://img.shields.io/badge/Python-3.7+-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

CHARMtools is a comprehensive toolkit designed for single-cell Hi-C and multi-omics data analysis, providing end-to-end solutions from data preprocessing to advanced analysis.

## 📋 Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Cell3D Object Usage (Example Data)](#-cell3d-object-usage-example-data)
- [Quick Start](#quick-start)
- [Usage Guide](#usage-guide)
- [Module Descriptions](#module-descriptions)
- [Contributing](#contributing)
- [License](#license)

## ✨ Features

- **Data Preprocessing**: Clean promiscuous legs, exon splicing, isolated contact points, and more
- **3D Structure Analysis**: 3D genome structure modeling and cleaning
- **Statistical Analysis**: Regression analysis, TSS enrichment, spatial statistics, and more
- **Visualization**: Rich plotting and visualization tools
- **Single-cell Analysis**: Specialized tools for single-cell Hi-C analysis
- **Command Line Interface**: Easy-to-use CLI tools

## 🚀 Installation

### Installation Steps

```bash
# Clone the repository
git clone https://github.com/skelviper/CHARMtools.git
cd CHARMtools
```

If you only use the preprocess workflow (`charm_preprocess_pipeline`), use the built-in environment `envs/charm.yml`; if you plan to run Cell3D and other downstream analyses, use `envs/charmtools.yml`.

```bash
# preprocess only
mamba env create -f envs/charm.yml
mamba activate charm

# Cell3D and downstream analysis
mamba env create -f envs/charmtools.yml
mamba activate analysis
```

## 🧪 Cell3D Object Usage (Example Data)

`example_with_data/example.ipynb` is a quick test notebook for part of the Cell3D workflow, with demo files in `example_with_data/data/`.

## 📁 Module Descriptions

### Core Modules

#### `charm_preprocess/` - Data Preprocessing
- `clean_leg.py` - Clean promiscuous legs that contact with multiple legs
- `clean_splicing.py` - Clean exon splicing from mRNA in contact files
- `clean_isolated.py` - Clean isolated contact points
- `sep_clean.py` - Separation cleaning tools
- `clean3.py` - 3D structure cleaning
- `tdg2pairs.py` - Convert 3DG files to pairs format

#### `obj/` - Core Objects
- `Cell3D.py` - 3D cell object, main 3D genome data structure
- `MultiCell3D.py` - Multi-cell 3D object

#### `analysis/` - Analysis and Visualization
- `regression.py` - Regression analysis models
- `tss_enrichment.py` - TSS enrichment analysis
- `spatialstat.py` - Spatial statistical analysis
- `simpleDiff.py` - Simple differential analysis
- `scGAD.py` - Single-cell GAD analysis
- `scAB.py` - Single-cell AB analysis
- `saddle.py` - Saddle point analysis
- `plot.py` - Visualization and plotting tools
- `peak_enrichment.py` - Peak enrichment analysis
- `loop.py` - Loop analysis
- `imputation.py` - Data imputation
- `coverage.py` - Coverage analysis
- `compartment.py` - Compartment analysis
- `cellcycle.py` - Cell cycle analysis
- `calc_escaping_score.py` - Escaping score calculation
- `TAD.py` - TAD analysis

#### `utils/` - Utility Tools
- `CHARMio.py` - Input/output tools, core IO functionality
- `pairs_manipulations.py` - Pairs file manipulation tools
- `helper.py` - Helper functions
- `generateColor2.py` - Color generation tools

### Other Directories
- `ref/` - Reference data
- `Rlibs/` - R language libraries
- `ref/mdso/` - MDSO (Multidimensional Scaling and Ordering) modules
- `archieve/` - Archive files

## 🤝 Contributing

We welcome all forms of contributions!

1. Fork this repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📞 Contact

- Project Homepage: [GitHub Repository](https://github.com/skelviper/CHARMtools)
- Issue Tracker: [Issues](https://github.com/skelviper/CHARMtools/issues)
- Email: skelviper@hotmail.com

---

**Note**: This is an actively developed project, and APIs may change. Please check the latest documentation for up-to-date information.
