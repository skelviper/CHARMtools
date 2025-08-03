# Quick Start Guide

Get up and running with CHARMtools in just a few minutes! This guide covers the essential steps to start analyzing 3D chromatin structure data.

## 🎯 Prerequisites

- CHARMtools installed ([Installation Guide](Installation-Guide))
- Python environment activated
- Basic familiarity with Python and pandas

## 📊 Your First Analysis

### Step 1: Import CHARMtools

```python
import charmtools
from charmtools.obj import Cell3D
from charmtools.utils import CHARMio
import pandas as pd
import numpy as np
```

### Step 2: Create a Cell3D Object

```python
# Create a new Cell3D object
cell = Cell3D(
    cellname="GM12878",
    resolution=40000,  # 40kb resolution
    on_disk=False      # Keep data in memory for small datasets
)

print(cell)
```

### Step 3: Load 3D Structure Data

#### Option A: Load from 3DG file

```python
# Load 3D coordinates from a .3dg file
cell = Cell3D(
    cellname="GM12878",
    resolution=40000,
    tdg_path="path/to/your/data.3dg"
)
```

#### Option B: Load from DataFrame

```python
# Create sample data
data = pd.DataFrame({
    'chrom': ['chr1'] * 100 + ['chr2'] * 100,
    'pos': list(range(0, 4000000, 40000)) * 2,
    'x': np.random.randn(200),
    'y': np.random.randn(200),
    'z': np.random.randn(200)
})

# Load from DataFrame
cell = Cell3D(
    cellname="sample",
    resolution=40000,
    tdg_path=data
)
```

### Step 4: Explore Your Data

```python
# Get basic information
print(f"Cell name: {cell.cellname}")
print(f"Resolution: {cell.resolution}")
print(f"Number of points: {len(cell)}")
print(f"Chromosomes: {cell.chromosomes()}")

# Get coordinate ranges
ranges = cell.coordinate_ranges()
print(f"Coordinate ranges: {ranges}")

# Get center of mass
com = cell.center_of_mass()
print(f"Center of mass: {com}")
```

### Step 5: Basic Data Operations

```python
# Get data for a specific chromosome
chr1_data = cell.get_data("chr1")
print(f"Chr1 has {len(chr1_data)} points")

# Get data for a specific region
region_data = cell.get_data("chr1:1000000-5000000")
print(f"Region has {len(region_data)} points")

# Create a subset
chr1_cell = cell.subset("chr1")
print(f"Chr1 subset: {chr1_cell.cellname}")
```

## 🧬 Adding Genomic Features

### Load Feature Data

```python
# Add ChIP-seq peaks (BED format)
cell.add_bed_feature(
    bed_path="path/to/H3K4me3_peaks.bed",
    feature_name="H3K4me3",
    binary=True  # Convert to binary (peak/no peak)
)

# Add expression data (bedGraph format)
cell.add_bedgraph_feature(
    bedgraph_path="path/to/expression.bedgraph",
    feature_name="RNA_expression"
)

# Add custom feature from DataFrame
feature_data = pd.DataFrame({
    'chrom': ['chr1', 'chr1', 'chr2'],
    'start': [1000000, 2000000, 1000000],
    'end': [1040000, 2040000, 1040000],
    'value': [1.5, 2.3, 0.8]
})

cell.add_feature_from_dataframe(
    df=feature_data,
    feature_name="custom_feature"
)

print(f"Features: {cell.features}")
```

## 📍 Spatial Analysis

### Distance Analysis

```python
# Calculate pairwise distances
distances = cell.calculate_distances(metric="euclidean")
print(f"Distance matrix shape: {distances.shape}")

# Calculate distances to nuclear center
radial_distances = cell.radial_distances()
print(f"Mean radial distance: {radial_distances.mean():.2f}")
```

### Clustering Analysis

```python
# Perform K-means clustering
cluster_labels = cell.cluster_points(
    n_clusters=5,
    method="kmeans"
)

# Add cluster labels as a feature
cell.add_cluster_feature(cluster_labels, "spatial_clusters")

# Hierarchical clustering
hier_labels = cell.cluster_points(
    n_clusters=5,
    method="hierarchical",
    linkage="ward"
)
```

### Neighborhood Analysis

```python
# Find K-nearest neighbors
neighbors = cell.find_knn(k=10)
print(f"Neighbors shape: {neighbors.shape}")

# Calculate local density
density = cell.calculate_density(k=10)
print(f"Mean density: {density.mean():.3f}")
```

## 🎨 Visualization

### 3D Structure Plots

```python
# Basic 3D plot
fig = cell.plot_3d(
    color_by="chromosome",
    size=3.0,
    title="3D Chromatin Structure"
)
fig.show()

# Color by feature
if "H3K4me3" in cell.features:
    fig = cell.plot_3d(
        color_by="H3K4me3",
        colorscale="viridis",
        size=2.5
    )
    fig.show()

# Plot specific chromosome
fig = cell.plot_3d(
    genome_coord="chr1",
    color_by="position",
    size=4.0
)
fig.show()
```

### 2D Analysis Plots

```python
# Distance matrix heatmap
fig = cell.plot_distance_matrix(
    genome_coord="chr1:1000000-10000000"
)
fig.show()

# Feature distribution
if "RNA_expression" in cell.features:
    fig = cell.plot_feature_distribution("RNA_expression")
    fig.show()

# Radial position analysis
fig = cell.plot_radial_distribution()
fig.show()
```

## 🔬 Structure Analysis

### Basic Structure Parameters

```python
# Calculate radius of gyration
rg = cell.radius_of_gyration()
print(f"Radius of gyration: {rg:.2f}")

# Calculate nuclear volume
volume = cell.nuclear_volume()
print(f"Nuclear volume: {volume:.2f}")

# Calculate asphericity
asphericity = cell.asphericity()
print(f"Asphericity: {asphericity:.3f}")
```

### Chromosome-specific Analysis

```python
# Analyze individual chromosomes
for chrom in cell.chromosomes():
    rg_chrom = cell.radius_of_gyration(chrom=chrom)
    print(f"{chrom} radius of gyration: {rg_chrom:.2f}")

# Calculate chromosome territories
territories = cell.chromosome_territories()
print("Chromosome territories calculated")
```

## 💾 Data Export

### Export Structure Data

```python
# Export to CIF format (for molecular viewers)
cell.to_cif(
    output_path="structure.cif",
    chain_by="chromosome",
    include_features=True
)

# Export to PDB format
cell.to_pdb(
    output_path="structure.pdb",
    chain_by="chromosome"
)

# Export coordinates to CSV
coords = cell.get_data()
coords.to_csv("coordinates.csv", index=False)
```

### Export Analysis Results

```python
# Export distance matrix
dist_matrix = cell.calculate_distances()
np.savetxt("distance_matrix.txt", dist_matrix)

# Export features
if cell.features:
    feature_data = cell.get_feature_data()
    feature_data.to_csv("features.csv", index=False)
```

## 🔄 Working with Multiple Cells

```python
from charmtools.obj import MultiCell3D

# Create multiple cells
cells = {
    "cell1": Cell3D(cellname="cell1", resolution=40000, tdg_path="cell1.3dg"),
    "cell2": Cell3D(cellname="cell2", resolution=40000, tdg_path="cell2.3dg"),
    "cell3": Cell3D(cellname="cell3", resolution=40000, tdg_path="cell3.3dg")
}

# Create MultiCell3D object
multi_cell = MultiCell3D(cells)

# Compare structures
comparison = multi_cell.compare_structures()
print("Structure comparison completed")

# Batch analysis
batch_results = multi_cell.batch_analysis([
    "radius_of_gyration",
    "nuclear_volume",
    "asphericity"
])
print(batch_results)
```

## 💡 Performance Tips

### Memory Management

```python
# For large datasets, use disk storage
large_cell = Cell3D(
    cellname="large_dataset",
    resolution=10000,
    tdg_path="large_data.3dg",
    on_disk=True,
    on_disk_path="./cache/"
)

# Load data to memory when needed
large_cell.to_memory()

# Move back to disk to free memory
large_cell.to_disk()
```

### Efficient Data Access

```python
# Use specific regions instead of whole genome
region_data = cell.get_data("chr1:1000000-5000000")

# Use dense format for complete regions
dense_data = cell.get_data("chr1", if_dense=True)

# Cache frequently used results
cell.build_kdtree()  # Build spatial index for fast queries
```

## 🚨 Common Pitfalls

### 1. Resolution Mismatch
```python
# Make sure resolution matches your data
# If your data has 40kb bins, use resolution=40000
cell = Cell3D(cellname="test", resolution=40000)  # Correct
# cell = Cell3D(cellname="test", resolution=10000)  # Wrong!
```

### 2. Memory Issues
```python
# For large datasets, always use on_disk=True
large_cell = Cell3D(
    cellname="large",
    resolution=10000,
    on_disk=True  # Essential for large data
)
```

### 3. Feature Coordinate Systems
```python
# Make sure feature coordinates match your resolution
# Features should be aligned to resolution boundaries
feature_pos = (feature_pos // resolution) * resolution
```

## 🎯 Next Steps

Now that you've learned the basics, explore more advanced topics:

- **[Cell3D Objects](Cell3D-Objects)** - Deep dive into 3D structure analysis
- **[Data Preprocessing](Data-Preprocessing)** - Clean and prepare your data
- **[Analysis Tools](Analysis-Tools)** - Advanced statistical methods
- **[Examples](Examples)** - Real-world case studies
- **[API Reference](API-Reference)** - Complete function documentation

## 📞 Need Help?

If you run into issues:

1. Check the [API Reference](API-Reference) for detailed function documentation
2. Browse [Examples](Examples) for similar use cases
3. Search [GitHub Issues](https://github.com/skelviper/CHARMtools/issues)
4. Ask questions in GitHub Discussions

---

**Ready for more?** Continue to [Cell3D Objects](Cell3D-Objects) for detailed analysis workflows!