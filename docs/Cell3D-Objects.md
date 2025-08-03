# Cell3D Objects

The Cell3D object is the core data structure in CHARMtools for representing and analyzing 3D chromatin structures. This guide covers everything you need to know about working with Cell3D objects.

## 📋 Table of Contents

- [Overview](#overview)
- [Creating Cell3D Objects](#creating-cell3d-objects)
- [Data Management](#data-management)
- [Genomic Features](#genomic-features)
- [Spatial Analysis](#spatial-analysis)
- [Structure Analysis](#structure-analysis)
- [Visualization](#visualization)
- [Data Export](#data-export)
- [Advanced Usage](#advanced-usage)
- [MultiCell3D Objects](#multicell3d-objects)

## 🔬 Overview

### What is a Cell3D Object?

A Cell3D object represents the 3D chromatin structure of a single cell, containing:

- **3D coordinates** of genomic loci
- **Genomic features** (ChIP-seq, RNA-seq, etc.)
- **Spatial analysis** results
- **Metadata** and annotations

### Architecture

Cell3D uses a modular architecture with specialized components:

```
Cell3D
├── Cell3DCore      # Basic data management
├── Cell3DData      # Data loading and processing
├── Cell3DFeatures  # Genomic feature handling
├── Cell3DSpatial   # Spatial analysis methods
├── Cell3DStructure # Structure analysis
├── Cell3DVisualization # Plotting and visualization
├── Cell3DIO        # Input/output operations
├── Cell3DAnalysis  # Advanced analysis methods
└── Cell3DUtils     # Utility functions
```

## 🏗️ Creating Cell3D Objects

### Basic Creation

```python
from charmtools.obj import Cell3D

# Create empty Cell3D object
cell = Cell3D(
    cellname="GM12878",
    resolution=40000
)
```

### From 3DG File

```python
# Load from 3DG file
cell = Cell3D(
    cellname="GM12878",
    resolution=40000,
    tdg_path="data/GM12878.3dg"
)
```

### From DataFrame

```python
import pandas as pd
import numpy as np

# Create sample data
data = pd.DataFrame({
    'chrom': ['chr1'] * 1000,
    'pos': range(0, 40000000, 40000),
    'x': np.random.randn(1000),
    'y': np.random.randn(1000),
    'z': np.random.randn(1000)
})

# Create Cell3D from DataFrame
cell = Cell3D(
    cellname="sample",
    resolution=40000,
    tdg_path=data
)
```

### Memory vs Disk Storage

```python
# In-memory storage (default, faster access)
cell = Cell3D(
    cellname="small_cell",
    resolution=40000,
    tdg_path="data.3dg",
    on_disk=False
)

# Disk storage (memory efficient for large datasets)
cell = Cell3D(
    cellname="large_cell",
    resolution=10000,
    tdg_path="large_data.3dg",
    on_disk=True,
    on_disk_path="./cache/"
)
```

## 💾 Data Management

### Basic Properties

```python
# Basic information
print(f"Cell name: {cell.cellname}")
print(f"Resolution: {cell.resolution}")
print(f"Number of points: {len(cell)}")
print(f"Storage mode: {'Disk' if cell.on_disk else 'Memory'}")

# Chromosomes and ranges
chromosomes = cell.chromosomes()
print(f"Chromosomes: {chromosomes}")

ranges = cell.coordinate_ranges()
print(f"Coordinate ranges: {ranges}")

# Center of mass
com = cell.center_of_mass()
print(f"Center of mass: {com}")
```

### Data Access

```python
# Get all data
all_data = cell.get_data()
print(f"All data shape: {all_data.shape}")

# Get chromosome data
chr1_data = cell.get_data("chr1")
print(f"Chr1 data shape: {chr1_data.shape}")

# Get region data
region_data = cell.get_data("chr1:1000000-5000000")
print(f"Region data shape: {region_data.shape}")

# Get dense data (fill missing positions)
dense_data = cell.get_data("chr1", if_dense=True)
print(f"Dense data shape: {dense_data.shape}")
```

### Data Transformation

```python
# Rotate structure
rotated_data = cell.get_data(
    rotate=True,
    rotate_x_angle=np.pi/4,
    rotate_y_angle=np.pi/6,
    rotate_z_angle=np.pi/3
)

# Random rotation
random_rotated = cell.get_data(rotate=True)
```

### Subsetting

```python
# Create subset (new object)
chr1_cell = cell.subset("chr1")
print(f"Subset cell: {chr1_cell.cellname}")

# Subset with query
active_regions = cell.subset(
    genome_coord="chr1",
    query="H3K4me3 > 0.5"  # Requires H3K4me3 feature
)

# In-place subsetting
cell.subset("chr1:1000000-10000000", in_place=True)
```

### Memory Management

```python
# Move to disk
cell.to_disk("./cache/")
print(f"On disk: {cell.on_disk}")

# Move to memory
cell.to_memory()
print(f"On disk: {cell.on_disk}")

# Copy object
cell_copy = cell.copy()
```

## 🧬 Genomic Features

### Adding Features from Files

```python
# Add BED file (peaks, regions)
cell.add_bed_feature(
    bed_path="data/H3K4me3_peaks.bed",
    feature_name="H3K4me3_peaks",
    binary=True,  # Convert to binary (0/1)
    value_column=None  # Use presence/absence
)

# Add bedGraph file (continuous values)
cell.add_bedgraph_feature(
    bedgraph_path="data/H3K27ac.bedgraph",
    feature_name="H3K27ac_signal",
    value_column=4
)

# Add RNA-seq data
cell.add_rna_feature(
    rna_path="data/expression.txt",
    feature_name="RNA_expression",
    gene_column="gene_name",
    value_column="FPKM"
)
```

### Adding Features from DataFrames

```python
# Create feature data
feature_df = pd.DataFrame({
    'chrom': ['chr1', 'chr1', 'chr2'],
    'start': [1000000, 2000000, 1000000],
    'end': [1040000, 2040000, 1040000],
    'value': [2.5, 1.8, 3.2]
})

# Add feature
cell.add_feature_from_dataframe(
    df=feature_df,
    feature_name="custom_signal",
    value_column="value"
)
```

### Feature Management

```python
# List features
print(f"Available features: {cell.features}")

# Get feature data
feature_data = cell.get_feature_data()
print(feature_data.head())

# Get specific feature
h3k4me3_data = cell.get_feature_data("H3K4me3_peaks")

# Remove feature
cell.remove_feature("old_feature")

# Rename feature
cell.rename_feature("old_name", "new_name")
```

### Feature Analysis

```python
# Feature statistics
stats = cell.feature_statistics("H3K27ac_signal")
print(f"Feature stats: {stats}")

# Feature correlation
corr_matrix = cell.feature_correlation()
print(f"Feature correlation:\n{corr_matrix}")

# Feature enrichment in clusters
enrichment = cell.feature_enrichment_in_clusters(
    feature_name="H3K4me3_peaks",
    cluster_feature="spatial_clusters"
)
```

## 📍 Spatial Analysis

### Distance Calculations

```python
# Pairwise distances
distances = cell.calculate_distances(metric="euclidean")
print(f"Distance matrix shape: {distances.shape}")

# Distances to center
radial_dist = cell.radial_distances()
print(f"Mean radial distance: {radial_dist.mean():.2f}")

# Distances between chromosomes
inter_chrom_dist = cell.inter_chromosome_distances()
print(f"Inter-chromosome distances calculated")
```

### Clustering

```python
# K-means clustering
kmeans_labels = cell.cluster_points(
    n_clusters=5,
    method="kmeans",
    features=None  # Use coordinates only
)

# Hierarchical clustering
hier_labels = cell.cluster_points(
    n_clusters=5,
    method="hierarchical",
    linkage="ward",
    metric="euclidean"
)

# Feature-based clustering
feature_labels = cell.cluster_points(
    n_clusters=3,
    method="kmeans",
    features=["H3K4me3_peaks", "H3K27ac_signal"]
)

# Add cluster labels as features
cell.add_cluster_feature(kmeans_labels, "spatial_clusters")
cell.add_cluster_feature(hier_labels, "hierarchical_clusters")
```

### Neighborhood Analysis

```python
# Build spatial index
cell.build_kdtree()

# Find K-nearest neighbors
neighbors = cell.find_knn(k=10)
print(f"Neighbors shape: {neighbors.shape}")

# Calculate local density
density = cell.calculate_density(k=10, method="knn")
print(f"Mean density: {density.mean():.3f}")

# Find neighbors within radius
radius_neighbors = cell.find_neighbors_radius(radius=2.0)
```

### Spatial Statistics

```python
# Spatial autocorrelation
autocorr = cell.spatial_autocorrelation("H3K4me3_peaks")
print(f"Spatial autocorrelation: {autocorr:.3f}")

# Ripley's K function
ripley_k = cell.ripleys_k_function(distances=np.linspace(0.5, 5.0, 10))

# Nearest neighbor distances
nn_distances = cell.nearest_neighbor_distances()
print(f"Mean NN distance: {nn_distances.mean():.2f}")
```

## 🏗️ Structure Analysis

### Basic Structure Parameters

```python
# Radius of gyration
rg_total = cell.radius_of_gyration()
print(f"Total radius of gyration: {rg_total:.2f}")

# Per chromosome
for chrom in cell.chromosomes():
    rg_chrom = cell.radius_of_gyration(chrom=chrom)
    print(f"{chrom} Rg: {rg_chrom:.2f}")

# Nuclear volume
volume = cell.nuclear_volume(method="convex_hull")
print(f"Nuclear volume: {volume:.2f}")

# Asphericity
asphericity = cell.asphericity()
print(f"Asphericity: {asphericity:.3f}")

# Compactness
compactness = cell.compactness()
print(f"Compactness: {compactness:.3f}")
```

### Advanced Structure Analysis

```python
# Chromosome territories
territories = cell.chromosome_territories(method="convex_hull")
print(f"Territory volumes: {territories}")

# Chromosome positioning
positions = cell.chromosome_positioning()
print(f"Chromosome positions: {positions}")

# Structure comparison
other_cell = Cell3D(cellname="other", resolution=40000, tdg_path="other.3dg")
comparison = cell.compare_structure(other_cell)
print(f"Structure similarity: {comparison['similarity']:.3f}")

# Fractal dimension
fractal_dim = cell.fractal_dimension()
print(f"Fractal dimension: {fractal_dim:.3f}")
```

### Contact Analysis

```python
# Contact probability vs distance
contact_prob = cell.contact_probability_vs_distance(
    distances=np.logspace(0, 2, 20)
)

# Compartment analysis
compartments = cell.compartment_analysis(
    feature_name="H3K4me3_peaks",
    n_compartments=2
)

# TAD-like domain detection
domains = cell.detect_domains(
    method="density_based",
    min_size=5,
    density_threshold=0.1
)
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
fig = cell.plot_3d(
    color_by="H3K4me3_peaks",
    colorscale="viridis",
    size=2.5,
    opacity=0.8
)
fig.show()

# Plot specific region
fig = cell.plot_3d(
    genome_coord="chr1:1000000-10000000",
    color_by="position",
    size=4.0,
    show_connections=True
)
fig.show()

# Interactive plot with hover info
fig = cell.plot_3d(
    color_by="RNA_expression",
    hover_data=["chrom", "pos", "H3K4me3_peaks"],
    size=3.0
)
fig.show()
```

### 2D Analysis Plots

```python
# Distance matrix heatmap
fig = cell.plot_distance_matrix(
    genome_coord="chr1",
    colorscale="viridis"
)
fig.show()

# Feature distribution
fig = cell.plot_feature_distribution(
    feature_name="H3K27ac_signal",
    plot_type="histogram"
)
fig.show()

# Radial distribution
fig = cell.plot_radial_distribution(
    feature_name="H3K4me3_peaks",
    bins=20
)
fig.show()

# Cluster visualization
fig = cell.plot_clusters(
    cluster_feature="spatial_clusters",
    plot_type="3d"
)
fig.show()
```

### Comparative Plots

```python
# Feature correlation heatmap
fig = cell.plot_feature_correlation()
fig.show()

# Structure parameters radar plot
fig = cell.plot_structure_parameters()
fig.show()

# Contact probability plot
fig = cell.plot_contact_probability()
fig.show()
```

## 💾 Data Export

### Structure Export

```python
# Export to CIF format
cell.to_cif(
    output_path="structure.cif",
    chain_by="chromosome",
    include_features=True,
    feature_as_bfactor="H3K4me3_peaks"
)

# Export to PDB format
cell.to_pdb(
    output_path="structure.pdb",
    chain_by="chromosome",
    atom_type="CA"
)

# Export to XYZ format
cell.to_xyz(
    output_path="structure.xyz",
    include_features=False
)
```

### Data Export

```python
# Export coordinates
coords = cell.get_data()
coords.to_csv("coordinates.csv", index=False)

# Export features
features = cell.get_feature_data()
features.to_csv("features.csv", index=False)

# Export analysis results
results = {
    "radius_of_gyration": cell.radius_of_gyration(),
    "nuclear_volume": cell.nuclear_volume(),
    "asphericity": cell.asphericity()
}

import json
with open("analysis_results.json", "w") as f:
    json.dump(results, f, indent=2)
```

### Custom Export

```python
# Export specific regions
region_data = cell.get_data("chr1:1000000-5000000")
region_data.to_csv("chr1_region.csv", index=False)

# Export distance matrix
dist_matrix = cell.calculate_distances()
np.savetxt("distance_matrix.txt", dist_matrix)

# Export cluster assignments
if "spatial_clusters" in cell.features:
    cluster_data = cell.get_feature_data("spatial_clusters")
    cluster_data.to_csv("clusters.csv", index=False)
```

## 🚀 Advanced Usage

### Custom Analysis Pipelines

```python
def custom_analysis_pipeline(cell):
    """Custom analysis pipeline"""
    results = {}
    
    # Structure analysis
    results['structure'] = {
        'radius_of_gyration': cell.radius_of_gyration(),
        'nuclear_volume': cell.nuclear_volume(),
        'asphericity': cell.asphericity()
    }
    
    # Spatial analysis
    clusters = cell.cluster_points(n_clusters=5, method="kmeans")
    cell.add_cluster_feature(clusters, "analysis_clusters")
    
    results['spatial'] = {
        'n_clusters': 5,
        'cluster_sizes': np.bincount(clusters)
    }
    
    # Feature analysis
    if cell.features:
        results['features'] = {
            'correlation': cell.feature_correlation().to_dict(),
            'statistics': {feat: cell.feature_statistics(feat) 
                         for feat in cell.features}
        }
    
    return results

# Run custom pipeline
results = custom_analysis_pipeline(cell)
print(json.dumps(results, indent=2, default=str))
```

### Batch Processing

```python
def process_multiple_cells(cell_paths, resolution=40000):
    """Process multiple cells in batch"""
    results = {}
    
    for cell_name, path in cell_paths.items():
        print(f"Processing {cell_name}...")
        
        # Load cell
        cell = Cell3D(
            cellname=cell_name,
            resolution=resolution,
            tdg_path=path,
            on_disk=True  # Use disk for memory efficiency
        )
        
        # Run analysis
        results[cell_name] = custom_analysis_pipeline(cell)
        
        # Clean up
        del cell
    
    return results

# Example usage
cell_paths = {
    "GM12878": "data/GM12878.3dg",
    "K562": "data/K562.3dg",
    "HeLa": "data/HeLa.3dg"
}

batch_results = process_multiple_cells(cell_paths)
```

### Performance Optimization

```python
# Use spatial indexing for repeated queries
cell.build_kdtree()

# Cache frequently used calculations
@functools.lru_cache(maxsize=128)
def cached_distance_calculation(cell_id, metric):
    return cell.calculate_distances(metric=metric)

# Use numba for custom calculations
from numba import jit

@jit(nopython=True)
def fast_custom_metric(coords):
    """Fast custom distance metric"""
    n = coords.shape[0]
    distances = np.zeros((n, n))
    
    for i in range(n):
        for j in range(i+1, n):
            dist = np.sqrt(np.sum((coords[i] - coords[j])**2))
            distances[i, j] = distances[j, i] = dist
    
    return distances

# Use with Cell3D
coords = cell.get_data()[["x", "y", "z"]].values
fast_distances = fast_custom_metric(coords)
```

## 👥 MultiCell3D Objects

### Creating MultiCell3D

```python
from charmtools.obj import MultiCell3D

# Create individual cells
cells = {
    "GM12878": Cell3D(cellname="GM12878", resolution=40000, tdg_path="GM12878.3dg"),
    "K562": Cell3D(cellname="K562", resolution=40000, tdg_path="K562.3dg"),
    "HeLa": Cell3D(cellname="HeLa", resolution=40000, tdg_path="HeLa.3dg")
}

# Create MultiCell3D object
multi_cell = MultiCell3D(cells)
print(f"MultiCell3D with {len(multi_cell)} cells")
```

### Comparative Analysis

```python
# Compare structures
comparison = multi_cell.compare_structures()
print(f"Structure comparison: {comparison}")

# Batch analysis
batch_results = multi_cell.batch_analysis([
    "radius_of_gyration",
    "nuclear_volume",
    "asphericity"
])
print(batch_results)

# Feature comparison
if all("H3K4me3" in cell.features for cell in cells.values()):
    feature_comparison = multi_cell.compare_features("H3K4me3")
    print(f"Feature comparison: {feature_comparison}")
```

### Visualization

```python
# Plot multiple cells
fig = multi_cell.plot_multiple_3d(
    color_by="cell_type",
    layout="grid"
)
fig.show()

# Comparative plots
fig = multi_cell.plot_structure_comparison()
fig.show()

# Feature distribution comparison
fig = multi_cell.plot_feature_comparison("H3K4me3")
fig.show()
```

## 🔧 Troubleshooting

### Common Issues

#### Memory Issues
```python
# Use disk storage for large datasets
cell = Cell3D(cellname="large", resolution=10000, on_disk=True)

# Process data in chunks
for chrom in cell.chromosomes():
    chrom_data = cell.get_data(chrom)
    # Process chrom_data
```

#### Performance Issues
```python
# Build spatial index once
cell.build_kdtree()

# Use appropriate resolution
# Higher resolution = more data points = slower analysis
cell = Cell3D(cellname="test", resolution=40000)  # Good for most analyses
```

#### Feature Issues
```python
# Check feature alignment
print(f"Resolution: {cell.resolution}")
print(f"Feature coordinates should be multiples of {cell.resolution}")

# Validate features
cell.validate_features()
```

### Best Practices

1. **Use appropriate resolution** - 40kb for most analyses, 10kb for detailed studies
2. **Enable disk storage** for large datasets (>100MB)
3. **Build spatial indices** before repeated spatial queries
4. **Validate data** after loading and feature addition
5. **Use specific regions** instead of whole genome when possible
6. **Cache results** for expensive calculations
7. **Clean up objects** when done to free memory

---

**Next Steps:**
- [Data Preprocessing](Data-Preprocessing) - Clean and prepare your data
- [Analysis Tools](Analysis-Tools) - Advanced statistical methods
- [Visualization](Visualization) - Create publication-quality plots
- [Examples](Examples) - Real-world case studies