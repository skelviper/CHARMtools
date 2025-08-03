# Visualization

CHARMtools provides comprehensive visualization capabilities for 3D chromatin structure data, enabling researchers to create publication-ready figures and interactive visualizations for data exploration and presentation.

## Overview

The visualization system in CHARMtools supports:
- 3D structure visualization
- 2D contact matrix heatmaps
- Statistical plots and charts
- Interactive visualizations
- Comparative analysis plots
- Publication-ready figure generation

## Visualization Types

### 1. 3D Structure Visualization

**Purpose**: Visualize 3D chromatin structures in three-dimensional space.

#### Basic 3D Plots

```python
from CHARMtools.obj.Cell3D import Cell3D

# Load Cell3D object
cell = Cell3D(tdg_path="structure.3dg")

# Basic 3D structure plot
cell.plot_3d_structure(
    chromosomes=['chr1', 'chr2'],  # Specific chromosomes
    color_by='chromosome',         # Color scheme
    size=2,                       # Point size
    alpha=0.7,                    # Transparency
    show_connections=True         # Show backbone connections
)
```

#### Advanced 3D Visualization

```python
# Feature-based coloring
cell.add_feature_from_bed("enhancers.bed", feature_name="enhancers")
cell.plot_3d_structure(
    color_by='enhancers',         # Color by genomic features
    colormap='viridis',           # Custom colormap
    background_color='white',     # Background color
    camera_position='top'         # Camera angle
)

# Interactive 3D plot
cell.plot_3d_interactive(
    features=['enhancers', 'promoters'],
    enable_selection=True,
    show_labels=True
)
```

#### 3D Structure Comparison

```python
from CHARMtools.obj.MultiCell3D import MultiCell3D

# Compare multiple structures
multi_cell = MultiCell3D([cell1, cell2, cell3])
multi_cell.plot_comparison_3d(
    layout='grid',                # Layout arrangement
    sync_camera=True,            # Synchronized camera
    color_scheme='condition'     # Color by condition
)
```

### 2. Contact Matrix Visualization

**Purpose**: Display contact frequency matrices as heatmaps.

#### Basic Contact Heatmaps

```python
# Generate contact matrix heatmap
cell.plot_contact_matrix(
    chromosome='chr1',
    resolution=40000,
    region='chr1:1000000-5000000',
    colormap='Reds',
    log_scale=True,
    show_colorbar=True
)
```

#### Advanced Matrix Visualization

```python
# Multi-resolution heatmap
cell.plot_multi_resolution_matrix(
    resolutions=[10000, 40000, 100000],
    region='chr1:1000000-3000000',
    layout='triangular',
    normalize=True
)

# Comparative contact matrices
cell.plot_differential_matrix(
    other_cell=cell2,
    region='chr1:1000000-5000000',
    method='log2_ratio',
    significance_threshold=0.05
)
```

### 3. Genomic Feature Visualization

**Purpose**: Visualize genomic annotations and features in relation to 3D structure.

#### Feature Distribution Plots

```python
# Feature enrichment around structures
cell.plot_feature_enrichment(
    feature_name='enhancers',
    distance_range=(-50000, 50000),
    bin_size=5000,
    normalize=True
)

# Feature co-localization
cell.plot_feature_colocalization(
    features=['enhancers', 'promoters'],
    distance_threshold=10000,
    plot_type='scatter'
)
```

#### Genomic Track Visualization

```python
# Integrated genomic tracks
cell.plot_genomic_tracks(
    region='chr1:1000000-2000000',
    tracks=[
        {'type': 'structure', 'data': cell},
        {'type': 'features', 'data': 'enhancers.bed'},
        {'type': 'expression', 'data': 'expression.bw'}
    ],
    layout='vertical'
)
```

### 4. Spatial Analysis Visualization

**Purpose**: Visualize results from spatial analyses.

#### Distance Analysis Plots

```python
# Distance distribution plots
cell.plot_distance_distribution(
    feature_pairs=[('enhancers', 'promoters')],
    distance_type='euclidean',
    bins=50,
    show_statistics=True
)

# Spatial clustering visualization
cell.plot_spatial_clusters(
    clustering_method='kmeans',
    n_clusters=5,
    color_by_cluster=True,
    show_centroids=True
)
```

#### Neighborhood Analysis

```python
# Neighborhood enrichment heatmap
cell.plot_neighborhood_enrichment(
    features=['enhancers', 'promoters', 'insulators'],
    radius=50000,
    method='fisher_exact'
)
```

### 5. Comparative Visualization

**Purpose**: Compare structures and features across conditions or time points.

#### Structure Comparison

```python
# Side-by-side structure comparison
cell.plot_structure_comparison(
    other_cells=[cell2, cell3],
    alignment_method='procrustes',
    highlight_differences=True,
    difference_threshold=5.0
)

# Temporal dynamics visualization
multi_cell.plot_temporal_dynamics(
    time_points=[0, 2, 4, 6, 8],
    feature='compartments',
    animation=True
)
```

#### Statistical Comparison Plots

```python
# Box plots for feature comparisons
multi_cell.plot_feature_boxplots(
    features=['contact_density', 'clustering_coefficient'],
    group_by='condition',
    statistical_test='wilcoxon'
)

# Correlation heatmaps
multi_cell.plot_correlation_heatmap(
    features=['all'],
    method='pearson',
    cluster_features=True
)
```

### 6. Quality Control Visualization

**Purpose**: Visualize data quality metrics and validation results.

#### Data Quality Plots

```python
# Contact coverage visualization
cell.plot_coverage_metrics(
    resolution=40000,
    plot_type='histogram',
    show_statistics=True
)

# Validation plots
cell.plot_validation_metrics(
    metrics=['contact_probability', 'distance_decay'],
    reference_data='reference.cool',
    show_confidence_intervals=True
)
```

## Customization Options

### Color Schemes and Palettes

```python
# Custom color palettes
from CHARMtools.utils.generateColor2 import generate_palette

# Generate custom colors
custom_colors = generate_palette(
    n_colors=23,
    palette_type='categorical',
    colorblind_friendly=True
)

# Apply to visualization
cell.plot_3d_structure(
    color_by='chromosome',
    custom_colors=custom_colors
)
```

### Plot Styling

```python
# Set global plot style
import matplotlib.pyplot as plt
plt.style.use('seaborn-v0_8-whitegrid')

# Custom styling options
plot_config = {
    'figure_size': (12, 8),
    'dpi': 300,
    'font_family': 'Arial',
    'font_size': 12,
    'line_width': 2,
    'marker_size': 4
}

cell.set_plot_config(plot_config)
```

### Export Options

```python
# High-resolution export
cell.plot_3d_structure(
    save_path='structure_3d.png',
    dpi=300,
    format='png',
    transparent_background=True
)

# Vector format export
cell.plot_contact_matrix(
    save_path='contact_matrix.svg',
    format='svg',
    optimize_size=True
)

# Interactive HTML export
cell.plot_3d_interactive(
    save_path='interactive_structure.html',
    include_controls=True
)
```

## Interactive Visualizations

### Web-based Interactive Plots

```python
# Interactive 3D structure explorer
cell.create_interactive_explorer(
    features=['enhancers', 'promoters'],
    enable_selection=True,
    show_feature_info=True,
    export_selections=True
)

# Interactive contact matrix browser
cell.create_matrix_browser(
    multi_resolution=True,
    enable_zoom=True,
    show_annotations=True
)
```

### Jupyter Notebook Integration

```python
# Jupyter widget integration
from CHARMtools.visualization.widgets import StructureWidget

# Interactive widget
widget = StructureWidget(cell)
widget.add_controls(['chromosome', 'resolution', 'features'])
widget.display()
```

## Batch Visualization

### Automated Report Generation

```python
# Generate comprehensive visualization report
from CHARMtools.visualization.reports import generate_analysis_report

report = generate_analysis_report(
    cells=[cell1, cell2, cell3],
    output_dir='visualization_report',
    include_sections=[
        'structure_overview',
        'contact_analysis',
        'feature_analysis',
        'comparative_analysis'
    ],
    format='html'
)
```

### Batch Processing

```python
# Process multiple cells
from CHARMtools.visualization.batch import batch_visualize

batch_visualize(
    input_dir='cell_structures/',
    output_dir='visualizations/',
    plot_types=['3d_structure', 'contact_matrix'],
    parallel=True,
    n_jobs=4
)
```

## Performance Optimization

### Large Dataset Handling

```python
# Optimize for large datasets
cell.set_visualization_mode('performance')

# Use data subsampling
cell.plot_3d_structure(
    subsample_rate=0.1,  # Use 10% of data points
    level_of_detail=True  # Adaptive detail levels
)

# Progressive loading
cell.plot_contact_matrix(
    progressive_loading=True,
    chunk_size=1000000
)
```

### Memory Management

```python
# Memory-efficient plotting
cell.plot_large_structure(
    memory_limit='4GB',
    use_disk_cache=True,
    compression='lz4'
)
```

## Best Practices

### 1. Publication-Ready Figures

- Use high DPI (300+) for print publications
- Choose colorblind-friendly palettes
- Ensure proper font sizes and line weights
- Include scale bars and legends
- Use vector formats (SVG, PDF) when possible

### 2. Data Exploration

- Start with overview plots before detailed analysis
- Use interactive visualizations for exploration
- Apply appropriate data transformations (log scale, normalization)
- Include quality control visualizations

### 3. Comparative Analysis

- Ensure consistent scales across comparisons
- Use appropriate statistical visualizations
- Highlight significant differences
- Provide context with reference data

### 4. Performance

- Subsample large datasets for initial exploration
- Use appropriate plot types for data size
- Cache intermediate results
- Consider parallel processing for batch operations

## Troubleshooting

### Common Issues

1. **Memory Errors**
   - Reduce data size or use subsampling
   - Increase system memory
   - Use disk-based caching

2. **Slow Rendering**
   - Optimize data structures
   - Use appropriate plot types
   - Enable hardware acceleration

3. **Display Issues**
   - Check display backend compatibility
   - Verify color space settings
   - Update graphics drivers

### Getting Help

For visualization issues:
- Check the [Troubleshooting Guide](Troubleshooting.md)
- Review example notebooks
- Contact support for specific rendering problems

## Integration with External Tools

### Export to Other Software

```python
# Export for external visualization
cell.export_for_chimera('structure.pdb')  # ChimeraX
cell.export_for_pymol('structure.pse')    # PyMOL
cell.export_for_vmd('structure.tcl')      # VMD
```

### Web Integration

```python
# Generate web-compatible formats
cell.export_web_viewer(
    output_dir='web_viewer/',
    include_controls=True,
    responsive_design=True
)
```

## Next Steps

After creating visualizations:
1. Integrate with [Analysis Tools](Analysis-Tools.md) for comprehensive workflows
2. Use [Cell3D Objects](Cell3D-Objects.md) for data management
3. Export results using [Utilities](Utilities.md) functions

---

*For more information about CHARMtools, visit the [Home](Home.md) page or check the [Installation Guide](Installation-Guide.md).*