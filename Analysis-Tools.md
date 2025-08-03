# Analysis Tools

The `analysis` module provides a comprehensive suite of tools for analyzing 3D chromatin structure data. These tools enable various types of genomic and spatial analyses commonly used in chromatin biology research.

## Overview

The analysis module includes tools for:
- Topological domain analysis
- Compartment identification
- Loop detection and analysis
- Spatial statistics
- Single-cell specific analyses
- Comparative genomics

## Available Analysis Tools

### 1. Topological Associated Domains (TAD.py)

**Purpose**: Identify and analyze topological associated domains in 3D chromatin structure.

**Key Functions**:
- TAD boundary detection
- Domain size analysis
- Boundary strength calculation
- Cross-cell TAD comparison

**Input**: Contact matrices or Cell3D objects
**Output**: TAD boundaries, domain annotations, statistical summaries

### 2. Compartment Analysis (compartment.py)

**Purpose**: Identify A/B compartments and analyze compartmentalization patterns.

**Key Functions**:
- Principal component analysis for compartment identification
- Compartment strength calculation
- Switching analysis between conditions
- Visualization of compartment patterns

**Input**: Contact matrices, genomic features
**Output**: Compartment assignments, strength scores, switching regions

### 3. Loop Detection (loop.py)

**Purpose**: Detect and analyze chromatin loops from contact data.

**Key Functions**:
- Loop calling algorithms
- Loop strength quantification
- Loop clustering and classification
- Differential loop analysis

**Input**: Contact pairs, Cell3D objects
**Output**: Loop coordinates, strength scores, classifications

### 4. Spatial Statistics (spatialstat.py)

**Purpose**: Perform spatial statistical analyses on 3D chromatin structures.

**Key Functions**:
- Spatial clustering analysis
- Distance distribution analysis
- Neighborhood enrichment
- Spatial correlation metrics

**Input**: 3D coordinates, genomic annotations
**Output**: Statistical measures, clustering results, correlation matrices

### 5. Single-Cell A/B Compartments (scAB.py)

**Purpose**: Single-cell specific compartment analysis.

**Key Functions**:
- Single-cell compartment calling
- Cell-to-cell variability analysis
- Population-level compartment consensus
- Compartment dynamics analysis

**Input**: Single-cell contact data
**Output**: Cell-specific compartment calls, variability metrics

### 6. Single-Cell Gene Activity Domains (scGAD.py)

**Purpose**: Analyze gene activity domains in single cells.

**Key Functions**:
- Gene activity domain identification
- Expression-structure correlation
- Domain boundary analysis
- Cell type specific patterns

**Input**: Single-cell data, gene expression
**Output**: Activity domains, correlation scores

### 7. Coverage Analysis (coverage.py)

**Purpose**: Analyze contact coverage and data quality metrics.

**Key Functions**:
- Contact density calculation
- Coverage uniformity assessment
- Quality control metrics
- Bias detection and correction

**Input**: Contact data
**Output**: Coverage maps, quality metrics

### 8. Peak Enrichment (peak_enrichment.py)

**Purpose**: Analyze enrichment of genomic features at specific loci.

**Key Functions**:
- Peak calling and annotation
- Enrichment score calculation
- Statistical significance testing
- Comparative enrichment analysis

**Input**: Contact data, genomic features
**Output**: Peak locations, enrichment scores, p-values

### 9. TSS Enrichment (tss_enrichment.py)

**Purpose**: Analyze transcription start site (TSS) enrichment patterns.

**Key Functions**:
- TSS contact enrichment
- Promoter interaction analysis
- Gene expression correlation
- Regulatory network inference

**Input**: Contact data, TSS annotations
**Output**: Enrichment profiles, interaction networks

### 10. Saddle Plot Analysis (saddle.py)

**Purpose**: Generate saddle plots for compartment strength visualization.

**Key Functions**:
- Saddle plot generation
- Compartment strength quantification
- Comparative saddle analysis
- Statistical testing

**Input**: Contact matrices, compartment data
**Output**: Saddle plots, strength metrics

### 11. Regression Analysis (regression.py)

**Purpose**: Perform regression analyses on chromatin structure data.

**Key Functions**:
- Linear and non-linear regression
- Feature importance analysis
- Predictive modeling
- Cross-validation

**Input**: Structural features, genomic annotations
**Output**: Regression models, predictions, feature importance

### 12. Differential Analysis (simpleDiff.py)

**Purpose**: Compare chromatin structures between conditions.

**Key Functions**:
- Differential contact analysis
- Statistical testing
- Effect size calculation
- Multiple testing correction

**Input**: Multiple Cell3D objects or contact matrices
**Output**: Differential regions, statistical results

### 13. Imputation (imputation.py)

**Purpose**: Impute missing contacts in sparse data.

**Key Functions**:
- Matrix completion algorithms
- Missing data imputation
- Quality assessment
- Cross-validation

**Input**: Sparse contact matrices
**Output**: Imputed matrices, quality metrics

### 14. Cell Cycle Analysis (cellcycle.py)

**Purpose**: Analyze cell cycle-related chromatin changes.

**Key Functions**:
- Cell cycle phase classification
- Phase-specific structure analysis
- Temporal dynamics
- Progression markers

**Input**: Time-series or cell cycle data
**Output**: Phase classifications, dynamic patterns

### 15. Escaping Score (calc_escaping_score.py)

**Purpose**: Calculate compartment escaping scores.

**Key Functions**:
- Escaping score calculation
- Boundary crossing analysis
- Compartment stability
- Dynamic range assessment

**Input**: Compartment data, contact matrices
**Output**: Escaping scores, stability metrics

### 16. Visualization (plot.py)

**Purpose**: Generate publication-ready plots for analysis results.

**Key Functions**:
- Contact matrix heatmaps
- 3D structure visualization
- Statistical plots
- Comparative visualizations

**Input**: Analysis results, Cell3D objects
**Output**: Publication-ready figures

## Usage Examples

### Basic Analysis Workflow

```python
from CHARMtools.obj.Cell3D import Cell3D
from CHARMtools.analysis import TAD, compartment, spatialstat

# Load data
cell = Cell3D(tdg_path="structure.3dg")

# TAD analysis
tad_results = TAD.detect_domains(cell)

# Compartment analysis
compartments = compartment.identify_compartments(cell)

# Spatial statistics
spatial_stats = spatialstat.calculate_spatial_metrics(cell)
```

### Comparative Analysis

```python
from CHARMtools.analysis import simpleDiff

# Compare two conditions
cell1 = Cell3D(tdg_path="condition1.3dg")
cell2 = Cell3D(tdg_path="condition2.3dg")

# Differential analysis
diff_results = simpleDiff.compare_structures(cell1, cell2)
```

## Integration with Core Objects

All analysis tools are designed to work seamlessly with:
- **Cell3D Objects**: Primary input for most analyses
- **MultiCell3D Objects**: Population-level analyses
- **Standard Data Formats**: Direct matrix and coordinate input

## Performance Considerations

- Most analyses support parallel processing
- Memory usage scales with data size
- Some algorithms may require significant computational resources
- Consider data subsampling for exploratory analysis

## Best Practices

1. **Data Quality**: Ensure proper preprocessing before analysis
2. **Parameter Selection**: Use appropriate parameters for your data type
3. **Statistical Testing**: Apply multiple testing corrections where appropriate
4. **Validation**: Cross-validate results using independent datasets
5. **Documentation**: Keep detailed records of analysis parameters

## Getting Help

For analysis-specific questions:
- Check individual module documentation
- Review example workflows in [Quick Start](Quick-Start.md)
- Consult the [Troubleshooting Guide](Troubleshooting.md)
- Contact the development team for advanced usage

## Next Steps

After running analyses:
1. Use [Visualization](Visualization.md) tools to explore results
2. Export data using [Cell3D Objects](Cell3D-Objects.md) methods
3. Integrate results with external tools and databases

---

*For more information about CHARMtools, visit the [Home](Home.md) page or check the [Installation Guide](Installation-Guide.md).*