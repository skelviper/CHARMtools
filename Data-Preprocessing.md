# Data Preprocessing

The `charm_preprocess` module provides essential tools for cleaning and preprocessing 3D chromatin structure data before analysis. These tools ensure data quality and prepare datasets for downstream analysis with CHARMtools.

## Overview

Data preprocessing is a critical step in 3D chromatin analysis that involves:
- Removing artifacts and noise from contact data
- Filtering promiscuous contacts
- Converting between different data formats
- Quality control and validation

## Available Tools

### 1. Promiscuous Contact Cleaning (`clean_leg.py`)

**Purpose**: Removes promiscuous contacts using Tan's phantom leg method to eliminate artifacts from over-represented genomic loci.

**Input**:
- `.pairs` file in 4DN standard format
- Maximum distance threshold (default: varies)
- Maximum count threshold (default: varies)
- Number of threads for parallel processing

**Output**:
- Cleaned `.pairs` file with promiscuous contacts removed

**Usage**:
```bash
python -m CHARMtools.charm_preprocess.clean_leg \
    --input input.pairs.gz \
    --output cleaned.pairs.gz \
    --max-distance 1000 \
    --max-count 10 \
    --threads 8
```

**Key Features**:
- Multi-threaded processing for large datasets
- Configurable distance and count thresholds
- Memory-efficient processing
- Preserves original data format

### 2. Splicing Artifact Removal (`clean_splicing.py`)

**Purpose**: Removes contacts that may result from RNA splicing artifacts by filtering contacts within exonic regions.

**Input**:
- `.pairs` file in 4DN standard format
- GTF annotation file (GENCODE/ENSEMBL format)
- Number of threads for parallel processing

**Output**:
- Cleaned `.pairs` file with splicing artifacts removed

**Usage**:
```bash
python -m CHARMtools.charm_preprocess.clean_splicing \
    --input input.pairs.gz \
    --gtf annotations.gtf \
    --output cleaned.pairs.gz \
    --threads 8
```

**Key Features**:
- GTF-based exon annotation parsing
- Interval-based efficient searching
- Support for multiple genome assemblies
- Special handling for rabbit genome

### 3. Format Conversion (`tdg2pairs.py`)

**Purpose**: Converts 3D genome structure files (.3dg) to contact pairs format (.pairs) using spatial proximity.

**Input**:
- `.3dg` file containing 3D coordinates
- Distance threshold for contact definition

**Output**:
- `.pairs` file with contacts derived from 3D proximity

**Usage**:
```bash
python -m CHARMtools.charm_preprocess.tdg2pairs \
    --input structure.3dg \
    --output contacts.pairs \
    --distance 3
```

**Key Features**:
- KD-tree based efficient spatial querying
- Configurable distance thresholds
- Phase information preservation
- Automatic coordinate sorting

### 4. Format Conversion (`tdg2cif.py`)

**Purpose**: Converts 3D genome structure files to CIF format for structural visualization.

**Input**:
- `.3dg` file containing 3D coordinates

**Output**:
- `.cif` file compatible with molecular visualization software

**Usage**:
```bash
python -m CHARMtools.charm_preprocess.tdg2cif \
    --input structure.3dg \
    --output structure.cif
```

### 5. Additional Cleaning Tools

#### Isolated Contact Removal (`clean_isolated.py`)
**Purpose**: Removes isolated contacts that may represent noise or artifacts.

#### General Cleaning (`clean3.py`)
**Purpose**: Comprehensive cleaning pipeline combining multiple filtering steps.

#### Separation-based Cleaning (`sep_clean.py`)
**Purpose**: Filters contacts based on genomic separation criteria.

## Data Formats

### Input Formats

#### .pairs Format (4DN Standard)
```
# Standard 4DN pairs format
# columns: readID chr1 pos1 chr2 pos2 strand1 strand2 [phase0 phase1]
read1    chr1    1000000    chr1    1500000    +    -    0    1
read2    chr2    2000000    chr3    3000000    +    +    1    0
```

#### .3dg Format
```
# 3D genome structure format
# columns: chr pos x y z
chr1_pat    1000000    10.5    20.3    15.7
chr1_mat    1000000    12.1    18.9    16.2
```

#### GTF Format
```
# Gene annotation format
chr1    HAVANA    exon    1000    2000    .    +    .    gene_id "ENSG00000000001.1";
```

### Output Formats

All preprocessing tools maintain compatibility with downstream CHARMtools analysis modules by preserving standard formats and metadata.

## Best Practices

### 1. Processing Order
1. **Format Conversion**: Convert raw data to standard formats first
2. **Quality Control**: Apply basic filtering and validation
3. **Artifact Removal**: Remove promiscuous contacts and splicing artifacts
4. **Final Validation**: Verify data integrity before analysis

### 2. Parameter Selection
- **Distance Thresholds**: Adjust based on resolution and biological context
- **Count Thresholds**: Balance between noise removal and data retention
- **Thread Count**: Optimize based on available computational resources

### 3. Quality Control
- Monitor filtering statistics to ensure reasonable data retention
- Validate output formats before downstream analysis
- Keep logs of all preprocessing steps for reproducibility

## Performance Considerations

### Memory Usage
- Large datasets may require chunked processing
- Monitor memory usage during multi-threaded operations
- Consider using disk-based processing for very large files

### Computational Requirements
- Multi-threading significantly improves performance
- KD-tree operations are memory-intensive but fast
- GTF parsing and indexing require substantial initial overhead

## Troubleshooting

### Common Issues

1. **Memory Errors**
   - Reduce thread count
   - Process data in smaller chunks
   - Increase available system memory

2. **Format Errors**
   - Verify input file formats
   - Check for proper compression (gzip)
   - Validate column structures

3. **Performance Issues**
   - Optimize thread count for your system
   - Use SSD storage for temporary files
   - Monitor system resources during processing

### Getting Help

For additional support:
- Check the [Troubleshooting Guide](Troubleshooting.md)
- Review example workflows in [Quick Start](Quick-Start.md)
- Contact the development team for specific issues

## Next Steps

After preprocessing your data:
1. Load cleaned data into [Cell3D Objects](Cell3D-Objects.md)
2. Explore [Analysis Tools](Analysis-Tools.md) for downstream analysis
3. Use [Visualization](Visualization.md) tools for data exploration

---

*For more information about CHARMtools, visit the [Home](Home.md) page or check the [Installation Guide](Installation-Guide.md).*