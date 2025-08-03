# Command Line Interface (CLI)

CHARMtools provides a comprehensive command-line interface for data preprocessing and analysis tasks. The CLI is designed for batch processing, pipeline integration, and automated workflows.

## Overview

The CHARMtools CLI provides access to:
- Data preprocessing tools
- Format conversion utilities
- Quality control functions
- Batch processing capabilities

## Basic Usage

```bash
# General syntax
charmtools <command> [options] <input_files>

# Get help
charmtools --help
charmtools <command> --help
```

## Available Commands

### Data Cleaning Commands

#### 1. clean_leg - Remove Promiscuous Contacts

**Purpose**: Remove promiscuous contacts that interact with multiple genomic loci using Tan's phantom leg method.

**Syntax**:
```bash
charmtools clean_leg [options] INPUT_FILE
```

**Options**:
- `-t, --thread <N>`: Number of threads (default: 4)
- `-d, --distance <DIST>`: Maximum distance to calculate adjacent legs (default: 1000)
- `-n, --count <COUNT>`: Number threshold of adjacent legs (default: 15)
- `-o, --output <FILE>`: Output file name (required)

**Example**:
```bash
charmtools clean_leg \
    -t 8 \
    -d 1000 \
    -n 15 \
    -o cleaned_contacts.pairs.gz \
    raw_contacts.pairs.gz
```

**Input**: `.pairs` file in 4DN standard format
**Output**: Cleaned `.pairs` file with promiscuous contacts removed

#### 2. clean_splicing - Remove Splicing Artifacts

**Purpose**: Remove contacts that may result from RNA splicing artifacts by filtering contacts within exonic regions.

**Syntax**:
```bash
charmtools clean_splicing [options] INPUT_FILE
```

**Options**:
- `-r, --reference <GTF>`: Annotation GTF file (required)
- `-t, --thread <N>`: Number of threads (default: 4)
- `-o, --output <FILE>`: Output file name (required)

**Example**:
```bash
charmtools clean_splicing \
    -r gencode.v38.annotation.gtf \
    -t 8 \
    -o splicing_cleaned.pairs.gz \
    input_contacts.pairs.gz
```

**Input**: `.pairs` file and GTF annotation file
**Output**: Cleaned `.pairs` file with splicing artifacts removed

#### 3. clean3 - Clean 3D Structure Particles

**Purpose**: Remove 3D genome particles that are poorly supported by contact data.

**Syntax**:
```bash
charmtools clean3 [options]
```

**Options**:
- `-i, --input <FILE>`: Structure file to clean (.3dg/.xyz format) (required)
- `-r, --reference <FILE>`: Contact file (.pairs format) (required)
- `-o, --output <FILE>`: Output cleaned structure file (required)
- `-q, --quantile <Q>`: Quantile of particles to remove (default: 0.06)
- `-d, --distance <DIST>`: Max distance (bp) from contact leg to 3D particle (default: 500000)

**Example**:
```bash
charmtools clean3 \
    -i structure.3dg \
    -r contacts.pairs.gz \
    -o cleaned_structure.3dg \
    -q 0.05 \
    -d 500000
```

**Input**: 3D structure file and contact pairs file
**Output**: Cleaned 3D structure file

#### 4. clean_isolated - Remove Isolated Contacts

**Purpose**: Remove isolated contacts according to L-0.5 distance criteria.

**Syntax**:
```bash
charmtools clean_isolated [options] INPUT_FILE
```

**Options**:
- `-t, --thread <N>`: Number of threads (default: 4)
- `-o, --output <FILE>`: Output file name (required)

**Example**:
```bash
charmtools clean_isolated \
    -t 6 \
    -o isolated_cleaned.pairs.gz \
    input_contacts.pairs.gz
```

### Format Conversion Commands

#### 5. tdg2pairs - Convert 3D Structure to Contact Pairs

**Purpose**: Convert 3D genome structure files (.3dg) to contact pairs format (.pairs) using spatial proximity.

**Syntax**:
```bash
charmtools tdg2pairs [options]
```

**Options**:
- `-i, --input <FILE>`: Input .3dg file (required)
- `-o, --output <FILE>`: Output .pairs file (required)
- `-d, --distance <DIST>`: Distance threshold for contact definition (default: 3)

**Example**:
```bash
charmtools tdg2pairs \
    -i structure.3dg \
    -o derived_contacts.pairs \
    -d 3
```

**Input**: `.3dg` file containing 3D coordinates
**Output**: `.pairs` file with contacts derived from 3D proximity

#### 6. tdg2cif - Convert 3D Structure to CIF Format

**Purpose**: Convert 3D genome structure files to CIF format for structural visualization.

**Syntax**:
```bash
charmtools tdg2cif [options]
```

**Options**:
- `-i, --input <FILE>`: Input .3dg file (required)
- `-o, --output <FILE>`: Output .cif file (required)

**Example**:
```bash
charmtools tdg2cif \
    -i structure.3dg \
    -o structure.cif
```

**Input**: `.3dg` file containing 3D coordinates
**Output**: `.cif` file compatible with molecular visualization software

#### 7. sep_clean - Separation-based Cleaning

**Purpose**: Filter contacts based on genomic separation criteria.

**Syntax**:
```bash
charmtools sep_clean [options] INPUT_FILE
```

**Options**:
- `-s, --separation <DIST>`: Minimum genomic separation (default: 1000)
- `-o, --output <FILE>`: Output file name (required)

**Example**:
```bash
charmtools sep_clean \
    -s 5000 \
    -o separation_filtered.pairs.gz \
    input_contacts.pairs.gz
```

## Workflow Examples

### Complete Preprocessing Pipeline

```bash
#!/bin/bash
# Complete preprocessing workflow

# Step 1: Remove promiscuous contacts
charmtools clean_leg \
    -t 8 -d 1000 -n 15 \
    -o step1_clean_leg.pairs.gz \
    raw_contacts.pairs.gz

# Step 2: Remove splicing artifacts
charmtools clean_splicing \
    -r gencode.v38.annotation.gtf \
    -t 8 \
    -o step2_clean_splicing.pairs.gz \
    step1_clean_leg.pairs.gz

# Step 3: Remove isolated contacts
charmtools clean_isolated \
    -t 8 \
    -o step3_clean_isolated.pairs.gz \
    step2_clean_splicing.pairs.gz

# Step 4: Apply separation filter
charmtools sep_clean \
    -s 1000 \
    -o final_cleaned.pairs.gz \
    step3_clean_isolated.pairs.gz

echo "Preprocessing complete: final_cleaned.pairs.gz"
```

### Batch Processing Multiple Files

```bash
#!/bin/bash
# Process multiple files in batch

for file in *.pairs.gz; do
    base=$(basename "$file" .pairs.gz)
    echo "Processing $file..."
    
    charmtools clean_leg \
        -t 4 -d 1000 -n 15 \
        -o "${base}_cleaned.pairs.gz" \
        "$file"
done

echo "Batch processing complete"
```

### Structure Analysis Pipeline

```bash
#!/bin/bash
# Structure analysis workflow

# Convert 3D structure to contacts
charmtools tdg2pairs \
    -i cell_structure.3dg \
    -o derived_contacts.pairs \
    -d 3

# Clean the derived contacts
charmtools clean_leg \
    -t 8 -d 1000 -n 15 \
    -o cleaned_derived.pairs.gz \
    derived_contacts.pairs

# Convert structure for visualization
charmtools tdg2cif \
    -i cell_structure.3dg \
    -o structure_for_viz.cif

echo "Structure analysis pipeline complete"
```

## Integration with Other Tools

### Using with Make

```makefile
# Makefile for CHARMtools preprocessing

# Variables
THREADS = 8
DISTANCE = 1000
COUNT = 15
GTF = gencode.v38.annotation.gtf

# Rules
%.clean_leg.pairs.gz: %.pairs.gz
	charmtools clean_leg -t $(THREADS) -d $(DISTANCE) -n $(COUNT) -o $@ $<

%.clean_splicing.pairs.gz: %.clean_leg.pairs.gz
	charmtools clean_splicing -r $(GTF) -t $(THREADS) -o $@ $<

%.final.pairs.gz: %.clean_splicing.pairs.gz
	charmtools clean_isolated -t $(THREADS) -o $@ $<

# Default target
all: sample1.final.pairs.gz sample2.final.pairs.gz

clean:
	rm -f *.clean_*.pairs.gz *.final.pairs.gz
```

### Using with Snakemake

```python
# Snakefile for CHARMtools preprocessing

SAMPLES = ["sample1", "sample2", "sample3"]
GTF = "gencode.v38.annotation.gtf"

rule all:
    input:
        expand("{sample}.final.pairs.gz", sample=SAMPLES)

rule clean_leg:
    input:
        "{sample}.pairs.gz"
    output:
        "{sample}.clean_leg.pairs.gz"
    threads: 8
    shell:
        "charmtools clean_leg -t {threads} -d 1000 -n 15 -o {output} {input}"

rule clean_splicing:
    input:
        pairs="{sample}.clean_leg.pairs.gz",
        gtf=GTF
    output:
        "{sample}.clean_splicing.pairs.gz"
    threads: 8
    shell:
        "charmtools clean_splicing -r {input.gtf} -t {threads} -o {output} {input.pairs}"

rule clean_isolated:
    input:
        "{sample}.clean_splicing.pairs.gz"
    output:
        "{sample}.final.pairs.gz"
    threads: 8
    shell:
        "charmtools clean_isolated -t {threads} -o {output} {input}"
```

## Performance Optimization

### Thread Usage

```bash
# Optimize thread usage based on system
NTHREADS=$(nproc)  # Use all available cores
# or
NTHREADS=$(($(nproc) - 2))  # Leave 2 cores free

charmtools clean_leg -t $NTHREADS -o output.pairs.gz input.pairs.gz
```

### Memory Management

```bash
# For large files, process in chunks or use streaming
# Monitor memory usage
charmtools clean_leg \
    -t 4 \
    -o large_output.pairs.gz \
    large_input.pairs.gz &

# Monitor process
watch -n 5 'ps aux | grep charmtools'
```

## Error Handling and Logging

### Basic Error Handling

```bash
#!/bin/bash
set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

# Function for error handling
handle_error() {
    echo "Error occurred in script at line: $1"
    exit 1
}
trap 'handle_error $LINENO' ERR

# Run CHARMtools commands
charmtools clean_leg -t 8 -o output.pairs.gz input.pairs.gz
echo "Processing completed successfully"
```

### Logging

```bash
#!/bin/bash
# Create log file
LOGFILE="charmtools_$(date +%Y%m%d_%H%M%S).log"

# Function to log messages
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOGFILE"
}

log_message "Starting CHARMtools preprocessing"

# Run commands with logging
charmtools clean_leg -t 8 -o output.pairs.gz input.pairs.gz 2>&1 | tee -a "$LOGFILE"

log_message "Processing completed"
```

## Troubleshooting CLI Issues

### Common Problems

1. **Command not found**:
```bash
# Check installation
which charmtools
# or
charmtools --version
```

2. **Permission errors**:
```bash
# Check file permissions
ls -la input_file.pairs.gz
# Fix permissions if needed
chmod +r input_file.pairs.gz
```

3. **Memory errors**:
```bash
# Reduce thread count
charmtools clean_leg -t 2 -o output.pairs.gz input.pairs.gz
```

4. **File format errors**:
```bash
# Validate input file format
zcat input.pairs.gz | head -n 10
# Check for proper pairs format
```

### Getting Help

```bash
# General help
charmtools --help

# Command-specific help
charmtools clean_leg --help
charmtools clean_splicing --help

# Version information
charmtools --version
```

## Best Practices

1. **File Organization**:
   - Use descriptive filenames
   - Keep raw and processed data separate
   - Use version control for scripts

2. **Resource Management**:
   - Monitor system resources during processing
   - Use appropriate thread counts
   - Process large datasets in chunks

3. **Quality Control**:
   - Validate input files before processing
   - Check output file integrity
   - Keep processing logs

4. **Reproducibility**:
   - Document all parameters used
   - Use version-controlled scripts
   - Keep track of software versions

## Next Steps

After using the CLI tools:
1. Load processed data into [Cell3D Objects](Cell3D-Objects.md)
2. Perform analysis using [Analysis Tools](Analysis-Tools.md)
3. Create visualizations with [Visualization](Visualization.md) tools
4. Use the [API](API_DOCUMENTATION.md) for programmatic access

---

*For more information about CHARMtools, visit the [Home](Home.md) page or check the [Installation Guide](Installation-Guide.md).*