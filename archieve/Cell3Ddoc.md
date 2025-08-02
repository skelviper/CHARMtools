# Cell3D class documentation

## Overview
The `Cell3D` and `MultiCell3D` class is intended for handling and analyzing three-dimensional genome structure cellular data. 

## Classes and Functions

### Class: `Cell3D`
- **Initialization (`__init__`):**
  - Parameters:
    - `cellname`: Name of the cell.
    - `tdg_path`: Path to the 'tdg' file containing cellular data.
    - `resolution`: Resolution of the data.
  - Description: Initializes the `Cell3D` object with basic properties, loads the 'tdg' data, and sets up a k-d tree for spatial data management.

- **Method: `add_density_in_radius`**
  - Parameters:
    - `radius`: The radius within which to calculate density.
  - Description: Calculates the density of points within a specified radius for each point in the cellular data. It adds this density as a new feature to the data.

- **Method: `add_feature_in_radius`**
  - Parameters:
    - `feature`: The feature name to be smoothed.
    - `radius`: Radius of the sphere for smoothing.
    - `type`: Either "mean" or "sum" to define the smoothing type.
  - Description: Smooths a given feature over a sphere of a specified radius by averaging or summing the feature values.

- **Method: `add_chrom_length`**
  - Parameters:
    - `chrom_length_path`: Path to the file containing chromosome lengths.
  - Description: Adds chromosome length data to the `Cell3D` object from a given file.

- **Method: `calc_distance_matrix`**
  - Parameters:
    - `chrom`: Chromosome name.
  - Description: Calculates the distance matrix for a given chromosome. It requires chromosome length data to be available.

- **Method: `plot_feature_distribution`**
  - Parameters:
    - `feature`: Feature to be plotted.
    - `bins`: Number of bins for the histogram.
  - Description: Plots the distribution of a specified feature using a histogram.

- **Method: `calc_feature_distances`**
  - Parameters:
    - `feature`: The feature based on which distances are calculated.
    - `radius`: The radius within which to consider the points.
    - `points_num_per_chrom`: Number of points per chromosome for calculations.
    - `points_num_other_chrom`: Number of points from other chromosomes.
    - `random_seed`: Seed for random number generation.
  - Description: Calculates distances based on a specific feature within a given radius. It includes statistical tests to compare distances of feature points against random points.

- **Method: `calc_distance_matrix`**
  - Description: This method is designed to return a distance matrix of a specific region or band within a chromosome. It's important for analyzing spatial relationships within cellular structures.


### Class: `MultiCell3D`

- **Initialization (`__init__`):**
  - Parameters:
    - `cells`: List of `Cell3D` objects.

get_info
get_cell
get_data
subset
calc_distance_matrix
calc_3dproximity_matrix

calc_scABC_pred_gene
calc_feature_matrix
calc_feature_proximity_matrix
calc_radial_position_matrix
get_feature_vec
FindMarkers
zoomify_matrix
simple_diff

