# Convert 3dg file to cif file using Cell3D object
# @Date 2024
# @Author Updated implementation based on Cell3D

import sys
import pandas as pd
import argparse
from utils.CHARMio import parse_3dg
from obj.Cell3D import Cell3D

def cli(args):
    """
    CLI function for tdg2cif command
    """
    input_file, output_file, factorBpath, maxGap, expansion = \
        args.input_file, args.output_file, args.factorBpath, args.maxGap, args.expansion
    
    print("input_file: ", input_file)
    print("output_file: ", output_file)
    
    # Convert using Cell3D object
    tdg_to_cif_cell3d(input_file, output_file, factorBpath, maxGap, expansion)

def tdg_to_cif_cell3d(tdg_path: str, output_cif_path: str, factorBpath: str = None, maxGap: int = 1000000, expansion: float = 1.0):
    """
    Convert 3dg file to cif format using Cell3D object
    
    Parameters:
    -----------
    tdg_path : str
        Path to input 3dg file
    output_cif_path : str
        Path to output cif file
    factorBpath : str, optional
        Path to factor B file
    maxGap : int, optional
        Maximum gap between bins (not used in Cell3D implementation)
    expansion : float, optional
        Expansion factor for coordinates
    """
    
    # Load 3dg data
    tdg_data = parse_3dg(tdg_path)
    
    # Apply expansion factor if specified
    if expansion != 1.0:
        tdg_data[['x', 'y', 'z']] = tdg_data[['x', 'y', 'z']] * expansion
    
    # Load factor B data if provided
    if factorBpath:
        factor_b_data = pd.read_csv(factorBpath, sep='\t', header=None)
        # Assume factor B file has same order as 3dg file
        if len(factor_b_data) == len(tdg_data):
            tdg_data['factor_b'] = factor_b_data.iloc[:, -1]  # Use last column as factor B
        else:
            print("Warning: Factor B file length doesn't match 3dg file. Using default values.")
            tdg_data['factor_b'] = 50.0  # Default factor B value
    else:
        tdg_data['factor_b'] = 50.0  # Default factor B value
    
    # Determine resolution from data (estimate from position differences)
    if len(tdg_data) > 1:
        # Get resolution from first chromosome
        first_chrom_data = tdg_data[tdg_data['chrom'] == tdg_data['chrom'].iloc[0]]
        if len(first_chrom_data) > 1:
            resolution = int(first_chrom_data['pos'].iloc[1] - first_chrom_data['pos'].iloc[0])
        else:
            resolution = 40000  # Default resolution
    else:
        resolution = 40000  # Default resolution
    
    # Create Cell3D object
    # Extract cellname from input file path
    cellname = tdg_path.split('/')[-1].split('.')[0]
    
    # Create Cell3D object with the loaded data
    cell3d = Cell3D(cellname=cellname, resolution=resolution)
    cell3d.tdg = tdg_data
    
    # Write to CIF format using Cell3D's write_cif method
    cell3d.write_cif(factor_b='factor_b', outputpath=output_cif_path)
    
    print(f"Successfully converted {tdg_path} to {output_cif_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert 3dg file to cif file using Cell3D object')
    parser.add_argument('-i', '--input', dest='input_file', required=True, help='Input 3dg file path')
    parser.add_argument('-o', '--output', dest='output_file', required=True, help='Output cif file path')
    parser.add_argument('-f', '--factorB', dest='factorBpath', help='Factor B file path')
    parser.add_argument('-g', '--maxGap', dest='maxGap', type=int, default=1000000, help='Max gap between bins')
    parser.add_argument('-e', '--expansion', dest='expansion', type=float, default=1.0, help='Expansion factor')
    
    args = parser.parse_args()
    cli(args)