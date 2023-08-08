import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
import glob
import tqdm

from mpl_toolkits import mplot3d
from sklearn.preprocessing import LabelEncoder
import argparse


def charm_get_3dplot_data(CpG_path,tdg_path,atac_path=None,ct_path=None,resolution=1000000):
    """
    load CHARM 3dg data and merge with CpG, ATAC, and CT data
    Required:
        CpG_path: path to CpG file
        tdg_path: path to tdg file
    Optional:
        atac_path: path to ATAC file
        ct_path: path to CT file
        resolution: resolution of ATAC and CT data, default is 1Mb
    Return:
        tdg: pd.DataFrame containing 3dg data merged with CpG, ATAC, and CT data
    """
    # load cpg
    CpG = pd.read_csv(CpG_path,header=None,sep="\t")
    CpG.columns = ["chrom","pos","CpG"]
    # if chrom column not start with "chr", add "chr"
    if CpG["chrom"].str.startswith("chr").sum() == 0:
        CpG["chrom"] = "chr"+CpG["chrom"]
    CpG["pos"] = CpG["pos"].astype(int)

    CpG_by_allele = pd.concat([
        CpG.assign(chrom = lambda x: x["chrom"]+"a"),
        CpG.assign(chrom = lambda x: x["chrom"]+"b")
    ])
    # load tdg
    tdg = pd.read_csv(tdg_path,sep="\t",header=None,comment="#")
    tdg.columns = ["chrom","pos","x","y","z"]
    # if chrom column contains "(mat)", replace with "a"
    tdg["chrom"] = tdg["chrom"].str.replace("\(mat\)","a",regex=True)
    # if chrom column contains "(pat)", replace with "b"
    tdg["chrom"] = tdg["chrom"].str.replace("\(pat\)","b",regex=True)
    # merge with CpG
    tdg = pd.merge(tdg,CpG_by_allele,on=["chrom","pos"],how="left").dropna()
    tdg.chrom = pd.Categorical(tdg.chrom)
    LE = LabelEncoder()
    tdg['chrom_code'] = LE.fit_transform(tdg['chrom'])
    # merge with atac
    if atac_path is not None:
        fragments = pd.read_csv(atac_path,sep="\t",header=None)
        fragments.columns = ["chrom","start","end","allele","score","strand"]
        fragments = fragments.query("chrom.str.contains('chr')").query('allele != "."')
        fragments = fragments.assign(chrom = np.where(fragments["allele"] == "0",fragments["chrom"]+"a",fragments["chrom"]+"b"))
        fragments["pos"] = ((fragments["start"] + fragments["end"])/2 + (resolution / 2)) // resolution * resolution
        fragments["pos"] = fragments["pos"].astype(int)
        fragments = fragments.groupby(["chrom","pos"]).size().reset_index().rename(columns={0:"count"})
        tdg = pd.merge(tdg,fragments,on=["chrom","pos"],how="left")
        tdg["count_atac"] = tdg["count"].fillna(0)
        tdg = tdg.drop(columns=["count"],axis=1)
    # merge with ct
    if ct_path is not None:
        fragments = pd.read_csv(ct_path,sep="\t",header=None)
        fragments.columns = ["chrom","start","end","allele","score","strand"]
        fragments = fragments.query("chrom.str.contains('chr')").query('allele != "."')
        fragments = fragments.assign(chrom = np.where(fragments["allele"] == "0",fragments["chrom"]+"a",fragments["chrom"]+"b"))
        fragments["pos"] = ((fragments["start"] + fragments["end"])/2 + (resolution / 2)) // resolution * resolution
        fragments["pos"] = fragments["pos"].astype(int)
        fragments = fragments.groupby(["chrom","pos"]).size().reset_index().rename(columns={0:"count"})
        tdg = pd.merge(tdg,fragments,on=["chrom","pos"],how="left")
        tdg["count_ct"] = tdg["count"].fillna(0)
        tdg = tdg.drop(columns=["count"],axis=1)
    return tdg

def sort_chromosomes(chrom):
    num_part = ''.join(filter(str.isdigit, chrom)) 
    if 'X' in chrom:
        num = 100  
    elif 'Y' in chrom:
        num = 101 
    else:
        num = int(num_part)
    return (chrom[-1], num)  
# Return a tuple where the first element is 'a' or 'b' and the second element is the number

def tdg2cif(cellname,tdg,outputpath = None ,resolution=200000):
    """Convert a DataFrame of 3D coordinates to a CIF file."""

    if outputpath is None:
        outputpath = "./" +cellname+".cif"

    file_head_name ="data_" + cellname + "_res" + str(int(resolution/1000)) + "k"
    # Create CIF format string block 1 
    cif_str = "#\nloop_\n_entity_poly.entity_id\n_entity_poly.type\n_entity_poly.nstd_linkage\n_entity_poly.nstd_monomer\n_entity_poly.pdbx_seq_one_letter_code\n_entity_poly.pdbx_seq_one_letter_code_can\n_entity_poly.pdbx_strand_id\n_entity_poly.pdbx_target_identifier\n"
    cif_str = file_head_name + "\n" + cif_str
    # Get unique chroms
    unique_chroms = tdg['chrom'].unique()
    # Sort the array
    unique_chroms = sorted(unique_chroms, key=sort_chromosomes)
    tdg["chrom"] = pd.Categorical(tdg["chrom"], categories=unique_chroms)
    # Add each chrom as a new line in the CIF block
    for i, chrom in enumerate(unique_chroms, start=1):
        cif_str += f"{i} 'Chromatin' no no ? ? ? ?\n"


    # Create CIF format string for the second block
    cif_str2 = "#\nloop_\n_entity.id\n_entity.type\n_entity.src_method\n_entity.pdbx_description\n_entity.formula_weight\n_entity.pdbx_number_of_molecules\n_entity.pdbx_ec\n_entity.pdbx_mutation\n_entity.pdbx_fragment\n_entity.details\n"

    # Add each chrom as a new line in the CIF block
    for i, chrom in enumerate(unique_chroms, start=1):
        # Determine if the chrom is maternal or paternal
        chrom_type = 'maternal' if 'a' in chrom else 'paternal'
        # Get the chrom number
        chrom_num = ''.join(filter(str.isdigit, chrom))
        if 'X' in chrom:
            chrom_num = 'X'
        elif 'Y' in chrom:
            chrom_num = 'Y'
        cif_str2 += f"{i} polymer man 'Chromosome{chrom_num} ({chrom_type})' ? ? ? ? ? ?\n"

    #print(cif_str2)

    # Continue from previous code

    # Create CIF format string for the third block
    cif_str3 = "#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.B_iso_or_equiv\n"

    # Create a dictionary to store the current index for each chrom
    chrom_indices = {chrom: 0 for chrom in unique_chroms}

    # Add each row of the DataFrame as a new line in the CIF block
    for i, row in tdg.iterrows():
        # Get the current index for this chrom
        chrom_index = chrom_indices[row['chrom']] + 1
        # Update the index for this chrom
        chrom_indices[row['chrom']] = chrom_index
        # Get the entity id for this chrom
        entity_id = unique_chroms.index(row['chrom']) + 1
        cif_str3 += f"ATOM {i+1} C CA . GLY {row['chrom']} {entity_id} {chrom_index} ? {row['x']} {row['y']} {row['z']} {row['CpG']}\n"

    #print(cif_str3)

    # Open the file in write mode
    with open(outputpath, 'w') as f:
        # Write the three blocks to the file
        f.write(cif_str)
        f.write(cif_str2)
        f.write(cif_str3)

    print("Done " + cellname)
    return None

# Example api usage

# import xxx
# resolution = 200000
# cellname = "d4B78"

# tdg = charm_get_3dplot_data(#CpG_path = "/share/Data/public/ref_genome/human_ref/GRCh37d5/raw_data/cpg/hg19.cpg.1m.txt",
#                             CpG_path = "/share/home/zliu/share/Data/public/ref_genome/mouse_ref/M23/CpG/dipc_cpg/mm10.CpG.200000.txt",
#                             tdg_path= "/shareb/zliu/analysis/hires_xci/esc_qualify/3d_info/dip_3dg/d4B78/clean.200k.0.3dg",
#                             resolution=resolution)
# tdg2cif(cellname,tdg,outputpath = "./"+cellname+".cif",resolution=200000)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert 3d coordinates to CIF file')
    parser.add_argument('--cellname', type=str, help='cellname',default="Cell001")
    parser.add_argument('--tdg_path', type=str, help='path to 3d coordinates',required=True)
    parser.add_argument('--output_path', type=str, help='path to output cif file',required=False)
    parser.add_argument('--resolution', type=int, help='resolution of 3d coordinates',default=1000000)
    parser.add_argument('--cpg_path', type=str, help='path to CpG file',required=True)
    parser.add_argument('--atac_path', type=str, help='path to ATAC file',required=False)
    parser.add_argument('--ct_path', type=str, help='path to CT file',required=False)               

    args = parser.parse_args()
    tdg = charm_get_3dplot_data(CpG_path = args.cpg_path,
                                tdg_path= args.tdg_path,
                                atac_path=args.atac_path,
                                ct_path=args.ct_path,
                                resolution=args.resolution)
    tdg2cif(args.cellname,tdg,outputpath = args.output_path ,resolution=args.resolution)

# Example cli usage
#python tdg2cif_v2.py --cellname d4B22 --tdg_path |\
# /shareb/zliu/analysis/hires_xci/esc_qualify/3d_info/dip_3dg/d4B22/clean.200k.0.3dg |\
# --cpg_path /share/home/zliu/share/Data/public/ref_genome/mouse_ref/M23/CpG/dipc_cpg/mm10.CpG.200000.txt

# visualize with UCSF chimeraX
# 1. open the cif file in chimeraX
# 2. Click Graphics-Lighting & Effects-Soft
# 3. type ribbon in the command line
# 4. type color bychain in the command line
