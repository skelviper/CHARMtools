# Cell3D Output Module - Output format conversion
import warnings
import pandas as pd
import numpy as np
from io import StringIO
import sys
import time
from scipy import spatial
from utils.CHARMio import write_pairs
from utils.pairs_manipulations import sortPairs
import subprocess
import shlex
import io

class Cell3DIO:
    """Output format conversion for Cell3D objects"""
    
    def write_cif(self,factor_b=None,output_path = None):
        """
        Convert a DataFrame of 3D coordinates to a CIF file.
        Parameters:
            cellname : str
            tdg : pandas.DataFrame containing at least for 'chrom', 'pos', 'x', 'y', 'z', 'CpG'
            outputpath : str
            resolution : int
        Returns:
            None and a output file to outputpath
        """

        cellname = self.cellname
        tdg = self.tdg
        resolution = self.resolution

        if factor_b is None:
            if 'CpG' in tdg.columns:
                factor_b = 'CpG'
            else:
                factor_b = 'pos'  # Default to 'pos' if 'CpG' is not available

        if outputpath is None:
            outputpath = "./" +cellname+".cif"

        file_head_name ="data_" + cellname + "_res" + str(int(resolution/1000)) + "k"
        # Create CIF format string block 1 
        cif_str = "#\nloop_\n_entity_poly.entity_id\n_entity_poly.type\n_entity_poly.nstd_linkage\n_entity_poly.nstd_monomer\n_entity_poly.pdbx_seq_one_letter_code\n_entity_poly.pdbx_seq_one_letter_code_can\n_entity_poly.pdbx_strand_id\n_entity_poly.pdbx_target_identifier\n"
        cif_str = file_head_name + "\n" + cif_str
        # Get unique chroms
        unique_chroms = tdg['chrom'].unique()
        # Sort the array
        def sort_chromosomes(chrom):
            num_part = ''.join(filter(str.isdigit, chrom)) 
            if 'X' in chrom:
                num = 100  
            elif 'Y' in chrom:
                num = 101 
            else:
                num = int(num_part)
            return (chrom[-1], num)  
        
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
            cif_str3 += f"ATOM {i+1} C CA . GLY {row['chrom']} {entity_id} {chrom_index} ? {row['x']} {row['y']} {row['z']} {row[factor_b]}\n"

        #print(cif_str3)

        # Open the file in write mode
        with open(output_path, 'w') as f:
            # Write the three blocks to the file
            f.write(cif_str)
            f.write(cif_str2)
            f.write(cif_str3)

        print("Done " + cellname)
        return None

    def write_tdg2pairs(self, output_path, distance=2, mark_pat_mat=True):
        """
        Generate pairs from reconstructed 3D structures.
        
        Parameters:
        -----------
        output_path : str
            Path to save the output pairs file.
        distance : int, optional
            Maximum distance between pairs in base pairs. Default is 2.
        """
        if self.chrom_length is None:
                raise ValueError("Chrom length is not available, please run add_chrom_length first")

        df3d = self.get_data()[["chrom", "pos", "x", "y", "z"]]

        if self.kdtree is None:
            self.kdtree = spatial.KDTree(df3d[["x","y","z"]].values)
        kdtree = self.kdtree
        res = np.array(list(kdtree.query_pairs(r=distance)))
        if res.size == 0:
            pairs_df = pd.DataFrame(columns=["readID","chrom1","pos1","chrom2","pos2","phase0","phase1"])
        else:
            chrom = df3d["chrom"].to_numpy()
            pos = df3d["pos"].to_numpy()

            def base_and_phase(c):
                if isinstance(c, str) and len(c) > 0:
                    last = c[-1]
                    if last in ("a", "b"):
                        return c[:-1], 0 if last == "a" else 1
                return c, 0

            rows = []
            for i, j in res:
                c1, ph1 = base_and_phase(chrom[i])
                c2, ph2 = base_and_phase(chrom[j])
                rows.append([".", c1, int(pos[i]), c2, int(pos[j]), ph1, ph2])
            pairs_df = pd.DataFrame(rows, columns=["readID","chrom1","pos1","chrom2","pos2","phase0","phase1"])

        pairs_df['strand1'] = '+'
        pairs_df['strand2'] = '+'
        pairs_df = pairs_df[["readID","chrom1","pos1","chrom2","pos2","strand1","strand2","phase0","phase1"]]
        pairs_df = sortPairs(pairs_df)

        #add phasing probability
        pairs_df["phase_prob00"] = np.where((pairs_df["phase0"]==0) & (pairs_df["phase1"]==0) ,float(1),float(0))
        pairs_df["phase_prob01"] = np.where((pairs_df["phase0"]==0) & (pairs_df["phase1"]==1) ,float(1),float(0))
        pairs_df["phase_prob10"] = np.where((pairs_df["phase0"]==1) & (pairs_df["phase1"]==0) ,float(1),float(0))
        pairs_df["phase_prob11"] = np.where((pairs_df["phase0"]==1) & (pairs_df["phase1"]==1) ,float(1),float(0))

        if mark_pat_mat:
            # rename chrom
            pairs_df['chrom1'] = pairs_df['chrom1'] + np.where(pairs_df['phase0']==0, '(pat)', '(mat)')
            pairs_df['chrom2'] = pairs_df['chrom2'] + np.where(pairs_df['phase1']==0, '(pat)', '(mat)')

        # add header to pairs_df
        chrom_df = self.chrom_length.copy()
        chrom_df['chrom'] = chrom_df['chrom'].str[:-1]
        chrom_df = chrom_df.drop_duplicates().reset_index(drop=True)
        chrom_header = list('#chromosome: ' + chrom_df['chrom'].astype(str) + ' ' + chrom_df['size'].astype(str))
        header = [
            '## pairs format v1.0',
            '#sorted: chr1-chr2-pos1-pos2',
            '#shape: upper triangle'
        ] + chrom_header + ['#columns: readID\tchr1\tpos1\tchr2\tpos2\tstrand1\tstrand2\tphase0\tphase1\tphase_prob00\tphase_prob01\tphase_prob10\tphase_prob11']

        if isinstance(output_path, io.StringIO):  # if output_path is a StringIO buffer
            output_path.write("\n".join(header) + "\n")
            pairs_df.to_csv(output_path, index=False, header=False, sep="\t")
        else:
            if output_path.endswith('.gz'):
                import gzip
                with gzip.open(output_path, "wt") as f:
                    for line in header:
                        f.write(line + "\n")
                    pairs_df.to_csv(f, index=False, header=False, sep="\t")
            else:
                with open(output_path, "w") as f:
                    for line in header:
                        f.write(line + "\n")
                    pairs_df.to_csv(f, index=False, header=False, sep="\t")

        # write_pairs(pairs_df, output_path)

        return None

    
    def write_csv(self, output_path=None, cellname=None, tdg=None, include_features=True):
        """Write Cell3D data to CSV format"""
        if cellname is None:
            cellname = self.cellname
        if tdg is None:
            if self.on_disk:
                self.to_memory()
            tdg = self.tdg
        
        # Sort by chromosome and position
        tdg_sorted = tdg.sort_values(['chrom', 'pos']).reset_index(drop=True)
        
        if not include_features:
            # Only include basic columns
            columns_to_include = ['chrom', 'pos', 'x', 'y', 'z']
            tdg_output = tdg_sorted[columns_to_include]
        else:
            tdg_output = tdg_sorted
        
        if output_path:
            tdg_output.to_csv(output_path, index=False)
            return None
        else:
            return tdg_output.to_csv(index=False)
    
    def export_summary(self, output_path=None):
        """Export summary statistics of the Cell3D object"""
        if self.on_disk:
            self.to_memory()
        
        summary = {
            'cellname': self.cellname,
            'resolution': self.resolution,
            'n_points': len(self.tdg),
            'chromosomes': list(self.tdg['chrom'].unique()),
            'n_chromosomes': self.tdg['chrom'].nunique(),
            'features': self.features,
            'n_features': len(self.features),
            'coordinate_ranges': {
                'x_min': self.tdg['x'].min(),
                'x_max': self.tdg['x'].max(),
                'y_min': self.tdg['y'].min(),
                'y_max': self.tdg['y'].max(),
                'z_min': self.tdg['z'].min(),
                'z_max': self.tdg['z'].max()
            }
        }
        
        if hasattr(self, 'metadata') and self.metadata:
            summary['metadata'] = self.metadata
        
        if output_path:
            import json
            with open(output_path, 'w') as f:
                json.dump(summary, f, indent=2, default=str)
            return None
        else:
            return summary


    ## untested code below
    def write_pdb(self, output_path=None, cellname=None, tdg=None):
        """Write Cell3D data to PDB format"""
        if cellname is None:
            cellname = self.cellname
        if tdg is None:
            if self.on_disk:
                self.to_memory()
            tdg = self.tdg
        
        # Sort by chromosome and position
        tdg_sorted = tdg.sort_values(['chrom', 'pos']).reset_index(drop=True)
        
        # Create PDB content
        pdb_content = StringIO()
        
        # Write header
        pdb_content.write(f"HEADER    CHROMATIN STRUCTURE                    {cellname}\n")
        pdb_content.write(f"TITLE     3D CHROMATIN STRUCTURE FOR CELL {cellname}\n")
        pdb_content.write("REMARK   1 GENERATED BY CHARM TOOLS\n")
        
        # Write atom records
        for idx, row in tdg_sorted.iterrows():
            chrom = row['chrom']
            pos = row['pos']
            x, y, z = row['x'], row['y'], row['z']
            
            # Extract chromosome information
            if chrom.endswith('a') or chrom.endswith('b'):
                chrom_base = chrom[:-1]
                allele = chrom[-1]
            else:
                chrom_base = chrom
                allele = 'A'
            
            # Calculate residue number
            res_num = int(pos // self.resolution) + 1
            
            # Write ATOM record
            pdb_content.write(f"ATOM  {idx+1:>5}  CA  BIN {allele.upper()}{res_num:>4}    ")
            pdb_content.write(f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00 20.00           C  \n")
        
        pdb_content.write("END\n")
        
        # Get the final content
        pdb_string = pdb_content.getvalue()
        pdb_content.close()
        
        if output_path:
            with open(output_path, 'w') as f:
                f.write(pdb_string)
            return None
        else:
            return pdb_string
    
    def write_xyz(self, output_path=None, cellname=None, tdg=None):
        """Write Cell3D data to XYZ format"""
        if cellname is None:
            cellname = self.cellname
        if tdg is None:
            if self.on_disk:
                self.to_memory()
            tdg = self.tdg
        
        # Sort by chromosome and position
        tdg_sorted = tdg.sort_values(['chrom', 'pos']).reset_index(drop=True)
        
        # Create XYZ content
        xyz_content = StringIO()
        
        # Write number of atoms
        xyz_content.write(f"{len(tdg_sorted)}\n")
        
        # Write comment line
        xyz_content.write(f"3D chromatin structure for cell {cellname}\n")
        
        # Write coordinates
        for idx, row in tdg_sorted.iterrows():
            chrom = row['chrom']
            x, y, z = row['x'], row['y'], row['z']
            
            # Use chromosome as atom type
            atom_type = chrom.replace('chr', 'C')
            
            xyz_content.write(f"{atom_type:<4} {x:>12.6f} {y:>12.6f} {z:>12.6f}\n")
        
        # Get the final content
        xyz_string = xyz_content.getvalue()
        xyz_content.close()
        
        if output_path:
            with open(output_path, 'w') as f:
                f.write(xyz_string)
            return None
        else:
            return xyz_string
    
    def reconstruct_3dg(self,hickit_command,distance=2):
        """
        Reconstruct 3D genome from pairs under a distance threshold.
        hickit_command: str, command line arguments for hickit, do not include "-i" and input file.
        distance: float, distance threshold to define pairs.
        clean3: bool, whether to clean up incorrect structure.
        """
        # Create an in-memory buffer (instead of writing to a file)
        pairs_buffer = io.StringIO()
        self.write_tdg2pairs(output_path=pairs_buffer, distance=distance, mark_pat_mat=False)
        pairs_buffer.seek(0)

        # run hickit command
        command = ["/shared/jliu/unnamed/structureReconGPU/hickit/hickit", 
               "-i", "-"]  # "-" indicates reading from stdin
        parsed_command = shlex.split(hickit_command)
        command.extend(parsed_command)

        try:
            subprocess.run(command, input=pairs_buffer.read(), text=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error occurred while running hickit: {e}")

        pairs_buffer.close()
        return None
        
    
    def clean3(self,tdg_file_path, pair_path, clean_3dg_path,distance=2):
        """
        Clean up incorrect 3D structure based on distance consistency
        """
        self.write_tdg2pairs(output_path=pair_path, distance=distance, mark_pat_mat=True)
        
        import subprocess
        command = [
            'python', '/shared/jliu/unnamed/github/CHARMtools/CHARMtools.py', 'clean3',
            '-r', pair_path, '-i', tdg_file_path, '-o', clean_3dg_path
        ]
        
        try:
            subprocess.run(command, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error occurred while running clean3: {e}")
        
        return None