import pandas as pd
import numpy as np
from sklearn.preprocessing import LabelEncoder
import warnings
from scipy.spatial import cKDTree
import rmsd


class Cell3D:
    def __init__(self, cellname,tdg_path, resolution):
        self.cellname = cellname
        self.tdg = self._load_tdg(tdg_path)
        self.resolution = resolution
        self.features = []
        self.kdtree = cKDTree(self.tdg[["x", "y", "z"]].values)

    def __repr__(self):
        self.get_info()
        return ""

    def get_info(self):
        print("CHARMtools Cell3D object v0.1")
        print("Cell name: {}".format(self.cellname))
        print("Resolution: {}".format(self.resolution))
        print("Features: {}".format(self.features))
        print(f"<Chromatin3DStructure with {len(self.tdg)} records>")

    # Construct a Cell3D object
    def _load_tdg(self, tdg_path):
        tdg = pd.read_csv(tdg_path, sep="\t", header=None, comment="#")
        tdg.columns = ["chrom", "pos", "x", "y", "z"]
        tdg["chrom"] = tdg["chrom"].str.replace("\(pat\)", "a", regex=True)
        tdg["chrom"] = tdg["chrom"].str.replace("\(mat\)", "b", regex=True)
        LE = LabelEncoder()
        tdg.chrom = pd.Categorical(tdg.chrom)
        tdg['chrom_code'] = LE.fit_transform(tdg['chrom'])
        return tdg

    def _load_bed_fragments(path, resolution, type = "allelic_resolved"):
        fragments = pd.read_csv(path, sep="\t", header=None)
        fragments.columns = ["chrom", "start", "end", "allele", "score", "strand"][:len(fragments.columns)]
        if type == "allelic_resolved":
            fragments = fragments.query("chrom.str.contains('chr')").query('allele != "."')
            fragments = fragments.assign(chrom=np.where(fragments["allele"] == "0", fragments["chrom"] + "a", fragments["chrom"] + "b"))
        else:
            fragments = pd.concat([
                fragments.assign(chrom=lambda x: x["chrom"] + "a"),
                fragments.assign(chrom=lambda x: x["chrom"] + "b")
            ])
        fragments["pos"] = ((fragments["start"] + fragments["end"]) / 2 + (resolution / 2)) // resolution * resolution
        fragments["pos"] = fragments["pos"].astype(int)
        return fragments.groupby(["chrom", "pos"]).size().reset_index().rename(columns={0: "count"})

    def _load_CpG(CpG_path):
        CpG = pd.read_csv(CpG_path, header=None, sep="\t")
        CpG.columns = ["chrom", "pos", "CpG"]
        if CpG["chrom"].str.startswith("chr").sum() == 0:
            CpG["chrom"] = "chr" + CpG["chrom"]
        CpG["pos"] = CpG["pos"].astype(int)
        return pd.concat([
            CpG.assign(chrom=lambda x: x["chrom"] + "a"),
            CpG.assign(chrom=lambda x: x["chrom"] + "b")
        ])    

    def add_bed_data(self, path, column_name, resolution=None,type="allelic_resolved"):
        if resolution is None:
            resolution = self.resolution
        if column_name in self.tdg.columns:
            warnings.warn("Column {} already exists, will be overwritten".format(column_name))
        fragments = Cell3D._load_bed_fragments(path, resolution,type)
        self.tdg = pd.merge(self.tdg, fragments, on=["chrom", "pos"], how="left")
        self.tdg[column_name] = self.tdg["count"].fillna(0)
        self.tdg = self.tdg.drop(columns=["count"], axis=1)
        self.features.append(column_name)

    def add_CpG_data(self, path):
        CpG = Cell3D._load_CpG(path)
        self.features.append("CpG")
        self.tdg = pd.merge(self.tdg, CpG, on=["chrom", "pos"], how="left")#.dropna()

    def add_density(self, radius):
        self.features.append("density_"+str(radius))
        densities = []

        for point in self.tdg[["x","y","z"]].values:
            count = self.kdtree.query_ball_point(point,radius)
            density = len(count) / ((4/3) * np.pi * radius **3)
            densities.append(density)

        self.tdg["density_"+str(radius)] = np.array(densities)
        return None

    def add_feature_in_radius(self, feature, radius,type = "mean"):
        """
        feature: str
        radius: float
        type: str
            "mean" or "sum"

        Takes about 20 seconds for 250k points
        """
        indices_list = self.kdtree.query_ball_tree(self.kdtree, r=radius)
        if type =="mean":
            avgs = [self.tdg[feature].iloc[indices].mean() for indices in indices_list]
        elif type =="sum":
            avgs = [self.tdg[feature].iloc[indices].sum() for indices in indices_list]
        else:
            raise ValueError("type should be mean or sum")
        self.tdg[feature + "_" + type + "_in_radius_" + str(radius)] = avgs
        self.features.append(feature + "_avg_in_radius" + str(radius))
        return None

    # mutate the Cell3D object
    def _point_cloud_rotation(point_cloud, x_angle=None,y_angle=None,z_angle=None):
        """
        point_cloud: numpy array of shape (n,3)
        [xyz]_angle: float in radian

        Rotate point cloud by x,y,z angle and return the rotated point cloud in numpy array of shape (n,3)
        """
        if x_angle:
            rotation_matrix = np.array([[1,0,0],[0,np.cos(x_angle),-np.sin(x_angle)],[0,np.sin(x_angle),np.cos(x_angle)]])
            point_cloud = np.dot(point_cloud,rotation_matrix)
        if y_angle:
            rotation_matrix = np.array([[np.cos(y_angle),0,np.sin(y_angle)],[0,1,0],[-np.sin(y_angle),0,np.cos(y_angle)]])
            point_cloud = np.dot(point_cloud,rotation_matrix)
        if z_angle:
            rotation_matrix = np.array([[np.cos(z_angle),-np.sin(z_angle),0],[np.sin(z_angle),np.cos(z_angle),0],[0,0,1]])
            point_cloud = np.dot(point_cloud,rotation_matrix)

        return point_cloud

    # o
    def get_data(self,query = "",slice=False,random_slice=False,slice_width = 3,rotate=False,
                 rotate_x_angle=None,rotate_y_angle=None,rotate_z_angle=None):
        if query != "":
            tdg_temp = self.tdg.query(query)
        else:
            tdg_temp = self.tdg
        if rotate:
            rotated = Cell3D._point_cloud_rotation(tdg_temp[["x","y","z"]].values,
                                        x_angle=(rotate_x_angle is None and np.random.uniform(0, 2*np.pi) or rotate_x_angle),
                                        y_angle=(rotate_y_angle is None and np.random.uniform(0, 2*np.pi) or rotate_y_angle),
                                        z_angle=(rotate_z_angle is None and np.random.uniform(0, 2*np.pi) or rotate_z_angle)  
                                        )
            tdg_temp = tdg_temp.assign(x=rotated[:,0],y=rotated[:,1],z=rotated[:,2])

        if slice:
            if random_slice:
                xrange = (tdg_temp["x"].min(),tdg_temp["x"].max())
                slice_upper = np.random.uniform(xrange[0]+slice_width,xrange[1])
                slice_lower = slice_upper - slice_width
                tdg_temp = tdg_temp.query(f'x > @slice_lower & x < @slice_upper')
            else:
                slice_width = slice_width / 2 
                tdg_temp = tdg_temp.query("x > -@slice_width & x < @slice_width")

        return tdg_temp.copy()
    
    # analysis

    # data visualize
    def fast_plot3D(ax,x,y,z,c,view=(0,0,0),**kwargs):
        ax.scatter(x,y,z,c,**kwargs)
        ax.view_init(*view)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        return ax

    # Visualize the Cell3D object
    # 1. matplotlib

    # 2. output to cif
    def write_cif(self,outputpath = None):
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


# Analysis Functions

def calculate_RMSD(*dataframes):
        """
        Calculate RMSD for a list of dataframes.

        Parameters:
            *dataframes : variable length dataframe arguments
                Each dataframe should have columns 'chrom', 'pos', 'x', 'y', 'z'

        Returns:
            float, float : RMS RMSD and median RMSD of all pair combinations
        """
        
        num_structures = len(dataframes)
        if num_structures < 2:
            raise ValueError("At least 2 structures are required.")
        
        # Find common points across all structures
        common_points = dataframes[0][['chrom', 'pos']]
        for df in dataframes[1:]:
            common_points = pd.merge(common_points, df[['chrom', 'pos']], on=['chrom', 'pos'],how='inner')
        
        all_rmsds = []

        for i in range(num_structures):
            for j in range(i + 1, num_structures):
                coords_i = pd.merge(common_points, dataframes[i], on=['chrom', 'pos'])[['x', 'y', 'z']].values
                coords_j = pd.merge(common_points, dataframes[j], on=['chrom', 'pos'])[['x', 'y', 'z']].values
                
                # Subtract centroid
                centroid_i = rmsd.centroid(coords_i)
                centroid_j = rmsd.centroid(coords_j)
                coords_i -= centroid_i
                coords_j -= centroid_j
                
                # Calculate pairwise RMSD
                if rmsd.kabsch_rmsd(coords_i, coords_j) > rmsd.kabsch_rmsd(coords_i, -1.0 * coords_j):
                    coords_j *= -1.0            
                rmsd_value = rmsd.kabsch_rmsd(coords_i, coords_j)
                all_rmsds.append(rmsd_value)

        # Calculate RMS RMSD
        rms_rmsd = np.sqrt(np.mean(np.array(all_rmsds)**2))
        
        # Calculate median RMSD
        median_rmsd = np.median(all_rmsds)
        
        return rms_rmsd, median_rmsd

def feature_radial_distribution(cell, feature,random = False,random_seed = 42):
    """
    Calculate the radial distribution of a feature in a Cell3D object.

    Parameters:
        cell : Cell3D object
        feature : str
            Name of the feature to calculate radial distribution for
    Returns:
        pd.DataFrame : DataFrame containing the radial distribution
    """
    # 1. get center of mass
    tdg = cell.get_data()
    center = tdg[["x", "y", "z"]].mean()
    # 2. get distance to center
    tdg["radial_distance"] = np.sqrt((tdg["x"] - center["x"]) ** 2 + (tdg["y"] - center["y"]) ** 2 + (tdg["z"] - center["z"]) ** 2) 
    tdg["radial_distance"] = tdg["radial_distance"] / tdg["radial_distance"].mean()
    # 3. get radial distribution
    if random:
        np.random.seed(random_seed)
        tdg[feature] = tdg[feature].values[np.random.permutation(len(tdg[feature]))]
        #print("This is a randomization test.")

    tdg["feature_radial_distribution"] = tdg[feature] / tdg[feature].mean() * tdg["radial_distance"]
    return tdg