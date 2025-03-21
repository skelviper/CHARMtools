import pandas as pd
import numpy as np
from sklearn.preprocessing import LabelEncoder
import warnings
from scipy.spatial import cKDTree
import rmsd
import pybedtools
from scipy import stats
import re
import warnings

# TODO：
# 1. 3D inputation for feature
# 2. escape score calculation in cell 3D
# 3. "var explained method"
# 4. high resolution support


#import dask.dataframe as dd

def dev_only(func):
    import functools
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        if __name__ == "__main__":
            return func(*args, **kwargs)
        else:
            print(f"Skipping {func.__name__} as it's being imported")
    
    return wrapper


class Cell3D:
    #def __init__(self, cellname,resolution,tdg_path=None,pairs_path=None,type="3dg",on_disk=False,cache_dir=None):
    def __init__(self, cellname, resolution, tdg_path=None, pairs_path=None, type="3dg"):
        self.cellname = cellname
        if type == "3dg":
            self.tdg = self._load_tdg(tdg_path)
        elif type == "pairs":
            self.tdg = self._fdg_from_pairs(pairs_path)
        #self.tdg = self._load_tdg(tdg_path)
        self.resolution = resolution
        self.features = []
        self.kdtree = cKDTree(self.tdg[["x", "y", "z"]].values)
        self.chrom_length = None
        self.metadata = {}
        # self.on_disk = on_disk

        # if on_disk:
        #     path = cache_dir+"/"+cellname+".h5"
        #     self._save_to_disk(path)
        #     self.tdg = path

    def __repr__(self):
        self.get_info()
        return ""

    def get_info(self):
        print("CHARMtools Cell3D object v0.1")
        print("Cell name: {}".format(self.cellname))
        print("Resolution: {}".format(self.resolution))
        print("Features: {}".format(self.features))
        print(f"<Chromatin3DStructure with {len(self.tdg)} records>")
        print("Metadata:")
        for key, value in self.metadata.items():
            print(f"  {key}: {value}")

    # Construct a Cell3D object
    @dev_only
    def _fdg_read_pairs(self,pairs_path):
        """
        Read pairs file and return a pandas.DataFrame
        """
        data = CHARMio.parse_pairs("./data/mESCP2H3K27me3073.impute.pairs.gz")
        data = data[data[['phase_prob00','phase_prob01','phase_prob10','phase_prob11']].max(axis=1)>0.9]
        data["type"] = data[['phase_prob00','phase_prob01','phase_prob10','phase_prob11']].idxmax(axis=1)
        data["type"] = data["type"].apply(lambda x: x.replace("0","a").replace("1","b").replace("phase_prob",""))
        data["chr1"] = data["chr1"] + data["type"].apply(lambda x: x[0])
        data["chr2"] = data["chr2"] + data["type"].apply(lambda x: x[1])
        data = data[['chr1','pos1','chr2','pos2']]
        data = data.groupby(['chr1','pos1','chr2','pos2']).size().reset_index().rename(columns={0:"count"})
        data["count"] = np.log2(data["count"] / data["count"].sum() * 1000 +1)

        return data

    @dev_only
    def _calc_forces(points, contacts, k_backbone=1, k_repulsion=1, k_contact=1, rep_dist = 1.5,kd_tree="scipy"):
        num_points = points.shape[0]
        forces = np.zeros((num_points, 3))
        
        # Vectorized backbone forces calculation
        diffs = points[1:] - points[:-1]
        dists = np.linalg.norm(diffs, axis=1)
        diffs /= dists[:, None]  # Normalize
        force_magnitudes = k_backbone * (dists - 1)
        forces[:-1] += (force_magnitudes[:, None] * diffs)
        forces[1:] += (force_magnitudes[:, None] * diffs)

        if kd_tree == "scipy":
            tree = cKDTree(points)
        else:
            tree = build_kdtree(points)
        for i in range(num_points):
            if kd_tree == "scipy":
                idx = tree.query_ball_point(points[i], r=rep_dist)
            else:
                idx = query_kdtree(tree, points[i], rep_dist)
            for j in idx:
                if i < j:  
                    diff = points[i] - points[j]
                    dist = np.linalg.norm(diff)
                    if dist > 0:
                        force = k_repulsion * diff / (dist * dist * dist)  
                        forces[i] += force
                        forces[j] -= force

        # Calculate contact forces
        for c in contacts.itertuples():  # Assuming contacts is a DataFrame
            ind1, ind2, count = c.pos1, c.pos2, c.count
            diff = points[ind1] - points[ind2]
            dist = np.linalg.norm(diff)
            diff = diff / dist
            if dist > 0:
                forces[ind1] -= k_contact * diff * dist * count
                forces[ind2] += k_contact * diff * dist * count
        
        return forces

    @dev_only
    def _fdg(points, contacts, k_backbone=1, k_repulsion=1, k_contact=1, rep_dist=2,num_iter=500, lr=1):
        points = points.copy()
        num_points = points.shape[0]
        for i in tqdm.tqdm(range(num_iter)):
            forces = calc_forces(points, contacts, k_backbone, k_repulsion, k_contact, rep_dist)
            points += lr * forces
            #print(points)
            # print RMS force
            if i % 100 == 0:
                print(np.sqrt(np.mean(forces**2)))
        return points

    @dev_only
    def _fdg_from_pairs(self,pairs_path):
        """
        Construct a Cell3D object from pairs file.
        """
        #init_points()

        contacts = self._fdg_read_pairs(pairs_path)

        points = _fdg(points, contacts, k_backbone=1, k_repulsion=1, k_contact=1, rep_dist=2,num_iter=500, lr=1)
        return points

    def _load_tdg(self, tdg_path):
        if isinstance(tdg_path, pd.DataFrame):
            tdg = tdg_path.copy()
        elif isinstance(tdg_path, str):
            tdg = pd.read_csv(tdg_path, sep="\t", header=None, comment="#")
        else:
            raise ValueError("tdg_path should be a pandas.DataFrame or a string")
        
        tdg.columns = ["chrom", "pos", "x", "y", "z"]
        tdg["chrom"] = tdg["chrom"].str.replace("\(pat\)", "a", regex=True)
        tdg["chrom"] = tdg["chrom"].str.replace("\(mat\)", "b", regex=True)
        LE = LabelEncoder()
        tdg.chrom = pd.Categorical(tdg.chrom)
        tdg['chrom_code'] = LE.fit_transform(tdg['chrom'])
        return tdg

    def _load_bed_fragments(path, resolution, type = "allelic_resolved",peaks = None,flank = 0):
        """
        path: file path of a tsv file containing chrom, start, end, allele, score, strand
        peaks: a pd.DataFrame containing chrom, start, end
        """
        FRIP = None
        fragments = pd.read_csv(path, sep="\t", header=None)
        fragments.columns = ["chrom", "start", "end", "allele", "score", "strand"][:len(fragments.columns)]
        if peaks is not None:
            peaks = peaks.copy()
            peaks.columns = ["chrom", "start", "end"]
            #peaks["start"] = min(peaks["start"] - flank,0)
            peaks["start"] = peaks["start"] - flank
            peaks["start"] = peaks["start"].clip(lower=0)
            peaks["end"] = peaks["end"] + flank
            fragments_bed = pybedtools.BedTool.from_dataframe(fragments)
            peak_bed = pybedtools.BedTool.from_dataframe(peaks)
            intersect = peak_bed.intersect(fragments_bed,wa=True,wb=True)
            intersect = intersect.to_dataframe()
            intersect = intersect.iloc[:,3:9]
            intersect.columns =  ["chrom", "start", "end", "allele", "score", "strand"]
            FRIP = intersect.shape[0] / fragments.shape[0]
            fragments = intersect
            
        if type == "allelic_resolved":
            fragments = fragments.query("chrom.str.contains('chr')").query('allele != "."')
            fragments = fragments.assign(chrom=np.where(fragments["allele"] == "0", fragments["chrom"] + "a", fragments["chrom"] + "b"))
        else:
            fragments = pd.concat([
                fragments.assign(chrom=lambda x: x["chrom"] + "a"),
                fragments.assign(chrom=lambda x: x["chrom"] + "b")
            ])
        fragments["pos"] = ((fragments["start"] + fragments["end"]) / 2) // resolution * resolution
        fragments["pos"] = fragments["pos"].astype(int)
        return fragments.groupby(["chrom", "pos"]).size().reset_index().rename(columns={0: "count"}),FRIP

    def _load_CpG(CpG_path):
        """
        deprecated , use add_bedGraph_data instead
        """
        CpG = pd.read_csv(CpG_path, header=None, sep="\t")
        CpG.columns = ["chrom", "pos", "CpG"]
        if CpG["chrom"].str.startswith("chr").sum() == 0:
            CpG["chrom"] = "chr" + CpG["chrom"]
        CpG["pos"] = CpG["pos"].astype(int)
        return pd.concat([
            CpG.assign(chrom=lambda x: x["chrom"] + "a"),
            CpG.assign(chrom=lambda x: x["chrom"] + "b")
        ])    

    def add_bed_data(self, path, column_name, resolution=None,type="allelic_resolved",peaks = None,flank=0):

        if resolution is None:
            resolution = self.resolution
        if column_name in self.tdg.columns:
            warnings.warn("Column {} already exists, will be overwritten".format(column_name))
        
        fragments,FRIP = Cell3D._load_bed_fragments(path, resolution,type,peaks=peaks,flank=flank)
        self.tdg = pd.merge(self.tdg, fragments, on=["chrom", "pos"], how="left")
        self.tdg[column_name] = self.tdg["count"].fillna(0)
        self.tdg = self.tdg.drop(columns=["count"], axis=1)
        self.features.append(column_name)
        if FRIP is not None:
            self.metadata[column_name + "_FRIP"] = FRIP

    def add_bedGraph_data(self,path,column_name, resolution=None,type="allelic_resolved"):
        """
        Add features from bedGraph file to the Cell3D object.

        Parameters:
            path : str
                Path to the bedGraph file
            column_name : str
                Name of the column to be added
            resolution : int
                Resolution of the bedGraph file, default is the resolution of the Cell3D object
            type : str
                "allelic_resolved" or "non_allelic_resolved" not implemented yet

        Returns:
            None
        """
        if resolution is None:
            resolution = self.resolution
        positions = self.tdg[["chrom","pos"]].copy()
        positions["start"] = positions["pos"]
        positions["end"] = positions["pos"] + resolution
        positions = positions[["chrom","start","end"]]

        if type == "allelic_resolved":
            bedgraph = pd.read_csv(path,sep="\t",header=None)
            bedgraph.columns = ["chrom","start","end","allele",column_name]
            bedgraph = bedgraph.query("chrom.str.contains('chr')").query('allele != "."')
            bedgraph = bedgraph.assign(chrom=np.where(bedgraph["allele"] == 0, bedgraph["chrom"] + "a", bedgraph["chrom"] + "b"))
            bedgraph = bedgraph[["chrom","start","end",column_name]]
            #print(bedgraph.query('chrom.str.contains("b") & RNA > 0'))

        else: 
            bedgraph = pd.read_csv(path,sep="\t",header=None)
            bedgraph.columns = ["chrom","start","end",column_name]
            bedgraph = pd.concat([
                bedgraph.assign(chrom=lambda x: x["chrom"] + "a"),
                bedgraph.assign(chrom=lambda x: x["chrom"] + "b")
            ])

        
        # merge positions and bedgraph
        positions_bed = pybedtools.BedTool.from_dataframe(positions)
        bedgraph_bed = pybedtools.BedTool.from_dataframe(bedgraph)
        merged_bed = positions_bed.intersect(bedgraph_bed,wa=True,wb=True)
        newcol = merged_bed.to_dataframe()[["chrom","start","end","thickStart"]].groupby(["chrom","start","end"]).sum().reset_index()
        newcol.columns = ["chrom","start","end",column_name]
        newcol["pos"] = newcol["start"]
        #self.tdg =
        temp = pd.merge(self.tdg,newcol[["chrom","pos",column_name]],on=["chrom","pos"],how="left")
        temp[column_name] = temp[column_name].fillna(0)
        self.tdg = temp
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
    
    def add_knn_density(self, k):
        """
        use mean distance of k nearest neighbors as density
        """
        self.features.append("knn_density_"+str(k))
        densities = []
        for point in self.tdg[["x","y","z"]].values:
            dist,ind = self.kdtree.query(point,k=k+1)
            density = dist.mean()
            densities.append(density)
        self.tdg["knn_density_"+str(k)] = np.array(densities)
        return None

    def add_feature_in_radius(self, feature, radius,type = "mean",if_self = False,if_rank = False):
        """
        Smooth given feature by averaging or summing over a sphere of given radius.

        Parameters:
            feature: str
                Name of the feature to be smoothed
            radius: float
                Radius of the sphere
            type: str
                "mean" or "sum"
        Additional information:
            Takes about 20 seconds for 250k points
        
        Returns:
            None
        """
        from scipy.stats import rankdata
        indices_list = self.kdtree.query_ball_tree(self.kdtree, r=radius)
        if not if_self:
            indices_list = [indices[1:] for indices in indices_list]
        if type =="mean":
            avgs = [self.tdg[feature].iloc[indices].mean() for indices in indices_list]
        elif type =="sum":
            avgs = [self.tdg[feature].iloc[indices].sum() for indices in indices_list]
        else:
            raise ValueError("type should be mean or sum")
        if if_rank:
            avgs = rankdata(avgs,nan_policy="omit") / sum(np.isfinite(avgs))

        self.tdg[feature + "_" + type + "_in_radius_" + str(radius)] = avgs
        self.features.append(feature + "_avg_in_radius" + str(radius))
        return None

    def add_chrom_length(self,chrom_length_path):
        chrom_length = pd.read_csv(chrom_length_path,sep="\t",header=None)
        chrom_length.columns = ["chrom","size"]
        chrom_length["chrom"] = chrom_length["chrom"].str.replace("\(pat\)", "a", regex=True)
        chrom_length["chrom"] = chrom_length["chrom"].str.replace("\(mat\)", "b", regex=True)
        chrom_length["chrom"] = chrom_length["chrom"].str.replace("pat", "a", regex=True)
        chrom_length["chrom"] = chrom_length["chrom"].str.replace("mat", "b", regex=True)

        self.features.append("chrom_length")
        self.chrom_length = chrom_length
    
    def calc_distance_matrix(self,genome_coord):
        """
        INPUT:
            genome_coord: str, format like chrom:start-end or list/tuple of chrom,start,end. \|
                          whole chromosome is also acceptable. e.g. "chr1a:10000-20000" or ["chr1a",10000,20000] or "chr1a
        OUTPUT:
            distance_matrix: symettrical distance matrix of the given region, np.array of shape (n,n)
        """
        chrom,start,end = _auto_genome_coord(genome_coord)

        if start is None and self.chrom_length is None:
            raise ValueError("Running whole chromosome calculation with chrom_length is not available, |\
                                please run add_chrom_length first")
        
        if start is None:        
            matsize = self.chrom_length.query("chrom == @chrom")["size"].values[0] // self.resolution + 1
            reconstruct_df = self.tdg.query("chrom == @chrom").copy()
        else:
            matsize = (end - start - 1) // self.resolution + 1
            reconstruct_df = self.tdg.query("chrom == @chrom & pos >= @start & pos < @end").copy()
            reconstruct_df["pos"] = reconstruct_df["pos"] - start

        reconstruct_df["pos"] = reconstruct_df["pos"] // self.resolution
        coordinates = reconstruct_df.iloc[:, 2:5].values
        diff = coordinates[:, np.newaxis, :] - coordinates[np.newaxis, :, :]
        dist_matrix_reconstruct = np.linalg.norm(diff, axis=-1)
        full_dist_matrix = np.full((matsize, matsize), np.nan)
        for i, pos_i in enumerate(reconstruct_df['pos']):
            for j, pos_j in enumerate(reconstruct_df['pos']):
                full_dist_matrix[pos_i, pos_j] = dist_matrix_reconstruct[i, j]
    
        return full_dist_matrix
    
    def calc_feature_matrix(self,genome_coord,feature):
        """

        !!! Danger Zone!!!

        INPUT:
            genome_coord: str, format like chrom:start-end or list/tuple of chrom,start,end. \|
                          whole chromosome is also acceptable. e.g. "chr1a:10000-20000" or ["chr1a",10000,20000] or "chr1a
            feature: str, feature to calculate distance matrix
        OUTPUT:
            feature_matrix: symettrical feature matrix of the given region, np.array of shape (n,n)
            feature_vec: np.array of shape (n,) indicating whether the feature is missing
        """
        chrom,start,end = _auto_genome_coord(genome_coord)
        vec_size = (end-start-1) // self.resolution + 1
        feature_vec = np.zeros(vec_size)
        data = self.get_data(genome_coord)
        for row in data.iterrows():
            row = row[1]
            pos = row["pos"] 
            feature_vec[(pos-start)//self.resolution] = row[feature]
        feature_vec = feature_vec == 0
        mat = 1/(self.calc_distance_matrix(genome_coord) + 1)
        feature_mat = mat.copy()
        # it‘s kind of tricky here since 0 leads to CNV and
        #  NaN leads to overestimate of noise distance/difficult to normalize
        # feature_mat[feature_vec,:] = np.nan
        # feature_mat[:,feature_vec] = np.nan
        feature_mat[feature_vec,:] = 0
        feature_mat[:,feature_vec] = 0
        return feature_mat,feature_vec
    
    def calc_feature_proximity_matrix(self,genome_coord,feature,distance_threshold=3):
        """
        Calculate the proximity matrix of a feature in a given region.

        Parameters:
            genome_coord : str
                Genome coordinates in the format of "chrom:start-end"
            feature : str
                Name of the feature
            distance_threshold : float
                Distance threshold to consider two points as proximal

        Returns:
            np.array : Proximity matrix
        """
        chrom,start,end = _auto_genome_coord(genome_coord)
        vec_size = (end-start-1) // self.resolution + 1
        feature_vec = np.zeros(vec_size)
        data = self.get_data(genome_coord)
        for row in data.iterrows():
            row = row[1]
            pos = row["pos"] 
            feature_vec[(pos-start)//self.resolution] = row[feature]
        feature_vec = feature_vec == 0

        #mat = 1/(self.calc_distance_matrix(genome_coord) + 1)
        mat = self.calc_distance_matrix(genome_coord) < distance_threshold
        feature_mat = mat.copy()

        feature_mat[feature_vec,:] = 0
        feature_mat[:,feature_vec] = 0
        return feature_mat,feature_vec


    # mutate the Cell3D object
    def _point_cloud_rotation(point_cloud, x_angle=None,y_angle=None,z_angle=None):
        """
        Rotate point cloud by x,y,z angle and return the rotated point cloud in numpy array of shape (n,3)

        Parameters:
            point_cloud: numpy array of shape (n,3)
            [xyz]_angle: float in radian

        Returns:
            point_cloud: numpy array of shape (n,3)
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
    def get_data(self,genome_coord = "",slice=False,random_slice=False,slice_width = 3,rotate=False,
                 rotate_x_angle=None,rotate_y_angle=None,rotate_z_angle=None):
        """
        Get the data of the Cell3D object.

        Parameters:
            query: str
                Query string to filter the data
            slice: bool
                Whether to slice the data
            random_slice: bool
                Whether to slice randomly
            slice_width: float
                Width of the slice
            rotate: bool
                Whether to rotate the data
            rotate_[xyz]_angle: float
                Rotation angle in radian
        Returns:
            pandas.DataFrame
        """
        if genome_coord != "":
            chrom,start,end = _auto_genome_coord(genome_coord)
            if start is None:
                tdg_temp = self.tdg.query("chrom == @chrom").copy()
            else:
                start = start // self.resolution * self.resolution
                resolution = self.resolution
                tdg_temp = self.tdg.query("chrom == @chrom & pos >= @start & pos < @end").copy()
        else:
            tdg_temp = self.tdg.copy()
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

        return tdg_temp
    
    # 2. output to cif
    def write_cif(self,factor_b,outputpath = None):
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
            cif_str3 += f"ATOM {i+1} C CA . GLY {row['chrom']} {entity_id} {chrom_index} ? {row['x']} {row['y']} {row['z']} {row[factor_b]}\n"

        #print(cif_str3)

        # Open the file in write mode
        with open(outputpath, 'w') as f:
            # Write the three blocks to the file
            f.write(cif_str)
            f.write(cif_str2)
            f.write(cif_str3)

        print("Done " + cellname)
        return None
    
    # data visualize
    def plot3D(self,query=None,**kwargs):
        """
        Plot the 3D structure of the Cell3D object.
        """
        if query is None:
            plot3D(cell = self.tdg, **kwargs)
        else:
            plot3D(cell = self.tdg.query(query), **kwargs)

    def _roated_hic_mat(matrix, h):
        matrix = np.nan_to_num(matrix)
        rotated =  rotate(matrix, 45, reshape=True)
        return rotated[rotated.shape[0]//2-h:rotated.shape[0]//2, :]

    def plotDistanceMatrix(self,genome_coord,type=None,h=None,**kwargs):
        """
        Plot the distance matrix of a given region.
        Parameters:
            genome_coord : str
                Genome coordinates in the format of "chrom:start-end"
            type : str
                "normal" or "rotate"
        """
        import matplotlib.pyplot as plt
        chrom,start,end = _auto_genome_coord(genome_coord)
        mat = self.calc_distance_matrix(genome_coord)
        # extra params
        cmap = kwargs.get("cmap","YlOrRd")
        vmin = kwargs.get("vmin",None)
        vmax = kwargs.get("vmax",None)

        if type == "normal":
            if vmin and vmax:
                plt.imshow(mat,cmap=cmap,vmin=vmin,vmax=vmax,extent=[start,end,end,start])
            else:
                plt.imshow(mat,cmap=cmap,extent=[start,end,end,start])
        elif type == "rotate":
            rotated = _roated_hic_mat(mat,h)
            if vmin and vmax:
                plt.imshow(rotated,cmap=cmap,vmin=vmin,vmax=vmax)
            else:
                plt.imshow(rotated,cmap=cmap)

        plt.colorbar()
        plt.show()
        return None


    # analysis
    def calc_radial_position(self,key="radial_position",if_rank = False,if_norm_max = False, if_norm_mean = True):
        """
        Calculate the radial position of each point.
        """
        data = self.tdg
        center = data[["x", "y", "z"]].mean()
        data[key] = np.sqrt((data[["x", "y", "z"]] - center).pow(2).sum(axis=1))
        if if_norm_max:
            data[key] = data[key] / data[key].max()
        if if_rank:
            data[key] = data[key].rank(method='first')
            data[key] = data[key] / data[key].max()
        if if_norm_mean:
            data[key] = data[key] / data[key].mean()
        self.tdg = data
        self.features.append(key)
        return None
    
    def feature_radial_distribution(cell, feature,random = False,random_seed = 42,if_normalize_avg = False,if_rank = False):
        """
        TODO: use calc_radial_position if present to avoid recalculation 
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
        if if_normalize_avg:
            tdg["radial_distance"] = tdg["radial_distance"] / tdg["radial_distance"].mean()

        if if_rank:
            tdg["radial_distance"] = tdg["radial_distance"].rank(method='first')

        # 3. get radial distribution
        if random:
            np.random.seed(random_seed)
            tdg[feature] = tdg[feature].values[np.random.permutation(len(tdg[feature]))]
            #print("This is a randomization test.")

        tdg["feature_radial_distribution"] = tdg[feature] / tdg[feature].sum() * tdg["radial_distance"]
        return tdg

    def calc_radius_gyration(self,genome_coord):
        """
        Calculate the radius of gyration for a given region.
        """
        data = self.get_data(genome_coord)
        # center data
        data_centered = data[["x", "y", "z"]] - data[["x", "y", "z"]].mean()
        # calculate radius of gyration
        radius_gyration = np.sqrt((data_centered ** 2).sum(axis=1).mean())
        return radius_gyration
    
    def calc_feature_distances(self,feature_key,quantile,random_seed=0,k=10):
        """
        INPUT:
            df: charm tdg dataframe from charm_get_3dplot_data
            feature_key: feature to calculate distance
            threshold: threshold to divide points into two groups
        OUTPUT:
            distances: np.array of inter-chromosomal distances of sampled points with feature > threshold and random distances
        """
        tdg = self.tdg
        selected_points = tdg.sort_values(feature_key,ascending=False).iloc[:int(tdg.shape[0]*quantile)]
        kdtree_seleted_points = cKDTree(selected_points[['x', 'y', 'z']])

        random_points = tdg.sample(n=selected_points.shape[0], random_state=random_seed)
        kdtree_random_points = cKDTree(random_points[['x', 'y', 'z']])

        distances_selected = kdtree_seleted_points.query(selected_points[['x', 'y', 'z']], k=k)[0].mean(axis=1)
        distances_random = kdtree_random_points.query(random_points[['x', 'y', 'z']], k=k)[0].mean(axis=1)
        return np.array([distances_selected, distances_random])

    def calc_singlecell_compartment(self,method="cooltools"):
        """
        Calculate compartment score for each cell,use CpG to correcct the sign
        INPUT:
        OUTPUT:
            add PC1~3 to self.tdg by default
        """
        import cooltools
        from sklearn.decomposition import PCA
        from . import compartment

        chroms = self.tdg.chrom.unique()
        tdg_per_chroms = []
        for i in chroms:

            mat = self.calc_distance_matrix(i)
            mat = 1/(mat + 1)
            if method == "cooltools":
                eigval,eigvec = cooltools.api.eigdecomp.cis_eig(mat)
            else:
                mat = np.nan_to_num(mat)
                mat_corr = np.corrcoef(compartment.getOEMatrix(mat))
                mat_corr= np.nan_to_num(mat_corr)
                pca=PCA(n_components=3)
                pca.fit(mat_corr)
                eigvec = pca.components_
                eigval,eigvec = cooltools.api.eigdecomp.cis_eig(mat_corr)

            pca_df = pd.DataFrame(eigvec.T)
            pca_df.columns = ["PC1", "PC2", "PC3"]
            pca_df["pos"] = np.arange(0, len(pca_df)) * self.resolution 

            cpg = self.tdg.query('chrom == @i')[["chrom","pos","CpG"]]
            temp = pd.merge(cpg,pca_df)

            # correct the sign of PC1~3
            for i in range(3):
                i = i + 1
                if mat_cor_with_na(temp["PC" + str(i)].values,temp["CpG"].values)[1] < 0:
                    temp["PC" + str(i)] = -temp["PC" + str(i)]

            temp['PC1_sign'] = temp['PC1'].apply(lambda x: 'positive' if x >= 0 else 'negative')
            temp['CompID'] = (temp['PC1_sign'] != temp['PC1_sign'].shift()).cumsum()

            tdg_per_chroms.append(temp)

            # merge all chroms and add to self.tdg
        temp = pd.concat(tdg_per_chroms)

        self.tdg = pd.merge(self.tdg,temp[["chrom","pos","PC1","PC2","PC3","PC1_sign","CompID"]])
            

#Example
# from scipy.stats import ks_2samp
# distances_cell = cell.calc_feature_distances("Random", 0.01, random_seed=1,k=20)
# #distances_cell_random = cell.calc_feature_distances("Random", 0.05, random_seed=1,k=20)
# #distances = distances_cell[0] 
# #random_distances = distances_cell_random[0]
# distances, random_distances = distances_cell[0], distances_cell[1]
# fig, axes = plt.subplots(figsize=(6, 3), ncols=2, dpi=120)

# # histogram of vector, use fraction as y axis
# ax1 = axes[0]
# ax1.hist(random_distances, bins=20, color="grey", alpha=0.5, label="random")
# ax1.hist(distances, bins=20, color="red", alpha=0.5, label="H3K9me3")
# ax1.set_xlabel("Distance")
# ax1.set_ylabel("Frequency")
# #ax1.set_xlim(0,4)
# ax1.legend()

# # 计算累积分布
# hist_random, edges_random = np.histogram(random_distances, bins=20, density=True)
# hist_open, edges_open = np.histogram(distances, bins=20, density=True)
# cdf_random = np.cumsum(hist_random)/np.sum(hist_random)
# cdf_open = np.cumsum(hist_open)/np.sum(hist_open)

# # 绘制累积分布曲线
# ax2 = axes[1]
# ax2.plot(edges_random[:-1], cdf_random, color="grey", alpha=0.5, label="random")
# ax2.plot(edges_open[:-1], cdf_open, color="red", alpha=0.5, label="H3K9me3")
# ax2.set_xlabel("Distance")
# ax2.set_ylabel("Cumulative frequency")
# #ax2.set_xlim(0,4)
# ax2.legend()

# # 进行Kolmogorov-Smirnov检验
# ks_stat, ks_pval = ks_2samp(random_distances, distances)
# print(f"KS statistic: {ks_stat}")
# print(f"KS test p-value: {ks_pval}")
# # t检验
# from scipy.stats import ttest_ind
# t_stat, t_pval = ttest_ind(random_distances, distances)
# print(f"t statistic: {t_stat}")
# print(f"t test p-value: {t_pval}")

# plt.tight_layout()
# plt.show()
    
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

def calc_feature_distances_v1(cell, feature_key, threshold, points_num_per_chrom=50, points_num_other_chrom=100, random_seed=0):
    """
    This is similar to LongCai's method. 
    INPUT: 
        df: charm tdg dataframe from charm_get_3dplot_data
        feature_key: feature to calculate distance
        threshold: threshold to divide points into two groups
    OUTPUT:
        distances: np.array of inter-chromosomal distances of sampled points with feature > threshold and random distances

    """
    from scipy.spatial import distance_matrix
    from scipy.stats import mannwhitneyu
    import matplotlib.pyplot as plt
    # generate randomized inter-distance
    np.random.seed(random_seed)
    df = cell.get_data()
    random_distances = []
    for chrom in set(df.chrom.values):
        chrom_points = df.query("chrom == @chrom").sample(n=min(len(df.query("chrom == @chrom")), points_num_per_chrom), random_state=random_seed)[["x", "y", "z"]].values
        else_chrom_points = df.query("chrom != @chrom").sample(n=min(len(df.query("chrom != @chrom")), points_num_other_chrom), random_state=random_seed)[["x", "y", "z"]].values
        random_distances.extend(distance_matrix(chrom_points, else_chrom_points).flatten())
        #random_distances.extend(distance_matrix(chrom_points, chrom_points).flatten())
    random_distances = np.array(random_distances)

    df['status'] = np.where(df[feature_key] > threshold, 'open', 'close')
    df = df.query("status == 'open'")

    distances = []
    for chrom in set(df.chrom.values):
        chrom_points = df.query("chrom == @chrom").sample(n=min(len(df.query("chrom == @chrom")), points_num_per_chrom), random_state=random_seed)[["x", "y", "z"]].values
        else_chrom_points = df.query("chrom != @chrom").sample(n=min(len(df.query("chrom != @chrom")), points_num_other_chrom), random_state=random_seed)[["x", "y", "z"]].values
        distances.extend(distance_matrix(chrom_points, else_chrom_points).flatten())
        #distances.extend(distance_matrix(chrom_points, chrom_points).flatten())
    distances = np.array(distances)
    print("Mean distance of random points: %.2f" % np.mean(random_distances))
    print("Mean distance of feature points: %.2f" % np.mean(distances))
    # wilcox test of random vs open
    # Perform the Mann-Whitney U test
    U_stat, p_val = mannwhitneyu(random_distances, distances)

    print(f"Mann-Whitney U statistic: {U_stat}")
    print(f"Mann-Whitney test p-value: {p_val}")

    return np.array([distances, np.random.choice(random_distances, size=distances.shape[0], replace=False)])


def calc_feature_distances_v2():
    """
    This is a simple method 
    """

### Eample usage ###
# distances_cell = calc_feature_distances(tdg,"count_ct", 1, points_num_per_chrom=50, points_num_other_chrom=100, random_seed=0)
# distances,random_distances = distances_cell
# fig, axes = plt.subplots(figsize=(6, 3), ncols=2, dpi=120)

# histogram of vector, use fraction as y axis
# ax1 = axes[0]
# ax1.hist(random_distances, bins=20, color="grey", alpha=0.5, label="random")
# ax1.hist(distances, bins=20, color="red", alpha=0.5, label="H3K9me3")
# ax1.set_xlabel("Distance")
# ax1.set_ylabel("Frequency")
# ax1.legend()

# # 计算累积分布
# hist_random, edges_random = np.histogram(random_distances, bins=20, density=True)
# hist_open, edges_open = np.histogram(distances, bins=20, density=True)
# cdf_random = np.cumsum(hist_random)/np.sum(hist_random)
# cdf_open = np.cumsum(hist_open)/np.sum(hist_open)

# # 绘制累积分布曲线
# ax2 = axes[1]
# ax2.plot(edges_random[:-1], cdf_random, color="grey", alpha=0.5, label="random")
# ax2.plot(edges_open[:-1], cdf_open, color="red", alpha=0.5, label="H3K9me3")
# ax2.set_xlabel("Distance")
# ax2.set_ylabel("Cumulative frequency")
# ax2.legend()

# # 进行Kolmogorov-Smirnov检验
# ks_stat, ks_pval = ks_2samp(random_distances, distances)
# print(f"KS statistic: {ks_stat}")
# print(f"KS test p-value: {ks_pval}")

# plt.tight_layout()
# plt.show()

def mat_cor_with_na(mat1,mat2):
    # Calculate distance matrices
    distance_matrix_1 = mat1.flatten()
    distance_matrix_2 = mat2.flatten()

    # Replace inf values with nan
    distance_matrix_1 = np.where(np.isinf(distance_matrix_1), np.nan, distance_matrix_1)
    distance_matrix_2 = np.where(np.isinf(distance_matrix_2), np.nan, distance_matrix_2)

    # Remove any NaN values from both arrays (only where both have NaNs in the same position)
    mask = ~np.isnan(distance_matrix_1) & ~np.isnan(distance_matrix_2)
    distance_matrix_1 = distance_matrix_1[mask]
    distance_matrix_2 = distance_matrix_2[mask]

    # Check if there are any remaining NaNs or infs
    if not np.isfinite(distance_matrix_1).all() or not np.isfinite(distance_matrix_2).all():
        raise ValueError("The input arrays contain infs or NaNs after preprocessing.")

    # Now you can safely call pearsonr
    pearsonr_value,_ = stats.pearsonr(distance_matrix_1, distance_matrix_2)
    spearmanr_value,_ = stats.spearmanr(distance_matrix_1, distance_matrix_2)

    return [pearsonr_value,spearmanr_value]

def plot3D(cell, color_by, resolution=200000, smooth=True, smoothness=2,
            spline_degree=3,cmap="viridis",title='3D Chromatin Structure Visualization',
            vmax=None,vmin=None,width = 5):
    """
    Function to plot a 3D chromatin structure using Plotly. Supports both smooth and straight line
    representations. The plot is colored based on a specified column and allows interactive manipulation.

    :param dataframe: DataFrame containing the chromatin data.
    :param color_by: Column name in the dataframe to color the plot by.
    :param resolution: Resolution for determining continuity. Points closer than this
                    in the 'pos' column are considered continuous.
    :param smooth: Boolean to determine if the lines should be smooth or straight.
    :param smoothness: Number of points to interpolate between two consecutive points.
    :param spline_degree: Degree of the spline curve.
    :param cmap: Colormap to use for coloring the plot.
    :param title: Title for the plot.
    """
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    from scipy.interpolate import splprep, splev
    import plotly.express as px

    # Creating a Plotly figure
    fig = make_subplots(rows=1, cols=1, specs=[[{'type': 'scatter3d'}]])
    # if cell is not a pd.DataFrame
    if not isinstance(cell, pd.DataFrame):
        try:
            cell = cell.get_data()
        except:
            raise ValueError("cell should be a pandas.DataFrame or a Cell3D object")

    dataframe = cell.copy()
    # Defining a colormap based on the color_by column
    # color_scale = 'Viridis' if dataframe[color_by].dtype.kind in 'biufc' else 'Plotly3'
    if dataframe[color_by].dtype.kind in 'biufc':
        # Use a continuous colorscale for numeric data
        color_scale = cmap
        if vmax is None:
            vmax = dataframe[color_by].quantile(0.95)
        if vmin is None:
            vmin = dataframe[color_by].quantile(0.05)
        dataframe[color_by] = np.clip(dataframe[color_by], vmin, vmax)
    else:
        # For categorical data, create a discrete color map
        unique_categories = dataframe[color_by].unique()
        # str to sth like plotly_utils.colors.qualitative.Plotly
        color_scale = px.colors.sequential.Rainbow
        color_map = {category: color_scale[i % len(color_scale)] for i, category in enumerate(unique_categories)}


        # Grouping data by chromosome
    grouped_data = dataframe.groupby('chrom')
    show_colorbar=0
    for chrom, group in grouped_data:
        group = group.sort_values('pos')
        pos_diff = group['pos'].diff().fillna(0)
        break_indices = np.where(pos_diff > resolution)[0]
        segments = np.split(group, break_indices)

        # 初始化存储合并后的坐标点
        x_all, y_all, z_all, color_values_all = [], [], [], []

        for segment in segments:
            if len(segment) < 2:
                continue

            if smooth:
                if len(segment) <= spline_degree:
                    continue
                x, y, z = segment['x'].values, segment['y'].values, segment['z'].values
                tck, _ = splprep([x, y, z], s=smoothness, k=spline_degree)
                u_new = np.linspace(0, 1, (len(segment) - 1) * 20 + 1)
                x_new, y_new, z_new = splev(u_new, tck)
                
                x_all.extend(x_new)
                y_all.extend(y_new)
                z_all.extend(z_new)

                # Adjust color values for interpolated points
                if dataframe[color_by].dtype.kind in 'biufc':
                    color_values = np.repeat(segment[color_by].values[:-1], 20)
                else:
                    color_values = [color_map[val] for val in segment[color_by].values[:-1]]
                    color_values = np.repeat(color_values, 20, axis=0)

                color_values_all.extend(color_values)
                #print(len(color_values_all)

            else:
                # 直接使用原始点
                x_all.extend(segment['x'])
                y_all.extend(segment['y'])
                z_all.extend(segment['z'])
                if dataframe[color_by].dtype.kind in 'biufc':
                    color_values_all.extend(segment[color_by].values)
                else:
                    color_values_all.extend([color_map[val] for val in segment[color_by].values])

        fig.add_trace(go.Scatter3d(x=x_all, y=y_all, z=z_all,
                                mode='lines',
                                line=dict(color=color_values_all, colorscale=color_scale, width=width,showscale=True if show_colorbar == 0 else False,
                                          colorbar=dict(title=color_by, x=1.25)),
                                name=chrom))
        show_colorbar += 1
    fig.update_layout(title=title,
                    scene=dict(xaxis_title='X', yaxis_title='Y', zaxis_title='Z'),
                    width=800, height=600)
    fig.update_scenes(aspectmode="data")
    fig.show()


def calc_volume(point_cloud,method="convexhull",alpha=0.2):
    """
    Calculate the volume of a point cloud.

    Parameters:
        point_cloud: numpy array of shape (n,3)
        method: str
            "alphashape" or "convexhull"
        alpha: float
            alpha value for alphashape method
    Output:
        volume: float
    """
    
    if method =="alphashape":
        import alphashape
        alpha_shape = alphashape.alphashape(point_cloud, alpha=alpha)
        return alpha_shape.volume
    
    if method =="convexhull":
        from scipy.spatial import ConvexHull
        hull = ConvexHull(point_cloud)
        return hull.volume

    else:
        raise ValueError("method should be alphashape or convexhull")

def _auto_genome_coord(genome_coord):
    """
    Automatically convert genome_coord to chrom,start,end format
    INPUT:

    OUTPUT:
    """
    # determine the genome_coord format
    if isinstance(genome_coord,str):
        if ":" in genome_coord:
            chrom,start,end = re.split(":|-",genome_coord)
            start,end = int(start),int(end)
            mat_type = "region"
        else:
            chrom,start,end = genome_coord,None,None
            mat_type = "chrom"
    elif isinstance(genome_coord,(list,tuple)):
        chrom,start,end = genome_coord
        mat_type = "region"
    else:
        raise ValueError('Genome_coord should be str or list/tuple. e.g. "chr1a:10000-20000" or ["chr1a",10000,20000] or "chr1a"')
    
    return chrom,start,end
