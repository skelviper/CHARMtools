# A re-implantation of Cell3D object with hdf5 support
from scipy.spatial import cKDTree
import warnings
import pandas as pd
import numpy as np
from scipy import stats
import scipy
import re
import pybedtools
import os
import tqdm
import copy

from sklearn.cluster import DBSCAN
from sklearn.preprocessing import LabelEncoder
import warnings

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
    # Initialization
    def __init__(self, cellname,resolution,tdg_path=None,pairs_path=None,type="3dg",on_disk=False,on_disk_path=None):
        self.cellname = cellname
        if type == "3dg":
            self.tdg = self._load_tdg(tdg_path)
        elif type == "pairs":
            self.tdg = self._fdg_from_pairs(pairs_path)
        self.record_size = len(self.tdg)
        self.resolution = resolution
        self.features = []
        #self.kdtree = cKDTree(self.tdg[["x", "y", "z"]].values)
        self.kdtree = None
        self.chrom_length = None
        self.on_disk = on_disk
        self.on_disk_path = on_disk_path
        self.metadata = {}
        self.extra = {}

        if on_disk:
            self.to_disk(on_disk_path)

    def __repr__(self):
        self.get_info()
        return ""

    def dense_tdg(self):
        if self.on_disk:
            self.to_memory()
        if self.chrom_length is None:
            raise ValueError("Chrom length is not available, please run add_chrom_length first")
        result_df = pd.DataFrame(columns = ['chrom','pos'])
        for index,row in self.chrom_length.iterrows():
            chrom = row["chrom"]
            length = row["size"]
            positions = np.arange(0,length,self.resolution).astype(int)
            temp_df = pd.DataFrame({'chrom':chrom,'pos':positions}) 
            result_df = pd.concat([result_df,temp_df],ignore_index=True)
        new_tdg = pd.merge(result_df,self.tdg,on=['chrom','pos'],how='left').copy()
        new_tdg.sort_values(by=['chrom','pos'],inplace=True)
        if self.kdtree is not None:
            self.kdtree = cKDTree(new_tdg[["x", "y", "z"]].values)
        self.tdg = new_tdg
        return None
    
    def build_kdtree(self):
        if self.on_disk:
            self.to_memory()
        self.kdtree = cKDTree(self.tdg[["x", "y", "z"]].values)
    
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
            forces = _calc_forces(points, contacts, k_backbone, k_repulsion, k_contact, rep_dist)
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

    def _load_tdg(self, tdg):
        if isinstance(tdg, pd.DataFrame):
            tdg = tdg.copy()
        elif isinstance(tdg, str):
            tdg = pd.read_csv(tdg, sep="\t", header=None, comment="#")
        else:
            raise ValueError("tdg should be a pandas.DataFrame or a string")
        
        tdg.columns = ["chrom", "pos", "x", "y", "z"]
        tdg["chrom"] = tdg["chrom"].astype("category")
        tdg["pos"] = tdg["pos"].astype(int)
        tdg["x"] = tdg["x"].astype(np.float32)
        tdg["y"] = tdg["y"].astype(np.float32)
        tdg["z"] = tdg["z"].astype(np.float32)

        tdg["chrom"] = tdg["chrom"].str.replace("\(pat\)", "a", regex=True)
        tdg["chrom"] = tdg["chrom"].str.replace("\(mat\)", "b", regex=True)
        LE = LabelEncoder()
        tdg.chrom = pd.Categorical(tdg.chrom)
        tdg['chrom_code'] = LE.fit_transform(tdg['chrom'])
        return tdg

    def to_disk(self, on_disk_path = None):
        if on_disk_path == None and self.on_disk_path == None:
            self.on_disk_path = f"{self.cellname}.h5"
            warnings.warn("No path provided, saving to default path ./")
        elif on_disk_path != None:
            self.on_disk_path = on_disk_path + "/" + f"{self.cellname}.h5"
        if self.kdtree is not None:
            self.kdtree = None
        self.on_disk = True
        # create path if folder does not exist
        if not os.path.exists(os.path.dirname(self.on_disk_path)):
            os.makedirs(os.path.dirname(self.on_disk_path))
        # remove the file if it already exists
        if os.path.exists(self.on_disk_path):
            os.remove(self.on_disk_path)
        with pd.HDFStore(self.on_disk_path) as diskfile:
            diskfile.put("tdg", self.tdg.set_index(['chrom','pos']).sort_index(),format="table",data_columns=True)

        self.tdg = None
    
    def to_memory(self):
        if self.on_disk:
            self.tdg = pd.read_hdf(self.on_disk_path, 'tdg').reset_index().copy()
            self.on_disk = False
        else:
            warnings.warn("Object is already in memory")

    def get_info(self):
        print("CHARMtools Cell3D object v0.2")
        print("Cell name: {}".format(self.cellname))
        print("Resolution: {}".format(self.resolution))
        print("Features: {}".format(self.features))
        print(f"Chromatin3DStructure with {self.record_size} records")
        print(f"Object is stored on disk: {self.on_disk}")  
        print("Metadata:")
        for key, value in self.metadata.items():
            print(f"  {key}: {value}") 

    def _sparse_to_dense(self,genome_coord,sparse_df):
        chrom,start,end = _auto_genome_coord(genome_coord)
        positions = np.arange(start,end,self.resolution)
        dense_df = pd.DataFrame(positions,columns=['pos'])
        sparse_df['pos'] = sparse_df['pos'].astype(int)
        dense_df = pd.merge(dense_df, sparse_df, on='pos', how='left')
        dense_df['chrom'] = chrom

        columns_order = ['chrom', 'pos'] + [col for col in dense_df.columns if col not in ['chrom', 'pos']]
        dense_df = dense_df[columns_order]
        return dense_df

    def get_data(self,genome_coord = "",if_dense=False,rotate=False,rotate_x_angle=None,rotate_y_angle=None,rotate_z_angle=None):
        """
        Get the data of the Cell3D object.

        Parameters:
            genome_coord : str
                Genome coordinates in the format of "chrom:start-end" / "chrom" 
                if blank, return the whole data
            if_dense : bool
                Whether to convert sparse data to dense data
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
                if self.on_disk:
                    conditions = f"(chrom == {chrom})"
                    tdg_temp = pd.read_hdf(self.on_disk_path, 'tdg', where=conditions).reset_index().copy()
                else:
                    tdg_temp = self.tdg.query("chrom == @chrom").copy()
            else:
                start = int(start // self.resolution * self.resolution)
                end = int((end - 1) // self.resolution * self.resolution + self.resolution)
                
                if self.on_disk:
                    conditions = f"(chrom == {chrom}) & (pos >= {start}) & (pos < {end})"
                    tdg_temp = pd.read_hdf(self.on_disk_path, 'tdg', where=conditions).reset_index().copy()
                else:
                    tdg_temp = self.tdg.query("chrom == @chrom & pos >= @start & pos < @end").copy()
        else:
            if self.on_disk:
                tdg_temp = pd.read_hdf(self.on_disk_path, 'tdg').reset_index().copy()
            else:
                tdg_temp = self.tdg.copy()

        if rotate:
            rotated = Cell3D._point_cloud_rotation(tdg_temp[["x","y","z"]].values,
                                        x_angle=(rotate_x_angle is None and np.random.uniform(0, 2*np.pi) or rotate_x_angle),
                                        y_angle=(rotate_y_angle is None and np.random.uniform(0, 2*np.pi) or rotate_y_angle),
                                        z_angle=(rotate_z_angle is None and np.random.uniform(0, 2*np.pi) or rotate_z_angle)  
                                        )
            tdg_temp = tdg_temp.assign(x=rotated[:,0],y=rotated[:,1],z=rotated[:,2])

        if if_dense:
            if genome_coord == "":
                # if object without chrom_length, raise error
                if self.chrom_length is None:
                    raise ValueError("Running whole chromosome calculation with chrom_length is not available, |\
                                     please run add_chrom_length first")
                result_df = pd.DataFrame(columns = ['chrom','pos'])
                for index,row in self.chrom_length.iterrows():
                    chrom = row["chrom"]
                    length = row["size"]
                    positions = np.arange(0,length,self.resolution).astype(int)
                    temp_df = pd.DataFrame({'chrom':chrom,'pos':positions}) 
                    result_df = pd.concat([result_df,temp_df],ignore_index=True) 
                tdg_temp = pd.merge(result_df,tdg_temp,on=['chrom','pos'],how='left')    
            else:
                tdg_temp = self._sparse_to_dense(genome_coord,tdg_temp)
        
        return tdg_temp


    def get_data_slice(self,genome_coord="",random_slice_x = False,if_full = True,
                       slice_width = 3,if_rotate=False,
                       rotate_x_angle=None,rotate_y_angle=None,rotate_z_angle=None):
        """
        Method for insilico-GAM
        """
        tdg_temp = self.get_data(genome_coord=genome_coord,
                                 rotate = if_rotate,
                                 rotate_x_angle = rotate_x_angle,
                                 rotate_y_angle = rotate_y_angle,
                                 rotate_z_angle = rotate_z_angle,
                                 if_dense=True
                        )
        if if_full:
            tdg_temp["in_slice"] = tdg_temp["x"].apply(lambda x: -slice_width < x < slice_width)
            return tdg_temp
        else:
            if random_slice_x:
                xrange = (tdg_temp["x"].min(),tdg_temp["x"].max())
                slice_upper = np.random.uniform(xrange[0]+slice_width,xrange[1])
                slice_lower = slice_upper - slice_width
                tdg_temp = tdg_temp.query(f'x > @slice_lower & x < @slice_upper')
            else:
                slice_width = slice_width / 2 
                tdg_temp = tdg_temp.query("x > -@slice_width & x < @slice_width")
            
            return tdg_temp
        

    def get_data_cluster(self,eps=1.2,min_samples=5,cluster_name = "cluster",type ="normal",random_state=42):
        """
        Method for insilico-SPRITE(Cluster)
        """
        if cluster_name in self.tdg.columns:
            warnings.warn(f"Column {cluster_name} already exists, will be overwritten")

        self.calc_3D_cluster(eps=eps,min_samples=min_samples,random_state=random_state)

        return self.tdg[['chrom','pos','x','y','z',cluster_name]]
    
    def get_data_sphere(self,genome_coord="",radius=3,sample_frac=0.1):
        if genome_coord == "":
            # warning not fully implemented
            warnings.warn("Whole genome calculation is not fully implemented")
            if self.kdtree is None:
                raise ValueError("kdtree is not built, please run build_kdtree first")
            random_point = self.tdg.sample(1)[["x","y","z"]].values[0]
            indices = self.kdtree.query_ball_point(random_point,r=radius)
            points_in_ball = self.tdg.iloc[indices]
            return points_in_ball.sample(frac=sample_frac)
        else:
            tdg_temp = self.get_data(genome_coord=genome_coord,if_dense=True)
            kdtree_temp = cKDTree(tdg_temp[["x","y","z"]].values)
            random_point = tdg_temp.sample(1)[["x","y","z"]].values[0]
            indices = kdtree_temp.query_ball_point(random_point,r=radius)
            indices_sample = np.random.choice(indices,int(sample_frac*len(indices)),replace=False)
            tdg_temp["in_ball"] = tdg_temp.index.isin(indices_sample).astype(int)
            return tdg_temp
        
    def get_feature_vec(self, genome_coord, column_name):
        """
        INPUT:
            genome_coord: str, format like chrom:start-end or list/tuple of chrom,start,end. \|
                        whole chromosome is also acceptable. e.g. "chr1a:10000-20000" or ["chr1a",10000,20000] or "chr1a"
            column_name: str, the name of the column to extract the feature vector from. Or list of str.
        OUTPUT:
            feature_vec: np.array of the feature vector for the given region.
        """
        chrom, start, end = _auto_genome_coord(genome_coord)

        if start is None and self.chrom_length is None:
            raise ValueError("Running whole chromosome calculation with chrom_length is not available, |\
                                please run add_chrom_length first")

        if start is None:
            matsize = self.chrom_length.query("chrom == @chrom")["size"].values[0] // self.resolution + 1
            reconstruct_df = self.get_data(genome_coord)
        else:
            matsize = (end - start - 1) // self.resolution + 1
            reconstruct_df = self.get_data(genome_coord)
            reconstruct_df["pos"] = reconstruct_df["pos"] - start

        reconstruct_df["pos"] = reconstruct_df["pos"] // self.resolution

        if isinstance(column_name, str):
            if column_name not in reconstruct_df.columns:
                raise ValueError(f"Column '{column_name}' not found in the data")
            feature_vec = reconstruct_df[column_name].values
            for pos, value in zip(reconstruct_df["pos"], reconstruct_df[column_name]):
                feature_vec[pos] = value
            return feature_vec

        elif isinstance(column_name, list):
            for i in column_name:
                if i not in reconstruct_df.columns:
                    raise ValueError(f"Column '{i}' not found in the data")

            feature_vec = np.full((len(column_name), matsize), np.nan)
            for pos, values in zip(reconstruct_df["pos"], reconstruct_df[column_name].values):
                feature_vec[:, pos] = values
            return feature_vec
        else:
            raise ValueError("column_name should be a string or a list of strings")

    
    def subset(self,genome_coord="",query=None,in_place=False):
        """
        Subset the Cell3D object to a given genome coordinate.

        Parameters:
            genome_coord : str
                Genome coordinates in the format of "chrom:start-end"
        """
        if in_place:
            if query is not None:
                self.tdg = self.get_data(genome_coord).query(query).reset_index(drop=True)
            else:
                self.tdg = self.get_data(genome_coord).reset_index(drop=True)
            self.record_size = len(self.tdg)
            return None
        else:
            new_cell = copy.deepcopy(self)
            new_cell.cellname = self.cellname + "_subset"
            if query is not None:
                new_cell.tdg = self.get_data(genome_coord).query(query).reset_index(drop=True)
            else:
                new_cell.tdg = new_cell.get_data(genome_coord).reset_index(drop=True)

            new_cell.record_size = len(new_cell.tdg)
            return new_cell


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
    
    # Functions to load charm features
    def _load_bed_fragments(path, resolution, type = "allelic_resolved",peaks = None,flank = 0,keep_3prime = True):
        """
        path: file path of a tsv file containing chrom, start, end, allele, score, strand
        peaks: a pd.DataFrame containing chrom, start, end

        keep_3prime: bool, whether to trim the fragments to 3' end only
        """
        FRIP = None
        fragments = pd.read_csv(path, sep="\t", header=None)
        if fragments.shape[1] > 6:
            fragments = fragments.iloc[:, :6]   
        fragments.columns = ["chrom", "start", "end", "allele", "score", "strand"][:len(fragments.columns)]
        if keep_3prime:
            fragments = fragments.assign(start=lambda x: np.where(x["strand"] == "+", x["end"] - 1, x["start"]))
            fragments = fragments.assign(end=lambda x: np.where(x["strand"] == "+", x["end"], x["start"] + 1))
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
        elif type == "allelic_resolved_rev":
            fragments = fragments.query("chrom.str.contains('chr')").query('allele != "."')
            fragments = fragments.assign(chrom=np.where(fragments["allele"] == "1", fragments["chrom"] + "a", fragments["chrom"] + "b"))
        else:
            fragments = pd.concat([
                fragments.assign(chrom=lambda x: x["chrom"] + "a"),
                fragments.assign(chrom=lambda x: x["chrom"] + "b")
            ])
        fragments["pos"] = ((fragments["start"] + fragments["end"]) / 2) // resolution * resolution
        fragments["pos"] = fragments["pos"].astype(int)
        return fragments.groupby(["chrom", "pos"]).size().reset_index().rename(columns={0: "count"}),FRIP

    # def _load_CpG(CpG_path):
    #     """
    #     deprecated , use add_bedGraph_data instead
    #     """
    #     CpG = pd.read_csv(CpG_path, header=None, sep="\t")
    #     CpG.columns = ["chrom", "pos", "CpG"]
    #     if CpG["chrom"].str.startswith("chr").sum() == 0:
    #         CpG["chrom"] = "chr" + CpG["chrom"]
    #     CpG["pos"] = CpG["pos"].astype(int)
    #     return pd.concat([
    #         CpG.assign(chrom=lambda x: x["chrom"] + "a"),
    #         CpG.assign(chrom=lambda x: x["chrom"] + "b")
    #     ])    

    def add_bed_data(self, path, column_name, resolution=None,type="allelic_resolved",peaks = None,flank=0,keep_3prime=True):

        if resolution is None:
            resolution = self.resolution
        if column_name in self.tdg.columns:
            warnings.warn("Column {} already exists, will be overwritten".format(column_name))

        if self.on_disk:
            self.to_memory()
        
        fragments,FRIP = Cell3D._load_bed_fragments(path, resolution,type,peaks=peaks,flank=flank,keep_3prime=keep_3prime)
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
        
        if self.on_disk:
            self.to_memory()

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
        merged_bed = positions_bed.intersect(bedgraph_bed,wa=True,wb=True,nonamecheck=True)
        newcol = merged_bed.to_dataframe()[["chrom","start","end","thickStart"]].groupby(["chrom","start","end"]).sum().reset_index()
        newcol.columns = ["chrom","start","end",column_name]
        newcol["pos"] = newcol["start"]
        #self.tdg =
        temp = pd.merge(self.tdg,newcol[["chrom","pos",column_name]],on=["chrom","pos"],how="left")
        temp[column_name] = temp[column_name].fillna(0)
        self.tdg = temp
        self.features.append(column_name)
        
    def add_RNA_data(self, rnag1, rnag2, genes,cellname = None, column_name=None,type="tss"):
        """
        Add RNA data to the Cell3D object.

        Parameters:
            genome1_rnamatrix : np.array
                RNA matrix of genome 1 (chr[]b)
            genome2_rnamatrix : np.array
                RNA matrix of genome 2 (chr[]a)
            column_name : str
                Name of the column to be added

        Returns:
            None
        """
        if self.on_disk:
            self.to_memory()

        if cellname is None:
            cellname = self.cellname
        if column_name is None:
            column_name = "UMI"

        resolution = self.resolution
        genes = genes.copy()

        if type == "tss":
            genes.columns = ['chrom','start','end','id','gene','strand']
            genes['tss'] = genes.apply(lambda x: x['start'] if x['strand'] == '+' else x['end'], axis=1)
            genes['tss'] = genes['tss']//resolution*resolution
            genes = genes[["chrom","tss","gene"]]

            rnag1m = rnag1[['gene',cellname]]
            rnag2m = rnag2[['gene',cellname]]

            rnag1m = rnag1m.merge(genes, left_on='gene', right_on='gene', how='left')[['chrom','tss',cellname]]
            rnag2m = rnag2m.merge(genes, left_on='gene', right_on='gene', how='left')[['chrom','tss',cellname]]
            rnag1m['chrom'] = rnag1m['chrom'].apply(lambda x: x + "b")
            rnag2m['chrom'] = rnag2m['chrom'].apply(lambda x: x + "a")

            mergedf = pd.concat([rnag1m,rnag2m])
            mergedf.columns = ['chrom','pos',column_name]
            mergedf = mergedf.groupby(['chrom','pos']).sum().reset_index()

            df = pd.merge(self.tdg,mergedf,on=['chrom','pos'],how='left')

        if type == "gene":
            # the gene options assume that gene table is already in binned foramt, you can use create_bins_genetable to create the binned gene table
            df = self.tdg.copy()
            rnag1m = rnag1[["gene",cellname]].copy()
            rnag2m = rnag2[["gene",cellname]].copy()

            rnag1m.columns = ["gene","UMIs"]
            rnag2m.columns = ["gene","UMIs"]
            rnag1m = genes[["chrom","start","end","gene"]].merge(rnag1m.query('UMIs > 0'),on="gene",how="inner").drop("gene",axis=1).groupby(["chrom","start","end"]).sum().reset_index()
            rnag2m = genes[["chrom","start","end","gene"]].merge(rnag2m.query('UMIs > 0'),on="gene",how="inner").drop("gene",axis=1).groupby(["chrom","start","end"]).sum().reset_index()
            
            rnag1m['chrom'] = rnag1m['chrom'] + 'b'
            rnag2m['chrom'] = rnag2m['chrom'] + 'a'
            rna = pd.concat([rnag1m,rnag2m]).drop('end',axis=1)
            rna.columns = ['chrom','pos',column_name]
            df = df.merge(rna,on=["chrom","pos"],how="left")
        df[column_name] = df[column_name].fillna(0)

        self.tdg = df
        self.features.append(column_name)
        return None
    
    def add_knn_density(self, k):
        """
        use mean distance of k nearest neighbors as density
        """
        if self.on_disk:
            self.to_memory()
        if self.kdtree is None:
            self.build_kdtree()

        self.features.append("knn_density_"+str(k))
        densities = []
        for point in self.tdg[["x","y","z"]].values:
            dist,ind = self.kdtree.query(point,k=k+1)
            density = dist.mean()
            densities.append(density)
        self.tdg["knn_density_"+str(k)] = np.array(densities)
        return None

    def add_feature_in_radius(self, feature, radius,type = "mean",if_self = True,if_rank = False,add_self=False):
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
        if self.on_disk:
            self.to_memory()
        if self.kdtree is None:
            self.build_kdtree()
        from scipy.stats import rankdata
        indices_list = self.kdtree.query_ball_tree(self.kdtree, r=radius)
        if not if_self:
            #indices_list = [indices[1:] for indices in indices_list]
            raise ValueError("if_self parameter is not implemented yet")
        if type =="mean":
            avgs = [self.tdg[feature].iloc[indices].mean() for indices in indices_list]
        elif type =="sum":
            avgs = [self.tdg[feature].iloc[indices].sum() for indices in indices_list]
        elif type == "countall": 
            avgs = [len(indices) for indices in indices_list]
        else:
            raise ValueError("type should be mean or sum or countall")
        if if_rank:
            avgs = rankdata(avgs,nan_policy="omit") / sum(np.isfinite(avgs))
        
        # avgs nan to 0
        avgs = np.array(avgs)
        avgs[np.isnan(avgs)] = 0

        self.tdg[feature + "_" + type + "_in_radius_" + str(radius)] = avgs
        self.features.append(feature + "_" + type + "_in_radius" + str(radius))
        if add_self:
            self.tdg[feature + "_" + type + "_in_radius_" + str(radius)] = self.tdg[feature + "_" + type + "_in_radius_" + str(radius)] + self.tdg[feature]
        return None

    def calc_intermingling(self,radius=3):
        """
        Calculate the intermingling score of the given region.
        intermingling_ratio = number of different chromosome points in radius / number of points in radius
        intermingling_index = −∑pilnpi, where pi denoted the fraction of nearby particles from chromosome i)

        Parameters:
            radius : float
        """
        if self.on_disk:
            self.to_memory()
        if self.kdtree is None:
            self.build_kdtree()
        indices_list = self.kdtree.query_ball_tree(self.kdtree, r=radius)
        intermingling_ratios = []
        intermingling_indexes = []
        for indices in indices_list:
            temp_df = self.tdg.iloc[indices]
            # for intermingling ratio
            self_chrom = temp_df["chrom"].iloc[0]
            intermingling_ratios.append((temp_df["chrom"] != self_chrom).sum() / len(temp_df))
            # for intermingling index
            tempindxes = []
            for chroms in temp_df["chrom"].unique():
                pi = (temp_df["chrom"] == chroms).sum() / len(temp_df)
                tempindxes.append(pi * np.log(pi))
            intermingling_indexes.append(-sum(tempindxes))
        self.tdg["intermingling_ratio"] = intermingling_ratios
        self.tdg["intermingling_index"] = intermingling_indexes
        self.features.append("intermingling_ratio") 
        self.features.append("intermingling_index")
        return None



    def add_chrom_length(self,chrom_length_path):
        chrom_length = pd.read_csv(chrom_length_path,sep="\t",header=None)
        chrom_length.columns = ["chrom","size"]
        chrom_length["chrom"] = chrom_length["chrom"].str.replace("\(pat\)", "a", regex=True)
        chrom_length["chrom"] = chrom_length["chrom"].str.replace("\(mat\)", "b", regex=True)
        chrom_length["chrom"] = chrom_length["chrom"].str.replace("pat", "a", regex=True)
        chrom_length["chrom"] = chrom_length["chrom"].str.replace("mat", "b", regex=True)

        valid_chroms = self.tdg["chrom"].unique()
        chrom_length = chrom_length.query("chrom in @valid_chroms")

        self.features.append("chrom_length")
        self.chrom_length = chrom_length

    # Functions for output sub-region matrix
    
    def calc_distance_matrix(self,genome_coord,obsexp=False):
        """
        Calculate the distance matrix of a given region.
        INPUT:
            genome_coord: str, format like chrom:start-end or list/tuple of chrom,start,end. \|
                          whole chromosome is also acceptable. e.g. "chr1a:10000-20000" or ["chr1a",10000,20000] or "chr1a
        OUTPUT:
            distance_matrix: symettrical distance matrix of the given region, np.array of shape (n,n)
        """
        temp_df = self.get_data(genome_coord,if_dense=True)
        mat = scipy.spatial.distance.squareform(
            scipy.spatial.distance.pdist(temp_df[["x","y","z"]].values)
        )
        if obsexp:

            def _obs_exp(mat, expected):
                mat_norm = np.zeros_like(mat)
                n = mat.shape[0]
                for k in range(1, n):
                    mat_norm[np.arange(k, n), np.arange(0, n - k)] = mat[np.arange(k, n), np.arange(0, n - k)] / expected[k]
                mat_norm = mat_norm + mat_norm.T - np.diag(np.diag(mat_norm))
                return mat_norm
            
            if self.expected is None:
                raise ValueError("Expected vector is not available, please run calc_expected first")
            chrom = genome_coord.split(":")[0]
            mat = _obs_exp(mat, self.expected[chrom])

        return mat

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

        mat = self.calc_distance_matrix(genome_coord) < distance_threshold
        feature_mat = mat.copy()

        feature_mat[feature_vec,:] = 0
        feature_mat[:,feature_vec] = 0
        return feature_mat,feature_vec

    def calc_3D_cluster(self, query="",eps=1.2,min_samples=5,cluster_name="cluster",type="normal",random_seed=42):
        """
        Clustering given points in 3D space using DBSCAN.
        """
        if self.on_disk:
            self.to_memory()
        if query != "":
            data = self.tdg.query(query).copy()
        else:
            data = self.tdg.copy()
        if type == "random":
            data = self.tdg.sample(data.shape[0],random_state=random_seed).copy()

        cluster = DBSCAN(eps=eps, min_samples=min_samples).fit(data[["x","y","z"]])

        # if cluster_name already in data
        if cluster_name in data.columns:
            print(f"Warning: {cluster_name} already in data, remove the old one")
            self.tdg.drop(columns=[cluster_name],inplace=True)

        data.loc[:,cluster_name] = cluster.labels_
        self.tdg = pd.concat([self.tdg,data.loc[:,cluster_name]],axis=1)
        # cluster name column na to -1
        self.tdg[cluster_name] = self.tdg[cluster_name].fillna(-1)
        self.tdg[cluster_name] = self.tdg[cluster_name].astype(int).astype(str)
        self.features.append(cluster_name)
        return None

    def point_to_3Dcluster(self,feature,new_column_name):
        """
        For each cluster, calculate the bounding ball and assign points inside the ball to the cluster.
        """
        from miniball import get_bounding_ball

        if self.on_disk:
            self.to_memory()
        if self.kdtree is None:
            self.build_kdtree()

        data = self.tdg.copy()
        if new_column_name in data.columns:
            warnings.warn(f"Column {new_column_name} already exists, will be overwritten")
        if feature not in data.columns:
            raise ValueError(f"Feature {feature} not found in the data")
        
        data[new_column_name] = "-1"
        valid_clusters = data[data[feature] != "-1"][feature].unique()

        radius_list = []
        for cluster in valid_clusters:
            points = data[data[feature] == cluster][['x','y','z']].values
            center,radius_squared = get_bounding_ball(points)
            radius = np.sqrt(radius_squared)
            data.loc[self.kdtree.query_ball_point(center, radius),new_column_name] = cluster
            radius_list.append(radius)
        self.extra[new_column_name + "_radius"] = radius_list
        self.tdg = data
        return None


    def calc_scABC_pred_gene(self,tss_genome_coord,flank = 2000000,expression_key = "UMIs_tss",
                             activity_keys = ["atac_sum_in_radius_2","ct_sum_in_radius_2"],distance_type="3d"):
        """
        Calculate the predicted gene expression by scABC model.
        """
        resolution = self.resolution
        chrom,start,end = _auto_genome_coord(tss_genome_coord)
        region = f"{chrom}:{start-flank}-{end+flank}"
        temp_tdg = self.get_data(region,if_dense=True)
        tss_row = temp_tdg.loc[flank//resolution,:]
        temp_tdg["activity"] = temp_tdg[activity_keys].prod(axis=1)
        tss_location = tss_row[['x','y','z']].values

        if distance_type == "3d":
            temp_tdg["distance"] = np.sqrt((temp_tdg['x'] - tss_location[0])**2 + (temp_tdg['y'] - tss_location[1])**2 + (temp_tdg['z'] - tss_location[2])**2)
        elif distance_type == "2d":
            temp_tdg["distance"] = np.abs(temp_tdg['pos'] - tss_row['pos'])
        temp_tdg['abc'] = temp_tdg['activity'] * (1/temp_tdg['distance']) 

        return tss_row[expression_key], temp_tdg.abc.values
    
    def calc_expected(self,n_diag=None):
        """
        calculate expected distance for each chromosome, for calculation efficiency, you can only calculate the first n_diag diagonals
        """
        if self.on_disk:
            self.to_memory()
        if self.chrom_length is None:
            raise ValueError("Chromosome length is needed for calculating expected, please run add_chrom_length first")
        self.expected = {}
        for chrom,length in self.chrom_length.values:
            genome_coord = chrom + ":0-" + str(length)
            mat = self.calc_distance_matrix(genome_coord)
            means = []
            if n_diag is None:
                n_diag = mat.shape[0]
            for i in range(n_diag):
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    means.append(np.nanmean(np.diag(mat,i)))
                    self.expected[chrom] = np.array(means)

        self.features.append("expected")
        return None
        

    # data visualize
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
        if self.on_disk:
            self.to_memory()
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
    
    def plot3D(self, color_by,genome_coord=None, smooth=True, smoothness=2,
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
        resolution = self.resolution
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        from scipy.interpolate import splprep, splev
        import plotly.express as px

        # Creating a Plotly figure
        fig = make_subplots(rows=1, cols=1, specs=[[{'type': 'scatter3d'}]])
        # if cell is not a pd.DataFrame
        dataframe = self.get_data(genome_coord)
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

    def _roated_hic_mat(matrix, h):
        from scipy.ndimage import rotate
        matrix = np.nan_to_num(matrix)
        rotated =  rotate(matrix, 45, reshape=True)
        return rotated[rotated.shape[0]//2-h:rotated.shape[0]//2, :]

    def plotDistanceMatrix(self,genome_coord,type="normal",h=None,**kwargs):
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
        cmap = kwargs.get("cmap","YlOrRd_r")
        vmin = kwargs.get("vmin",None)
        vmax = kwargs.get("vmax",None)

        if type == "normal":
            if vmin and vmax:
                plt.imshow(mat,cmap=cmap,vmin=vmin,vmax=vmax,extent=[start,end,end,start])
            else:
                plt.imshow(mat,cmap=cmap,extent=[start,end,end,start])
        elif type == "rotate":
            rotated = Cell3D._roated_hic_mat(mat,h)
            if vmin and vmax:
                plt.imshow(rotated,cmap=cmap,vmin=vmin,vmax=vmax)
            else:
                plt.imshow(rotated,cmap=cmap)

        plt.colorbar()
        plt.show()
        return None
    
    def expandStructure(self, expand_factor=3, inplace=False):
        """
        Center the structure ,calculate center for each chromosome and expand the structure by multiple dist(chromcenter to center)
        Parameters:
            expand_factor: int, default 3
        """
        if self.on_disk:
            self.to_memory()
        tdg = self.tdg.copy()
        center = tdg[["x", "y", "z"]].mean()
        tdg[["x", "y", "z"]] = (tdg[["x", "y", "z"]] - center)
        
        chrom_centers = tdg.groupby("chrom")[["x", "y", "z"]].transform("mean")
        tdg["chrom_center_x"] = chrom_centers["x"]
        tdg["chrom_center_y"] = chrom_centers["y"]
        tdg["chrom_center_z"] = chrom_centers["z"]
        
        tdg[["x", "y", "z"]] = tdg[["x", "y", "z"]] + chrom_centers * expand_factor
        # drop chrom_center_x,y,z columns
        tdg = tdg.drop(columns=["chrom_center_x", "chrom_center_y", "chrom_center_z"], axis=1)
        
        if inplace:
            self.tdg = tdg
        else:
            cellnew = copy.deepcopy(self)
            cellnew.tdg = tdg
            cellnew.cellname = self.cellname + "_expanded"
            return cellnew




    # analysis
    def calc_radial_position(self,key="radial_position",if_rank = False,if_norm_max = False, if_norm_mean = True):
        """
        Calculate the radial position of each point.
        """
        if self.on_disk:
            self.to_memory()
        data = self.get_data()
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
    
    def feature_radial_distribution(self, feature,random = False,random_seed = 42,if_normalize_avg = False,if_rank = False):
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
        tdg = self.get_data()
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

        if self.on_disk:
            self.to_memory()
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

# analysis & helper functions
        
def calculate_RMSD(*dataframes):
    import rmsd
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
    pass

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

def create_bins_genetable(df,resolution=5000):
    # Prepare output columns
    df = df.copy()
    df['start'] = df['start'].astype(int) // resolution * resolution
    df['end'] = df['end'].astype(int) // resolution * resolution + resolution
    chroms, ids, names, strands, starts_all, ends_all = [], [], [], [], [], []

    # Process each row
    for _, row in df.iterrows():
        starts = np.arange(row['start'], row['end'], resolution)
        ends = np.clip(starts + resolution, None, row['end'])
        
        # Append results
        chroms.extend([row['chrom']] * len(starts))
        ids.extend([row['id']] * len(starts))
        names.extend([row['gene']] * len(starts))
        strands.extend([row['strand']] * len(starts))
        starts_all.extend(starts)
        ends_all.extend(ends)
    
    # Create a new DataFrame
    new_df = pd.DataFrame({
        'chrom': chroms,
        'start': starts_all,
        'end': ends_all,
        'id': ids,
        'gene': names,
        'strand': strands
    })

    return new_df
