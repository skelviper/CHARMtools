import pandas as pd
import numpy as np
from sklearn.preprocessing import LabelEncoder
import warnings
from scipy.spatial import KDTree
import rmsd


class Cell3D:
    def __init__(self, cellname,tdg_path, resolution):
        self.cellname = cellname
        self.tdg = self._load_tdg(tdg_path)
        self.resolution = resolution
        self.features = []
        self.kdtree = KDTree(self.tdg[["x", "y", "z"]].values)

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
        tdg["chrom"] = tdg["chrom"].str.replace("\(mat\)", "a", regex=True)
        tdg["chrom"] = tdg["chrom"].str.replace("\(pat\)", "b", regex=True)
        LE = LabelEncoder()
        tdg.chrom = pd.Categorical(tdg.chrom)
        tdg['chrom_code'] = LE.fit_transform(tdg['chrom'])
        return tdg

    def _load_bed_fragments(path, resolution):
        fragments = pd.read_csv(path, sep="\t", header=None)
        fragments.columns = ["chrom", "start", "end", "allele", "score", "strand"]
        fragments = fragments.query("chrom.str.contains('chr')").query('allele != "."')
        fragments = fragments.assign(chrom=np.where(fragments["allele"] == "0", fragments["chrom"] + "a", fragments["chrom"] + "b"))
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

    def add_bed_data(self, path, column_name, resolution=None):
        if resolution is None:
            resolution = self.resolution
        if column_name in self.tdg.columns:
            warnings.warn("Column {} already exists, will be overwritten".format(column_name))
        fragments = Cell3D._load_bed_fragments(path, resolution)
        self.tdg = pd.merge(self.tdg, fragments, on=["chrom", "pos"], how="left")
        self.tdg[column_name] = self.tdg["count"].fillna(0)
        self.tdg = self.tdg.drop(columns=["count"], axis=1)
        self.features.append(column_name)

    def add_CpG_data(self, path):
        CpG = Cell3D._load_CpG(path)
        self.features.append("CpG")
        self.tdg = pd.merge(self.tdg, CpG, on=["chrom", "pos"], how="left").dropna()

    def add_density(self, radius):
        self.features.append("density_"+str(radius))
        densities = []

        for point in self.tdg[["x","y","z"]].values:
            count = self.kdtree.query_ball_point(point,radius)
            density = len(count) / ((4/3) * np.pi * radius **3)
            densities.append(density)

        self.tdg["density_"+str(radius)] = np.array(densities)
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

        return tdg_temp
    
    # analysis
    #calc_rmsd = lambda x,y: np.sqrt(np.mean((x-y)**2))
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

