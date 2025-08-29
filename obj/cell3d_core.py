# Cell3D Core Module - Basic functionality and data management
import warnings
import pandas as pd
import numpy as np
import copy
import os
from sklearn.preprocessing import LabelEncoder
from ..utils.helper import auto_genome_coord

class Cell3DCore:
    """Core functionality for Cell3D object including initialization and basic data operations"""
    
    def __init__(self, cellname, resolution, tdg_path=None, pairs_path=None, type="3dg", on_disk=False, on_disk_path=None):
        self.cellname = cellname
        if type == "3dg":
            self.tdg = self._load_tdg(tdg_path)
        elif type == "pairs":
            self.tdg = self._fdg_from_pairs(pairs_path)
        self.record_size = self.tdg.shape[0]
        self.resolution = resolution
        self.features = []
        self.kdtree = None
        self.chrom_length = None
        self.on_disk = on_disk
        self.on_disk_path = on_disk_path
        self.metadata = {}
        self.extra = {}
        self.expected = None

        if on_disk:
            self.to_disk(on_disk_path)

    def __repr__(self):
        self.get_info()
        return ""

    def copy(self):
        return copy.deepcopy(self)

    def _load_tdg(self, tdg):
        if isinstance(tdg, pd.DataFrame):
            tdg = tdg.copy()
        elif isinstance(tdg, str):
            import os
            if not os.path.exists(tdg):
                raise FileNotFoundError(f"TDG file not found: {tdg}")
            try:
                tdg = pd.read_csv(tdg, sep="\t", header=None, comment="#")
                if tdg.empty:
                    raise ValueError(f"TDG file is empty: {tdg}")
            except Exception as e:
                raise ValueError(f"Error reading TDG file {tdg}: {str(e)}")
        else:
            raise ValueError("tdg should be a pandas.DataFrame or a string")
        
        tdg.columns = ["chrom", "pos", "x", "y", "z"]
        tdg["chrom"] = tdg["chrom"].astype("category")
        tdg["pos"] = tdg["pos"].astype(int)
        tdg["x"] = tdg["x"].astype(np.float32)
        tdg["y"] = tdg["y"].astype(np.float32)
        tdg["z"] = tdg["z"].astype(np.float32)

        tdg["chrom"] = tdg["chrom"].str.replace(r"\(pat\)", "a", regex=True)  
        tdg["chrom"] = tdg["chrom"].str.replace(r"\(mat\)", "b", regex=True)  
        LE = LabelEncoder()
        tdg.chrom = pd.Categorical(tdg.chrom)
        tdg['chrom_code'] = LE.fit_transform(tdg['chrom'])
        
        # Update record_size after loading data
        self.record_size = len(tdg)
        
        return tdg

    def to_disk(self, on_disk_path=None):
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

        return ""

    def get_data(self, genome_coord="", if_dense=False, rotate=False, rotate_x_angle=None, rotate_y_angle=None, rotate_z_angle=None):
        """Get the data of the Cell3D object."""
        if genome_coord != "":
            chrom, start, end = auto_genome_coord(genome_coord)
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
            from .cell3d_data import Cell3DData
            rotated = Cell3DData._point_cloud_rotation(tdg_temp[["x","y","z"]].values,
                                        x_angle=(rotate_x_angle is None and np.random.uniform(0, 2*np.pi) or rotate_x_angle),
                                        y_angle=(rotate_y_angle is None and np.random.uniform(0, 2*np.pi) or rotate_y_angle),
                                        z_angle=(rotate_z_angle is None and np.random.uniform(0, 2*np.pi) or rotate_z_angle)  
                                        )
            tdg_temp = tdg_temp.assign(x=rotated[:,0],y=rotated[:,1],z=rotated[:,2])

        if if_dense:
            if genome_coord == "":
                if self.chrom_length is None:
                    raise ValueError("Running whole chromosome calculation with chrom_length is not available, please run add_chrom_length first")
                result_df = pd.DataFrame(columns = ['chrom','pos'])
                for index,row in self.chrom_length.iterrows():
                    chrom = row["chrom"]
                    length = row["size"]
                    positions = np.arange(0,length,self.resolution).astype(int)
                    temp_df = pd.DataFrame({'chrom':chrom,'pos':positions}) 
                    result_df = pd.concat([result_df,temp_df],ignore_index=True) 
                tdg_temp = pd.merge(result_df,tdg_temp,on=['chrom','pos'],how='left')    
            else:
                from .cell3d_data import Cell3DData
                tdg_temp = Cell3DData._sparse_to_dense(self, genome_coord, tdg_temp)
        
        return tdg_temp

    def subset(self, genome_coord="", query=None, in_place=False):
        """Subset the Cell3D object to a given genome coordinate."""
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

    def _fdg_from_pairs(self, pairs_path):
        """Placeholder for FDG from pairs - to be implemented in dev module"""
        raise NotImplementedError("FDG from pairs functionality moved to dev module")