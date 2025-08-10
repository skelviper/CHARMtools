# Cell3D Features Module - Feature addition and management
import warnings
import pandas as pd
import numpy as np
import pybedtools

class Cell3DFeatures:
    """
    Feature addition and management for Cell3D objects
    给结构数据添加来自其他组学的信息。
    """
    
    @staticmethod
    def _load_bed_fragments(path, resolution, type="allelic_resolved", peaks=None, flank=0, keep_3prime=True, merge_overlap=True):
        """Load bed fragments with various processing options"""
        FRIP = None
        fragments = pd.read_csv(path, sep="\t", header=None)
        if fragments.shape[1] > 6:
            fragments = fragments.iloc[:, :6]   
        fragments.columns = ["chrom", "start", "end", "allele", "score", "strand"][:len(fragments.columns)]
        
        if merge_overlap:
            def merge_intervals(group):
                sorted_group = group.sort_values('start')
                intervals = sorted_group[['start', 'end']].values
                
                if not intervals.size:
                    return pd.DataFrame()
                
                merged = []
                current_start, current_end = intervals[0]
                
                for start, end in intervals[1:]:
                    if start <= current_end: 
                        current_end = max(current_end, end)
                    else:
                        merged.append((current_start, current_end))
                        current_start, current_end = start, end
                
                merged.append((current_start, current_end))
                result_df = pd.DataFrame(merged, columns=['start', 'end'])

                result_df['chrom'] = group.name[0]
                result_df['allele'] = group.name[1]
                result_df['strand'] = group.name[2]
                
                return result_df[['chrom', 'start', 'end', 'allele', 'strand']]

            if fragments.shape[1] > 3:
                merged_df = (
                    fragments.groupby(['chrom', 'allele', 'strand'], group_keys=False)
                    .apply(merge_intervals)
                    .reset_index(drop=True)
                    .sort_values(['chrom', 'start'])
                )
                merged_df["score"] = "."
                fragments = merged_df[["chrom", "start", "end", "allele", "score", "strand"]]
            else:
                fragments = fragments
            
        if keep_3prime:
            fragments = fragments.assign(start=lambda x: np.where(x["strand"] == "+", x["end"] - 1, x["start"]))
            fragments = fragments.assign(end=lambda x: np.where(x["strand"] == "+", x["end"], x["start"] + 1))
            
        if peaks is not None:
            peaks = peaks.copy()
            peaks.columns = ["chrom", "start", "end"] 
            peaks["start"] = peaks["start"] - flank
            peaks["start"] = peaks["start"].clip(lower=0)
            peaks["end"] = peaks["end"] + flank
            fragments_bed = pybedtools.BedTool.from_dataframe(fragments)
            peak_bed = pybedtools.BedTool.from_dataframe(peaks)
            intersect = peak_bed.intersect(fragments_bed,wa=True,wb=True,nonamecheck=True)
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
        return fragments.groupby(["chrom", "pos"]).size().reset_index().rename(columns={0: "count"}), FRIP

    def add_bed_data(self, path, column_name, resolution=None, type="allelic_resolved", peaks=None, flank=0, keep_3prime=True):
        """Add features from bed file to the Cell3D object"""
        if resolution is None:
            resolution = self.resolution
        if column_name in self.tdg.columns:
            warnings.warn("Column {} already exists, will be overwritten".format(column_name))

        if self.on_disk:
            self.to_memory()
        
        fragments, FRIP = self._load_bed_fragments(path, resolution, type, peaks=peaks, flank=flank, keep_3prime=keep_3prime)
        self.tdg = pd.merge(self.tdg, fragments, on=["chrom", "pos"], how="left")
        self.tdg[column_name] = self.tdg["count"].fillna(0)
        self.tdg = self.tdg.drop(columns=["count"], axis=1)
        self.features.append(column_name)
        if FRIP is not None:
            self.metadata[column_name + "_FRIP"] = FRIP

    def add_bedGraph_data(self, path, column_name, resolution=None, type="allelic_resolved"):
        """Add features from bedGraph file to the Cell3D object"""
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
        temp = pd.merge(self.tdg,newcol[["chrom","pos",column_name]],on=["chrom","pos"],how="left")
        temp[column_name] = temp[column_name].fillna(0)
        self.tdg = temp
        self.features.append(column_name)
        
    def add_RNA_data(self, rnag1, rnag2, genes, cellname=None, column_name=None, type="tss"):
        """Add RNA data to the Cell3D object"""
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

    def add_chrom_length(self, chrom_length_path):
        """Add chromosome length information"""
        chrom_length = pd.read_csv(chrom_length_path,sep="\t",header=None)
        chrom_length.columns = ["chrom","size"]
        chrom_length["chrom"] = chrom_length["chrom"].str.replace(r"\(pat\)", "a", regex=True)  
        chrom_length["chrom"] = chrom_length["chrom"].str.replace(r"\(mat\)", "b", regex=True)  
        chrom_length["chrom"] = chrom_length["chrom"].str.replace("pat", "a", regex=True)
        chrom_length["chrom"] = chrom_length["chrom"].str.replace("mat", "b", regex=True)
        
        self.features.append("chrom_length")
        self.chrom_length = chrom_length

    def calc_expected(self, n_diag=None):
        """Calculate expected distance for each chromosome"""
        if self.on_disk:
            self.to_memory()
        if self.chrom_length is None:
            raise ValueError("Chromosome length is needed for calculating expected, please run add_chrom_length first")
        self.expected = {}
        
        for chrom, length in self.chrom_length.values:
            temp_df = self.get_data(genome_coord=chrom, if_dense=True).copy()
            means = []
            if n_diag is None:
                n_diag = temp_df.shape[0]

            for n in range(n_diag):
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    means.append(np.nanmean(np.linalg.norm(temp_df[['x','y','z']].values - temp_df[['x','y','z']].shift(-n).values, axis=1)))
            self.expected[chrom] = np.array(means)

        self.features.append("expected")
        return None