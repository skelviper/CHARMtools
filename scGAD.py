import pybedtools
from concurrent.futures import ProcessPoolExecutor
import tqdm
import pandas as pd

from . import CHARMio


def _calc_scgad_raw_cell(pairs_path, genes_df):
    pairs = CHARMio.parse_pairs(pairs_path)
    pairs = pairs.query("chr1 == chr2")
    pairs = pairs[["chr1","pos1","pos2"]]
    pairs = pairs.rename(columns={"chr1":"chr","pos1":"start","pos2":"end"}).reset_index(drop=True)
    pairs_bed = pybedtools.BedTool.from_dataframe(pairs)
    genes_bed = pybedtools.BedTool.from_dataframe(genes_df)
    pairs_bed = pairs_bed.intersect(genes_bed,wb=True,f=1)
    pairs_bed = pairs_bed.to_dataframe()
    sc_raw_gad = pairs_bed.groupby("thickEnd").size().reset_index()
    sc_raw_gad.columns = ["gene","count"]
    sc_raw_gad = sc_raw_gad.set_index("gene")
    return sc_raw_gad

def calc_scgad(cellnames,pairs_paths, genes_df,num_cores=20):
    """
    Calculate scGAD for a list of cellnames and pairs_paths. Return a scGAD matrix for interested genes.
    Parameters:
    cellnames: list
        List of cellnames.
    pairs_paths: list
        List of paths to pairs files.
    genes_df: pd.DataFrame
        DataFrame of genes.contains columns: ["chr", "start", "end", "id","name","strand"]
    """
    with ProcessPoolExecutor(num_cores) as executor:
        res = []
        for i in tqdm.tqdm(executor.map(_calc_scgad_raw_cell, pairs_paths, [genes_df]*len(pairs_paths)), total=len(pairs_paths)):
            res.append(i)
        scGADraw = pd.concat(res, axis=1)
    scGADraw.columns = cellnames
    scGADraw.fillna(0, inplace=True)

    pybedtools.helpers.cleanup(remove_all=True)

    # R_ij~
    df_normalized = scGADraw.div(scGADraw.sum(axis=0), axis=1)
    # mean_i
    mean_i = df_normalized.mean(axis=1)
    # scGAD_global
    scGAD_global = (df_normalized.sub(mean_i, axis=0).pow(2))
    scGAD_global = scGAD_global / scGAD_global.sum(axis=0) * df_normalized.shape[1]
    # sqrt 
    scGAD_global = scGAD_global.pow(0.5)

    return scGAD_global