import numpy as np
import pandas as pd
import tqdm


def peak_enrichment(frag_df,peak_df,chrom_length,resolution=100,flank=2000,npeaks=2000):
    """
    Calculate the enrichment of fragments around peaks
    Parameters
        frag_df: a bed-like dataframe with at least 3 columns representing chrom, start, end
        peak_df: similar to frag_df
        chrom_length: a dataframe with 2 columns representing chrom, length
        resolution: resolution of the data
        flank: the flanking region around the peak
        npeaks: number of peaks to sample
    Returns
        res: a numpy array with shape (npeaks, (flank*2)//resolution+1)
    """
    temp_frags = frag_df.copy().iloc[:,0:3]
    temp_frags.columns = ["chrom","start","end"]
    temp_peaks = peak_df.copy().iloc[:,0:3]
    temp_peaks.columns = ["chrom","start","end"]
    chrom_length = chrom_length.copy()
    chrom_length.columns = ["chrom","length"]

    if temp_peaks.shape[0] >= npeaks:
        temp_peaks = temp_peaks.sample(npeaks)
        
    temp_peaks["start"] = ((temp_peaks["start"] + temp_peaks["end"]) / 2 // resolution * resolution).astype(int)
    temp_peaks["end"] = (temp_peaks["start"] + resolution).astype(int)

    temp_frags["start"] = (temp_frags["start"]//resolution*resolution).astype(int)
    temp_frags["end"] = temp_frags["end"]//resolution*resolution
    temp_frags = temp_frags.groupby(["chrom","start","end"]).size().reset_index().rename(columns={0:"count"})

    # join with chrom_length, dense the fragments
    chrom_bins = {chrom: np.arange(0, length, resolution) for chrom, length in zip(chrom_length["chrom"], chrom_length["length"])}
    result_list = []
    for chrom, positions in chrom_bins.items():
        temp_df = pd.DataFrame({"chrom": chrom, "start": positions, "end": positions + resolution})
        temp_merged = temp_df.merge(temp_frags[temp_frags["chrom"] == chrom], on=["chrom", "start", "end"], how="left").fillna(0)
        result_list.append(temp_merged)

    temp_frags = pd.concat(result_list, ignore_index=True)
    temp_frags.set_index(["chrom", "start", "end"], inplace=True)

    res = []
    for peak_row in tqdm.tqdm(temp_peaks.iterrows()):
        peak_row = peak_row[1].values
        try:
            # generate list of regions for each peak
            min_start = peak_row[1] - flank
            max_start = peak_row[2] + flank

            region_select = [np.array([peak_row[0] for i in range( ((flank // resolution) * 2 + 1))]) ,
                            np.arange(min_start,max_start,resolution),
                            np.arange(min_start+resolution,max_start+resolution,resolution)]
            region_select = np.array(region_select).T
            region_select = pd.DataFrame(region_select).rename(columns={0:"chrom",1:"start",2:"end"})
            region_select["start"] = region_select["start"].astype(int)
            region_select["end"] = region_select["end"].astype(int)
            index = region_select.set_index(["chrom","start","end"]).index

            res.append(temp_frags.loc[index,:]["count"].values)
        except:
            continue

    res = np.array(res)
    return res