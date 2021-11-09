#loop pileup analysis
#requirements
import cooler
import cooltools
import pandas as pd
import numpy as np
from cooltools import snipping
import multiprocess

#loop should first be called and save in a tsv-like file contain at least chrom-start-end for both loop anchors
def cooltoolsGetOEPileUp(clr:cooler.Cooler,flank:int,resolution:int,expected:pd.DataFrame,loopAnchor:pd.DataFrame,arms:pd.DataFrame,nthreads:int)->np.ndarray:
    """
    cooltools warpper for OE pileup
    """
    loopAnchor.loc[:, 'mid1'] = (loopAnchor['start1'] + loopAnchor['end1'])//2
    loopAnchor.loc[:, 'mid2'] = (loopAnchor['start2'] + loopAnchor['end2'])//2

    windows1 = snipping.make_bin_aligned_windows(
        resolution,
        loopAnchor['chrom1'],
        loopAnchor['mid1'],
        flank_bp=flank)

    windows2 = snipping.make_bin_aligned_windows(
        resolution,
        loopAnchor['chrom2'],
        loopAnchor['mid2'],
        flank_bp=flank)

    windows = pd.merge(windows1, windows2, left_index=True, right_index=True, suffixes=('1', '2'))
    windows = snipping.assign_regions(windows, arms)
    oe_snipper = cooltools.snipping.ObsExpSnipper(clr, expected, regions=arms)
    # create the stack of snips:
    with multiprocess.Pool(nthreads) as pool:
        stack = cooltools.snipping.pileup(
                windows,
                oe_snipper.select,
                oe_snipper.snip,
                map=pool.map
                )

    mtx = np.nanmean(stack, axis=2)


    return mtx,stack
