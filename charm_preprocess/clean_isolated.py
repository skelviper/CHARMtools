from concurrent import futures
import sys
import time
from functools import partial

import pandas as pd
import numpy as np

from ..utils.CHARMio import parse_pairs, write_pairs


def cli(args):
    filename, out_name, num_thread, up_dense, up_distance = \
        args.filename[0], args.output_file, args.thread, args.dense, args.distance
    pairs = parse_pairs(filename)
    write_pairs(clean_isolated(pairs, num_thread, up_dense, up_distance), out_name)


def clean_isolated(pairs, num_thread, up_dense, up_distance):
    t0 = time.time()
    input_data = [(key, value) for key, value in pairs.groupby(["chrom1", "chrom2"], sort=False)]
    working_func = partial(_clean_contacts_in_pair_fast, up_dense=up_dense, up_distance=up_distance)
    if num_thread > 1:
        with futures.ProcessPoolExecutor(num_thread) as executor:
            cleaned_groups = list(executor.map(working_func, input_data))
    else:
        cleaned_groups = [working_func(item) for item in input_data]
    cleaned = pd.concat(cleaned_groups, axis=0)
    cleaned.attrs = pairs.attrs # groupby don't keep attrs
    print("clean_isolated: %d contacts removed in %s" % (len(pairs)-len(cleaned), pairs.attrs["name"]))
    sys.stderr.write("clean_isolated: finished in %.2fs\n"%(time.time()-t0))
    return cleaned


def _clean_contacts_in_pair_fast(item, up_dense, up_distance) -> "dataframe":
    """
    NumPy-based replacement of the original per-contact apply loop.
    Keeps contacts that have at least ``up_dense + 1`` neighbours (including itself)
    within the L-0.5 distance threshold.
    """
    key, contacts = item
    if contacts.empty:
        return contacts
    sorted_contacts = contacts.sort_values(by="pos1", axis=0, kind="mergesort")
    pos1 = sorted_contacts["pos1"].to_numpy()
    pos2 = sorted_contacts["pos2"].to_numpy()
    n = pos1.shape[0]
    keep = np.zeros(n, dtype=bool)
    left = 0
    right = 0
    threshold = up_dense + 1
    for i in range(n):
        current_pos = pos1[i]
        while right < n and pos1[right] - current_pos <= up_distance:
            right += 1
        while current_pos - pos1[left] > up_distance:
            left += 1
        window_len = right - left
        if window_len < threshold:
            continue
        idx = slice(left, right)
        dx = np.abs(pos1[idx] - current_pos)
        dy = np.abs(pos2[idx] - pos2[i])
        proximity = np.count_nonzero(np.square(np.sqrt(dx) + np.sqrt(dy)) <= up_distance)
        keep[i] = proximity >= threshold
    kept_contacts = sorted_contacts.loc[keep]
    kept_contacts = kept_contacts.sort_index()
    sys.stderr.write("(%s, %s): %d --> %d\n" % (key[0], key[1], len(contacts), len(kept_contacts)))
    return kept_contacts

