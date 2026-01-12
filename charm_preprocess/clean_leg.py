import pandas as pd
import numpy as np
import sys
import time
from concurrent import futures
from functools import partial

from ..utils.CHARMio import parse_pairs, write_pairs


def _promiscuous_flags_for_group(item, sorted_positions_by_chrom, max_distance, max_count):
    chrom, pos_series = item
    sorted_all = sorted_positions_by_chrom.get(chrom)
    if sorted_all is None or sorted_all.size == 0:
        return pos_series.index.to_numpy(), np.zeros(pos_series.shape[0], dtype=bool)
    positions = pos_series.to_numpy()
    left_idx = np.searchsorted(sorted_all, positions - max_distance, side="left")
    candidate_idx = left_idx + max_count
    flags = np.zeros(positions.shape[0], dtype=bool)
    in_bounds = candidate_idx < sorted_all.shape[0]
    flags[in_bounds] = (sorted_all[candidate_idx[in_bounds]] - positions[in_bounds]) <= max_distance
    return pos_series.index.to_numpy(), flags


def _promiscuous_mask_for_column(pairs, chrom_col, pos_col, sorted_positions_by_chrom, max_distance, max_count, num_thread):
    groups = [(chrom, group[pos_col]) for chrom, group in pairs.groupby(chrom_col, sort=False)]
    mask = np.zeros(len(pairs), dtype=bool)
    worker = partial(
        _promiscuous_flags_for_group,
        sorted_positions_by_chrom=sorted_positions_by_chrom,
        max_distance=max_distance,
        max_count=max_count,
    )
    if num_thread > 1:
        with futures.ThreadPoolExecutor(max_workers=num_thread) as executor:
            for idxs, flags in executor.map(worker, groups):
                mask[idxs] = flags
    else:
        for item in groups:
            idxs, flags = worker(item)
            mask[idxs] = flags
    return pd.Series(mask, index=pairs.index, dtype=bool)


def cli(args):
    in_name, num_thread, out_name, max_distance, max_count = \
        args.filename[0], args.thread, args.out_name, args.max_distance, args.max_count
    pairs = parse_pairs(in_name)
    res = clean_leg(pairs, num_thread, max_distance, max_count)
    write_pairs(res, out_name)


def clean_leg(pairs, num_thread: int, max_distance: int, max_count: int):
    t0 = time.time()
    t_sort = time.time()
    left = pairs[["chrom1", "pos1"]].copy()
    right = pairs[["chrom2", "pos2"]].copy()
    left.columns, right.columns = ("chr", "pos"), ("chr", "pos")
    all_legs = pd.concat((left, right), axis=0, ignore_index=True)
    sorted_positions_by_chrom = {
        key: np.sort(value["pos"].to_numpy())
        for key, value in all_legs.groupby("chr", sort=False)
    }
    sys.stderr.write("clean_leg: group sort in %.2fs\n" % (time.time() - t_sort))
    left_mask = _promiscuous_mask_for_column(
        pairs, "chrom1", "pos1", sorted_positions_by_chrom, max_distance, max_count, num_thread
    )
    right_mask = _promiscuous_mask_for_column(
        pairs, "chrom2", "pos2", sorted_positions_by_chrom, max_distance, max_count, num_thread
    )
    mask = left_mask | right_mask
    result = pairs.loc[~mask].copy()
    result.attrs.update(pairs.attrs)
    print("clean_leg: remove %d contacts in %s\n" % (len(pairs) - len(result), pairs.attrs["name"]))
    sys.stderr.write("clean_leg: finished in %.2fs\n" % (time.time() - t0))
    return result

