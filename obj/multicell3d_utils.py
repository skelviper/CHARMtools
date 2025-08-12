import numpy as np
import pandas as pd
import concurrent.futures
import tqdm
import re

# Utility functions
def _concat_in_chunks(data_chunk):
    return pd.concat(data_chunk, axis=1)

def _parallel_concat(data, nproc=10):
    chunk_size = len(data) // nproc
    data_chunks = [data[i:i + chunk_size] for i in range(0, len(data), chunk_size)]

    with concurrent.futures.ProcessPoolExecutor(nproc) as executor:
        concatenated_chunks = list(tqdm.tqdm(executor.map(_concat_in_chunks, data_chunks), total=len(data_chunks)))
    
    final_result = pd.concat(concatenated_chunks, axis=1)
    
    return final_result

def _auto_genome_coord(genome_coord):
    """
    Automatically convert genome_coord to chrom,start,end format
    INPUT:

    OUTPUT:
    """
    # determine the genome_coord format
    if isinstance(genome_coord, str):
        if ":" in genome_coord:
            chrom, start, end = re.split(":|-", genome_coord)
            start, end = int(start), int(end)
            mat_type = "region"
        else:
            chrom, start, end = genome_coord, None, None
            mat_type = "chrom"
    elif isinstance(genome_coord, (list, tuple)):
        chrom, start, end = genome_coord
        mat_type = "region"
    else:
        raise ValueError('Genome_coord should be str or list/tuple. e.g. "chr1a:10000-20000" or ["chr1a",10000,20000] or "chr1a"')
    
    return chrom, start, end

def dev_only(func):
    """
    Decorator for development-only functions.
    """
    import functools
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        if __name__ == "__main__":
            return func(*args, **kwargs)
        else:
            print(f"Skipping {func.__name__} as it's being imported")
    
    return wrapper