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