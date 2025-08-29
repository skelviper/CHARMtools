# wrap some useful functions i can never remember their name

import re
from typing import Optional

def get_mem_size(obj):
    from pympler import asizeof
    return asizeof.asizeof(obj)

def auto_genome_coord(genome_coord):
    """
    Automatically convert genome coordinate to chrom, start, end format.
    
    This function parses various genome coordinate formats and returns standardized
    chromosome, start, and end values. It supports both region-based coordinates
    (with start and end positions) and chromosome-only specifications.
    
    Parameters:
    -----------
    genome_coord : str, list, or tuple
        Genome coordinate in one of the following formats:
        - String with colon/hyphen separator: "chr1:10000-20000" or "chr1a:10000-20000"
        - List or tuple: ["chr1a", 10000, 20000]
        - Chromosome name only: "chr1a"
    
    Returns:
    --------
    tuple
        A tuple containing (chrom, start, end) where:
        - chrom (str): Chromosome name
        - start (int or None): Start position (None for chromosome-only input)
        - end (int or None): End position (None for chromosome-only input)
    
    Raises:
    -------
    ValueError
        If genome_coord is not in a supported format
    
    Examples:
    ---------
    >>> auto_genome_coord("chr1:10000-20000")
    ('chr1', 10000, 20000)
    >>> auto_genome_coord(["chr1a", 10000, 20000])
    ('chr1a', 10000, 20000)
    >>> auto_genome_coord("chr1a")
    ('chr1a', None, None)
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

def natural_chr_key(ch: str):
    """染色体自然排序的键函数
    
    Parameters:
    -----------
    ch : str
        染色体名称，如 'chr1', 'chr2', 'chrX', 'chrY'
    
    Returns:
    --------
    tuple
        排序键，数字染色体按数值排序，非数字染色体按字母排序
    
    Examples:
    ---------
    >>> sorted(['chr10', 'chr1', 'chr2', 'chrX'], key=natural_chr_key)
    ['chr1', 'chr2', 'chr10', 'chrX']
    """
    # 去掉 'chr' 前缀
    ch_clean = ch.replace('chr', '')
    # 数字染色体按数值排序，非数字染色体按字母排序
    if ch_clean.isdigit():
        return (0, int(ch_clean))
    else:
        return (1, ch_clean)