# wrap some useful functions i can never remember their name

import re
from typing import Optional

def get_mem_size(obj):
    from pympler import asizeof
    return asizeof.asizeof(obj)

def parse_genome_coord(coord: Optional[str]):
    """解析基因组坐标字符串，返回 (chrom, start, end) 或 None
    
    Parameters:
    -----------
    coord : str, optional
        基因组坐标字符串，支持格式：
        - chr1:1000000-2000000 (区间)
        - chr1:1000000 (单点)
        - "" 或 None (返回None)
    
    Returns:
    --------
    tuple or None
        (chrom, start, end) 或 None
    
    Examples:
    ---------
    >>> parse_genome_coord("chr1:1000000-2000000")
    ('chr1', 1000000, 2000000)
    >>> parse_genome_coord("chr1:1000000")
    ('chr1', 1000000, 1000000)
    >>> parse_genome_coord("")
    None
    """
    if not coord or coord.strip() == "":
        return None
    # 支持格式：chr1:1000000-2000000 或 chr1:1000000
    match = re.match(r'^(chr[^:]+):(\d+)(?:-(\d+))?$', coord.strip())
    if match:
        chrom, start_str, end_str = match.groups()
        start = int(start_str)
        end = int(end_str) if end_str else start
        return chrom, start, end
    return None

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