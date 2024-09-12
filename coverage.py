def f(X, C, N):
    """
    Function representing the Lander-Waterman equation:
    C/X = 1 - exp(-N/X)
    Rearranged to:
    f(X) = C/X - (1 - exp(-N/X))
    """
    return C / X - (1 - math.exp(-N / X))

def estimate_library_size(read_pairs, unique_read_pairs):
    """
    Estimates the size of a library based on the number of paired end molecules observed
    and the number of unique pairs observed.

    :param read_pairs: total number of read pairs (N)
    :param unique_read_pairs: number of distinct fragments observed in read pairs (C)
    :return: estimated number of distinct molecules in the library (X) or None if invalid input
    """
    read_pair_duplicates = read_pairs - unique_read_pairs

    if read_pairs > 0 and read_pair_duplicates > 0:
        m = 1.0
        M = 100.0

        # Check initial condition for the bounds
        if unique_read_pairs >= read_pairs or f(m * unique_read_pairs, unique_read_pairs, read_pairs) < 0:
            raise ValueError("Invalid values for pairs and unique pairs: {}, {}".format(read_pairs, unique_read_pairs))

        # Find value of M, large enough to act as other side for bisection method
        while f(M * unique_read_pairs, unique_read_pairs, read_pairs) > 0:
            M *= 10.0

        # Use bisection method (no more than 40 times) to find solution
        for _ in range(40):
            r = (m + M) / 2.0
            u = f(r * unique_read_pairs, unique_read_pairs, read_pairs)
            if u == 0:
                break
            elif u > 0:
                m = r
            else:
                M = r

        return int(unique_read_pairs * (m + M) / 2.0)
    else:
        return None