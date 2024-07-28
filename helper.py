# wrap some useful functions i can never remember their name

def get_mem_size(obj):
    from pympler import asizeof
    return asizeof.asizeof(obj)