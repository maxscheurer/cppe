import itertools


def get_symm_indices(k):
    combinations = []
    for x in range(k, -1, -1):
        for y in range(k, -1, -1):
            for z in range(k, -1, -1):
                if (x + y + z) != k:
                    continue
                combinations.append((x, y, z))
    component_products = list(itertools.product((0, 1, 2), repeat=k))
    tuple_products = []
    for c in component_products:
        t = (c.count(0), c.count(1), c.count(2))
        tuple_products.append(t)

    sym_indices = []
    for tu in combinations:
        sym_indices.append(tuple_products.index(tu))
    return sym_indices
