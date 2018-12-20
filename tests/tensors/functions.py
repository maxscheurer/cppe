import numpy as np
from numba import jit
import itertools

"""Contains various helper functions"""
#if we want to go higher than this (unlikely) we have to change dtype from int to float (overflow)
TABLE_SIZE = 20

def fdouble_factorial(n):
    if n < 2:
        return 1
    else:
        return n * fdouble_factorial(n-2)


def ffactorial(n):
    if n < 2:
        return 1
    else:
        return n * ffactorial(n-1)

factorial = np.zeros(TABLE_SIZE, dtype=np.int64)
double_factorial = np.zeros(TABLE_SIZE, dtype=np.int64)
for i in range(TABLE_SIZE):
    factorial[i]        = ffactorial(i)
    double_factorial[i] = fdouble_factorial(i)


def fill(values, dim):
    out = np.zeros([3]*dim)
    counter = 0
    for c in itertools.combinations_with_replacement((0,1,2), dim):
        for index in list(set(itertools.permutations(c))):
            out[index] = values[counter]
        counter += 1
    return out
