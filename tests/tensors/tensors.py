from __future__ import division
from numba import jit
from numpy import sqrt
import numpy as np
@jit(nopython=True, cache=True)
def T0(x,y,z):
    return np.array(1/sqrt(x**2 + y**2 + z**2))
@jit(nopython=True, cache=True)
def Tx(x,y,z):
    return -x/(x**2 + y**2 + z**2)**(3/2)
@jit(nopython=True, cache=True)
def Ty(x,y,z):
    return -y/(x**2 + y**2 + z**2)**(3/2)
@jit(nopython=True, cache=True)
def Tz(x,y,z):
    return -z/(x**2 + y**2 + z**2)**(3/2)
@jit(nopython=True, cache=True)
def T1(x,y,z):
    arr = np.zeros((3,), dtype=np.float64)
    arr[(0,)] = Tx(x,y,z)
    arr[(1,)] = Ty(x,y,z)
    arr[(2,)] = Tz(x,y,z)
    return arr
@jit(nopython=True, cache=True)
def Txx(x,y,z):
    return (3*x**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(3/2)
@jit(nopython=True, cache=True)
def Txy(x,y,z):
    return 3*x*y/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Txz(x,y,z):
    return 3*x*z/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Tyx(x,y,z):
    return 3*x*y/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Tyy(x,y,z):
    return (3*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(3/2)
@jit(nopython=True, cache=True)
def Tyz(x,y,z):
    return 3*y*z/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Tzx(x,y,z):
    return 3*x*z/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Tzy(x,y,z):
    return 3*y*z/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Tzz(x,y,z):
    return (3*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(3/2)
@jit(nopython=True, cache=True)
def T2(x,y,z):
    arr = np.zeros((3, 3), dtype=np.float64)
    arr[(0, 0)] = Txx(x,y,z)
    arr[(0, 1)] = Txy(x,y,z)
    arr[(0, 2)] = Txz(x,y,z)
    arr[(1, 0)] = Tyx(x,y,z)
    arr[(1, 1)] = Tyy(x,y,z)
    arr[(1, 2)] = Tyz(x,y,z)
    arr[(2, 0)] = Tzx(x,y,z)
    arr[(2, 1)] = Tzy(x,y,z)
    arr[(2, 2)] = Tzz(x,y,z)
    return arr
@jit(nopython=True, cache=True)
def Txxx(x,y,z):
    return 3*x*(-5*x**2/(x**2 + y**2 + z**2) + 3)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Txxy(x,y,z):
    return 3*y*(-5*x**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Txxz(x,y,z):
    return 3*z*(-5*x**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Txyx(x,y,z):
    return 3*y*(-5*x**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Txyy(x,y,z):
    return 3*x*(-5*y**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Txyz(x,y,z):
    return -15*x*y*z/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txzx(x,y,z):
    return 3*z*(-5*x**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Txzy(x,y,z):
    return -15*x*y*z/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txzz(x,y,z):
    return 3*x*(-5*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Tyxx(x,y,z):
    return 3*y*(-5*x**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Tyxy(x,y,z):
    return 3*x*(-5*y**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Tyxz(x,y,z):
    return -15*x*y*z/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyyx(x,y,z):
    return 3*x*(-5*y**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Tyyy(x,y,z):
    return 3*y*(-5*y**2/(x**2 + y**2 + z**2) + 3)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Tyyz(x,y,z):
    return 3*z*(-5*y**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Tyzx(x,y,z):
    return -15*x*y*z/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyzy(x,y,z):
    return 3*z*(-5*y**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Tyzz(x,y,z):
    return 3*y*(-5*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Tzxx(x,y,z):
    return 3*z*(-5*x**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Tzxy(x,y,z):
    return -15*x*y*z/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzxz(x,y,z):
    return 3*x*(-5*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Tzyx(x,y,z):
    return -15*x*y*z/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzyy(x,y,z):
    return 3*z*(-5*y**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Tzyz(x,y,z):
    return 3*y*(-5*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Tzzx(x,y,z):
    return 3*x*(-5*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Tzzy(x,y,z):
    return 3*y*(-5*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Tzzz(x,y,z):
    return 3*z*(-5*z**2/(x**2 + y**2 + z**2) + 3)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def T3(x,y,z):
    arr = np.zeros((3, 3, 3), dtype=np.float64)
    arr[(0, 0, 0)] = Txxx(x,y,z)
    arr[(0, 0, 1)] = Txxy(x,y,z)
    arr[(0, 0, 2)] = Txxz(x,y,z)
    arr[(0, 1, 0)] = Txyx(x,y,z)
    arr[(0, 1, 1)] = Txyy(x,y,z)
    arr[(0, 1, 2)] = Txyz(x,y,z)
    arr[(0, 2, 0)] = Txzx(x,y,z)
    arr[(0, 2, 1)] = Txzy(x,y,z)
    arr[(0, 2, 2)] = Txzz(x,y,z)
    arr[(1, 0, 0)] = Tyxx(x,y,z)
    arr[(1, 0, 1)] = Tyxy(x,y,z)
    arr[(1, 0, 2)] = Tyxz(x,y,z)
    arr[(1, 1, 0)] = Tyyx(x,y,z)
    arr[(1, 1, 1)] = Tyyy(x,y,z)
    arr[(1, 1, 2)] = Tyyz(x,y,z)
    arr[(1, 2, 0)] = Tyzx(x,y,z)
    arr[(1, 2, 1)] = Tyzy(x,y,z)
    arr[(1, 2, 2)] = Tyzz(x,y,z)
    arr[(2, 0, 0)] = Tzxx(x,y,z)
    arr[(2, 0, 1)] = Tzxy(x,y,z)
    arr[(2, 0, 2)] = Tzxz(x,y,z)
    arr[(2, 1, 0)] = Tzyx(x,y,z)
    arr[(2, 1, 1)] = Tzyy(x,y,z)
    arr[(2, 1, 2)] = Tzyz(x,y,z)
    arr[(2, 2, 0)] = Tzzx(x,y,z)
    arr[(2, 2, 1)] = Tzzy(x,y,z)
    arr[(2, 2, 2)] = Tzzz(x,y,z)
    return arr
@jit(nopython=True, cache=True)
def Txxxx(x,y,z):
    return 3*(35*x**4/(x**2 + y**2 + z**2)**2 - 30*x**2/(x**2 + y**2 + z**2) + 3)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Txxxy(x,y,z):
    return 15*x*y*(7*x**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txxxz(x,y,z):
    return 15*x*z*(7*x**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txxyx(x,y,z):
    return 15*x*y*(7*x**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txxyy(x,y,z):
    return 3*(35*x**2*y**2/(x**2 + y**2 + z**2)**2 - 5*x**2/(x**2 + y**2 + z**2) - 5*y**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Txxyz(x,y,z):
    return 15*y*z*(7*x**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txxzx(x,y,z):
    return 15*x*z*(7*x**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txxzy(x,y,z):
    return 15*y*z*(7*x**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txxzz(x,y,z):
    return 3*(35*x**2*z**2/(x**2 + y**2 + z**2)**2 - 5*x**2/(x**2 + y**2 + z**2) - 5*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Txyxx(x,y,z):
    return 15*x*y*(7*x**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txyxy(x,y,z):
    return 3*(35*x**2*y**2/(x**2 + y**2 + z**2)**2 - 5*x**2/(x**2 + y**2 + z**2) - 5*y**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Txyxz(x,y,z):
    return 15*y*z*(7*x**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txyyx(x,y,z):
    return 3*(35*x**2*y**2/(x**2 + y**2 + z**2)**2 - 5*x**2/(x**2 + y**2 + z**2) - 5*y**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Txyyy(x,y,z):
    return 15*x*y*(7*y**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txyyz(x,y,z):
    return 15*x*z*(7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txyzx(x,y,z):
    return 15*y*z*(7*x**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txyzy(x,y,z):
    return 15*x*z*(7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txyzz(x,y,z):
    return 15*x*y*(7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txzxx(x,y,z):
    return 15*x*z*(7*x**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txzxy(x,y,z):
    return 15*y*z*(7*x**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txzxz(x,y,z):
    return 3*(35*x**2*z**2/(x**2 + y**2 + z**2)**2 - 5*x**2/(x**2 + y**2 + z**2) - 5*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Txzyx(x,y,z):
    return 15*y*z*(7*x**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txzyy(x,y,z):
    return 15*x*z*(7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txzyz(x,y,z):
    return 15*x*y*(7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txzzx(x,y,z):
    return 3*(35*x**2*z**2/(x**2 + y**2 + z**2)**2 - 5*x**2/(x**2 + y**2 + z**2) - 5*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Txzzy(x,y,z):
    return 15*x*y*(7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txzzz(x,y,z):
    return 15*x*z*(7*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyxxx(x,y,z):
    return 15*x*y*(7*x**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyxxy(x,y,z):
    return 3*(35*x**2*y**2/(x**2 + y**2 + z**2)**2 - 5*x**2/(x**2 + y**2 + z**2) - 5*y**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Tyxxz(x,y,z):
    return 15*y*z*(7*x**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyxyx(x,y,z):
    return 3*(35*x**2*y**2/(x**2 + y**2 + z**2)**2 - 5*x**2/(x**2 + y**2 + z**2) - 5*y**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Tyxyy(x,y,z):
    return 15*x*y*(7*y**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyxyz(x,y,z):
    return 15*x*z*(7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyxzx(x,y,z):
    return 15*y*z*(7*x**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyxzy(x,y,z):
    return 15*x*z*(7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyxzz(x,y,z):
    return 15*x*y*(7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyyxx(x,y,z):
    return 3*(35*x**2*y**2/(x**2 + y**2 + z**2)**2 - 5*x**2/(x**2 + y**2 + z**2) - 5*y**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Tyyxy(x,y,z):
    return 15*x*y*(7*y**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyyxz(x,y,z):
    return 15*x*z*(7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyyyx(x,y,z):
    return 15*x*y*(7*y**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyyyy(x,y,z):
    return 3*(35*y**4/(x**2 + y**2 + z**2)**2 - 30*y**2/(x**2 + y**2 + z**2) + 3)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Tyyyz(x,y,z):
    return 15*y*z*(7*y**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyyzx(x,y,z):
    return 15*x*z*(7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyyzy(x,y,z):
    return 15*y*z*(7*y**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyyzz(x,y,z):
    return 3*(35*y**2*z**2/(x**2 + y**2 + z**2)**2 - 5*y**2/(x**2 + y**2 + z**2) - 5*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Tyzxx(x,y,z):
    return 15*y*z*(7*x**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyzxy(x,y,z):
    return 15*x*z*(7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyzxz(x,y,z):
    return 15*x*y*(7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyzyx(x,y,z):
    return 15*x*z*(7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyzyy(x,y,z):
    return 15*y*z*(7*y**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyzyz(x,y,z):
    return 3*(35*y**2*z**2/(x**2 + y**2 + z**2)**2 - 5*y**2/(x**2 + y**2 + z**2) - 5*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Tyzzx(x,y,z):
    return 15*x*y*(7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyzzy(x,y,z):
    return 3*(35*y**2*z**2/(x**2 + y**2 + z**2)**2 - 5*y**2/(x**2 + y**2 + z**2) - 5*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Tyzzz(x,y,z):
    return 15*y*z*(7*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzxxx(x,y,z):
    return 15*x*z*(7*x**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzxxy(x,y,z):
    return 15*y*z*(7*x**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzxxz(x,y,z):
    return 3*(35*x**2*z**2/(x**2 + y**2 + z**2)**2 - 5*x**2/(x**2 + y**2 + z**2) - 5*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Tzxyx(x,y,z):
    return 15*y*z*(7*x**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzxyy(x,y,z):
    return 15*x*z*(7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzxyz(x,y,z):
    return 15*x*y*(7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzxzx(x,y,z):
    return 3*(35*x**2*z**2/(x**2 + y**2 + z**2)**2 - 5*x**2/(x**2 + y**2 + z**2) - 5*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Tzxzy(x,y,z):
    return 15*x*y*(7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzxzz(x,y,z):
    return 15*x*z*(7*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzyxx(x,y,z):
    return 15*y*z*(7*x**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzyxy(x,y,z):
    return 15*x*z*(7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzyxz(x,y,z):
    return 15*x*y*(7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzyyx(x,y,z):
    return 15*x*z*(7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzyyy(x,y,z):
    return 15*y*z*(7*y**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzyyz(x,y,z):
    return 3*(35*y**2*z**2/(x**2 + y**2 + z**2)**2 - 5*y**2/(x**2 + y**2 + z**2) - 5*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Tzyzx(x,y,z):
    return 15*x*y*(7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzyzy(x,y,z):
    return 3*(35*y**2*z**2/(x**2 + y**2 + z**2)**2 - 5*y**2/(x**2 + y**2 + z**2) - 5*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Tzyzz(x,y,z):
    return 15*y*z*(7*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzzxx(x,y,z):
    return 3*(35*x**2*z**2/(x**2 + y**2 + z**2)**2 - 5*x**2/(x**2 + y**2 + z**2) - 5*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Tzzxy(x,y,z):
    return 15*x*y*(7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzzxz(x,y,z):
    return 15*x*z*(7*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzzyx(x,y,z):
    return 15*x*y*(7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzzyy(x,y,z):
    return 3*(35*y**2*z**2/(x**2 + y**2 + z**2)**2 - 5*y**2/(x**2 + y**2 + z**2) - 5*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def Tzzyz(x,y,z):
    return 15*y*z*(7*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzzzx(x,y,z):
    return 15*x*z*(7*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzzzy(x,y,z):
    return 15*y*z*(7*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzzzz(x,y,z):
    return 3*(35*z**4/(x**2 + y**2 + z**2)**2 - 30*z**2/(x**2 + y**2 + z**2) + 3)/(x**2 + y**2 + z**2)**(5/2)
@jit(nopython=True, cache=True)
def T4(x,y,z):
    arr = np.zeros((3, 3, 3, 3), dtype=np.float64)
    arr[(0, 0, 0, 0)] = Txxxx(x,y,z)
    arr[(0, 0, 0, 1)] = Txxxy(x,y,z)
    arr[(0, 0, 0, 2)] = Txxxz(x,y,z)
    arr[(0, 0, 1, 0)] = Txxyx(x,y,z)
    arr[(0, 0, 1, 1)] = Txxyy(x,y,z)
    arr[(0, 0, 1, 2)] = Txxyz(x,y,z)
    arr[(0, 0, 2, 0)] = Txxzx(x,y,z)
    arr[(0, 0, 2, 1)] = Txxzy(x,y,z)
    arr[(0, 0, 2, 2)] = Txxzz(x,y,z)
    arr[(0, 1, 0, 0)] = Txyxx(x,y,z)
    arr[(0, 1, 0, 1)] = Txyxy(x,y,z)
    arr[(0, 1, 0, 2)] = Txyxz(x,y,z)
    arr[(0, 1, 1, 0)] = Txyyx(x,y,z)
    arr[(0, 1, 1, 1)] = Txyyy(x,y,z)
    arr[(0, 1, 1, 2)] = Txyyz(x,y,z)
    arr[(0, 1, 2, 0)] = Txyzx(x,y,z)
    arr[(0, 1, 2, 1)] = Txyzy(x,y,z)
    arr[(0, 1, 2, 2)] = Txyzz(x,y,z)
    arr[(0, 2, 0, 0)] = Txzxx(x,y,z)
    arr[(0, 2, 0, 1)] = Txzxy(x,y,z)
    arr[(0, 2, 0, 2)] = Txzxz(x,y,z)
    arr[(0, 2, 1, 0)] = Txzyx(x,y,z)
    arr[(0, 2, 1, 1)] = Txzyy(x,y,z)
    arr[(0, 2, 1, 2)] = Txzyz(x,y,z)
    arr[(0, 2, 2, 0)] = Txzzx(x,y,z)
    arr[(0, 2, 2, 1)] = Txzzy(x,y,z)
    arr[(0, 2, 2, 2)] = Txzzz(x,y,z)
    arr[(1, 0, 0, 0)] = Tyxxx(x,y,z)
    arr[(1, 0, 0, 1)] = Tyxxy(x,y,z)
    arr[(1, 0, 0, 2)] = Tyxxz(x,y,z)
    arr[(1, 0, 1, 0)] = Tyxyx(x,y,z)
    arr[(1, 0, 1, 1)] = Tyxyy(x,y,z)
    arr[(1, 0, 1, 2)] = Tyxyz(x,y,z)
    arr[(1, 0, 2, 0)] = Tyxzx(x,y,z)
    arr[(1, 0, 2, 1)] = Tyxzy(x,y,z)
    arr[(1, 0, 2, 2)] = Tyxzz(x,y,z)
    arr[(1, 1, 0, 0)] = Tyyxx(x,y,z)
    arr[(1, 1, 0, 1)] = Tyyxy(x,y,z)
    arr[(1, 1, 0, 2)] = Tyyxz(x,y,z)
    arr[(1, 1, 1, 0)] = Tyyyx(x,y,z)
    arr[(1, 1, 1, 1)] = Tyyyy(x,y,z)
    arr[(1, 1, 1, 2)] = Tyyyz(x,y,z)
    arr[(1, 1, 2, 0)] = Tyyzx(x,y,z)
    arr[(1, 1, 2, 1)] = Tyyzy(x,y,z)
    arr[(1, 1, 2, 2)] = Tyyzz(x,y,z)
    arr[(1, 2, 0, 0)] = Tyzxx(x,y,z)
    arr[(1, 2, 0, 1)] = Tyzxy(x,y,z)
    arr[(1, 2, 0, 2)] = Tyzxz(x,y,z)
    arr[(1, 2, 1, 0)] = Tyzyx(x,y,z)
    arr[(1, 2, 1, 1)] = Tyzyy(x,y,z)
    arr[(1, 2, 1, 2)] = Tyzyz(x,y,z)
    arr[(1, 2, 2, 0)] = Tyzzx(x,y,z)
    arr[(1, 2, 2, 1)] = Tyzzy(x,y,z)
    arr[(1, 2, 2, 2)] = Tyzzz(x,y,z)
    arr[(2, 0, 0, 0)] = Tzxxx(x,y,z)
    arr[(2, 0, 0, 1)] = Tzxxy(x,y,z)
    arr[(2, 0, 0, 2)] = Tzxxz(x,y,z)
    arr[(2, 0, 1, 0)] = Tzxyx(x,y,z)
    arr[(2, 0, 1, 1)] = Tzxyy(x,y,z)
    arr[(2, 0, 1, 2)] = Tzxyz(x,y,z)
    arr[(2, 0, 2, 0)] = Tzxzx(x,y,z)
    arr[(2, 0, 2, 1)] = Tzxzy(x,y,z)
    arr[(2, 0, 2, 2)] = Tzxzz(x,y,z)
    arr[(2, 1, 0, 0)] = Tzyxx(x,y,z)
    arr[(2, 1, 0, 1)] = Tzyxy(x,y,z)
    arr[(2, 1, 0, 2)] = Tzyxz(x,y,z)
    arr[(2, 1, 1, 0)] = Tzyyx(x,y,z)
    arr[(2, 1, 1, 1)] = Tzyyy(x,y,z)
    arr[(2, 1, 1, 2)] = Tzyyz(x,y,z)
    arr[(2, 1, 2, 0)] = Tzyzx(x,y,z)
    arr[(2, 1, 2, 1)] = Tzyzy(x,y,z)
    arr[(2, 1, 2, 2)] = Tzyzz(x,y,z)
    arr[(2, 2, 0, 0)] = Tzzxx(x,y,z)
    arr[(2, 2, 0, 1)] = Tzzxy(x,y,z)
    arr[(2, 2, 0, 2)] = Tzzxz(x,y,z)
    arr[(2, 2, 1, 0)] = Tzzyx(x,y,z)
    arr[(2, 2, 1, 1)] = Tzzyy(x,y,z)
    arr[(2, 2, 1, 2)] = Tzzyz(x,y,z)
    arr[(2, 2, 2, 0)] = Tzzzx(x,y,z)
    arr[(2, 2, 2, 1)] = Tzzzy(x,y,z)
    arr[(2, 2, 2, 2)] = Tzzzz(x,y,z)
    return arr
@jit(nopython=True, cache=True)
def Txxxxx(x,y,z):
    return 15*x*(-63*x**4/(x**2 + y**2 + z**2)**2 + 70*x**2/(x**2 + y**2 + z**2) - 15)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txxxxy(x,y,z):
    return 45*y*(-21*x**4/(x**2 + y**2 + z**2)**2 + 14*x**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txxxxz(x,y,z):
    return 45*z*(-21*x**4/(x**2 + y**2 + z**2)**2 + 14*x**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txxxyx(x,y,z):
    return 45*y*(-21*x**4/(x**2 + y**2 + z**2)**2 + 14*x**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txxxyy(x,y,z):
    return 15*x*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 21*y**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txxxyz(x,y,z):
    return 315*x*y*z*(-3*x**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Txxxzx(x,y,z):
    return 45*z*(-21*x**4/(x**2 + y**2 + z**2)**2 + 14*x**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txxxzy(x,y,z):
    return 315*x*y*z*(-3*x**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Txxxzz(x,y,z):
    return 15*x*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 21*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txxyxx(x,y,z):
    return 45*y*(-21*x**4/(x**2 + y**2 + z**2)**2 + 14*x**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txxyxy(x,y,z):
    return 15*x*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 21*y**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txxyxz(x,y,z):
    return 315*x*y*z*(-3*x**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Txxyyx(x,y,z):
    return 15*x*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 21*y**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txxyyy(x,y,z):
    return 15*y*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 21*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txxyyz(x,y,z):
    return 15*z*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txxyzx(x,y,z):
    return 315*x*y*z*(-3*x**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Txxyzy(x,y,z):
    return 15*z*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txxyzz(x,y,z):
    return 15*y*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txxzxx(x,y,z):
    return 45*z*(-21*x**4/(x**2 + y**2 + z**2)**2 + 14*x**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txxzxy(x,y,z):
    return 315*x*y*z*(-3*x**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Txxzxz(x,y,z):
    return 15*x*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 21*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txxzyx(x,y,z):
    return 315*x*y*z*(-3*x**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Txxzyy(x,y,z):
    return 15*z*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txxzyz(x,y,z):
    return 15*y*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txxzzx(x,y,z):
    return 15*x*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 21*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txxzzy(x,y,z):
    return 15*y*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txxzzz(x,y,z):
    return 15*z*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 21*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txyxxx(x,y,z):
    return 45*y*(-21*x**4/(x**2 + y**2 + z**2)**2 + 14*x**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txyxxy(x,y,z):
    return 15*x*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 21*y**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txyxxz(x,y,z):
    return 315*x*y*z*(-3*x**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Txyxyx(x,y,z):
    return 15*x*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 21*y**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txyxyy(x,y,z):
    return 15*y*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 21*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txyxyz(x,y,z):
    return 15*z*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txyxzx(x,y,z):
    return 315*x*y*z*(-3*x**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Txyxzy(x,y,z):
    return 15*z*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txyxzz(x,y,z):
    return 15*y*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txyyxx(x,y,z):
    return 15*x*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 21*y**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txyyxy(x,y,z):
    return 15*y*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 21*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txyyxz(x,y,z):
    return 15*z*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txyyyx(x,y,z):
    return 15*y*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 21*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txyyyy(x,y,z):
    return 45*x*(-21*y**4/(x**2 + y**2 + z**2)**2 + 14*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txyyyz(x,y,z):
    return 315*x*y*z*(-3*y**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Txyyzx(x,y,z):
    return 15*z*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txyyzy(x,y,z):
    return 315*x*y*z*(-3*y**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Txyyzz(x,y,z):
    return 15*x*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txyzxx(x,y,z):
    return 315*x*y*z*(-3*x**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Txyzxy(x,y,z):
    return 15*z*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txyzxz(x,y,z):
    return 15*y*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txyzyx(x,y,z):
    return 15*z*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txyzyy(x,y,z):
    return 315*x*y*z*(-3*y**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Txyzyz(x,y,z):
    return 15*x*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txyzzx(x,y,z):
    return 15*y*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txyzzy(x,y,z):
    return 15*x*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txyzzz(x,y,z):
    return 315*x*y*z*(-3*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Txzxxx(x,y,z):
    return 45*z*(-21*x**4/(x**2 + y**2 + z**2)**2 + 14*x**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txzxxy(x,y,z):
    return 315*x*y*z*(-3*x**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Txzxxz(x,y,z):
    return 15*x*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 21*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txzxyx(x,y,z):
    return 315*x*y*z*(-3*x**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Txzxyy(x,y,z):
    return 15*z*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txzxyz(x,y,z):
    return 15*y*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txzxzx(x,y,z):
    return 15*x*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 21*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txzxzy(x,y,z):
    return 15*y*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txzxzz(x,y,z):
    return 15*z*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 21*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txzyxx(x,y,z):
    return 315*x*y*z*(-3*x**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Txzyxy(x,y,z):
    return 15*z*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txzyxz(x,y,z):
    return 15*y*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txzyyx(x,y,z):
    return 15*z*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txzyyy(x,y,z):
    return 315*x*y*z*(-3*y**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Txzyyz(x,y,z):
    return 15*x*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txzyzx(x,y,z):
    return 15*y*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txzyzy(x,y,z):
    return 15*x*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txzyzz(x,y,z):
    return 315*x*y*z*(-3*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Txzzxx(x,y,z):
    return 15*x*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 21*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txzzxy(x,y,z):
    return 15*y*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txzzxz(x,y,z):
    return 15*z*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 21*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txzzyx(x,y,z):
    return 15*y*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txzzyy(x,y,z):
    return 15*x*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txzzyz(x,y,z):
    return 315*x*y*z*(-3*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Txzzzx(x,y,z):
    return 15*z*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 21*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Txzzzy(x,y,z):
    return 315*x*y*z*(-3*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Txzzzz(x,y,z):
    return 45*x*(-21*z**4/(x**2 + y**2 + z**2)**2 + 14*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyxxxx(x,y,z):
    return 45*y*(-21*x**4/(x**2 + y**2 + z**2)**2 + 14*x**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyxxxy(x,y,z):
    return 15*x*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 21*y**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyxxxz(x,y,z):
    return 315*x*y*z*(-3*x**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tyxxyx(x,y,z):
    return 15*x*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 21*y**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyxxyy(x,y,z):
    return 15*y*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 21*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyxxyz(x,y,z):
    return 15*z*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyxxzx(x,y,z):
    return 315*x*y*z*(-3*x**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tyxxzy(x,y,z):
    return 15*z*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyxxzz(x,y,z):
    return 15*y*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyxyxx(x,y,z):
    return 15*x*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 21*y**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyxyxy(x,y,z):
    return 15*y*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 21*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyxyxz(x,y,z):
    return 15*z*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyxyyx(x,y,z):
    return 15*y*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 21*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyxyyy(x,y,z):
    return 45*x*(-21*y**4/(x**2 + y**2 + z**2)**2 + 14*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyxyyz(x,y,z):
    return 315*x*y*z*(-3*y**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tyxyzx(x,y,z):
    return 15*z*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyxyzy(x,y,z):
    return 315*x*y*z*(-3*y**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tyxyzz(x,y,z):
    return 15*x*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyxzxx(x,y,z):
    return 315*x*y*z*(-3*x**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tyxzxy(x,y,z):
    return 15*z*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyxzxz(x,y,z):
    return 15*y*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyxzyx(x,y,z):
    return 15*z*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyxzyy(x,y,z):
    return 315*x*y*z*(-3*y**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tyxzyz(x,y,z):
    return 15*x*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyxzzx(x,y,z):
    return 15*y*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyxzzy(x,y,z):
    return 15*x*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyxzzz(x,y,z):
    return 315*x*y*z*(-3*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tyyxxx(x,y,z):
    return 15*x*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 21*y**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyyxxy(x,y,z):
    return 15*y*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 21*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyyxxz(x,y,z):
    return 15*z*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyyxyx(x,y,z):
    return 15*y*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 21*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyyxyy(x,y,z):
    return 45*x*(-21*y**4/(x**2 + y**2 + z**2)**2 + 14*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyyxyz(x,y,z):
    return 315*x*y*z*(-3*y**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tyyxzx(x,y,z):
    return 15*z*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyyxzy(x,y,z):
    return 315*x*y*z*(-3*y**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tyyxzz(x,y,z):
    return 15*x*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyyyxx(x,y,z):
    return 15*y*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 21*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyyyxy(x,y,z):
    return 45*x*(-21*y**4/(x**2 + y**2 + z**2)**2 + 14*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyyyxz(x,y,z):
    return 315*x*y*z*(-3*y**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tyyyyx(x,y,z):
    return 45*x*(-21*y**4/(x**2 + y**2 + z**2)**2 + 14*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyyyyy(x,y,z):
    return 15*y*(-63*y**4/(x**2 + y**2 + z**2)**2 + 70*y**2/(x**2 + y**2 + z**2) - 15)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyyyyz(x,y,z):
    return 45*z*(-21*y**4/(x**2 + y**2 + z**2)**2 + 14*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyyyzx(x,y,z):
    return 315*x*y*z*(-3*y**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tyyyzy(x,y,z):
    return 45*z*(-21*y**4/(x**2 + y**2 + z**2)**2 + 14*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyyyzz(x,y,z):
    return 15*y*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 21*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyyzxx(x,y,z):
    return 15*z*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyyzxy(x,y,z):
    return 315*x*y*z*(-3*y**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tyyzxz(x,y,z):
    return 15*x*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyyzyx(x,y,z):
    return 315*x*y*z*(-3*y**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tyyzyy(x,y,z):
    return 45*z*(-21*y**4/(x**2 + y**2 + z**2)**2 + 14*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyyzyz(x,y,z):
    return 15*y*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 21*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyyzzx(x,y,z):
    return 15*x*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyyzzy(x,y,z):
    return 15*y*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 21*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyyzzz(x,y,z):
    return 15*z*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 21*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyzxxx(x,y,z):
    return 315*x*y*z*(-3*x**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tyzxxy(x,y,z):
    return 15*z*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyzxxz(x,y,z):
    return 15*y*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyzxyx(x,y,z):
    return 15*z*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyzxyy(x,y,z):
    return 315*x*y*z*(-3*y**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tyzxyz(x,y,z):
    return 15*x*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyzxzx(x,y,z):
    return 15*y*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyzxzy(x,y,z):
    return 15*x*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyzxzz(x,y,z):
    return 315*x*y*z*(-3*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tyzyxx(x,y,z):
    return 15*z*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyzyxy(x,y,z):
    return 315*x*y*z*(-3*y**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tyzyxz(x,y,z):
    return 15*x*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyzyyx(x,y,z):
    return 315*x*y*z*(-3*y**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tyzyyy(x,y,z):
    return 45*z*(-21*y**4/(x**2 + y**2 + z**2)**2 + 14*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyzyyz(x,y,z):
    return 15*y*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 21*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyzyzx(x,y,z):
    return 15*x*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyzyzy(x,y,z):
    return 15*y*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 21*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyzyzz(x,y,z):
    return 15*z*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 21*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyzzxx(x,y,z):
    return 15*y*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyzzxy(x,y,z):
    return 15*x*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyzzxz(x,y,z):
    return 315*x*y*z*(-3*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tyzzyx(x,y,z):
    return 15*x*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyzzyy(x,y,z):
    return 15*y*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 21*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyzzyz(x,y,z):
    return 15*z*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 21*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyzzzx(x,y,z):
    return 315*x*y*z*(-3*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tyzzzy(x,y,z):
    return 15*z*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 21*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tyzzzz(x,y,z):
    return 45*y*(-21*z**4/(x**2 + y**2 + z**2)**2 + 14*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzxxxx(x,y,z):
    return 45*z*(-21*x**4/(x**2 + y**2 + z**2)**2 + 14*x**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzxxxy(x,y,z):
    return 315*x*y*z*(-3*x**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tzxxxz(x,y,z):
    return 15*x*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 21*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzxxyx(x,y,z):
    return 315*x*y*z*(-3*x**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tzxxyy(x,y,z):
    return 15*z*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzxxyz(x,y,z):
    return 15*y*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzxxzx(x,y,z):
    return 15*x*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 21*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzxxzy(x,y,z):
    return 15*y*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzxxzz(x,y,z):
    return 15*z*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 21*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzxyxx(x,y,z):
    return 315*x*y*z*(-3*x**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tzxyxy(x,y,z):
    return 15*z*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzxyxz(x,y,z):
    return 15*y*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzxyyx(x,y,z):
    return 15*z*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzxyyy(x,y,z):
    return 315*x*y*z*(-3*y**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tzxyyz(x,y,z):
    return 15*x*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzxyzx(x,y,z):
    return 15*y*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzxyzy(x,y,z):
    return 15*x*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzxyzz(x,y,z):
    return 315*x*y*z*(-3*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tzxzxx(x,y,z):
    return 15*x*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 21*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzxzxy(x,y,z):
    return 15*y*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzxzxz(x,y,z):
    return 15*z*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 21*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzxzyx(x,y,z):
    return 15*y*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzxzyy(x,y,z):
    return 15*x*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzxzyz(x,y,z):
    return 315*x*y*z*(-3*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tzxzzx(x,y,z):
    return 15*z*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 21*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzxzzy(x,y,z):
    return 315*x*y*z*(-3*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tzxzzz(x,y,z):
    return 45*x*(-21*z**4/(x**2 + y**2 + z**2)**2 + 14*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzyxxx(x,y,z):
    return 315*x*y*z*(-3*x**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tzyxxy(x,y,z):
    return 15*z*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzyxxz(x,y,z):
    return 15*y*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzyxyx(x,y,z):
    return 15*z*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzyxyy(x,y,z):
    return 315*x*y*z*(-3*y**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tzyxyz(x,y,z):
    return 15*x*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzyxzx(x,y,z):
    return 15*y*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzyxzy(x,y,z):
    return 15*x*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzyxzz(x,y,z):
    return 315*x*y*z*(-3*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tzyyxx(x,y,z):
    return 15*z*(-63*x**2*y**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzyyxy(x,y,z):
    return 315*x*y*z*(-3*y**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tzyyxz(x,y,z):
    return 15*x*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzyyyx(x,y,z):
    return 315*x*y*z*(-3*y**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tzyyyy(x,y,z):
    return 45*z*(-21*y**4/(x**2 + y**2 + z**2)**2 + 14*y**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzyyyz(x,y,z):
    return 15*y*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 21*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzyyzx(x,y,z):
    return 15*x*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzyyzy(x,y,z):
    return 15*y*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 21*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzyyzz(x,y,z):
    return 15*z*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 21*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzyzxx(x,y,z):
    return 15*y*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzyzxy(x,y,z):
    return 15*x*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzyzxz(x,y,z):
    return 315*x*y*z*(-3*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tzyzyx(x,y,z):
    return 15*x*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzyzyy(x,y,z):
    return 15*y*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 21*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzyzyz(x,y,z):
    return 15*z*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 21*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzyzzx(x,y,z):
    return 315*x*y*z*(-3*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tzyzzy(x,y,z):
    return 15*z*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 21*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzyzzz(x,y,z):
    return 45*y*(-21*z**4/(x**2 + y**2 + z**2)**2 + 14*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzzxxx(x,y,z):
    return 15*x*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 21*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzzxxy(x,y,z):
    return 15*y*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzzxxz(x,y,z):
    return 15*z*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 21*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzzxyx(x,y,z):
    return 15*y*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzzxyy(x,y,z):
    return 15*x*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzzxyz(x,y,z):
    return 315*x*y*z*(-3*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tzzxzx(x,y,z):
    return 15*z*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 21*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzzxzy(x,y,z):
    return 315*x*y*z*(-3*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tzzxzz(x,y,z):
    return 45*x*(-21*z**4/(x**2 + y**2 + z**2)**2 + 14*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzzyxx(x,y,z):
    return 15*y*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 7*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzzyxy(x,y,z):
    return 15*x*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzzyxz(x,y,z):
    return 315*x*y*z*(-3*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tzzyyx(x,y,z):
    return 15*x*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzzyyy(x,y,z):
    return 15*y*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 7*y**2/(x**2 + y**2 + z**2) + 21*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzzyyz(x,y,z):
    return 15*z*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 21*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzzyzx(x,y,z):
    return 315*x*y*z*(-3*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tzzyzy(x,y,z):
    return 15*z*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 21*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzzyzz(x,y,z):
    return 45*y*(-21*z**4/(x**2 + y**2 + z**2)**2 + 14*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzzzxx(x,y,z):
    return 15*z*(-63*x**2*z**2/(x**2 + y**2 + z**2)**2 + 21*x**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzzzxy(x,y,z):
    return 315*x*y*z*(-3*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tzzzxz(x,y,z):
    return 45*x*(-21*z**4/(x**2 + y**2 + z**2)**2 + 14*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzzzyx(x,y,z):
    return 315*x*y*z*(-3*z**2/(x**2 + y**2 + z**2) + 1)/(x**2 + y**2 + z**2)**(9/2)
@jit(nopython=True, cache=True)
def Tzzzyy(x,y,z):
    return 15*z*(-63*y**2*z**2/(x**2 + y**2 + z**2)**2 + 21*y**2/(x**2 + y**2 + z**2) + 7*z**2/(x**2 + y**2 + z**2) - 3)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzzzyz(x,y,z):
    return 45*y*(-21*z**4/(x**2 + y**2 + z**2)**2 + 14*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzzzzx(x,y,z):
    return 45*x*(-21*z**4/(x**2 + y**2 + z**2)**2 + 14*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzzzzy(x,y,z):
    return 45*y*(-21*z**4/(x**2 + y**2 + z**2)**2 + 14*z**2/(x**2 + y**2 + z**2) - 1)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def Tzzzzz(x,y,z):
    return 15*z*(-63*z**4/(x**2 + y**2 + z**2)**2 + 70*z**2/(x**2 + y**2 + z**2) - 15)/(x**2 + y**2 + z**2)**(7/2)
@jit(nopython=True, cache=True)
def T5(x,y,z):
    arr = np.zeros((3, 3, 3, 3, 3), dtype=np.float64)
    arr[(0, 0, 0, 0, 0)] = Txxxxx(x,y,z)
    arr[(0, 0, 0, 0, 1)] = Txxxxy(x,y,z)
    arr[(0, 0, 0, 0, 2)] = Txxxxz(x,y,z)
    arr[(0, 0, 0, 1, 0)] = Txxxyx(x,y,z)
    arr[(0, 0, 0, 1, 1)] = Txxxyy(x,y,z)
    arr[(0, 0, 0, 1, 2)] = Txxxyz(x,y,z)
    arr[(0, 0, 0, 2, 0)] = Txxxzx(x,y,z)
    arr[(0, 0, 0, 2, 1)] = Txxxzy(x,y,z)
    arr[(0, 0, 0, 2, 2)] = Txxxzz(x,y,z)
    arr[(0, 0, 1, 0, 0)] = Txxyxx(x,y,z)
    arr[(0, 0, 1, 0, 1)] = Txxyxy(x,y,z)
    arr[(0, 0, 1, 0, 2)] = Txxyxz(x,y,z)
    arr[(0, 0, 1, 1, 0)] = Txxyyx(x,y,z)
    arr[(0, 0, 1, 1, 1)] = Txxyyy(x,y,z)
    arr[(0, 0, 1, 1, 2)] = Txxyyz(x,y,z)
    arr[(0, 0, 1, 2, 0)] = Txxyzx(x,y,z)
    arr[(0, 0, 1, 2, 1)] = Txxyzy(x,y,z)
    arr[(0, 0, 1, 2, 2)] = Txxyzz(x,y,z)
    arr[(0, 0, 2, 0, 0)] = Txxzxx(x,y,z)
    arr[(0, 0, 2, 0, 1)] = Txxzxy(x,y,z)
    arr[(0, 0, 2, 0, 2)] = Txxzxz(x,y,z)
    arr[(0, 0, 2, 1, 0)] = Txxzyx(x,y,z)
    arr[(0, 0, 2, 1, 1)] = Txxzyy(x,y,z)
    arr[(0, 0, 2, 1, 2)] = Txxzyz(x,y,z)
    arr[(0, 0, 2, 2, 0)] = Txxzzx(x,y,z)
    arr[(0, 0, 2, 2, 1)] = Txxzzy(x,y,z)
    arr[(0, 0, 2, 2, 2)] = Txxzzz(x,y,z)
    arr[(0, 1, 0, 0, 0)] = Txyxxx(x,y,z)
    arr[(0, 1, 0, 0, 1)] = Txyxxy(x,y,z)
    arr[(0, 1, 0, 0, 2)] = Txyxxz(x,y,z)
    arr[(0, 1, 0, 1, 0)] = Txyxyx(x,y,z)
    arr[(0, 1, 0, 1, 1)] = Txyxyy(x,y,z)
    arr[(0, 1, 0, 1, 2)] = Txyxyz(x,y,z)
    arr[(0, 1, 0, 2, 0)] = Txyxzx(x,y,z)
    arr[(0, 1, 0, 2, 1)] = Txyxzy(x,y,z)
    arr[(0, 1, 0, 2, 2)] = Txyxzz(x,y,z)
    arr[(0, 1, 1, 0, 0)] = Txyyxx(x,y,z)
    arr[(0, 1, 1, 0, 1)] = Txyyxy(x,y,z)
    arr[(0, 1, 1, 0, 2)] = Txyyxz(x,y,z)
    arr[(0, 1, 1, 1, 0)] = Txyyyx(x,y,z)
    arr[(0, 1, 1, 1, 1)] = Txyyyy(x,y,z)
    arr[(0, 1, 1, 1, 2)] = Txyyyz(x,y,z)
    arr[(0, 1, 1, 2, 0)] = Txyyzx(x,y,z)
    arr[(0, 1, 1, 2, 1)] = Txyyzy(x,y,z)
    arr[(0, 1, 1, 2, 2)] = Txyyzz(x,y,z)
    arr[(0, 1, 2, 0, 0)] = Txyzxx(x,y,z)
    arr[(0, 1, 2, 0, 1)] = Txyzxy(x,y,z)
    arr[(0, 1, 2, 0, 2)] = Txyzxz(x,y,z)
    arr[(0, 1, 2, 1, 0)] = Txyzyx(x,y,z)
    arr[(0, 1, 2, 1, 1)] = Txyzyy(x,y,z)
    arr[(0, 1, 2, 1, 2)] = Txyzyz(x,y,z)
    arr[(0, 1, 2, 2, 0)] = Txyzzx(x,y,z)
    arr[(0, 1, 2, 2, 1)] = Txyzzy(x,y,z)
    arr[(0, 1, 2, 2, 2)] = Txyzzz(x,y,z)
    arr[(0, 2, 0, 0, 0)] = Txzxxx(x,y,z)
    arr[(0, 2, 0, 0, 1)] = Txzxxy(x,y,z)
    arr[(0, 2, 0, 0, 2)] = Txzxxz(x,y,z)
    arr[(0, 2, 0, 1, 0)] = Txzxyx(x,y,z)
    arr[(0, 2, 0, 1, 1)] = Txzxyy(x,y,z)
    arr[(0, 2, 0, 1, 2)] = Txzxyz(x,y,z)
    arr[(0, 2, 0, 2, 0)] = Txzxzx(x,y,z)
    arr[(0, 2, 0, 2, 1)] = Txzxzy(x,y,z)
    arr[(0, 2, 0, 2, 2)] = Txzxzz(x,y,z)
    arr[(0, 2, 1, 0, 0)] = Txzyxx(x,y,z)
    arr[(0, 2, 1, 0, 1)] = Txzyxy(x,y,z)
    arr[(0, 2, 1, 0, 2)] = Txzyxz(x,y,z)
    arr[(0, 2, 1, 1, 0)] = Txzyyx(x,y,z)
    arr[(0, 2, 1, 1, 1)] = Txzyyy(x,y,z)
    arr[(0, 2, 1, 1, 2)] = Txzyyz(x,y,z)
    arr[(0, 2, 1, 2, 0)] = Txzyzx(x,y,z)
    arr[(0, 2, 1, 2, 1)] = Txzyzy(x,y,z)
    arr[(0, 2, 1, 2, 2)] = Txzyzz(x,y,z)
    arr[(0, 2, 2, 0, 0)] = Txzzxx(x,y,z)
    arr[(0, 2, 2, 0, 1)] = Txzzxy(x,y,z)
    arr[(0, 2, 2, 0, 2)] = Txzzxz(x,y,z)
    arr[(0, 2, 2, 1, 0)] = Txzzyx(x,y,z)
    arr[(0, 2, 2, 1, 1)] = Txzzyy(x,y,z)
    arr[(0, 2, 2, 1, 2)] = Txzzyz(x,y,z)
    arr[(0, 2, 2, 2, 0)] = Txzzzx(x,y,z)
    arr[(0, 2, 2, 2, 1)] = Txzzzy(x,y,z)
    arr[(0, 2, 2, 2, 2)] = Txzzzz(x,y,z)
    arr[(1, 0, 0, 0, 0)] = Tyxxxx(x,y,z)
    arr[(1, 0, 0, 0, 1)] = Tyxxxy(x,y,z)
    arr[(1, 0, 0, 0, 2)] = Tyxxxz(x,y,z)
    arr[(1, 0, 0, 1, 0)] = Tyxxyx(x,y,z)
    arr[(1, 0, 0, 1, 1)] = Tyxxyy(x,y,z)
    arr[(1, 0, 0, 1, 2)] = Tyxxyz(x,y,z)
    arr[(1, 0, 0, 2, 0)] = Tyxxzx(x,y,z)
    arr[(1, 0, 0, 2, 1)] = Tyxxzy(x,y,z)
    arr[(1, 0, 0, 2, 2)] = Tyxxzz(x,y,z)
    arr[(1, 0, 1, 0, 0)] = Tyxyxx(x,y,z)
    arr[(1, 0, 1, 0, 1)] = Tyxyxy(x,y,z)
    arr[(1, 0, 1, 0, 2)] = Tyxyxz(x,y,z)
    arr[(1, 0, 1, 1, 0)] = Tyxyyx(x,y,z)
    arr[(1, 0, 1, 1, 1)] = Tyxyyy(x,y,z)
    arr[(1, 0, 1, 1, 2)] = Tyxyyz(x,y,z)
    arr[(1, 0, 1, 2, 0)] = Tyxyzx(x,y,z)
    arr[(1, 0, 1, 2, 1)] = Tyxyzy(x,y,z)
    arr[(1, 0, 1, 2, 2)] = Tyxyzz(x,y,z)
    arr[(1, 0, 2, 0, 0)] = Tyxzxx(x,y,z)
    arr[(1, 0, 2, 0, 1)] = Tyxzxy(x,y,z)
    arr[(1, 0, 2, 0, 2)] = Tyxzxz(x,y,z)
    arr[(1, 0, 2, 1, 0)] = Tyxzyx(x,y,z)
    arr[(1, 0, 2, 1, 1)] = Tyxzyy(x,y,z)
    arr[(1, 0, 2, 1, 2)] = Tyxzyz(x,y,z)
    arr[(1, 0, 2, 2, 0)] = Tyxzzx(x,y,z)
    arr[(1, 0, 2, 2, 1)] = Tyxzzy(x,y,z)
    arr[(1, 0, 2, 2, 2)] = Tyxzzz(x,y,z)
    arr[(1, 1, 0, 0, 0)] = Tyyxxx(x,y,z)
    arr[(1, 1, 0, 0, 1)] = Tyyxxy(x,y,z)
    arr[(1, 1, 0, 0, 2)] = Tyyxxz(x,y,z)
    arr[(1, 1, 0, 1, 0)] = Tyyxyx(x,y,z)
    arr[(1, 1, 0, 1, 1)] = Tyyxyy(x,y,z)
    arr[(1, 1, 0, 1, 2)] = Tyyxyz(x,y,z)
    arr[(1, 1, 0, 2, 0)] = Tyyxzx(x,y,z)
    arr[(1, 1, 0, 2, 1)] = Tyyxzy(x,y,z)
    arr[(1, 1, 0, 2, 2)] = Tyyxzz(x,y,z)
    arr[(1, 1, 1, 0, 0)] = Tyyyxx(x,y,z)
    arr[(1, 1, 1, 0, 1)] = Tyyyxy(x,y,z)
    arr[(1, 1, 1, 0, 2)] = Tyyyxz(x,y,z)
    arr[(1, 1, 1, 1, 0)] = Tyyyyx(x,y,z)
    arr[(1, 1, 1, 1, 1)] = Tyyyyy(x,y,z)
    arr[(1, 1, 1, 1, 2)] = Tyyyyz(x,y,z)
    arr[(1, 1, 1, 2, 0)] = Tyyyzx(x,y,z)
    arr[(1, 1, 1, 2, 1)] = Tyyyzy(x,y,z)
    arr[(1, 1, 1, 2, 2)] = Tyyyzz(x,y,z)
    arr[(1, 1, 2, 0, 0)] = Tyyzxx(x,y,z)
    arr[(1, 1, 2, 0, 1)] = Tyyzxy(x,y,z)
    arr[(1, 1, 2, 0, 2)] = Tyyzxz(x,y,z)
    arr[(1, 1, 2, 1, 0)] = Tyyzyx(x,y,z)
    arr[(1, 1, 2, 1, 1)] = Tyyzyy(x,y,z)
    arr[(1, 1, 2, 1, 2)] = Tyyzyz(x,y,z)
    arr[(1, 1, 2, 2, 0)] = Tyyzzx(x,y,z)
    arr[(1, 1, 2, 2, 1)] = Tyyzzy(x,y,z)
    arr[(1, 1, 2, 2, 2)] = Tyyzzz(x,y,z)
    arr[(1, 2, 0, 0, 0)] = Tyzxxx(x,y,z)
    arr[(1, 2, 0, 0, 1)] = Tyzxxy(x,y,z)
    arr[(1, 2, 0, 0, 2)] = Tyzxxz(x,y,z)
    arr[(1, 2, 0, 1, 0)] = Tyzxyx(x,y,z)
    arr[(1, 2, 0, 1, 1)] = Tyzxyy(x,y,z)
    arr[(1, 2, 0, 1, 2)] = Tyzxyz(x,y,z)
    arr[(1, 2, 0, 2, 0)] = Tyzxzx(x,y,z)
    arr[(1, 2, 0, 2, 1)] = Tyzxzy(x,y,z)
    arr[(1, 2, 0, 2, 2)] = Tyzxzz(x,y,z)
    arr[(1, 2, 1, 0, 0)] = Tyzyxx(x,y,z)
    arr[(1, 2, 1, 0, 1)] = Tyzyxy(x,y,z)
    arr[(1, 2, 1, 0, 2)] = Tyzyxz(x,y,z)
    arr[(1, 2, 1, 1, 0)] = Tyzyyx(x,y,z)
    arr[(1, 2, 1, 1, 1)] = Tyzyyy(x,y,z)
    arr[(1, 2, 1, 1, 2)] = Tyzyyz(x,y,z)
    arr[(1, 2, 1, 2, 0)] = Tyzyzx(x,y,z)
    arr[(1, 2, 1, 2, 1)] = Tyzyzy(x,y,z)
    arr[(1, 2, 1, 2, 2)] = Tyzyzz(x,y,z)
    arr[(1, 2, 2, 0, 0)] = Tyzzxx(x,y,z)
    arr[(1, 2, 2, 0, 1)] = Tyzzxy(x,y,z)
    arr[(1, 2, 2, 0, 2)] = Tyzzxz(x,y,z)
    arr[(1, 2, 2, 1, 0)] = Tyzzyx(x,y,z)
    arr[(1, 2, 2, 1, 1)] = Tyzzyy(x,y,z)
    arr[(1, 2, 2, 1, 2)] = Tyzzyz(x,y,z)
    arr[(1, 2, 2, 2, 0)] = Tyzzzx(x,y,z)
    arr[(1, 2, 2, 2, 1)] = Tyzzzy(x,y,z)
    arr[(1, 2, 2, 2, 2)] = Tyzzzz(x,y,z)
    arr[(2, 0, 0, 0, 0)] = Tzxxxx(x,y,z)
    arr[(2, 0, 0, 0, 1)] = Tzxxxy(x,y,z)
    arr[(2, 0, 0, 0, 2)] = Tzxxxz(x,y,z)
    arr[(2, 0, 0, 1, 0)] = Tzxxyx(x,y,z)
    arr[(2, 0, 0, 1, 1)] = Tzxxyy(x,y,z)
    arr[(2, 0, 0, 1, 2)] = Tzxxyz(x,y,z)
    arr[(2, 0, 0, 2, 0)] = Tzxxzx(x,y,z)
    arr[(2, 0, 0, 2, 1)] = Tzxxzy(x,y,z)
    arr[(2, 0, 0, 2, 2)] = Tzxxzz(x,y,z)
    arr[(2, 0, 1, 0, 0)] = Tzxyxx(x,y,z)
    arr[(2, 0, 1, 0, 1)] = Tzxyxy(x,y,z)
    arr[(2, 0, 1, 0, 2)] = Tzxyxz(x,y,z)
    arr[(2, 0, 1, 1, 0)] = Tzxyyx(x,y,z)
    arr[(2, 0, 1, 1, 1)] = Tzxyyy(x,y,z)
    arr[(2, 0, 1, 1, 2)] = Tzxyyz(x,y,z)
    arr[(2, 0, 1, 2, 0)] = Tzxyzx(x,y,z)
    arr[(2, 0, 1, 2, 1)] = Tzxyzy(x,y,z)
    arr[(2, 0, 1, 2, 2)] = Tzxyzz(x,y,z)
    arr[(2, 0, 2, 0, 0)] = Tzxzxx(x,y,z)
    arr[(2, 0, 2, 0, 1)] = Tzxzxy(x,y,z)
    arr[(2, 0, 2, 0, 2)] = Tzxzxz(x,y,z)
    arr[(2, 0, 2, 1, 0)] = Tzxzyx(x,y,z)
    arr[(2, 0, 2, 1, 1)] = Tzxzyy(x,y,z)
    arr[(2, 0, 2, 1, 2)] = Tzxzyz(x,y,z)
    arr[(2, 0, 2, 2, 0)] = Tzxzzx(x,y,z)
    arr[(2, 0, 2, 2, 1)] = Tzxzzy(x,y,z)
    arr[(2, 0, 2, 2, 2)] = Tzxzzz(x,y,z)
    arr[(2, 1, 0, 0, 0)] = Tzyxxx(x,y,z)
    arr[(2, 1, 0, 0, 1)] = Tzyxxy(x,y,z)
    arr[(2, 1, 0, 0, 2)] = Tzyxxz(x,y,z)
    arr[(2, 1, 0, 1, 0)] = Tzyxyx(x,y,z)
    arr[(2, 1, 0, 1, 1)] = Tzyxyy(x,y,z)
    arr[(2, 1, 0, 1, 2)] = Tzyxyz(x,y,z)
    arr[(2, 1, 0, 2, 0)] = Tzyxzx(x,y,z)
    arr[(2, 1, 0, 2, 1)] = Tzyxzy(x,y,z)
    arr[(2, 1, 0, 2, 2)] = Tzyxzz(x,y,z)
    arr[(2, 1, 1, 0, 0)] = Tzyyxx(x,y,z)
    arr[(2, 1, 1, 0, 1)] = Tzyyxy(x,y,z)
    arr[(2, 1, 1, 0, 2)] = Tzyyxz(x,y,z)
    arr[(2, 1, 1, 1, 0)] = Tzyyyx(x,y,z)
    arr[(2, 1, 1, 1, 1)] = Tzyyyy(x,y,z)
    arr[(2, 1, 1, 1, 2)] = Tzyyyz(x,y,z)
    arr[(2, 1, 1, 2, 0)] = Tzyyzx(x,y,z)
    arr[(2, 1, 1, 2, 1)] = Tzyyzy(x,y,z)
    arr[(2, 1, 1, 2, 2)] = Tzyyzz(x,y,z)
    arr[(2, 1, 2, 0, 0)] = Tzyzxx(x,y,z)
    arr[(2, 1, 2, 0, 1)] = Tzyzxy(x,y,z)
    arr[(2, 1, 2, 0, 2)] = Tzyzxz(x,y,z)
    arr[(2, 1, 2, 1, 0)] = Tzyzyx(x,y,z)
    arr[(2, 1, 2, 1, 1)] = Tzyzyy(x,y,z)
    arr[(2, 1, 2, 1, 2)] = Tzyzyz(x,y,z)
    arr[(2, 1, 2, 2, 0)] = Tzyzzx(x,y,z)
    arr[(2, 1, 2, 2, 1)] = Tzyzzy(x,y,z)
    arr[(2, 1, 2, 2, 2)] = Tzyzzz(x,y,z)
    arr[(2, 2, 0, 0, 0)] = Tzzxxx(x,y,z)
    arr[(2, 2, 0, 0, 1)] = Tzzxxy(x,y,z)
    arr[(2, 2, 0, 0, 2)] = Tzzxxz(x,y,z)
    arr[(2, 2, 0, 1, 0)] = Tzzxyx(x,y,z)
    arr[(2, 2, 0, 1, 1)] = Tzzxyy(x,y,z)
    arr[(2, 2, 0, 1, 2)] = Tzzxyz(x,y,z)
    arr[(2, 2, 0, 2, 0)] = Tzzxzx(x,y,z)
    arr[(2, 2, 0, 2, 1)] = Tzzxzy(x,y,z)
    arr[(2, 2, 0, 2, 2)] = Tzzxzz(x,y,z)
    arr[(2, 2, 1, 0, 0)] = Tzzyxx(x,y,z)
    arr[(2, 2, 1, 0, 1)] = Tzzyxy(x,y,z)
    arr[(2, 2, 1, 0, 2)] = Tzzyxz(x,y,z)
    arr[(2, 2, 1, 1, 0)] = Tzzyyx(x,y,z)
    arr[(2, 2, 1, 1, 1)] = Tzzyyy(x,y,z)
    arr[(2, 2, 1, 1, 2)] = Tzzyyz(x,y,z)
    arr[(2, 2, 1, 2, 0)] = Tzzyzx(x,y,z)
    arr[(2, 2, 1, 2, 1)] = Tzzyzy(x,y,z)
    arr[(2, 2, 1, 2, 2)] = Tzzyzz(x,y,z)
    arr[(2, 2, 2, 0, 0)] = Tzzzxx(x,y,z)
    arr[(2, 2, 2, 0, 1)] = Tzzzxy(x,y,z)
    arr[(2, 2, 2, 0, 2)] = Tzzzxz(x,y,z)
    arr[(2, 2, 2, 1, 0)] = Tzzzyx(x,y,z)
    arr[(2, 2, 2, 1, 1)] = Tzzzyy(x,y,z)
    arr[(2, 2, 2, 1, 2)] = Tzzzyz(x,y,z)
    arr[(2, 2, 2, 2, 0)] = Tzzzzx(x,y,z)
    arr[(2, 2, 2, 2, 1)] = Tzzzzy(x,y,z)
    arr[(2, 2, 2, 2, 2)] = Tzzzzz(x,y,z)
    return arr
T = [T0, T1, T2, T3, T4, T5, ]
