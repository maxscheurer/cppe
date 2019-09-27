from numba import jit
from numpy import sqrt, exp
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

def T_damp_thole_0(x, y, z, a):
    shape = (3,)*0
    result = np.zeros(shape, dtype=np.float64)
    x0 = sqrt(x ** 2 + y ** 2 + z ** 2)
    x1 = a * x0
    result[...,] = -((x1 + 2) * exp(-x1) / 2 - 1) / x0
    return result


def T_damp_thole_1(x, y, z, a):
    shape = (3,)*1
    result = np.zeros(shape, dtype=np.float64)
    x0 = x ** 2 + y ** 2 + z ** 2
    x1 = a * sqrt(x0)
    x2 = -x1
    x3 = exp(x2)
    x4 = -a * x3 * (x2 - 1) / (2 * x0) + (x3 * (x1 + 2) - 2) / (2 * x0 ** (3 / 2))
    result[0] = x * x4
    result[1] = x4 * y
    result[2] = x4 * z
    return result


def T_damp_thole_2(x, y, z, a):
    shape = (3,)*2
    result = np.zeros(shape, dtype=np.float64)
    x0 = x ** 2
    x1 = y ** 2
    x2 = z ** 2
    x3 = x0 + x1 + x2
    x4 = sqrt(x3)
    x5 = a * x4
    x6 = exp(-x5)
    x7 = a * x6
    x8 = x7 * (x5 + 1) / x3 ** 2
    x9 = 1 / x3
    x10 = 3 * x9
    x11 = x3 ** (-3 / 2)
    x12 = x5 + 2
    x13 = x12 * x6 - 2
    x14 = x11 * x13 / 2
    x15 = x0 * x11
    x16 = a * x9
    x17 = 2 * x16
    x18 = x12 * x16
    x19 = 1 / x4
    x20 = -x12 * x19 + x19
    x21 = x19 * x7 / 2
    x22 = x11 * x12
    x23 = 3 * x13 / (2 * x3 ** (5 / 2)) + x21 * (-x11 - x17 + x18 + x22) + x8
    x24 = x * x23
    x25 = -x24 * y
    x26 = -x24 * z
    x27 = -x23 * y * z
    result[0, 0] = (
        -x0 * x8
        - x14 * (x0 * x10 - 1)
        - x21 * (-x0 * x17 + x0 * x18 + x12 * x15 - x15 + x20)
    )
    result[0, 1] = x25
    result[0, 2] = x26
    result[1, 0] = x25
    result[1, 1] = (
        -x1 * x8
        - x14 * (x1 * x10 - 1)
        - x21 * (-x1 * x11 - x1 * x17 + x1 * x18 + x1 * x22 + x20)
    )
    result[1, 2] = x27
    result[2, 0] = x26
    result[2, 1] = x27
    result[2, 2] = (
        -x14 * (x10 * x2 - 1)
        - x2 * x8
        - x21 * (-x11 * x2 - x17 * x2 + x18 * x2 + x2 * x22 + x20)
    )
    return result


def T_damp_thole_3(x, y, z, a):
    shape = (3,)*3
    result = np.zeros(shape, dtype=np.float64)
    x0 = x ** 2
    x1 = y ** 2
    x2 = z ** 2
    x3 = x0 + x1 + x2
    x4 = 1 / x3
    x5 = x0 * x4
    x6 = 3 * x5 - 1
    x7 = sqrt(x3)
    x8 = a * x7
    x9 = exp(-x8)
    x10 = x8 + 1
    x11 = a / x3 ** 2
    x12 = x10 * x11
    x13 = 3 * x12 * x9
    x14 = x8 + 2
    x15 = x14 * x9 - 2
    x16 = x3 ** (-5 / 2)
    x17 = 3 * x16
    x18 = x15 * x17
    x19 = x3 ** (-3 / 2)
    x20 = 3 * x19
    x21 = x0 * x19
    x22 = x14 * x21
    x23 = a * x4
    x24 = 2 * x23
    x25 = x14 * x23
    x26 = 1 / x7
    x27 = -x14 * x26 + x26
    x28 = -x0 * x24 + x0 * x25 - x21 + x22 + x27
    x29 = a * x9
    x30 = x28 * x29
    x31 = x0 * x17
    x32 = 6 * x11
    x33 = a ** 2
    x34 = x20 * x33
    x35 = x14 * x31
    x36 = x11 * x14
    x37 = 3 * x36
    x38 = -x14 * x20 + x20 + 6 * x23 - 3 * x25
    x39 = x26 * x29
    x40 = x / 2
    x41 = x10 * x29 / x3 ** 3
    x42 = 3 * x41
    x43 = x15 / x3 ** (7 / 2)
    x44 = 3 * x43
    x45 = x9 / 2
    x46 = x12 * x45
    x47 = 3 * x15 * x16 / 2
    x48 = x14 * x19
    x49 = -x19 + x48
    x50 = x29 * (-x24 + x25 + x49)
    x51 = x19 / 2
    x52 = x33 * x4 * x45
    x53 = 5 * x11
    x54 = 2 * x36
    x55 = -x23 + x49
    x56 = x39 / 2
    x57 = (
        x0 * x42
        + x0 * x44
        + x21 * x50
        + x28 * x52
        + x30 * x51
        + x46 * x6
        + x47 * x6
        - x56 * (x0 * x53 - x0 * x54 + x21 * x33 + x31 - x35 + x55)
    )
    x58 = x57 * y
    x59 = x57 * z
    x60 = 15 * x43 / 2
    x61 = x19 + x24 - x25 - x48
    x62 = x1 * x17
    x63 = x14 * x62
    x64 = x1 * x48
    x65 = -x1 * x32 - x1 * x34 + x1 * x37 + x33 * x64 - x62 + x63
    x66 = 9 * x41 / 2
    x67 = -x46 - x47
    x68 = x1 * x19
    x69 = -x1 * x24 + x1 * x25 + x27 + x64 - x68
    x70 = x29 * x51
    x71 = x50 * x68 + x69 * x70
    x72 = x * (x1 * x60 + x1 * x66 + x56 * (x61 + x65) + x67 + x71)
    x73 = x14 * x17
    x74 = x33 * x48
    x75 = (
        x40
        * y
        * z
        * (x20 * x50 + x39 * (-x17 - x32 - x34 + x37 + x73 + x74) + 9 * x41 + 15 * x43)
    )
    x76 = x19 * x2
    x77 = -x2 * x24 + x2 * x25 + x2 * x48 + x27 - x76
    x78 = -x17 * x2 - x2 * x32 - x2 * x34 + x2 * x37 + x2 * x73 + x2 * x74
    x79 = x2 * x60 + x2 * x66 + x50 * x76 + x56 * (x61 + x78) + x67 + x70 * x77
    x80 = x * x79
    x81 = x1 * x4
    x82 = 3 * x81 - 1
    x83 = x20 * x29
    x84 = z * (
        x1 * x42
        + x1 * x44
        + x46 * x82
        + x47 * x82
        + x52 * x69
        - x56 * (x1 * x53 - x1 * x54 + x33 * x68 + x55 + x62 - x63)
        + x71
    )
    x85 = x79 * y
    x86 = x2 * x4
    result[0, 0, 0] = x40 * (
        x13 * x6
        + x18 * (5 * x5 - 3)
        + x20 * x30
        + x39 * (-x0 * x32 - x0 * x34 + x0 * x37 + x22 * x33 - x31 + x35 + x38)
    )
    result[0, 0, 1] = x58
    result[0, 0, 2] = x59
    result[0, 1, 0] = x58
    result[0, 1, 1] = x72
    result[0, 1, 2] = x75
    result[0, 2, 0] = x59
    result[0, 2, 1] = x75
    result[0, 2, 2] = x80
    result[1, 0, 0] = x58
    result[1, 0, 1] = x72
    result[1, 0, 2] = x75
    result[1, 1, 0] = x72
    result[1, 1, 1] = (
        y * (x13 * x82 + x18 * (5 * x81 - 3) + x39 * (x38 + x65) + x69 * x83) / 2
    )
    result[1, 1, 2] = x84
    result[1, 2, 0] = x75
    result[1, 2, 1] = x84
    result[1, 2, 2] = x85
    result[2, 0, 0] = x59
    result[2, 0, 1] = x75
    result[2, 0, 2] = x80
    result[2, 1, 0] = x75
    result[2, 1, 1] = x84
    result[2, 1, 2] = x85
    result[2, 2, 0] = x80
    result[2, 2, 1] = x85
    result[2, 2, 2] = (
        z
        * (x13 * (3 * x86 - 1) + x18 * (5 * x86 - 3) + x39 * (x38 + x78) + x77 * x83)
        / 2
    )
    return result


def T_damp_thole_4(x, y, z, a):
    shape = (3,)*4
    result = np.zeros(shape, dtype=np.float64)
    x0 = x ** 2
    x1 = y ** 2
    x2 = z ** 2
    x3 = x0 + x1 + x2
    x4 = 1 / x3
    x5 = x0 * x4
    x6 = 5 * x5 - 3
    x7 = a / x3 ** 3
    x8 = x0 * x7
    x9 = sqrt(x3)
    x10 = a * x9
    x11 = exp(-x10)
    x12 = x10 + 1
    x13 = x11 * x12
    x14 = 6 * x13
    x15 = x ** 4
    x16 = x3 ** (-2)
    x17 = 35 * x16
    x18 = x10 + 2
    x19 = x11 * x18 - 2
    x20 = x3 ** (-5 / 2)
    x21 = 3 * x20 / 2
    x22 = x19 * x21
    x23 = 3 * x5 - 1
    x24 = x3 ** (-3 / 2)
    x25 = 3 * x24
    x26 = x0 * x24
    x27 = x18 * x26
    x28 = a * x4
    x29 = 2 * x28
    x30 = x18 * x28
    x31 = 1 / x9
    x32 = -x18 * x31 + x31
    x33 = -x0 * x29 + x0 * x30 - x26 + x27 + x32
    x34 = x11 * x33
    x35 = a * x34
    x36 = a * x16
    x37 = 6 * x36
    x38 = a ** 2
    x39 = x25 * x38
    x40 = x27 * x38
    x41 = 3 * x36
    x42 = x18 * x41
    x43 = 3 * x20
    x44 = x0 * x43
    x45 = x18 * x44
    x46 = -x44 + x45
    x47 = x18 * x25
    x48 = 6 * x28
    x49 = 3 * x30
    x50 = x25 - x47 + x48 - x49
    x51 = -x0 * x37 - x0 * x39 + x0 * x42 + x40 + x46 + x50
    x52 = a * x11
    x53 = x26 * x52
    x54 = x3 ** (-7 / 2)
    x55 = 15 * x54
    x56 = x15 * x55
    x57 = 18 * x20
    x58 = x0 * x57
    x59 = 30 * x7
    x60 = x38 * x57
    x61 = a ** 3
    x62 = x16 * x61
    x63 = x15 * x62
    x64 = x26 * x38
    x65 = x0 * x36
    x66 = 18 * x36
    x67 = x18 * x66
    x68 = x15 * x18
    x69 = x20 * x38
    x70 = 6 * x69
    x71 = 15 * x7
    x72 = -x25 + x47 - x48 + x49
    x73 = x11 / 2
    x74 = a * x31 * x73
    x75 = 9 * x12
    x76 = x3 ** (-4)
    x77 = x0 * x52
    x78 = x76 * x77
    x79 = x19 / x3 ** (9 / 2)
    x80 = x0 * x79
    x81 = x11 * x7
    x82 = x12 * x81
    x83 = 9 * x82 / 2
    x84 = 3 * x82 / 2
    x85 = x19 * x54
    x86 = 15 * x85 / 2
    x87 = x18 * x24
    x88 = -x24 + x87
    x89 = -x29 + x30 + x88
    x90 = a * x24
    x91 = x11 * x90
    x92 = 3 * x91 / 2
    x93 = x89 * x92
    x94 = 9 * x20 / 2
    x95 = 3 * x16 * x38 / 2
    x96 = x34 * x95
    x97 = 5 * x65
    x98 = 2 * x18 * x65
    x99 = -x28 + x88
    x100 = x11 * (x44 - x45 + x64 + x97 - x98 + x99)
    x101 = x100 * x90
    x102 = x91 / 2
    x103 = x38 * x4
    x104 = x103 * x73
    x105 = 12 * x0
    x106 = x18 * x7
    x107 = x0 * x55
    x108 = x107 * x18
    x109 = x107 - x108
    x110 = -x39
    x111 = 9 * x20
    x112 = x111 * x18
    x113 = x110 - x111 + x112 + x18 * x37 - 15 * x36
    x114 = x * (
        -3 * x101 / 2
        + x102 * x51
        + x104 * x51
        + x23 * x83
        + x23 * x93
        + x35 * x94
        + x6 * x84
        + x6 * x86
        - x74
        * (x0 * x62 - x105 * x106 + x105 * x69 + x109 + x113 - x38 * x45 + 27 * x8)
        + x75 * x78
        + 15 * x80
        + x96
    )
    x115 = -x114 * y
    x116 = -x114 * z
    x117 = x1 * x41
    x118 = 5 * x69
    x119 = x0 * x118
    x120 = 23 * x8
    x121 = 8 * x18
    x122 = x121 * x8
    x123 = x24 - x87
    x124 = x1 * x43
    x125 = x124 * x18
    x126 = -x124 + x125
    x127 = x123 + x126 + x28
    x128 = x46 - x64 - x97 + x98
    x129 = 30 * x80
    x130 = x23 * x86
    x131 = x1 * x24
    x132 = a * x100
    x133 = x1 * x4
    x134 = x100 * x38
    x135 = x123 + x29 - x30
    x136 = x1 * x87
    x137 = x136 * x38
    x138 = -x1 * x37 - x1 * x39 + x117 * x18 + x126 + x137
    x139 = x135 + x138
    x140 = x21 * x35
    x141 = -x1 * x29 + x1 * x30 - x131 + x136 + x32
    x142 = x102 * x141
    x143 = x61 * x73
    x144 = x143 * x33
    x145 = 21 * x12
    x146 = x145 * x78
    x147 = x20 * x89
    x148 = 6 * x147
    x149 = x148 * x77
    x150 = x1 * x7
    x151 = 3 * x13
    x152 = x150 * x151
    x153 = 3 * x85
    x154 = x0 * x153 + x102 * x33 + x104 * x33 + x151 * x8 + x22 * x23
    x155 = (
        -x1 * x129
        - x1 * x130
        - x1 * x140
        - x1 * x146
        - x1 * x149
        - x1 * x96
        + x131 * x132
        - x131 * x144
        + x133 * x134
        - x139 * x53
        - x142 * x23
        - x152 * x23
        + x154
        + x74
        * (
            x1 * x107
            - x1 * x108
            + x1 * x119
            + x1 * x120
            - x1 * x122
            - x117
            + x127
            + x128
        )
    )
    x156 = 3 * x82
    x157 = x156 * x23
    x158 = x102 * x89
    x159 = x38 * x87
    x160 = x18 * x43
    x161 = x160 - x43
    x162 = x110 + x159 + x161 - x37 + x42
    x163 = y * z
    x164 = -x163 * (
        -x100 * x103
        - x101
        + x129
        + x130
        + x140
        + x144 * x24
        + x146
        + x149
        + x157
        + x158 * x23
        + x162 * x53
        - x74 * (x109 + x119 + x120 - x122 + x161 - x41)
        + x96
    )
    x165 = x160 * x2 - x2 * x43
    x166 = x165 - x2 * x41
    x167 = x2 * x24
    x168 = x2 * x4
    x169 = x159 * x2
    x170 = x165 + x169 - x2 * x37 - x2 * x39 + x2 * x42
    x171 = x135 + x170
    x172 = -x167 - x2 * x29 + x2 * x30 + x2 * x87 + x32
    x173 = x102 * x172
    x174 = (
        -x129 * x2
        - x130 * x2
        + x132 * x167
        + x134 * x168
        - x140 * x2
        - x144 * x167
        - x146 * x2
        - x149 * x2
        + x154
        - x157 * x2
        - x171 * x53
        - x173 * x23
        - x2 * x96
        + x74
        * (
            x107 * x2
            - x108 * x2
            + x119 * x2
            + x120 * x2
            - x122 * x2
            + x123
            + x128
            + x166
            + x28
        )
    )
    x175 = x18 * x36
    x176 = x111 - x112 - 9 * x175 + 9 * x24 * x38 - x38 * x47 + x66
    x177 = x1 * x55
    x178 = x1 * x57
    x179 = x1 * x62
    x180 = x177 * x18
    x181 = x1 * x69
    x182 = 15 * x1
    x183 = (
        -x1 * x59
        + x106 * x182
        - x177
        - x178 * x38
        + x179 * x18
        - 4 * x179
        + 6 * x18 * x181
        + x180
    )
    x184 = x139 * x91
    x185 = x1 * x52
    x186 = x89 * x94
    x187 = 105 * x79 / 2
    x188 = x185 * x76
    x189 = 30 * x12
    x190 = x1 * x187 + x188 * x189
    x191 = x138 + x50
    x192 = x141 * x52
    x193 = x102 * x191 + x192 * x94
    x194 = -x75 * x81 - 45 * x85 / 2 - x93
    x195 = x * y
    x196 = -x195 * (
        3 * x184 / 2 + x185 * x186 + x190 + x193 + x194 + x74 * (x176 + x183)
    )
    x197 = -x159 - x160 + x37 + x39 - x42 + x43
    x198 = x192 * x21
    x199 = x162 * x52
    x200 = 15 * x147 / 2
    x201 = -x156 - x158 - x86
    x202 = x * z
    x203 = -x202 * (
        x131 * x199 + x184 + x185 * x200 + x190 + x198 + x201 + x74 * (x183 + x197)
    )
    x204 = x171 * x91
    x205 = x2 * x55
    x206 = x2 * x57
    x207 = x2 * x62
    x208 = x18 * x2
    x209 = (
        x18 * x205
        + x18 * x207
        - x2 * x59
        - x205
        - x206 * x38
        - 4 * x207
        + x208 * x70
        + x208 * x71
    )
    x210 = x172 * x52
    x211 = x2 * x52
    x212 = x187 * x2 + x189 * x211 * x76
    x213 = -x195 * (
        x167 * x199
        + x200 * x211
        + x201
        + x204
        + x21 * x210
        + x212
        + x74 * (x197 + x209)
    )
    x214 = x170 + x50
    x215 = (
        x102 * x214
        + x186 * x211
        + x194
        + 3 * x204 / 2
        + x210 * x94
        + x212
        + x74 * (x176 + x209)
    )
    x216 = -x202 * x215
    x217 = 5 * x133 - 3
    x218 = y ** 4
    x219 = 3 * x133 - 1
    x220 = x131 * x52
    x221 = x218 * x55
    x222 = x218 * x62
    x223 = x131 * x38
    x224 = x1 * x36
    x225 = x18 * x218
    x226 = 5 * x224
    x227 = 2 * x1 * x175
    x228 = x124 - x125 + x223 + x226 - x227 + x99
    x229 = x11 * x141 * x95
    x230 = -x163 * (
        x104 * x191
        + x182 * x79
        + x188 * x75
        + x193
        + x217 * x84
        + x217 * x86
        + x219 * x83
        + x219 * x93
        - x228 * x92
        + x229
        - x74
        * (
            x113
            - x125 * x38
            - 12 * x150 * x18
            + 27 * x150
            + x177
            + x179
            - x180
            + 12 * x181
        )
    )
    x231 = x1 * x2
    x232 = x231 * x52
    x233 = x2 * x219
    x234 = x167 * x52
    x235 = x150 * x2
    x236 = (
        x1 * x153
        + x104 * x141
        + x11 * x168 * x228 * x38
        - x141 * x143 * x167
        + x142
        - x145 * x232 * x76
        - x148 * x232
        + x152
        - x156 * x233
        - x171 * x220
        - x173 * x219
        - x198 * x2
        - x2 * x229
        + x219 * x22
        + x228 * x234
        - 30 * x231 * x79
        - x233 * x86
        + x74
        * (
            x118 * x231
            - x121 * x235
            + x127
            + x166
            + x177 * x2
            - x180 * x2
            - x223
            - x226
            + x227
            + 23 * x235
        )
    )
    x237 = -x163 * x215
    x238 = z ** 4
    x239 = x238 * x55
    x240 = x238 * x62
    x241 = x18 * x238
    result[0, 0, 0, 0] = (
        -x14 * x6 * x8
        - x22 * (x15 * x17 - 30 * x5 + 3)
        - x23 * x25 * x35
        - 2 * x51 * x53
        - x74
        * (
            -x0 * x67
            - x15 * x59
            - x15 * x60
            + x18 * x56
            - x18 * x58
            + x18 * x63
            - 6 * x40
            - x56
            + x58
            - 4 * x63
            + 18 * x64
            + 36 * x65
            + x68 * x70
            + x68 * x71
            + x72
        )
    )
    result[0, 0, 0, 1] = x115
    result[0, 0, 0, 2] = x116
    result[0, 0, 1, 0] = x115
    result[0, 0, 1, 1] = x155
    result[0, 0, 1, 2] = x164
    result[0, 0, 2, 0] = x116
    result[0, 0, 2, 1] = x164
    result[0, 0, 2, 2] = x174
    result[0, 1, 0, 0] = x115
    result[0, 1, 0, 1] = x155
    result[0, 1, 0, 2] = x164
    result[0, 1, 1, 0] = x155
    result[0, 1, 1, 1] = x196
    result[0, 1, 1, 2] = x203
    result[0, 1, 2, 0] = x164
    result[0, 1, 2, 1] = x203
    result[0, 1, 2, 2] = x213
    result[0, 2, 0, 0] = x116
    result[0, 2, 0, 1] = x164
    result[0, 2, 0, 2] = x174
    result[0, 2, 1, 0] = x164
    result[0, 2, 1, 1] = x203
    result[0, 2, 1, 2] = x213
    result[0, 2, 2, 0] = x174
    result[0, 2, 2, 1] = x213
    result[0, 2, 2, 2] = x216
    result[1, 0, 0, 0] = x115
    result[1, 0, 0, 1] = x155
    result[1, 0, 0, 2] = x164
    result[1, 0, 1, 0] = x155
    result[1, 0, 1, 1] = x196
    result[1, 0, 1, 2] = x203
    result[1, 0, 2, 0] = x164
    result[1, 0, 2, 1] = x203
    result[1, 0, 2, 2] = x213
    result[1, 1, 0, 0] = x155
    result[1, 1, 0, 1] = x196
    result[1, 1, 0, 2] = x203
    result[1, 1, 1, 0] = x196
    result[1, 1, 1, 1] = (
        -x14 * x150 * x217
        - 2 * x191 * x220
        - x192 * x219 * x25
        - x22 * (-30 * x133 + x17 * x218 + 3)
        - x74
        * (
            -x1 * x67
            - 6 * x137
            - x178 * x18
            + x178
            + x18 * x221
            + x18 * x222
            - x218 * x59
            - x218 * x60
            - x221
            - 4 * x222
            + 18 * x223
            + 36 * x224
            + x225 * x70
            + x225 * x71
            + x72
        )
    )
    result[1, 1, 1, 2] = x230
    result[1, 1, 2, 0] = x203
    result[1, 1, 2, 1] = x230
    result[1, 1, 2, 2] = x236
    result[1, 2, 0, 0] = x164
    result[1, 2, 0, 1] = x203
    result[1, 2, 0, 2] = x213
    result[1, 2, 1, 0] = x203
    result[1, 2, 1, 1] = x230
    result[1, 2, 1, 2] = x236
    result[1, 2, 2, 0] = x213
    result[1, 2, 2, 1] = x236
    result[1, 2, 2, 2] = x237
    result[2, 0, 0, 0] = x116
    result[2, 0, 0, 1] = x164
    result[2, 0, 0, 2] = x174
    result[2, 0, 1, 0] = x164
    result[2, 0, 1, 1] = x203
    result[2, 0, 1, 2] = x213
    result[2, 0, 2, 0] = x174
    result[2, 0, 2, 1] = x213
    result[2, 0, 2, 2] = x216
    result[2, 1, 0, 0] = x164
    result[2, 1, 0, 1] = x203
    result[2, 1, 0, 2] = x213
    result[2, 1, 1, 0] = x203
    result[2, 1, 1, 1] = x230
    result[2, 1, 1, 2] = x236
    result[2, 1, 2, 0] = x213
    result[2, 1, 2, 1] = x236
    result[2, 1, 2, 2] = x237
    result[2, 2, 0, 0] = x174
    result[2, 2, 0, 1] = x213
    result[2, 2, 0, 2] = x216
    result[2, 2, 1, 0] = x213
    result[2, 2, 1, 1] = x236
    result[2, 2, 1, 2] = x237
    result[2, 2, 2, 0] = x216
    result[2, 2, 2, 1] = x237
    result[2, 2, 2, 2] = (
        -6 * x2 * x82 * (5 * x168 - 3)
        - x210 * x25 * (3 * x168 - 1)
        - 2 * x214 * x234
        - x22 * (-30 * x168 + x17 * x238 + 3)
        - x74
        * (
            18 * x167 * x38
            - 6 * x169
            - x18 * x206
            + x18 * x239
            + x18 * x240
            + 36 * x2 * x36
            - x2 * x67
            + x206
            - x238 * x59
            - x238 * x60
            - x239
            - 4 * x240
            + x241 * x70
            + x241 * x71
            + x72
        )
    )
    return result


def T_damp_thole_5(x, y, z, a):
    shape = (3,)*5
    result = np.zeros(shape, dtype=np.float64)
    x0 = x ** 2
    x1 = y ** 2
    x2 = z ** 2
    x3 = x0 + x1 + x2
    x4 = 1 / x3
    x5 = x0 * x4
    x6 = x ** 4
    x7 = x3 ** (-2)
    x8 = x6 * x7
    x9 = -30 * x5 + 35 * x8 + 3
    x10 = x3 ** (-3)
    x11 = a * x10
    x12 = sqrt(x3)
    x13 = a * x12
    x14 = exp(-x13)
    x15 = x13 + 1
    x16 = x14 * x15
    x17 = x11 * x16
    x18 = 15 * x17 / 2
    x19 = x3 ** (-7 / 2)
    x20 = x13 + 2
    x21 = x14 * x20 - 2
    x22 = x19 * x21
    x23 = 15 * x22 / 2
    x24 = x3 ** (-3 / 2)
    x25 = x0 * x24
    x26 = x20 * x25
    x27 = a * x4
    x28 = 2 * x27
    x29 = x20 * x27
    x30 = 1 / x12
    x31 = -x20 * x30 + x30
    x32 = -x0 * x28 + x0 * x29 - x25 + x26 + x31
    x33 = x14 * x32
    x34 = x3 ** (-5 / 2)
    x35 = a * x34
    x36 = 5 * x5 - 3
    x37 = 15 * x36
    x38 = 3 * x5 - 1
    x39 = a * x7
    x40 = 6 * x39
    x41 = a ** 2
    x42 = 3 * x24
    x43 = x41 * x42
    x44 = x26 * x41
    x45 = 3 * x39
    x46 = x20 * x45
    x47 = 3 * x34
    x48 = x0 * x47
    x49 = x20 * x48
    x50 = -x48 + x49
    x51 = x20 * x42
    x52 = 6 * x27
    x53 = 3 * x29
    x54 = x42 - x51 + x52 - x53
    x55 = -x0 * x40 - x0 * x43 + x0 * x46 + x44 + x50 + x54
    x56 = a * x24
    x57 = x14 * x56
    x58 = 5 * x57
    x59 = 15 * x19
    x60 = x59 * x6
    x61 = 18 * x34
    x62 = x0 * x61
    x63 = 30 * x11
    x64 = x41 * x6
    x65 = x20 * x62
    x66 = a ** 3
    x67 = x66 * x8
    x68 = x25 * x41
    x69 = x0 * x39
    x70 = 18 * x39
    x71 = x20 * x70
    x72 = x34 * x41
    x73 = 6 * x20
    x74 = x72 * x73
    x75 = x11 * x20
    x76 = 15 * x75
    x77 = -x42 + x51 - x52 + x53
    x78 = (
        -x0 * x71
        + x20 * x60
        + x20 * x67
        - 6 * x44
        - x6 * x63
        + x6 * x74
        + x6 * x76
        - x60
        - x61 * x64
        + x62
        - x65
        - 4 * x67
        + 18 * x68
        + 36 * x69
        + x77
    )
    x79 = x57 * x78
    x80 = x3 ** (-9 / 2)
    x81 = 105 * x80
    x82 = x6 * x81
    x83 = x0 * x19
    x84 = 150 * x83
    x85 = a / x3 ** 4
    x86 = x6 * x85
    x87 = x19 * x64
    x88 = x10 * x66
    x89 = x6 * x88
    x90 = a ** 4
    x91 = x34 * x90
    x92 = x6 * x91
    x93 = x66 * x7
    x94 = x0 * x93
    x95 = x20 * x82
    x96 = x0 * x72
    x97 = x0 * x11
    x98 = x20 * x97
    x99 = 60 * x20
    x100 = 10 * x20
    x101 = 45 * x19
    x102 = x101 * x20
    x103 = x20 * x86
    x104 = 45 * x34
    x105 = x24 * x41
    x106 = x20 * x24
    x107 = x106 * x41
    x108 = x20 * x39
    x109 = x104 * x20 - x104 - 45 * x105 + 15 * x107 + 45 * x108 - 90 * x39
    x110 = x14 / 2
    x111 = a * x30
    x112 = x110 * x111
    x113 = a * x14
    x114 = 60 * x113
    x115 = x15 / x3 ** 5
    x116 = x114 * x115
    x117 = x0 * x85
    x118 = x117 * x16
    x119 = x21 * x80
    x120 = x0 * x119
    x121 = 30 * x120
    x122 = 3 * x14 / 2
    x123 = x11 * x122 * x15
    x124 = x106 - x24
    x125 = x124 - x28 + x29
    x126 = x14 * x35
    x127 = x125 * x126
    x128 = 6 * x127
    x129 = a * x33
    x130 = 18 * x83
    x131 = 9 * x34
    x132 = a * x131
    x133 = x33 * x38
    x134 = x41 * x7
    x135 = 3 * x134
    x136 = 5 * x69
    x137 = 2 * x108
    x138 = x0 * x137
    x139 = x124 - x27
    x140 = x136 - x138 + x139 + x48 - x49 + x68
    x141 = x14 * x140
    x142 = x0 * x55
    x143 = 6 * x126
    x144 = x134 * x14
    x145 = 2 * x144
    x146 = 12 * x96
    x147 = 27 * x97
    x148 = 12 * x98
    x149 = x41 * x49
    x150 = 15 * x83
    x151 = x150 * x20
    x152 = x150 - x151
    x153 = -x43
    x154 = 15 * x39
    x155 = x20 * x40
    x156 = x131 * x20
    x157 = -x131 + x156
    x158 = x153 - x154 + x155 + x157
    x159 = x14 * (x146 + x147 - x148 - x149 + x152 + x158 + x94)
    x160 = a * x159
    x161 = x4 * x41
    x162 = x110 * x161
    x163 = 90 * x83
    x164 = 30 * x20
    x165 = 4 * x20
    x166 = x131 - x156
    x167 = x154 - x155 + x43
    x168 = x166 + x167
    x169 = (
        -a * x141 * x38 * x42
        + x0 * x128 * x36
        - x112
        * (
            -90 * x103
            + x163 * x20
            - x163
            - x164 * x87
            - x165 * x89
            + x168
            + x41 * x65
            + x82
            + 195 * x86
            + 105 * x87
            + 22 * x89
            + x92
            - 6 * x94
            - x95
            - 72 * x96
            - 162 * x97
            + 72 * x98
        )
        + x116 * x6
        + 30 * x118 * x36
        + x121 * (7 * x5 - 3)
        + x123 * x9
        + x129 * x130
        + x132 * x133
        + x133 * x135
        + x142 * x143
        + x142 * x145
        - 2 * x160 * x25
        + x162 * x78
        + x23 * x9
        + x79 / 2
    )
    x170 = x169 * y
    x171 = x169 * z
    x172 = x21 / x3 ** (11 / 2)
    x173 = x0 * x172
    x174 = 210 * x173
    x175 = x1 * x101
    x176 = 69 * x11
    x177 = 15 * x72
    x178 = x175 * x20
    x179 = 7 * x88
    x180 = x0 * x179
    x181 = 24 * x75
    x182 = x41 * x83
    x183 = 72 * x182
    x184 = 177 * x117
    x185 = x117 * x20
    x186 = 72 * x185
    x187 = x1 * x151
    x188 = x0 * x81
    x189 = x1 * x188
    x190 = -x189 * x20 + x189
    x191 = -x150 + x151
    x192 = -x146 - x147 + x148 + x149 + x168 + x191 - x94
    x193 = 105 * x119 / 2
    x194 = x193 * x36
    x195 = x1 * x24
    x196 = x195 * x66
    x197 = x110 * x55
    x198 = x1 * x4
    x199 = x159 * x41
    x200 = x132 * x141
    x201 = x135 * x141
    x202 = x122 * x55
    x203 = x202 * x35
    x204 = x1 * x106
    x205 = -x1 * x28 + x1 * x29 - x195 + x204 + x31
    x206 = x122 * x35
    x207 = x205 * x206
    x208 = -x106 + x24
    x209 = x208 + x28 - x29
    x210 = x1 * x40
    x211 = x1 * x43
    x212 = x204 * x41
    x213 = x1 * x45
    x214 = x20 * x213
    x215 = x1 * x47
    x216 = x20 * x215
    x217 = -x215 + x216
    x218 = -x210 - x211 + x212 + x214 + x217
    x219 = x209 + x218
    x220 = 3 * x57 / 2
    x221 = x220 * x38
    x222 = x134 * x202
    x223 = 3 * x33 / 2
    x224 = x34 * x66
    x225 = x223 * x224
    x226 = x10 * x41
    x227 = x226 * x33
    x228 = 21 * x227 / 2
    x229 = 45 * x19 / 2
    x230 = x129 * x229
    x231 = x132 * x14
    x232 = x125 * x231
    x233 = x232 * x38
    x234 = x16 * x85
    x235 = x234 * x37
    x236 = x113 * x125
    x237 = x130 * x236
    x238 = x0 * x1
    x239 = x113 * x115
    x240 = 120 * x239
    x241 = -x134 * x223
    x242 = 45 * x234 / 2
    x243 = x242 * x38
    x244 = x1 * x243
    x245 = x241 + x244
    x246 = 5 * x96
    x247 = 23 * x97
    x248 = 8 * x98
    x249 = x208 + x217 + x27
    x250 = -x136 + x138 + x50 - x68
    x251 = x1 * x150 + x1 * x246 + x1 * x247 - x1 * x248 - x187 - x213 + x249 + x250
    x252 = 9 * x33 / 2
    x253 = 9 * x17 / 2
    x254 = -x252 * x35 - x253 * x38
    x255 = -x220 * x251 + x254
    x256 = x57 / 2
    x257 = -9 * x118 - 15 * x120 - x162 * x55 - x23 * x36 - x256 * x55
    x258 = x * (
        x1 * x174
        + x1 * x194
        - x1 * x200
        - x1 * x201
        + x1 * x203
        + x1 * x222
        + x1 * x225
        + x1 * x228
        + x1 * x230
        + x1 * x233
        + x1 * x235
        + x1 * x237
        - x112
        * (
            -x1 * x176
            - x1 * x177
            + x1 * x180
            + x1 * x181
            + x1 * x183
            + x1 * x184
            - x1 * x186
            - x175
            + x178
            - x187 * x41
            + x190
            + x192
        )
        - x160 * x195
        + x196 * x197
        - x198 * x199
        + x207 * x36
        + x219 * x221
        + x238 * x240
        + x245
        + x255
        + x257
    )
    x259 = x0 * x240
    x260 = x20 * x47
    x261 = x260 - x47
    x262 = x107 + x153 + x261 - x40 + x46
    x263 = x220 * x262
    x264 = x152 + x246 + x247 - x248 + x261 - x45
    x265 = x24 * x66
    x266 = x110 * x265
    x267 = x188 * x20
    x268 = x151 * x41
    x269 = x * y * z
    x270 = x269 * (
        -x112
        * (
            -x101
            + x102
            - x176
            - x177
            + x180
            + x181
            + x183
            + x184
            - x186
            + x188
            - x267
            - x268
        )
        + x125 * x206 * x36
        - x159 * x161
        - x159 * x56
        + x174
        + x194
        - x200
        - x201
        + x203
        - x220 * x264
        + x222
        + x225
        + x228
        + x230
        + x233
        + x235
        + x237
        + x243
        + x259
        + x263 * x38
        + x266 * x55
    )
    x271 = x188 * x2 - x2 * x267
    x272 = x101 * x2
    x273 = x102 * x2
    x274 = -x176 * x2 - x177 * x2 + x181 * x2 - x272 + x273
    x275 = x2 * x24
    x276 = x275 * x66
    x277 = x2 * x4
    x278 = x106 * x2 - x2 * x28 + x2 * x29 - x275 + x31
    x279 = x206 * x278
    x280 = x2 * x40
    x281 = x2 * x43
    x282 = x107 * x2
    x283 = x2 * x46
    x284 = x2 * x47
    x285 = x2 * x260
    x286 = -x284 + x285
    x287 = -x280 - x281 + x282 + x283 + x286
    x288 = x209 + x287
    x289 = x2 * x243
    x290 = x241 + x289
    x291 = -x2 * x45 + x286
    x292 = (
        x150 * x2
        - x151 * x2
        + x2 * x246
        + x2 * x247
        - x2 * x248
        + x208
        + x250
        + x27
        + x291
    )
    x293 = -x220 * x292 + x254
    x294 = x * (
        -x112
        * (
            x180 * x2
            + x183 * x2
            + x184 * x2
            - x186 * x2
            + x192
            - x2 * x268
            + x271
            + x274
        )
        - x160 * x275
        + x174 * x2
        + x194 * x2
        + x197 * x276
        - x199 * x277
        - x2 * x200
        - x2 * x201
        + x2 * x203
        + x2 * x222
        + x2 * x225
        + x2 * x228
        + x2 * x230
        + x2 * x233
        + x2 * x235
        + x2 * x237
        + x2 * x259
        + x221 * x288
        + x257
        + x279 * x36
        + x290
        + x293
    )
    x295 = 5 * x19
    x296 = x1 * x295
    x297 = x1 * x11
    x298 = 35 * x80
    x299 = x238 * x298
    x300 = x1 * x182
    x301 = x1 * x117
    x302 = x1 * x185
    x303 = -x260 + x47
    x304 = x191 - x246 - x247 + x248 + x45
    x305 = x303 + x304
    x306 = x111 * x122
    x307 = x122 * x161
    x308 = x1 * x59
    x309 = x20 * x308
    x310 = -x308 + x309
    x311 = x166 + x310
    x312 = 9 * x105
    x313 = 9 * x108
    x314 = x41 * x51
    x315 = x312 - x313 - x314 + x70
    x316 = x1 * x63
    x317 = x1 * x61
    x318 = x317 * x41
    x319 = x1 * x93
    x320 = 4 * x319
    x321 = x20 * x319
    x322 = x1 * x72
    x323 = x322 * x73
    x324 = x1 * x76
    x325 = -x316 - x318 - x320 + x321 + x323 + x324
    x326 = x311 + x315 + x325
    x327 = x113 * x25
    x328 = x218 + x54
    x329 = x256 * x328
    x330 = x113 * x83
    x331 = x205 * x330
    x332 = x219 * x231
    x333 = 9 * x205 / 2
    x334 = x126 * x333
    x335 = x0 * x236
    x336 = 45 * x22 / 2
    x337 = 3 * x141 / 2
    x338 = (
        -x0 * x232
        - 54 * x118
        - 90 * x120
        - x134 * x252
        + x140 * x220
        + x161 * x337
        - x223 * x265
        - x336 * x38
    )
    x339 = 315 * x173
    x340 = x1 * x193
    x341 = x110 * x7 * x90
    x342 = x32 * x341
    x343 = x33 * x66
    x344 = 9 * x141 / 2
    x345 = x1 * x344
    x346 = 15 * x1 / 2
    x347 = x129 * x19
    x348 = 195 * x239
    x349 = (
        x1 * x339
        + x1 * x342
        - x134 * x345
        - x196 * x337
        + x215 * x343
        + x227 * x346
        + x238 * x348
        + x340 * x38
        - x345 * x35
        + x346 * x347
    )
    x350 = y * (
        x0 * x332
        + x175 * x335
        + x244
        - x251 * x307
        + x255
        - x306
        * (
            x20 * x296
            - x20 * x299
            - x296
            - 5 * x297
            + x299
            + 11 * x300
            + 51 * x301
            - 16 * x302
            + x305
        )
        + x326 * x327
        + x329 * x38
        + 9 * x331
        + x334 * x38
        + x338
        + x349
    )
    x351 = x303 + x310
    x352 = -x107 + x40 + x43 - x46
    x353 = x325 + x351 + x352
    x354 = x219 * x256
    x355 = x113 * x264
    x356 = x14 * x41
    x357 = x264 * x356
    x358 = x113 * x48
    x359 = x215 * x236
    x360 = x143 * x262
    x361 = 51 * x236 * x83
    x362 = (
        -18 * x118
        - x121
        - x123 * x38
        + x140 * x162
        + x140 * x256
        - x223 * x35
        - x23 * x38
        - x236 * x48
        - x266 * x32
    )
    x363 = z * (
        x1 * x361
        - x112 * (x190 - 15 * x297 + 33 * x300 + 153 * x301 - 48 * x302 + x304 + x351)
        - x162 * x251
        - x195 * x355
        - x198 * x357
        + x207 * x38
        + x219 * x358
        + x238 * x360
        + x245
        - x251 * x256
        + x327 * x353
        + 3 * x331
        + x349
        + x354 * x38
        + x359 * x38
        + x362
    )
    x364 = x11 * x2
    x365 = x182 * x2
    x366 = x117 * x2
    x367 = x185 * x2
    x368 = x2 * x59
    x369 = -x368
    x370 = x20 * x368
    x371 = x303 + x369 + x370
    x372 = x2 * x63
    x373 = x2 * x61
    x374 = x373 * x41
    x375 = x2 * x93
    x376 = 4 * x375
    x377 = x20 * x375
    x378 = x2 * x72
    x379 = x378 * x73
    x380 = x2 * x76
    x381 = -x372 - x374 - x376 + x377 + x379 + x380
    x382 = x352 + x371 + x381
    x383 = x256 * x288
    x384 = x278 * x330
    x385 = x236 * x284
    x386 = x143 * x2
    x387 = x193 * x2
    x388 = x2 * x344
    x389 = 15 * x2 / 2
    x390 = (
        x0 * x2 * x348
        - x134 * x388
        + x2 * x339
        + x2 * x342
        + x227 * x389
        - x276 * x337
        + x284 * x343
        + x347 * x389
        - x35 * x388
        + x38 * x387
    )
    x391 = y * (
        x0 * x262 * x386
        - x112 * (x271 + x304 - 15 * x364 + 33 * x365 + 153 * x366 - 48 * x367 + x371)
        - x162 * x292
        + x2 * x361
        - x256 * x292
        - x275 * x355
        - x277 * x357
        + x279 * x38
        + x288 * x358
        + x290
        + x327 * x382
        + x362
        + x38 * x383
        + x38 * x385
        + 3 * x384
        + x390
    )
    x392 = x2 * x298
    x393 = x0 * x392
    x394 = x2 * x295
    x395 = x20 * x394 - 5 * x364 - x394
    x396 = x166 + x315 + x369 + x370 + x381
    x397 = x287 + x54
    x398 = x256 * x397
    x399 = x231 * x288
    x400 = 9 * x126 / 2
    x401 = x278 * x400
    x402 = z * (
        x0 * x399
        + x272 * x335
        + x289
        - x292 * x307
        + x293
        - x306 * (-x20 * x393 + x305 + 11 * x365 + 51 * x366 - 16 * x367 + x393 + x395)
        + x327 * x396
        + x338
        + x38 * x398
        + x38 * x401
        + 9 * x384
        + x390
    )
    x403 = 315 * x119
    x404 = y ** 4
    x405 = 945 * x172 / 2
    x406 = 90 * x19
    x407 = x1 * x406
    x408 = x20 * x407
    x409 = 108 * x72
    x410 = 90 * x11
    x411 = x1 * x20
    x412 = 36 * x20
    x413 = x157 - x312 + x313 + x314 - x70
    x414 = x404 * x81
    x415 = 210 * x85
    x416 = x404 * x41
    x417 = 135 * x19
    x418 = 40 * x88
    x419 = x404 * x91
    x420 = x20 * x414
    x421 = x404 * x88
    x422 = x404 * x85
    x423 = 105 * x20
    x424 = (
        x100 * x421
        + x102 * x416
        + x20 * x419
        - x404 * x415
        - x404 * x418
        - x414
        - x416 * x417
        - 5 * x419
        + x420
        + x422 * x423
    )
    x425 = x205 * x231
    x426 = x113 * x42
    x427 = 135 * x234
    x428 = x113 * x195
    x429 = 2 * x428
    x430 = x19 * x236
    x431 = 30 * x430
    x432 = x113 * x205
    x433 = 525 * x239 / 2
    x434 = x253 + x336
    x435 = x404 * x59
    x436 = x404 * x93
    x437 = x195 * x41
    x438 = 36 * x39
    x439 = (
        x1 * x438
        - x1 * x71
        - x20 * x317
        + x20 * x435
        + x20 * x436
        - 6 * x212
        + x317
        - x404 * x63
        + x404 * x74
        + x404 * x76
        - x416 * x61
        - x435
        - 4 * x436
        + 18 * x437
        + x77
    )
    x440 = x1 * x328
    x441 = x143 * x440 + x256 * x439
    x442 = x * (
        x1 * x332
        - x1 * x403
        - x1 * x427
        + x112
        * (
            x1 * x409
            + 180 * x297
            + 24 * x319
            - 6 * x321
            - x322 * x412
            + x407
            - x408
            - x410 * x411
            + x413
            + x424
        )
        + x175 * x432
        - x219 * x426
        - x236 * x317
        + x326 * x429
        + x404 * x405
        + x404 * x431
        + x404 * x433
        - x425
        + x434
        + x441
    )
    x443 = x1 * x405
    x444 = x1 * x81
    x445 = x1 * x415
    x446 = x41 * x417
    x447 = x1 * x446
    x448 = x1 * x418
    x449 = x1 * x91
    x450 = 5 * x449
    x451 = x20 * x444
    x452 = x20 * x449
    x453 = x100 * x88
    x454 = x1 * x453
    x455 = x178 * x41
    x456 = x411 * x85
    x457 = 105 * x456
    x458 = x41 * x61
    x459 = x20 * x93
    x460 = x101 - x102 - x20 * x458 + x410 - 3 * x459 + 54 * x72 - 45 * x75 + 12 * x93
    x461 = x206 * x328
    x462 = x113 * x229
    x463 = x205 * x462
    x464 = x1 * x400
    x465 = 105 * x430 / 2
    x466 = -315 * x119 / 2 - 27 * x127 / 2 - 135 * x234 / 2 - x263
    x467 = x269 * (
        x1 * x433
        + x1 * x465
        + x112
        * (-x444 - x445 - x447 - x448 - x450 + x451 + x452 + x454 + x455 + x457 + x460)
        + x220 * x353
        + x262 * x464
        + x326 * x57
        + x332
        + x443
        + x461
        + x463
        + x466
    )
    x468 = x2 * x242
    x469 = x1 * x2
    x470 = x19 * x469
    x471 = x113 * x278
    x472 = x19 * x471
    x473 = x19 * x432
    x474 = x389 * x473
    x475 = x2 * x400
    x476 = x113 * x275
    x477 = x2 * x308
    x478 = x2 * x309
    x479 = x215 - x216
    x480 = x2 * x444
    x481 = x2 * x451
    x482 = x308 - x309
    x483 = x * (
        -x1 * x242
        + x112
        * (
            -x2 * x445
            - x2 * x447
            - x2 * x448
            - x2 * x450
            + x2 * x452
            + x2 * x454
            + x2 * x455
            + x2 * x457
            + x262
            + x316
            + x318
            + x320
            - x321
            - x323
            - x324
            + x368
            - x370
            + x372
            + x374
            + x376
            - x377
            - x379
            - x380
            - x480
            + x481
            + x482
        )
        + x114 * x125 * x470
        + x123
        + x2 * x443
        - x207
        + x219 * x475
        + x23
        + x256
        * (
            x125
            - x2 * x316
            - x2 * x318
            - x2 * x320
            + x2 * x321
            + x2 * x323
            + x2 * x324
            + x210
            + x211
            - x212
            - x214
            + x280
            + x281
            - x282
            - x283
            + x284
            - x285
            - x477
            + x478
            + x479
        )
        - x279
        + x288 * x464
        - x340
        + x346 * x472
        + x353 * x476
        - x354
        - x359
        + x360 * x469
        + x382 * x428
        - x383
        - x385
        - x387
        + x433 * x469
        - x468
        + x474
    )
    x484 = x2 * x81
    x485 = x2 * x91
    x486 = x423 * x85
    x487 = x269 * (
        x112
        * (
            -x2 * x415
            - x2 * x418
            - x2 * x446
            + x2 * x453
            + x2 * x486
            + x20 * x484
            + x20 * x485
            + x273 * x41
            + x460
            - x484
            - 5 * x485
        )
        + x2 * x405
        + x2 * x433
        + x2 * x465
        + x206 * x397
        + x220 * x382
        + x262 * x475
        + x278 * x462
        + x396 * x57
        + x399
        + x466
    )
    x488 = z ** 4
    x489 = x488 * x59
    x490 = (
        x2 * x438
        - x2 * x71
        - x20 * x373
        + x20 * x489
        + 18 * x275 * x41
        - 6 * x282
        + x373
        - x458 * x488
        + x459 * x488
        - x488 * x63
        + x488 * x74
        + x488 * x76
        - 4 * x488 * x93
        - x489
        + x77
    )
    x491 = x2 * x406
    x492 = x488 * x81
    x493 = x488 * x91
    x494 = (
        x102 * x41 * x488
        + x20 * x492
        + x20 * x493
        - x415 * x488
        - x418 * x488
        - x446 * x488
        + x453 * x488
        + x486 * x488
        - x492
        - 5 * x493
    )
    x495 = (
        x112
        * (
            -x2 * x20 * x410
            + x2 * x409
            - x20 * x491
            + 180 * x364
            + 24 * x375
            - 6 * x377
            - x378 * x412
            + x413
            + x491
            + x494
        )
        + x2 * x399
        - x2 * x403
        - x2 * x427
        - x231 * x278
        - x236 * x373
        + x256 * x490
        + x272 * x471
        - x288 * x426
        + x386 * x397
        + 2 * x396 * x476
        + x405 * x488
        + x431 * x488
        + x433 * x488
        + x434
    )
    x496 = x * x495
    x497 = x404 * x7
    x498 = -30 * x198 + 35 * x497 + 3
    x499 = 5 * x198 - 3
    x500 = 15 * x126
    x501 = 3 * x198 - 1
    x502 = 5 * x57 / 2
    x503 = 150 * x19
    x504 = x1 * x503
    x505 = x1 * x75
    x506 = x1 * x119
    x507 = x19 * x416
    x508 = 5 * x1 * x39
    x509 = x1 * x137
    x510 = x139 + x437 + x479 + x508 - x509
    x511 = 12 * x322
    x512 = 27 * x297
    x513 = 12 * x505
    x514 = x216 * x41
    x515 = x158 + x319 + x482 + x511 + x512 - x513 - x514
    x516 = x135 * x14
    x517 = x1 * x234
    x518 = z * (
        x1 * x128 * x499
        + 18 * x1 * x473
        - x112
        * (
            -x164 * x507
            - x165 * x421
            + x168
            + x20 * x318
            - 90 * x20 * x422
            - 162 * x297
            - 6 * x319
            - 72 * x322
            - x407
            + x408
            + x414
            + x419
            - x420
            + 22 * x421
            + 195 * x422
            + 72 * x505
            + 105 * x507
        )
        + x116 * x404
        + x123 * x498
        + x145 * x440
        + x162 * x439
        + x205 * x501 * x516
        + x23 * x498
        + x425 * x501
        - x426 * x501 * x510
        - x429 * x515
        + x441
        + 30 * x499 * x517
        + 30 * x506 * (7 * x198 - 3)
    )
    x519 = x172 * x469
    x520 = x122 * x205
    x521 = x41 * x470
    x522 = x469 * x85
    x523 = x2 * x456
    x524 = x2 * x510
    x525 = x2 * x205
    x526 = x14 * x226
    x527 = 5 * x322
    x528 = 23 * x297
    x529 = 8 * x505
    x530 = (
        x2 * x527
        + x2 * x528
        - x2 * x529
        + x249
        + x291
        - x437
        + x477
        - x478
        - x508
        + x509
    )
    x531 = -x220 * x530 - x253 * x501 - x334 + x468 * x501
    x532 = y * (
        x110 * x276 * x328
        - x112
        * (
            x167
            + x179 * x469
            + x274
            + x311
            - x319
            - x41 * x478
            + x480
            - x481
            - x511
            - x512
            + x513
            + x514
            + 72 * x521
            + 177 * x522
            - 72 * x523
        )
        + x122 * x134 * x2 * x328
        - x134 * x520
        - x162 * x328
        + x2 * x224 * x520
        + x2 * x232 * x501
        + 15 * x2 * x234 * x499
        + x2 * x461
        + x2 * x463
        + x220 * x288 * x501
        - x23 * x499
        - x231 * x524
        + 18 * x236 * x470
        + x240 * x469
        - x277 * x356 * x515
        + x279 * x499
        - x329
        + x387 * x499
        - x476 * x515
        - 15 * x506
        - x516 * x524
        - 9 * x517
        + 210 * x519
        + 21 * x525 * x526 / 2
        + x531
    )
    x533 = x1 * x392
    x534 = z * (
        -x1 * x232
        + x1 * x399
        + 9 * x1 * x472
        - x122 * x276 * x510
        + x14 * x205 * x284 * x66
        - x144 * x333
        - 9 * x144 * x524 / 2
        + x175 * x2 * x236
        + x205 * x389 * x526
        + x220 * x510
        - x265 * x520
        - x306
        * (
            -x20 * x533
            + x351
            + x395
            + x45
            + 11 * x521
            + 51 * x522
            - 16 * x523
            - x527
            - x528
            + x529
            + x533
        )
        + x307 * x510
        - x307 * x530
        - x336 * x501
        + x341 * x525
        + x348 * x469
        + x387 * x501
        + x396 * x428
        + x398 * x501
        + x401 * x501
        + x474
        - x475 * x510
        - 90 * x506
        - 54 * x517
        + 315 * x519
        + x531
    )
    x535 = x495 * y
    x536 = x488 * x7
    x537 = x2 * x503
    result[0, 0, 0, 0, 0] = x * (
        x112
        * (
            x100 * x89
            - x100 * x94
            + x102 * x64
            + 105 * x103
            + x109
            - x20 * x84
            + x20 * x92
            - x82
            + x84
            - 210 * x86
            - 135 * x87
            - 40 * x89
            - 5 * x92
            + 40 * x94
            + x95
            - x96 * x99
            + 180 * x96
            + 300 * x97
            - 150 * x98
        )
        + x18 * x9
        + x23 * (-70 * x5 + 63 * x8 + 15)
        + x33 * x35 * x37
        + x38 * x55 * x58
        + 5 * x79 / 2
    )
    result[0, 0, 0, 0, 1] = x170
    result[0, 0, 0, 0, 2] = x171
    result[0, 0, 0, 1, 0] = x170
    result[0, 0, 0, 1, 1] = x258
    result[0, 0, 0, 1, 2] = x270
    result[0, 0, 0, 2, 0] = x171
    result[0, 0, 0, 2, 1] = x270
    result[0, 0, 0, 2, 2] = x294
    result[0, 0, 1, 0, 0] = x170
    result[0, 0, 1, 0, 1] = x258
    result[0, 0, 1, 0, 2] = x270
    result[0, 0, 1, 1, 0] = x258
    result[0, 0, 1, 1, 1] = x350
    result[0, 0, 1, 1, 2] = x363
    result[0, 0, 1, 2, 0] = x270
    result[0, 0, 1, 2, 1] = x363
    result[0, 0, 1, 2, 2] = x391
    result[0, 0, 2, 0, 0] = x171
    result[0, 0, 2, 0, 1] = x270
    result[0, 0, 2, 0, 2] = x294
    result[0, 0, 2, 1, 0] = x270
    result[0, 0, 2, 1, 1] = x363
    result[0, 0, 2, 1, 2] = x391
    result[0, 0, 2, 2, 0] = x294
    result[0, 0, 2, 2, 1] = x391
    result[0, 0, 2, 2, 2] = x402
    result[0, 1, 0, 0, 0] = x170
    result[0, 1, 0, 0, 1] = x258
    result[0, 1, 0, 0, 2] = x270
    result[0, 1, 0, 1, 0] = x258
    result[0, 1, 0, 1, 1] = x350
    result[0, 1, 0, 1, 2] = x363
    result[0, 1, 0, 2, 0] = x270
    result[0, 1, 0, 2, 1] = x363
    result[0, 1, 0, 2, 2] = x391
    result[0, 1, 1, 0, 0] = x258
    result[0, 1, 1, 0, 1] = x350
    result[0, 1, 1, 0, 2] = x363
    result[0, 1, 1, 1, 0] = x350
    result[0, 1, 1, 1, 1] = x442
    result[0, 1, 1, 1, 2] = x467
    result[0, 1, 1, 2, 0] = x363
    result[0, 1, 1, 2, 1] = x467
    result[0, 1, 1, 2, 2] = x483
    result[0, 1, 2, 0, 0] = x270
    result[0, 1, 2, 0, 1] = x363
    result[0, 1, 2, 0, 2] = x391
    result[0, 1, 2, 1, 0] = x363
    result[0, 1, 2, 1, 1] = x467
    result[0, 1, 2, 1, 2] = x483
    result[0, 1, 2, 2, 0] = x391
    result[0, 1, 2, 2, 1] = x483
    result[0, 1, 2, 2, 2] = x487
    result[0, 2, 0, 0, 0] = x171
    result[0, 2, 0, 0, 1] = x270
    result[0, 2, 0, 0, 2] = x294
    result[0, 2, 0, 1, 0] = x270
    result[0, 2, 0, 1, 1] = x363
    result[0, 2, 0, 1, 2] = x391
    result[0, 2, 0, 2, 0] = x294
    result[0, 2, 0, 2, 1] = x391
    result[0, 2, 0, 2, 2] = x402
    result[0, 2, 1, 0, 0] = x270
    result[0, 2, 1, 0, 1] = x363
    result[0, 2, 1, 0, 2] = x391
    result[0, 2, 1, 1, 0] = x363
    result[0, 2, 1, 1, 1] = x467
    result[0, 2, 1, 1, 2] = x483
    result[0, 2, 1, 2, 0] = x391
    result[0, 2, 1, 2, 1] = x483
    result[0, 2, 1, 2, 2] = x487
    result[0, 2, 2, 0, 0] = x294
    result[0, 2, 2, 0, 1] = x391
    result[0, 2, 2, 0, 2] = x402
    result[0, 2, 2, 1, 0] = x391
    result[0, 2, 2, 1, 1] = x483
    result[0, 2, 2, 1, 2] = x487
    result[0, 2, 2, 2, 0] = x402
    result[0, 2, 2, 2, 1] = x487
    result[0, 2, 2, 2, 2] = x496
    result[1, 0, 0, 0, 0] = x170
    result[1, 0, 0, 0, 1] = x258
    result[1, 0, 0, 0, 2] = x270
    result[1, 0, 0, 1, 0] = x258
    result[1, 0, 0, 1, 1] = x350
    result[1, 0, 0, 1, 2] = x363
    result[1, 0, 0, 2, 0] = x270
    result[1, 0, 0, 2, 1] = x363
    result[1, 0, 0, 2, 2] = x391
    result[1, 0, 1, 0, 0] = x258
    result[1, 0, 1, 0, 1] = x350
    result[1, 0, 1, 0, 2] = x363
    result[1, 0, 1, 1, 0] = x350
    result[1, 0, 1, 1, 1] = x442
    result[1, 0, 1, 1, 2] = x467
    result[1, 0, 1, 2, 0] = x363
    result[1, 0, 1, 2, 1] = x467
    result[1, 0, 1, 2, 2] = x483
    result[1, 0, 2, 0, 0] = x270
    result[1, 0, 2, 0, 1] = x363
    result[1, 0, 2, 0, 2] = x391
    result[1, 0, 2, 1, 0] = x363
    result[1, 0, 2, 1, 1] = x467
    result[1, 0, 2, 1, 2] = x483
    result[1, 0, 2, 2, 0] = x391
    result[1, 0, 2, 2, 1] = x483
    result[1, 0, 2, 2, 2] = x487
    result[1, 1, 0, 0, 0] = x258
    result[1, 1, 0, 0, 1] = x350
    result[1, 1, 0, 0, 2] = x363
    result[1, 1, 0, 1, 0] = x350
    result[1, 1, 0, 1, 1] = x442
    result[1, 1, 0, 1, 2] = x467
    result[1, 1, 0, 2, 0] = x363
    result[1, 1, 0, 2, 1] = x467
    result[1, 1, 0, 2, 2] = x483
    result[1, 1, 1, 0, 0] = x350
    result[1, 1, 1, 0, 1] = x442
    result[1, 1, 1, 0, 2] = x467
    result[1, 1, 1, 1, 0] = x442
    result[1, 1, 1, 1, 1] = y * (
        x112
        * (
            x109
            - x20 * x504
            + 300 * x297
            + 40 * x319
            - 10 * x321
            - x322 * x99
            + 180 * x322
            + x424
            + x504
            - 150 * x505
        )
        + x18 * x498
        + x205 * x499 * x500
        + x23 * (-70 * x198 + 63 * x497 + 15)
        + x328 * x501 * x58
        + x439 * x502
    )
    result[1, 1, 1, 1, 2] = x518
    result[1, 1, 1, 2, 0] = x467
    result[1, 1, 1, 2, 1] = x518
    result[1, 1, 1, 2, 2] = x532
    result[1, 1, 2, 0, 0] = x363
    result[1, 1, 2, 0, 1] = x467
    result[1, 1, 2, 0, 2] = x483
    result[1, 1, 2, 1, 0] = x467
    result[1, 1, 2, 1, 1] = x518
    result[1, 1, 2, 1, 2] = x532
    result[1, 1, 2, 2, 0] = x483
    result[1, 1, 2, 2, 1] = x532
    result[1, 1, 2, 2, 2] = x534
    result[1, 2, 0, 0, 0] = x270
    result[1, 2, 0, 0, 1] = x363
    result[1, 2, 0, 0, 2] = x391
    result[1, 2, 0, 1, 0] = x363
    result[1, 2, 0, 1, 1] = x467
    result[1, 2, 0, 1, 2] = x483
    result[1, 2, 0, 2, 0] = x391
    result[1, 2, 0, 2, 1] = x483
    result[1, 2, 0, 2, 2] = x487
    result[1, 2, 1, 0, 0] = x363
    result[1, 2, 1, 0, 1] = x467
    result[1, 2, 1, 0, 2] = x483
    result[1, 2, 1, 1, 0] = x467
    result[1, 2, 1, 1, 1] = x518
    result[1, 2, 1, 1, 2] = x532
    result[1, 2, 1, 2, 0] = x483
    result[1, 2, 1, 2, 1] = x532
    result[1, 2, 1, 2, 2] = x534
    result[1, 2, 2, 0, 0] = x391
    result[1, 2, 2, 0, 1] = x483
    result[1, 2, 2, 0, 2] = x487
    result[1, 2, 2, 1, 0] = x483
    result[1, 2, 2, 1, 1] = x532
    result[1, 2, 2, 1, 2] = x534
    result[1, 2, 2, 2, 0] = x487
    result[1, 2, 2, 2, 1] = x534
    result[1, 2, 2, 2, 2] = x535
    result[2, 0, 0, 0, 0] = x171
    result[2, 0, 0, 0, 1] = x270
    result[2, 0, 0, 0, 2] = x294
    result[2, 0, 0, 1, 0] = x270
    result[2, 0, 0, 1, 1] = x363
    result[2, 0, 0, 1, 2] = x391
    result[2, 0, 0, 2, 0] = x294
    result[2, 0, 0, 2, 1] = x391
    result[2, 0, 0, 2, 2] = x402
    result[2, 0, 1, 0, 0] = x270
    result[2, 0, 1, 0, 1] = x363
    result[2, 0, 1, 0, 2] = x391
    result[2, 0, 1, 1, 0] = x363
    result[2, 0, 1, 1, 1] = x467
    result[2, 0, 1, 1, 2] = x483
    result[2, 0, 1, 2, 0] = x391
    result[2, 0, 1, 2, 1] = x483
    result[2, 0, 1, 2, 2] = x487
    result[2, 0, 2, 0, 0] = x294
    result[2, 0, 2, 0, 1] = x391
    result[2, 0, 2, 0, 2] = x402
    result[2, 0, 2, 1, 0] = x391
    result[2, 0, 2, 1, 1] = x483
    result[2, 0, 2, 1, 2] = x487
    result[2, 0, 2, 2, 0] = x402
    result[2, 0, 2, 2, 1] = x487
    result[2, 0, 2, 2, 2] = x496
    result[2, 1, 0, 0, 0] = x270
    result[2, 1, 0, 0, 1] = x363
    result[2, 1, 0, 0, 2] = x391
    result[2, 1, 0, 1, 0] = x363
    result[2, 1, 0, 1, 1] = x467
    result[2, 1, 0, 1, 2] = x483
    result[2, 1, 0, 2, 0] = x391
    result[2, 1, 0, 2, 1] = x483
    result[2, 1, 0, 2, 2] = x487
    result[2, 1, 1, 0, 0] = x363
    result[2, 1, 1, 0, 1] = x467
    result[2, 1, 1, 0, 2] = x483
    result[2, 1, 1, 1, 0] = x467
    result[2, 1, 1, 1, 1] = x518
    result[2, 1, 1, 1, 2] = x532
    result[2, 1, 1, 2, 0] = x483
    result[2, 1, 1, 2, 1] = x532
    result[2, 1, 1, 2, 2] = x534
    result[2, 1, 2, 0, 0] = x391
    result[2, 1, 2, 0, 1] = x483
    result[2, 1, 2, 0, 2] = x487
    result[2, 1, 2, 1, 0] = x483
    result[2, 1, 2, 1, 1] = x532
    result[2, 1, 2, 1, 2] = x534
    result[2, 1, 2, 2, 0] = x487
    result[2, 1, 2, 2, 1] = x534
    result[2, 1, 2, 2, 2] = x535
    result[2, 2, 0, 0, 0] = x294
    result[2, 2, 0, 0, 1] = x391
    result[2, 2, 0, 0, 2] = x402
    result[2, 2, 0, 1, 0] = x391
    result[2, 2, 0, 1, 1] = x483
    result[2, 2, 0, 1, 2] = x487
    result[2, 2, 0, 2, 0] = x402
    result[2, 2, 0, 2, 1] = x487
    result[2, 2, 0, 2, 2] = x496
    result[2, 2, 1, 0, 0] = x391
    result[2, 2, 1, 0, 1] = x483
    result[2, 2, 1, 0, 2] = x487
    result[2, 2, 1, 1, 0] = x483
    result[2, 2, 1, 1, 1] = x532
    result[2, 2, 1, 1, 2] = x534
    result[2, 2, 1, 2, 0] = x487
    result[2, 2, 1, 2, 1] = x534
    result[2, 2, 1, 2, 2] = x535
    result[2, 2, 2, 0, 0] = x402
    result[2, 2, 2, 0, 1] = x487
    result[2, 2, 2, 0, 2] = x496
    result[2, 2, 2, 1, 0] = x487
    result[2, 2, 2, 1, 1] = x534
    result[2, 2, 2, 1, 2] = x535
    result[2, 2, 2, 2, 0] = x496
    result[2, 2, 2, 2, 1] = x535
    result[2, 2, 2, 2, 2] = z * (
        x112
        * (
            x109
            - 150 * x2 * x75
            - x20 * x537
            + 300 * x364
            + 40 * x375
            - 10 * x377
            - x378 * x99
            + 180 * x378
            + x494
            + x537
        )
        + x18 * (-30 * x277 + 35 * x536 + 3)
        + x23 * (-70 * x277 + 63 * x536 + 15)
        + x278 * x500 * (5 * x277 - 3)
        + x397 * x58 * (3 * x277 - 1)
        + x490 * x502
    )
    return result


T_damp_thole = [
    T_damp_thole_0,
    T_damp_thole_1,
    T_damp_thole_2,
    T_damp_thole_3,
    T_damp_thole_4,
    T_damp_thole_5,
]
