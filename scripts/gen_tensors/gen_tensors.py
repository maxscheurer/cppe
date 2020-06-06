#!/usr/bin/env python
import sympy
import sys
from sympy.codegen.rewriting import create_expand_pow_optimization


maxorder = int(sys.argv[1])
maxorder_damped = int(sys.argv[2])
maxpow = int(sys.argv[3])


expand_pow = create_expand_pow_optimization(maxpow)
print("// maximum order, maximum order for Thole-damped tensors, maximum expansion order of pow")
print(f"// maxorder = {maxorder}, maxorder_damped = {maxorder_damped}, maxpow = {maxpow}")


def T_cart(sequence):
    """
    T = 1/R
    Generates sympy cartesian T-tensor of rank len(sequence)
    """
    x, y, z = sympy.symbols("x y z")
    R = sympy.sqrt(x**2 + y**2 + z**2)
    T = 1 / R
    if sequence:
        return sympy.diff(T, *[i for i in sequence])
    else:
        return T


def T_damp_erf(sequence):
    """
    S = v*r = f(r)
    f(r) = r/(erf(sqrt(a)*r))
    T = 1/S
    """
    x, y, z, a = sympy.symbols("x y z a")
    R = sympy.sqrt(x**2 + y**2 + z**2)
    S = R / sympy.erf(sympy.sqrt(a) * R)
    T = 1 / S
    if sequence:
        return sympy.diff(T, *[i for i in sequence])
    else:
        return T


def T_damp_thole(sequence):
    """
    Swart, M., and P. Th van Duijnen. "DRF90: a polarizable force field." Molecular Simulation 32.6 (2006): 471-484.
    u = r/(a_i*a_j)**(1/6) 
    a_i: polarizability on atom i
    a_j: polarizability on atom j
    v = damp*u
    damp: screening factor
    v = damp*u = r/(a_i*a_j)**(1/6)*damp = r*a
    a = 1/(a_i*a_j)**(1/6)*damp
    """
    x, y, z, a = sympy.symbols("x y z a")
    R = sympy.sqrt(x**2 + y**2 + z**2)
    v = R * a
    fV = 1 - (sympy.Rational(1, 2) * v + 1) * sympy.exp(-v)
    T = fV / R
    if sequence:
        return sympy.diff(T, *[i for i in sequence])
    else:
        return T

def T_damp_thole_amoeba(sequence):
    x, y, z, a = sympy.symbols("x y z a")
    """
    From:
    The polarizable point dipoles method withelectrostatic damping: Implementation on amodel system
     J. Chem. Phys. 133, 234101 (2010);

    AMOEBA: exp(-au**3)  <- exp(-(ra)**3) =>
    a = 1/(alpha_i * alpha_j)**1/6 * (damp)**1/3.
    """
    R = sympy.sqrt(x**2 + y**2 + z**2)
    v = a * R
    s0 = 1 - sympy.exp(-(v)**3) + v * sympy.functions.special.gamma_functions.uppergamma(
        sympy.Rational(2, 3), (v)**3)
    s1 = 1 - sympy.exp(-(v)**3)
    T = s0 / R
    Tx = -s1 * x / R**3
    Ty = -s1 * y / R**3
    Tz = -s1 * z / R**3
    if sequence:
        if len(sequence) > 1:
            if sequence[0] == "x":
                return sympy.diff(Tx, *[i for i in sequence[1:]])
            elif sequence[0] == "y":
                return sympy.diff(Ty, *[i for i in sequence[1:]])
            elif sequence[0] == "z":
                return sympy.diff(Tz, *[i for i in sequence[1:]])
        else:
            if sequence[0] == "x":
                return Tx
            elif sequence[0] == "y":
                return Ty
            elif sequence[0] == "z":
                return Tz
    else:
        return T


print('#include \"tensors_autogen.hh\"')
print('namespace libcppe {')
print('namespace tensors {')


def print_tensors(k, T, postfix="", extra_args=""):
    seqs = [""]
    for order in range(k + 1):
        print(f"Eigen::VectorXd T{order}{postfix}(const Eigen::Vector3d& rij{extra_args})")
        print("{")
        print("    double x = rij(0); double y = rij(1); double z = rij(2);")
        print(f"    Eigen::VectorXd result({len(seqs)});")
        components = [T(seq) for seq in seqs]
        cse = sympy.cse(components, optimizations="basic")
        cse_intermediates = cse[0]
        cse_outputs = cse[1]

        for intermediate in cse_intermediates:
            icode = sympy.ccode(expand_pow(intermediate[1].evalf()))
            print(f"    double {intermediate[0].evalf()} = {icode};")
        for i, seq in enumerate(seqs):
            code = sympy.ccode(expand_pow(cse_outputs[i].evalf()))
            print(
                f"    result[{i}] = {code}; // {seqs[i]}"
            )
        print("    return result;")
        print("}")

        newseq = []
        for s in seqs:
            newseq.append(s + "x")
            newseq.append(s + "y")
            newseq.append(s + "z")
        newseq = [''.join(sorted(s)) for s in newseq]
        newseq = sorted(list(set(newseq)))
        seqs = newseq


print_tensors(maxorder, T_cart)
print_tensors(maxorder_damped, T_damp_thole, postfix="_damp_thole", extra_args=", double a")

print("} // namespace tensors")
print("} // namespace libcppe")
