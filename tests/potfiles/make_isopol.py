import numpy as np
from pepytools import Potential

pot = Potential.from_file("pna_6w.pot")

newpols = []
for pol in pot.getPolarizabilities():
    trace = pol[0] + pol[3] + pol[5]
    isopol = np.zeros_like(pol)
    isopol[0] = trace/3.0
    isopol[3] = trace/3.0
    isopol[5] = trace/3.0
    newpols.append(isopol)

pot.setPolarizabilities(newpols)

pot.save("pna_6w_isopol.pot")
