# import numpy as np
import elasticity_helper as eh
import elasticity as el
import plots3d

# S11, S12, S13, S33, S44, S66
cs = eh.StiffnessConstants('tetragonal_1', 29.10, -13.10, -1.81, 11.00, 213.00, 84.00)
C = cs.C
C11 = C[0, 0]
C12 = C[0, 1]
C13 = C[0, 2]
C33 = C[2, 2]
C44 = C[3, 3]
C66 = C[5, 5]

labels = ['Ex (GPa)', 'Ey (GPa)', 'Ez (GPa)']

cs = el.Elasticity('tetragonal_1', C11, C12, C13, C33, C44, C66)

N = 200
plots3d.YoungModulus3DPlot(cs, N, offscreen=False, labels=labels)
