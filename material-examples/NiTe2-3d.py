import elasticity as el
import plots3d

import numpy as np
import elasticity as el
import matplotlib.pyplot as plt
# import plots3d

# pbe = el.Elasticity('trigonal_1', 110.8, 38.2, 22.9, -5.0, 45.5, 10.2)
pbesoc = el.Elasticity('trigonal_1', 113.7, 36.6, 27.2, -6.5, 45.7, 11.2)
# pzsoc = el.Elasticity('trigonal_1', 145.5, 54.0, 43.3, -14.3, 76.6, 26.8)
vdw = el.Elasticity('trigonal_1', 127.4, 47.1, 26.8, -9.4, 75.7, 20.2)

labelsE = ['Ex', 'Ey', 'Ez']
labelsB = ['Bx', 'By', 'Bz']

N = 200

plots3d.YoungModulus3DPlot(pbesoc, N, m=100, d=700, filename='NiTe2-E-pbesoc.png', labels=labelsE)
plots3d.YoungModulus3DPlot(vdw, N, m=100, d=700, filename='NiTe2-E-vdw.png', labels=labelsE)

plots3d.LinearCompressibility3DPlot(pbesoc, N, m=240, d=1800, filename='NiTe2-B-pbesoc.png', labels=labelsB)
plots3d.LinearCompressibility3DPlot(vdw, N, m=240, d=1800, filename='NiTe2-B-vdw.png', labels=labelsB)

