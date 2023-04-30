import elasticity as el
import plots3d

# S11, S12, S13, S33, S44, S66
h2o = el.Elasticity('hexagonal', 17.1, 8.5, 7.13, 18.21, 3.62)

labels = ['Ex (GPa)', 'Ey (GPa)', 'Ez (GPa)']

N = 200
plots3d.YoungModulus3DPlot(h2o, N, m=15, d=101, offscreen=False, labels=labels)
