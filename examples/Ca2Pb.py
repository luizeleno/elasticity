import elasticity as el
import plots3d


# S11, S12, S13, S33, S44, S66
ca = el.ElasticityTheory('orthorhombic', 26.47, 29.77, 21.71, 74.54, 29.76, 70.52, 23.63, 21.40, 30.71)

N = 200
plots3d.mayaviplot(ca, N, m=60, d=450)
