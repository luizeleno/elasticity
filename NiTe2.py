import numpy as np
import elasticity as el
import matplotlib.pyplot as plt
# import plots3d

# pbe = el.Elasticity('trigonal_1', 110.8, 38.2, 22.9, -5.0, 45.5, 10.2)
pbesoc = el.Elasticity('trigonal_1', 113.7, 36.6, 27.2, -6.5, 45.7, 11.2)
# pzsoc = el.Elasticity('trigonal_1', 145.5, 54.0, 43.3, -14.3, 76.6, 26.8)
vdw = el.Elasticity('trigonal_1', 127.4, 47.1, 26.8, -9.4, 75.7, 20.2)

# pbe.VRH()
# print('B', pbe.BulkModulus)
# print('E', pbe.YoungModulus)
# print('G', pbe.ShearModulus)
# print('v', pbe.PoissonRatio)


N = 200
m = 100
d = 700
#
labels = ['Ex (GPa)', 'Ey (GPa)', 'Ez (GPa)']
#
# plots3d.YoungModulus3DPlot(pbesoc, N, m, d, filename='E-NiTe2-pbesoc.png', labels=labels)
# plots3d.YoungModulus3DPlot(pbe, N, m, d, filename='E-NiTe2-pbe.png', labels=labels, offscreen=False)
# plots3d.YoungModulus3DPlot(pzsoc, N, m, d, filename='E-NiTe2-pzsoc.png', labels=labels)
# plots3d.YoungModulus3DPlot(vdw, N, m, d, filename='E-NiTe2-vdw.png', labels=labels)

# labels = ['Bx (GPa)', 'By (GPa)', 'Bz (GPa)']

# m = 300
# d = 2050

# plots3d.LinearCompressibility3DPlot(pbe, N, m, d, filename='B-NiTe2-pbe.png', labels=labels)
# plots3d.LinearCompressibility3DPlot(pbesoc, N, m, d, filename='B-NiTe2-pbesoc.png', labels=labels)
# plots3d.LinearCompressibility3DPlot(pzsoc, N, m, d, filename='B-NiTe2-pzsoc.png', labels=labels)
# plots3d.LinearCompressibility3DPlot(vdw, N, m, d, filename='B-NiTe2-vdw.png', labels=labels)

N = 201
phi = np.pi / 4  # np.linspace(0, 2*np.pi, N)
theta = np.linspace(0, 2 * np.pi, N)

ldir = np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])

Epbesoc = pbesoc.YoungModulus(ldir)
Evdw = vdw.YoungModulus(ldir)

fig = plt.figure(figsize=(7,7))
ax = fig.add_subplot(111)
plt.polar(theta-np.pi/2, Epbesoc, label='PBE+SOC')
plt.polar(theta-np.pi/2, Evdw, label='optB86b-vdW')
plt.legend(loc=0)
plt.text(0, 1, r'$E (GPa)$', transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.5) )
plt.thetagrids(range(0,360,45), ('[110]', '', '[001]') )
plt.savefig('NiTe2-Young-110.png', bbox_inches='tight')
plt.show()
