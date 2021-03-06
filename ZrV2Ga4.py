'''
    Test for a known compound
'''

import numpy as np
# import matplotlib.pyplot as plt

# This import registers the 3D projection, but is otherwise unused.
# from mpl_toolkits.mplot3d import Axes3D

import plots3d

import elasticity as el

zr = el.Elasticity('tetragonal_1', 288.7, 133.6, 76.7, 283.0, 102.5, 157.9)

N = 200
phi = np.linspace(0, 2*np.pi, N)
theta = np.pi / 2.
ldir = np.array([np.cos(phi)*np.sin(theta), np.sin(phi)*np.sin(theta), np.cos(theta)])
E = zr.YoungModulus(ldir)

# fig = plt.figure(figsize=(5, 5))
# ax = fig.add_subplot(111)
#
# plt.polar(phi, E)
# plt.thetagrids(range(0, 360, 45), ('[010]', '', '[001]'))
# plt.rgrids(range(50, 351, 50))
# plt.text(0, 1, 'E[0kl]', transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.5))

# plt.show()

plots3d.mayaviplot(zr, N)

# fig = plt.figure()
# ax = fig.gca(projection='3d')
#
# # Plot the surface.
# surf = ax.plot_surface(x, y, z, cmap='jet', linewidth=2, antialiased=False)
#
# # surf = plt.imshow(E)
#
# # Add a color bar which maps values to colors.
# fig.colorbar(surf)
#
# plt.show()
