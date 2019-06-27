from mayavi import mlab
import numpy as np


def mayaviplot(case, N=151, m=.1, d=.7):

    N = int(N) * 1j
    theta, phi = np.mgrid[0:np.pi:N, 0:2*np.pi:N]
    ldir = np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])
    E = case.YoungModulus(ldir)

    x = E * ldir[0]
    y = E * ldir[1]
    z = E * ldir[2]

    mlab.options.offscreen = True
    mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(1050, 1000))
    mlab.clf()

    mlab.mesh(x, y, z, scalars=E, colormap='jet')  # , vmin=0, vmax=.15
    mlab.mesh(x, y, z, representation='wireframe', color=(.2, .2, .2), opacity=.1)

    # m = .1
    extent = [-m, m, -m, m, -m, m]
    ax = mlab.axes(extent=extent, xlabel='Ex (GPa)', ylabel='Ey (GPa)', zlabel='Ez (GPa)')
    # ax = mlab.axes(xlabel='Ex (GPa)', ylabel='Ey (GPa)', zlabel='Ez (GPa)')
    ax.axes.font_factor = 1
    mlab.outline(extent=extent)

    mlab.view(-120, 60, d, [0, 0, 0])

    # mlab.show()
    mlab.savefig('Young3D.png')
