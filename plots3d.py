from mayavi import mlab
import numpy as np


def mayaviplot(property, N=151, m=.1, d=.7, labels=['x (GPa)', 'y (GPa)', 'z (GPa)'], offscreen=True, filename='Young3D.png'):
    '''
        m : extent of axes
        d: POV distance
    '''

    N = int(N) * 1j
    theta, phi = np.mgrid[0:np.pi:N, 0:2*np.pi:N]
    ldir = np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])
    case = property(ldir)

    x = case * ldir[0]
    y = case * ldir[1]
    z = case * ldir[2]

    if offscreen:
        mlab.options.offscreen = True

    mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(1050, 1000))
    mlab.clf()

    mlab.mesh(x, y, z, scalars=case, colormap='jet')  # , vmin=0, vmax=.15
    mlab.mesh(x, y, z, representation='wireframe', color=(.2, .2, .2), opacity=.1)

    extent = [-m, m, -m, m, -m, m]
    ax = mlab.axes(extent=extent, xlabel=labels[0], ylabel=labels[1], zlabel=labels[2])
    ax.axes.font_factor = 1
    mlab.outline(extent=extent)

    mlab.view(-120, 60, d, [0, 0, 0])

    if offscreen:
        mlab.savefig(filename)
    else:
        mlab.show()

def YoungModulus3DPlot(material, *args, **kwargs):

    E = material.YoungModulus

    mayaviplot(E, *args, **kwargs)

def LinearCompressibility3DPlot(material, *args, **kwargs):

    B = material.LinearCompressibility

    mayaviplot(B, *args, **kwargs)
