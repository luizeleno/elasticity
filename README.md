# A new python library for micro- and macroscopic mechanical properties

## Luiz T. F. Eleno

### University of São Paulo, Lorena School of Engineering, Materials Eng. Dept.

The __elasticity__ python library helps to increase the workflow speed of ab-initio and experimentally determined mechanical properties investigations on crystalline systems.

The code has been implemented in python using the state-of-the-art numpy, scipy, matplotlib abd mayavi libraries. The __elasticity__ library can deal with any crystal structure, or even an isotropic continuous medium. As input, the user must provide the correct number of independent elastic constants for the system, taken from the literature, determined experimentally or using an ab-initio code. For instance, a cubic material has only three independent elastic constants, whereas a hexagonal system has five of them. 

Functions and routines of the library can then perform several different tasks, depending on the interests of the user.  For instance, there is a routine to calculate macroscopic properties such as Young’s, Shear and Bulk Moduli, as well as Poisson’s ratio, according to the well-known Voigt-Reuss-Hill approximation. The user can also determine the microscopic behavior of the previously mentioned properties, according to crystallographic directions, generating 3D plots in order to better investigate any possible anisotropy in the mechanical properties. It is also possible to determine wave velocities along high-symmetry crystallographic directions, with implications to the study of several phenomena, such as electron-phonon coupling effects in superconductors. Finally, the __elasticity__ library can determine the temperature-dependence of the Elastic Moduli mentioned above, besides the heat capacity and thermal expansion coefficient, following the Quasi-Harmonic (Debye-Grüneisen) approximation.

The __elasticity__ library is released and distributed on github as an open source package under the GNU general license.
