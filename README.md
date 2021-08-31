# A new python library for micro- and macroscopic mechanical properties

## Luiz T. F. Eleno

### University of São Paulo, Lorena School of Engineering, Materials Eng. Dept.

## Prerequisites
- mayavi
- numpy
- scipy
- matplotlib
- plot3d
- plotly
- pandas
- vtk

To install with pip or anaconda,

```shell
pip install mayavi numpy scipy matplotlib plot3d plotly pandas vtk
```

```shell
conda install mayavi numpy scipy matplotlib plot3d plotly pandas vtk
```

## Brief explanation about the library content
The __elasticity__ python library helps to increase the workflow speed of ab-initio and experimentally determined mechanical properties investigations on crystalline systems.

The code has been implemented in python using the state-of-the-art numpy, scipy, matplotlib abd mayavi libraries. The __elasticity__ library can deal with any crystal structure, or even an isotropic continuous medium. As input, the user must provide the correct number of independent elastic constants for the system, taken from the literature, determined experimentally or using an ab-initio code. For instance, a cubic material has only three independent elastic constants, whereas a hexagonal system has five of them. 

### What the library covers
Functions and routines of the library can then perform several different tasks, depending on the interests of the user.  For instance, there is a routine to calculate macroscopic properties such as __Young’s__, __Shear and Bulk Moduli__, as well as __Poisson’s ratio__, according to the well-known __Voigt-Reuss-Hill__ approximation. 

The user can also determine the microscopic behavior of the previously mentioned properties, according to crystallographic directions, generating 3D plots in order to better investigate any possible anisotropy in the mechanical properties. It is also possible to determine wave velocities along high-symmetry crystallographic directions, with implications to the study of several phenomena, such as electron-phonon coupling effects in superconductors. 

Finally, the __elasticity__ library can determine the temperature-dependence of the Elastic Moduli mentioned above, besides the heat capacity and thermal expansion coefficient, following the __Quasi-Harmonic (Debye-Grüneisen)__ approximation.

## How to use it

As an example, visualize Young's module heatmap of solid H20,

```python
python material-examples/H20.py 
```

## The directory structure
```shell
tree -d .
```
```shell
.
|-- __pycache__
|-- library
|-- material-examples
|-- outputs
`-- plot-libs-example
```
In the,
- `library` resides the implementation of the functions;
- `material-examples` python file with ready to run parameters;
- `output` saved plots from the `material-examples`;
- `plot-libs-examples` some examples of the plot framework used in our library;

## License
The __elasticity__ library is released and distributed on github as an open source package under the GNU general license.


