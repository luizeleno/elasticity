import numpy as np
from scipy import constants

#Constants
BOLTZMANN = constants.Boltzmann
PI = constants.pi
AVOGADRO = constants.Avogadro
R = constants.gas_constant
RY_J = constants.physical_constants[ "Rydberg constant times hc in J"][0]
RY_EV = constants.physical_constants["Rydberg constant times hc in eV"][0]
EH_J = constants.physical_constants["Hartree energy"][0]
AU = constants.physical_constants["Bohr radius"][0]
mass_unit = constants.physical_constants["atomic mass constant"][0]
ANG = constants.physical_constants["angstrom"][0]
CONV_FACTOR = 10**9*AU**3/RY_J

#Energy conversion - Joule to Rydberg
def j_to_ry(energy):
    return energy/RY_J

#Energy conversion - Rydberg to Joule
def ry_to_j(energy):
    return energy*RY_J

#Energy conversion - Eletron Volts to Rydberg
def ev_to_ry(energy):
    return energy/RY_EV

#Energy conversion - Rydberg to Eletron Volt
def ry_to_ev(energy):
    return energy*RY_EV

#Energy conversion - Hartree to Rydberg
def eh_to_ry(energy):
    return energy*2.

#Energy conversion - Rydberg to Hartree
def ry_to_eh(energy):
    return energy/2.

#Volume conversion - Angstrom^3 to Bohr Radius^3
def ang_to_au(volume):
    return volume*ANG**3/AU**3

#Volume conversion - Bohr Radius^3 to Angstrom^3 
def au_to_ang(volume):
    return volume*AU**3/ANG**3

#Pressure conversion - Giga Pascal to Rydberg/(Bohr Radius^3)
def gpa_to_ry_au3(pressure):
    return pressure/CONV_FACTOR

#Pressure conversion - Rydberg/(Bohr Radius^3) to Giga Pascal 
def ry_au3_to_gpa(pressure):
    return pressure*CONV_FACTOR

#Pressure conversion - Giga Pascal to Dyne/cm^2
def gpa_to_ry_au3(pressure):
    return pressure/10**(-10)

#Pressure conversion - Dyne/cm^2 to Giga Pascal 
def gpa_to_ry_au3(pressure):
    return pressure*10**(-10)