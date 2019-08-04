from elasticity import VRH
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar

#physical constants
from scipy import constants
BOLTZMANN = constants.Boltzmann
PI = constants.pi
AVOGADRO = constants.Avogadro
PLANCK = constants.Planck
REDUCED_PLANCK = constants.hbar


class StateEquation():

    '''
        Equations of state for solid materials:
        1. Murnaghan
        2. Birch-Mut=rnaghan
    '''
    
    def Murnaghan(self, *equilibrium_parameters):
        return E0 - ((B0*V0)/(dB0-1)) + ((B0*V)/dB0) * (1+((V0/V)**dB0)/(dB0-1)) 
        
    def BirchMurnaghan(self, *equilibrium_parameters):
        return E0 + ((9*V0*B0)/16) * \
            (dB0*((V0/V)**(2/3)-1)**3+(6-4*(V0/V)**(2/3))*((V0/V)**(2/3)-1)**2)

class DebyeGruneisen():

    '''
        Debye model applied using the Voigt-Reuss-Hill approximation
        to polycristalline elastic properties
    '''

    def __init__(self, GH, BH, density, atomic_mass, n_atoms, V0, temperature):
        self.GH = GH
        self.BH = BH
        self.density = density
        self.atomic_mass = atomic_mass
        self.n_atoms = n_atoms
        self.V0 = V0
        self.temperature = temperature

    def SoundVelocity(self, GH, BH, density): 
        #transversal sound velocity
        vs = np.sqrt(GH/density)  
        #longitudinal sound velocity
        vl = np.sqrt((BH+4/3*GH)/density) 
        #average sound velocity
        self.vm = ((2/vs**3+1/vl**3)/3) ** (-1/3) 
        return self.vm 

    def DebyeTemperature_0(self, density, atomic_mass, n_atoms):
        self.theta0 = (PLANCK/BOLTZMANN) * self.vm * \
            ((3*AVOGADRO*n_atoms*density)/(4*PI*atomic_mass)) ** (1/3)
        return self.theta0

    def GruneisenParameter(self, parameter_list):
        if self.temperature <= self.theta0:
            self.gruneisen = (1+dB0) * (1/2) - 1
        else:
            self.gruneisen = (1+dB0) * (1/2) - (2/3)
        return self.gruneisen

    def DebyeTemperature(self, V):
        self.theta = self.theta0 * (self.V0/V) ** (self.gruneisen)


class HelmholtzGruneisen(DebyeGruneisen):

    '''
        Calculation of the Helmholtz free energy and, consequently, the minimum volume and minimum
        energy at a given temperature. Therefore, calculation of the bulk modulus and the specific
        heats at constant volume and pressure.
    '''

    def integrand1(self, x):
        return (x**3)/(np.exp(x)-1)
    
    def DebyeFunction(self, V):
        I = quad(self.integrand1,0,self.DebyeTemperature(V)/self.temperature)
        self.debye_function = I[0] * 3 * (self.DebyeTemperature(V)/self.temperature) ** (-3)
        return self.debye_function

    def VibrationalEnergy(self, V):
        self.vibrational_energy = self.n_atoms * BOLTZMANN * self .temperature * \
            ((9*self.theta)/(8*self.temperature)+ \
                3*np.log(1-np.exp(-self.theta/self.temperature))-self.debye_function(V))
        return self.vibrational_energy

    def HelmholtzEnergy(self, V, method = "Murnaghan"):
        if method="Murnaghan":
            self.helmholtz_energy = self.VibrationalEnergy(V) + self.Murnaghan(V)
        elif method="Birch-Murnaghan":
            self.helmholtz_energy = self.VibrationalEnergy(V) + self.BirchMurnaghan(V)
        else:
            pass
        return self.helmholtz_energy

    def MinimumVolume(self):
        minimum_energy = minimize_scalar(self.HelmholtzEnergy, \
            bracket(0.9*self.V0,self.V0), method="brent")
        self.minimum_volume = minimum_energy.x
        return self.minimum_volume

    def BulkModulus(self):
        pass
    
    def Cv(self):
        pass

    def Cp(self):
        pass