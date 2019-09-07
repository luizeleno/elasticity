from elasticity import VRH
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
from scipy.optimize import curve_fit
import numpy as np

#physical constants
from scipy import constants
BOLTZMANN = constants.Boltzmann
PI = constants.pi
AVOGADRO = constants.Avogadro
PLANCK = constants.Planck
REDUCED_PLANCK = constants.hbar
R = constants.gas_constant

class EOS():
    '''
        Calculation of the equations of state for solids:
        - Murnaghan equation
        - Birch-Murnaghan equation
        It is necessary to provide the volumes and energies from experimental data.
    '''

    def __init__(self, energies, volumes, method = "Murnaghan"):
        self.energies = energies
        self.volumes = volumes
        self.method = method
        self.E0, self.B0, self.Bp, self.V0 = self.fit_eos()

    #Polynomial fit to guess the initial parameters for the equations of state
    def polynomial_fit(self):
        a, b, c = np.polyfit(self.volumes, self.energies, 2)
        V0 = -b/(2*a)
        B0 = 2*a*V0
        E0 = a*V0 + b*V0 + c
        Bp = 2. #assumption
        return E0, B0, Bp, V0

    #Equations of state
    def murnaghan(self, V, *parameters):
        E0, B0, Bp, V0 = parameters 
        return E0 + (B0*V0)/(Bp-1) + ((B0*V)/Bp) * (1+((V0/V)**Bp)/(Bp-1))

    def birch_murnaghan(self, V, *parameters):
        E0, B0, Bp, V0 = parameters
        return E0 + (9/16*B0*B0) * (Bp*((V0/V)**(2./3)-1)**3 + 
                    (6-4*(V0/V)**(2./3))*((V0/V)**(2./3)-1)**2)

    #Fit of the equations of state
    def fit_eos(self):
        if self.method == "Murnaghan":
            popt, pcov = curve_fit(self.murnaghan, self.volumes, self.energies, 
                         p0 = self.polynomial_fit(), method = 'trf')
        elif self.method == "Birch-Murnaghan":
            popt, pcov = curve_fit(self.birch_murnaghan, self.volumes, self.energies, 
                         p0 = self.polynomial_fit(), method = 'trf')
        else:
            popt, pcov = curve_fit(self.murnaghan, self.volumes, self.energies, 
                         p0 = self.polynomial_fit(), method = 'trf')
            print("Only the Murnaghan ad Birch-Murnaghan equations are implemented. \
                    We run with Murnaghan")
        return popt[0], popt[1], popt[2], popt[3]

class Debye():
    '''
        Implementation of the Debye Model, considering the VRH Theory.
        The shear modulus GH and the bulk modulus BH can be inserted 
        by the user or can be calculated from the VRH class, implemented in
        the "elasticity" module.
    '''

    def __init__(self, G, B, density, molar_volume):
        self.G = G
        self.B = B
        self.density = density
        self.molar_volume = molar_volume
        self.theta0 = self.temp_debye0()

    #The Debye sound velocity
    def sound_velocity(self):
        v_trans = np.sqrt(self.G/self.density)
        v_long = np.sqrt((3*self.B + 4*self.G)/(3*self.density))
        return (1/3 * (2/v_trans**3 + 1/v_long**3)) ** (-1./3)

    #The Debye temperature by the Debye Model
    def temp_debye0(self):
        return (PLANCK/BOLTZMANN) * self.sound_velocity() * \
                      ((6*PI*AVOGADRO)/self.molar_volume) ** (1./3)


class QuasiHarmonic():
    '''
        Thermal expansion of the material, using the quasi-harmonic approximation.
    '''

    def __init__(self, T, G, B, density, molar_volume, chi, energies, volumes, method = "Murnaghan"):
        self.T = T
        self.chi = chi
        self.method = method
        self.debye_class = Debye(G, B, density, molar_volume)
        self.theta0 = self.debye_class.theta0
        self.eos_class = EOS(energies, volumes)
        self.E0, self.B0, self.Bp, self.V0 = self.eos_class.fit_eos()
 
    #The Gruneisen parameter
    def gru_coef(self):
        return (1+self.Bp)/2 - self.chi

    #The Debye temperature in the Debye-Gruneisen model
    def theta(self, V):
        return self.theta0 * (self.V0/V)**self.gru_coef()

    #The Debye function
    def integrand_debye(self, x):
        return x**3 / (np.exp(x) - 1)

    def debye_func(self, V):
        I = quad(self.integrand_debye, 0, self.theta(V)/self.T)
        return I[0] * 3 * (self.T/self.theta(V))**3

    #The vibrational energy in the quasi-harmonic approximation
    def F_vib(self, V):
        return 9/8*BOLTZMANN*self.theta(V) + BOLTZMANN*self.T * \
                (3 * np.log(1 - np.exp(-self.theta(V)/self.T)) - self.debye_func(V))
    
    #Calculation of the Helmholtz free energy
    def F(self, V):
        if self.method == "Murnaghan":
            return self.F_vib(V) + self.eos_class.murnaghan(V, self.E0, self.B0, self.Bp, self.V0)
        elif self.method == "Birch-Murnaghan":
            return self.F_vib(V) + self.eos_class.birch_murnaghan(V, self.E0, self.B0, self.Bp, self.V0)
        else:
            print("Only the Murnaghan ad Birch-Murnaghan equations are implemented. \
                    We run with Murnaghan.")            
            return self.F_vib(V) + self.eos_class.murnaghan(V, self.E0, self.B0, self.Bp, self.V0)

    #Mininmization of the Helmholtz free energy
    def V_min(self):
        m = minimize_scalar(self.F,   
                            bracket=(0.9*self.V0, 1.1*self.V0),   
                            method="brent")
        return m.x

    #The Debye temperature at the new equilibrium position        
    def theta_eq(self):
        return self.theta(self.V_min())

    #Bulk modulus after the thermal expansion
    def bulk_modulus(self):
        #variables
        theta = self.theta_eq()
        gamma = self.gru_coef()
        D = self.debye_func(self.V_min)
        v = self.V0/self.V_min()
        if self.method == "Murnaghan":
            B = self.B0*v**self.Bp + 1/self.V_min() * \
                (9/8*BOLTZMANN*theta*gamma*(1+gamma)+3*(1-3*gamma)*BOLTZMANN*self.T*D + 
                (9*gamma**2*BOLTZMANN*theta)/(np.exp(theta*self.T)-1))
        elif self.method == "Birch-Murnaghan":
            B = 1/2*self.B0*(7*v**(7/3)-5*v**(5/3)) + 3/8*self.B0*(self.Bp-4) * \
                (9*v**3-14*v**(7/3) + 5*v**(5/3)) + 1/self.V_min() * \
                (9/8*BOLTZMANN*theta*gamma*(1+gamma)+3*(1-3*gamma)*BOLTZMANN*self.T*D + 
                (9*gamma**2*BOLTZMANN*theta)/(np.exp(theta*self.T)-1))
        else:
            print("Only the Murnaghan ad Birch-Murnaghan equations are implemented. \
                    We run with Murnaghan.")
            B = self.B0*v**self.Bp + 1/self.V_min() * \
                (9/8*BOLTZMANN*theta*gamma*(1+gamma)+3*(1-3*gamma)*BOLTZMANN*self.T*D + 
                (9*gamma**2*BOLTZMANN*theta)/(np.exp(theta*self.T)-1))
        return B

    #The specific heat at constant volume
    def integrand_heat(self, x):
        return ((x**4) * np.exp(x)) / ((np.exp(x) - 1)**2)

    def heat_debye(self):
        I = quad(self.integrand_heat,0,self.theta_eq()/self.T)
        return I[0] * (9*AVOGADRO*BOLTZMANN) * (self.T/self.theta_eq())**3

    #Thermal expansion coefficient at constant volume
    def thermal_coef(self):
        return self.gru_coef()*self.heat_debye()/self.V_min()/self.bulk_modulus()

    #The specific heat at constant pressure
    def heat_presure(self):
        return self.heat_debye() + self.thermal_coef()**2*self.bulk_modulus()*self.V_min()*self.T


