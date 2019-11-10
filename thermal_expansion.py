from elasticity import VRH
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
from scipy.optimize import curve_fit
import numpy as np

#Physical constants
from scipy import constants
BOLTZMANN = constants.Boltzmann
PI = constants.pi
AVOGADRO = constants.Avogadro
REDUCED_PLANCK = constants.hbar
R = constants.gas_constant
RY_J = constants.physical_constants[ "Rydberg constant times hc in J"][0]
AU = constants.physical_constants["Bohr radius"][0]
MASS_UNIT = constants.physical_constants["atomic mass constant"][0]
CONV_FACTOR = 10**9*AU**3/RY_J

'''
Default units:
 ----------
 [temperatures] = `Kelvin`
 [volume] = `cubic Bohr radius`
            r'1au$^3$' = 0.14818471r'A$^3$'
 [energy] = `Rydberg`
            1Ry = 2.179872325e-18J
 [elastic properties] = `giga Pascal`
            1GPa = 1e9Pa
 [specific heats] = `factor of R`
            1R = 8.3144598r'Jmol$^{-1}$K$^{-1}$'
 [thermal expansion coefficient] = `1e-6r'K$^{-1}$'`
'''


class EOS():
    '''
    Implemented equations of state (EOS) for solids:
    ----------
    >>> Murnaghan
    >>> Birch-Murnaghan

    Entrance data:
    ----------
    >>> energies = experimental or ab initio energies
    >>> volumes = experimental or ab initio volumes
    >>> method = chosen equation of state (eos):
                Murnaghan (default) or Birch-Murnaghan

    Calculated parameters from the eos fit:
    ----------
    >>> r'E$_0$': energy
    >>> r'B$_0$': bulk modulus
    >>> r'B$_{p0}$': first derivative of the bulk modulus
    >>> r'V$_0$': volume
    '''

    def __init__(self, energies, volumes, method = "Murnaghan"):
        self.energies = energies
        self.volumes = volumes
        self.method = method
        self.E0, self.B0, self.Bp, self.V0 = self.fit_eos()

    def polynomial_fit(self):
        '''
        Second order polynomial fit for initial guess of the eos parameters.
        '''
        a, b, c = np.polyfit(self.volumes, self.energies, 2)
        V0 = -b/(2*a)
        B0 = 2*a*V0
        E0 = a*V0 + b*V0 + c
        Bp = 2. #assumption
        return E0, B0, Bp, V0

    def murnaghan(self, V, *parameters):
        '''
        Murnaghan EOS.
        '''
        E0, B0, Bp, V0 = parameters 
        return E0 + (B0*V0)/(Bp-1) + ((B0*V)/Bp) * (1+((V0/V)**Bp)/(Bp-1))

    def birch_murnaghan(self, V, *parameters):
        '''
        Birch-Murnaghan EOS.
        '''
        E0, B0, Bp, V0 = parameters
        return E0 + (9/16*B0*V0) * (Bp*((V0/V)**(2./3)-1)**3 + 
                    (6-4*(V0/V)**(2./3))*((V0/V)**(2./3)-1)**2)

    #Fit of the EOSs
    def fit_eos(self):
        '''
        Data fit with the curve_fit method from Scipy.
        The fit is done using the chosen EOS and the initial guess for the
        parameters from the polynomial fit.
        '''
        if self.method == "Birch-Murnaghan":
            popt, pcov = curve_fit(self.birch_murnaghan, self.volumes, \
                        self.energies, p0 = self.polynomial_fit(), \
                        method = 'trf')
        else:
            popt, pcov = curve_fit(self.murnaghan, self.volumes, 
                        self.energies, p0 = self.polynomial_fit(), 
                        method = 'trf')
        return popt[0], popt[1], popt[2], popt[3]

class Debye():
    '''
        Implementation of the Debye Model, considering the elasticity theory.
        Entrance parameters:
        ----------
        >>> G = Shear Modulus
        >>> B = Bulk Modulus
        >>> r'V$_0$' = cell unit volume
        >>> Z =  number of formula units in the unit cell
        >>> atomic_mass = atomic mass
        Results:
        ----------
        >>> Speed of sound
        >>> Debye temperature
    '''

    def __init__(self, G, B, V0, Z, atomic_mass):
        self.G = G*1e9
        self.B = B*1e9
        self.molar_volume = (V0*AU**3) * AVOGADRO / Z
        self.density = atomic_mass * 1e-3 / self.molar_volume

    def sound_velocity(self):
        '''
        Debye speed of sound.
        Calculation using de the elastic properties B and G.
        '''
        v_trans = np.sqrt(self.G/self.density)
        v_long = np.sqrt((3.*self.B + 4.*self.G)/(3*self.density))
        return np.power(1./3.*(2./v_trans**3 + 1./v_long**3), -1./3.)

    def debye_temp(self):
        '''
        Debye temperature (Debye Model)
        '''
        return (REDUCED_PLANCK/BOLTZMANN) * self.sound_velocity() * \
            np.power((6*np.power(PI, 2)*AVOGADRO)/self.molar_volume, 1./3)

class QuasiHarmonic():

    def __init__(self, chi, poisson, debye_temp, energies, volumes,
                 method = "Murnaghan"):
        #In parameters
        self.method = method
        self.chi = chi
        self.poisson = poisson
        self.theta0 = debye_temp
        self.eos_class = EOS(energies, volumes)
        self.E0, self.B0, self.Bp, self.V0 = self.eos_class.fit_eos()
        #Used parameters
        self.gamma = self.gru_coef()

    def gru_coef(self):
        '''
        Gruneisen parameter.
        '''
        return (1+self.Bp) / 2 - self.chi

    def theta(self, V):
        '''
        Debye temperature (Debye-Gruneisen model).
        '''
        return self.theta0 * np.power(self.V0/V, self.gamma)

    def integrand_debye(self, x):
        '''
        Integrand of Debye function.
        '''
        return np.power(x,3) / (np.exp(x)-1)

    def debye_func(self, V, temp):
        '''
        Debye Function.
        '''
        I = quad(self.integrand_debye, 0, self.theta(V)/temp)
        return 3. * I[0] * np.power(temp/self.theta(V), 3)

    def F_vib(self, V, temp):
        '''
        Helmholtz vibrational energy (quasi-harmonic approximation).
        '''
        return 1. / RY_J * (9./8.*BOLTZMANN*self.theta(V)+BOLTZMANN*temp* \
            (3.*np.log(1-np.exp(-self.theta(V)/temp))-self.debye_func(V, temp)))

    def F(self, V, temp):
        '''
        Helmholtz free energy: vibrational energy + EOS energy.
        The model does not consider the electronic contribution.
        The EOS energy can be from Murnaghan or Birch-Murnaghan equations.
        '''
        if self.method == "Birch-Murnaghan":
            return self.F_vib(V, temp) + \
                self.eos_class.birch_murnaghan(V, self.E0, self.B0, self.Bp, self.V0)
        else:        
            return self.F_vib(V, temp) + \
                self.eos_class.murnaghan(V, self.E0, self.B0, self.Bp, self.V0)

    def V_min(self, temp):
        '''
        Minimization of the Helmholtz Free Energy (F) using the minimize_scalar
        method from Scipy.
        The result is the value of the volume that minimizes the energy F at a 
        given temperature T.
        '''
        m = minimize_scalar(lambda x: self.F(x, temp),   
                            bracket=(0.9*self.V0, 1.1*self.V0),   
                            method="brent")
        return m.x
     
    def bulk_modulus(self, vol, temp, theta, D):
        '''
            Bulk modulus calculated using the second derivative of the Helmholtz 
            free energy with volume.
            It is consider an variation of the bulk modulus with temperature.
        '''
        v = self.V0 / vol
        B0 = self.B0 * 1e9 / CONV_FACTOR
        if self.method == "Birch-Murnaghan":
            B = 1/2*B0*(7*v**(7/3)-5*v**(5/3)) + 3/8*B0*(self.Bp-4) * \
                (9*v**3-14*v**(7/3) + 5*v**(5/3)) + 1/vol * \
                (9/8*BOLTZMANN*theta*self.gamma*(1+self.gamma)+
                3*(1-3*self.gamma)*BOLTZMANN*temp*D + 
                (9*self.gamma**2*BOLTZMANN*theta)/(np.exp(theta/temp)-1))
        else:
            B = B0*v**self.Bp + 1/vol * \
                (9/8*BOLTZMANN*theta*self.gamma*(1+self.gamma)+
                3*(1-3*self.gamma)*BOLTZMANN*temp*D + 
                (9*self.gamma**2*BOLTZMANN*theta)/(np.exp(theta/temp)-1))
        return B / 1e9

    def young_modulus(self, poisson, B):
        '''
            Young Modulus, considering an constant Poisson coefficient.
        '''
        return 3 * B / (1-2*poisson)

    def integrand_heat(self, x):
        '''
            Specific heat integrand.
        '''
        return np.power(x, 4.) * np.exp(x) / np.power(np.exp(x)-1., 2.)

    def heat_volume(self, vol, temp, theta):
        '''
            Specific heat at constant volume.
        '''
        I = quad(self.integrand_heat, 0, theta/temp)
        return 9. * I[0] * AVOGADRO * BOLTZMANN * np.power(temp/theta, 3)

    def thermal_coef(self, vol, B, cv):
        '''
            Thermal expansion coefficient at constant volume.
        '''
        return 1/AVOGADRO * self.gamma * cv / (vol*(AU**3)*B*1e9)

    def heat_pressure(self, vol, temp, B, cv, alpha):
        '''
           Specific heat at constant pressure. 
        '''
        return cv + np.power(alpha, 2) * B * 1e9 * vol * (AU**3) * temp * AVOGADRO

class QuasiHarmonicVector(QuasiHarmonic):
    '''
        Thermal expansion of the material, using the quasi-harmonic 
        approximation.

        Entrance parameters:
        ----------
        >>> T = temperature
        >>> chi = constant to calculate de poisson coefficient
        >>> poisson = Poisson coefficient at 0K
        >>> debye_temp = Debye temperature at 0K
        >>> energies = energies at 0K
        >>> volumes = volumes at 0K
        >>> method = chosen equation of state (eos): 
                     Murnaghan or Birch-Murnaghan

        Dependency on others classes:
        ----------
        >>> thermal_expansion.EOS (required)
        >>> thermal_expansion.Debye (optional)
        >>> elasticity.VRH (optional)
        >>> thermal_expansion.QuasiHarmonic

        Results:
        ----------
        >>> Vmin
        >>> Fmin
        >>> B
        >>> E
        >>> Theta
        >>> alpha
        >>> cv
        >>> cp
    '''

    def __init__(self, T, chi, poisson, debye_temp, energies, volumes,
                 method = "Murnaghan"):
        #In parameters
        self.method = method
        self.chi = chi
        self.poisson = poisson
        self.theta0 = debye_temp
        self.eos_class = EOS(energies, volumes)
        self.E0, self.B0, self.Bp, self.V0 = self.eos_class.fit_eos()
        self.T = T
        #Out parameters
        self.gamma = self.gru_coef()
        self.Vmin = self.V_min_vec()
        self.F = self.F_vec()
        self.Theta = self.theta_vec()
        self.D = self.debye_func_vec()
        self.B = self.bulk_modulus_vec()
        self.E = self.young_modulus_vec()
        self.cv = self.heat_volume_vec()
        self.alpha = self.thermal_coef_vec()
        self.cp = self.heat_pressure_vec()

    def V_min_vec(self):
        '''
            Array of the minimum volume at the equilibrium points.
        '''
        return list(map(self.V_min, self.T))

    def F_vec(self):
        '''
            Array of the Helmholtz free energy at the equilibrium points.
        '''
        return list(map(self.F, self.Vmin, self.T))

    def theta_vec(self):
        '''
            Array of the Debye temperature at the equilibrium points.
        '''
        return list(map(lambda x: self.theta(x), self.Vmin))

    def debye_func_vec(self):
        '''
        Debye Function.
        '''
        t = self.Theta/self.T
        I = quad(self.integrand_debye, 0, t.all())
        return I[0] * 3. * np.power(self.T/self.Theta, 3)

    def bulk_modulus_vec(self):
        '''
            Array of the bulk modulus at the equilibrium points.
        '''
        return list(map(lambda x, y, w, z: self.bulk_modulus(x, y, w, z),
            self.Vmin, self.T, self.Theta, self.D))

    def young_modulus_vec(self):
        '''
            Array of the Young modulus at the equilibrium points.
        '''
        return list(map(lambda x: self.young_modulus(self.poisson, x), self.B))

    def heat_volume_vec(self):
        '''
            Array of the specific heat at constant volume
            at the equilibrium points.
        '''
        return list(map(lambda x, y, z: self.heat_volume(x, y, z),
            self.Vmin, self.T, self.Theta))

    def thermal_coef_vec(self):
        '''
            Array of the thermal expansion coefficient at constant volume
            at the equilibrium points.
        '''
        return list(map(lambda x, y, z: self.thermal_coef(x, y, z),
            self.Vmin, self.B, self.cv))

    def heat_pressure_vec(self):
        '''
            Array of the specific heat at constant pressure
            at the equilibrium points.
        '''
        return list(map(lambda v, x, y, w, z: self.heat_pressure(v, x, y, w, z),
            self.Vmin, self.T, self.B, self.cv, self.alpha))