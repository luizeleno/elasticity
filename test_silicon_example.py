import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import elasticity as el
import thermal_expansion as thermal
from scipy import constants
import time

#Font size for the graphics
plt.rcParams.update({'font.size': 14})

#Material properties  doi = "10.17188/1190959"
material = "silicon"
structure = "cubic"
Z = 2
atomic_mass, atomic_mass_unit = 28.0855, "u"
C11, C12, C44, elastic_constants_unit = 162.07, 63.51, 77.26, "GPa"
chi = 1.

#Experimental data for specific heat at constant pressure
temperatures_exp = np.array([298.15, 300, 400, 500, 600, 700, 800, 900, 1000, 
                1100, 1200,1300, 1400, 1500, 1600, 1687])

cp_exp = np.array([20.007,  20.006,  22.258,  23.588,  24.420,  25.050,  25.608,
                26.110, 26.569, 26.988,  27.360, 27.707, 28.045,  28.372, 28.674, 
                28.930])

#Ab initio data
'''
    Calculated by the author.
'''
energies = np.array([-19.19110706, -19.18747301, -19.18192993, -19.19300108, 
                -19.19331307, -19.19219081, -19.18977257])
volumes = np.array([262.146575223, 254.284242613, 246.573276899, 270.169326182, 
                278.354115558, 286.702563418, 295.224668962])
energies_unit = "Ry"
volumes_unit = "au3"  

#VRH Theory
VRH = el.VRH(structure, C11, C12, C44)
VRH.VRH()
poisson = VRH.PoissonRatio
print('Poisson ratio: ', VRH.PoissonRatio,
    '\n' 'GH: ', VRH.ShearModulus, 
    '\n' 'BH: ', VRH.BulkModulus)

#EOS fit and plot
eos_murnaghan = thermal.EOS(energies, volumes, method = "Murnaghan")
eos_birch = thermal.EOS(energies, volumes, method = "Birch-Murnaghan")
print("Murnaghan fit: ", eos_murnaghan.fit_eos())
print("Birch fit: ", eos_birch.fit_eos())

#Debye
V0 = eos_murnaghan.fit_eos()[3]
debye = thermal.Debye(VRH.ShearModulus, VRH.BulkModulus,V0 , Z, atomic_mass)
print("Debye Temp: ", debye.debye_temp(),
    '\n' "Density: ", debye.density,
    '\n' "Molar volume: ", debye.molar_volume,
    '\n' "Sound velocity: ", debye.sound_velocity())

#QuasiHarmonic
debye_temp = debye.debye_temp()
debye_temp_unit = 'K'
T = np.linspace(1, 2000, 2000)
T_unit = "K"
a = time.time()
#Comparison to experimental data
qha1_mur = thermal.QuasiHarmonicVector(temperatures_exp, chi, poisson, debye_temp, 
                            energies, volumes, method = 'Murnaghan')
data1_mur = {
        'Temperature': temperatures_exp, 
        'Cv Murnaghan': qha1_mur.cv,
        'Cp Murnaghan': qha1_mur.cp,
        'Cp experimental': cp_exp,
        }
df1_mur = pd.DataFrame(data1_mur, columns= ['Temperature', 'Cv Murnaghan', 
        'Cp Murnaghan', 'Cp experimental'])
df1_mur.to_csv(r'C:\Users\danie\OneDrive\Documentos\GitHub\elasticity\examples\Thermal\cp_data_murnaghan_' + material + '.txt', 
        index = None, header=True) 

qha1_birch = thermal.QuasiHarmonicVector(temperatures_exp, chi, poisson, debye_temp, 
                            energies, volumes, method = 'Murnaghan')
data1_birch = {
        'Temperature': temperatures_exp, 
        'Cv Birch-Murnaghan': qha1_birch.cv,
        'Cp Birch-Murnaghan': qha1_birch.cp,
        'Cp experimental': cp_exp,
        }
df1_birch = pd.DataFrame(data1_birch, columns= ['Temperature', 'Cv Birch-Murnaghan', 
            'Cp Birch-Murnaghan', 'Cp experimental'])
df1_birch.to_csv(r'C:\Users\danie\OneDrive\Documentos\GitHub\elasticity\examples\Thermal\cp_data_birch_' + material + '.txt', 
        index = None, header=True) 

qha2_mur = thermal.QuasiHarmonicVector(T, chi, poisson, debye_temp, 
                            energies, volumes, method = 'Murnaghan')
print('gamma: ', qha2_mur.gamma)

data2_mur = {
        'Temperature': T, 
        'Volume': qha2_mur.Vmin,
        'F': qha2_mur.F,
        'Debye temperature': qha2_mur.Theta,
        'B': qha2_mur.B,
        'E': qha2_mur.E,
        'alpha': qha2_mur.alpha,
        'Cv': qha2_mur.cv,
        'Cp': qha2_mur.cp,
        }

df2_mur = pd.DataFrame(data2_mur, columns= ['Temperature', 'Volume', 'F', 
        'Debye Temperature', 'B', 'E', 'alpha', 'Cv', 'Cp'])

df2_mur.to_csv(r'C:\Users\danie\OneDrive\Documentos\GitHub\elasticity\examples\Thermal\data_murnaghan_' + material + '.txt', 
        index = None, header=True) 

qha2_birch = thermal.QuasiHarmonicVector(T, chi, poisson, debye_temp, 
                            energies, volumes, method = 'Murnaghan')

data2_birch = {
        'Temperature': T, 
        'Volume': qha2_birch.Vmin,
        'F': qha2_birch.F,
        'Debye temperature': qha2_birch.Theta,
        'B': qha2_birch.B,
        'E': qha2_birch.E,
        'alpha': qha2_birch.alpha,
        'Cv': qha2_birch.cv,
        'Cp': qha2_birch.cp,
        }

df2_birch = pd.DataFrame(data2_birch, columns= ['Temperature', 'Volume', 'F',
         'Debye Temperature', 'B', 'E', 'alpha', 'Cv', 'Cp'])

df2_birch.to_csv(r'C:\Users\danie\OneDrive\Documentos\GitHub\elasticity\examples\Thermal\data_birch_' + material + '.txt', 
        index = None, header=True) 

plt.figure(1)
plt.autoscale(tight=False)
plt.ticklabel_format(useOffset=False)
plt.xlabel("Temperatura (K)")
plt.ylabel("Calor específico (J/mol.K)")
plt.plot(T, qha2_mur.cv, color = 'b', label = r"c$_V$ calculado")
plt.plot(T, qha2_mur.cp, color = 'g', label = r"c$_P$ calculado")
plt.plot(temperatures_exp, cp_exp, color = 'r', label = r"c$_P$ experimental")
plt.tight_layout()
plt.legend(loc = 4)
plt.savefig(r'C:\Users\danie\OneDrive\Documentos\GitHub\elasticity\examples\Thermal\comparative_specific_heat_murnaghan_' + material + '.pdf')

plt.figure(2)
plt.autoscale(tight=False)
plt.ticklabel_format(useOffset=False)
plt.xlabel("Temperatura (K)")
plt.ylabel("Calor específico (J/mol.K)")
plt.plot(T, qha2_birch.cv, color = 'b', label = r"c$_V$ calculado")
plt.plot(T, qha2_birch.cp, color = 'g', label = r"c$_P$ calculado")
plt.plot(temperatures_exp, cp_exp, color = 'r', label = r"c$_P$ experimental")
plt.tight_layout()
plt.legend(loc = 4)
plt.savefig(r'C:\Users\danie\OneDrive\Documentos\GitHub\elasticity\examples\Thermal\comparative_specific_heat_birch_' + material + '.pdf')

plt.figure(3)
plt.autoscale(tight=False)
plt.ticklabel_format(useOffset=False)
plt.xlabel("Temperatura (K)")
plt.ylabel(r"$\Theta _D$ (K)")
plt.plot(T, qha2_mur.Theta, color = 'b', label = r"Murnaghan EOS")
plt.plot(T, qha2_birch.Theta, color = 'g', label = r"Birch-Murnaghan EOS")
plt.tight_layout()
plt.legend(loc = 4)
plt.savefig(r'C:\Users\danie\OneDrive\Documentos\GitHub\elasticity\examples\Thermal\debye_temperature' + material + '.pdf')

print('time to run: ', time.time()-a)