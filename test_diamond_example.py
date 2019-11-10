import elasticity as el
import thermal_expansion as thermal
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import constants
R = constants.gas_constant
import time

#Materials properties
material = "diamond"
structure = "cubic"
Z = 2
atomic_mass, atomic_mass_unit = 12.0107, "u"
C11, C12, C44, elastic_constants_unit = 1079, 124, 578, "GPa"
chi = 1.

#Experimental data for specific heat at constant pressure
temperatures_exp = np.array([70.16, 75.37, 81.59, 88.65, 96.68, 105.1, 113.04, 
                        125.28, 134.29, 144.1, 153.71,162.76, 173.33, 181.96, 
                        191.44, 200.94, 211.84, 231.06, 241.09, 252.37, 264.31, 
                        276.61, 287.96, 300, 320, 340, 360, 380, 400, 420, 440, 
                        460, 480, 500, 550, 600,650, 700, 750, 800, 850, 900, 
                        950, 1000, 1050, 1100])

cp_exp = np.array([0.0938, 0.124, 0.152, 0.189, 0.248, 0.329, 0.415, 0.578, 0.729,
                 0.913, 1.118, 1.331,1.616, 1.859, 2.169, 2.491, 2.868, 3.592, 
                 3.843, 4.325, 4.760, 5.263, 5.673, 6.196, 7.030, 7.863, 8.679, 
                 9.475, 10.24, 10.97, 11.67, 12.34, 12.97, 13.57, 14.94, 16.13,
                 17.16, 18.05, 18.83, 19.51, 20.11, 20.65, 21.14, 21.61, 22.06, 
                 22.52])

temperatures_exp_unit = 'K'
cp_exp_unit = 'J/molK'

#Ab initio data
energies = np.array([-22.76745769, -22.76213240, -22.75431180, -22.77047178, 
                -22.77134749, -22.77024716, -22.76732281])
volumes = np.array([72.2765903507, 70.1081472678, 67.9835155729, 74.4892919379, 
                76.7466991459, 79.0492590911, 81.39741889])
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
T = np.linspace(1, 4000, 4000)
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

#Plots
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