import elasticity as el
import thermal_expansion as thermal
import numpy as np
import matplotlib.pyplot as plt
from pandas import DataFrame
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
temperatures_exp = np.array([70.16, 75.37, 81.59, 88.65, 96.68, 105.1, 113.04, 125.28, 134.29, 144.1, 153.71, 
                        162.76, 173.33, 181.96, 191.44, 200.94, 211.84, 231.06, 241.09, 252.37, 264.31, 
                        276.61, 287.96, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 550, 600,
                        650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100])

cp_exp = np.array([0.0938, 0.124, 0.152, 0.189, 0.248, 0.329, 0.415, 0.578, 0.729, 0.913, 1.118, 1.331,
                    1.616, 1.859, 2.169, 2.491, 2.868, 3.592, 3.843, 4.325, 4.760, 5.263, 5.673, 6.196, 
                    7.030, 7.863, 8.679, 9.475, 10.24, 10.97, 11.67, 12.34, 12.97, 13.57, 14.94, 16.13,
                    17.16, 18.05, 18.83, 19.51, 20.11, 20.65, 21.14, 21.61, 22.06, 22.52])
temperatures_exp_unit = 'K'
cp_exp_unit = 'J/molK'

#Ab initio data
energies = np.array([-22.76745769, -22.76213240, -22.75431180, -22.77047178, -22.77134749,
                    -22.77024716, -22.76732281])
volumes = np.array([72.2765903507, 70.1081472678, 67.9835155729, 74.4892919379, 76.7466991459, 
                    79.0492590911, 81.39741889])
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
    '\n' "Molar volume: ", debye.molar_volume)

#QuasiHarmonic
debye_temp = debye.debye_temp()
debye_temp_unit = 'K'
T = np.linspace(1, 1000, 10)
T_unit = "K"
a = time.time()
#Comparison to experimental data
qha1 = thermal.QuasiHarmonicVector(temperatures_exp, chi, poisson, debye_temp, 
                            energies, volumes, method = 'Murnaghan')
data1 = {'Temperature': temperatures_exp, 
        'Cv Murnaghan': qha1.cv,
        'Cp Murnaghan': qha1.cp,
        'Cp experimental': cp_exp,
        }
df1 = DataFrame(data1, columns= ['Temperature', 'Cv Murnaghan', 'Cp Murnaghan', 'Cp experimental'])
#export_excel = df.to_excel (r'C:\\Users\\danie\\OneDrive\Documentos\\GitHub\\elasticity\\cp_data.xlsx', index = None, header=True) 
print (df1)

plt.figure(1)
plt.autoscale(tight=False)
plt.ticklabel_format(useOffset=False)
plt.xlabel("Temperatura (K)")
plt.ylabel("Calor específico [J/mol.K]")
plt.plot(temperatures_exp, qha1.cp, color = 'g', label = r"c$_P$ calculado")
plt.plot(temperatures_exp, cp_exp, color = 'r', label = r"c$_P$ experimental")
plt.tight_layout()
plt.legend(loc = 4)
plt.show()

print('Temperature: ', temperatures_exp, 
    '\n' 'Vmin: ', qha1.Vmin,
    '\n' 'F: ', qha1.F,
    '\n' 'Theta: ', qha1.Theta,
    '\n' 'B: ', qha1.B,
    '\n' 'Alpha: ', qha1.alpha,
    '\n' 'Cv Murnaghan: ', qha1.cv,
    '\n' 'Cp Murnaghan: ', qha1.cp)

#Data
#qha2 = thermal.QuasiHarmonicVector(T, chi, poisson, debye_temp, 
#                            energies, volumes, method = 'Murnaghan')

#data2 = {'Temperature': T, 
#        'Vmin': qha2.Vmin,
#        'F': qha2.F,
#        'Theta': qha2.theta,
#        'B': qha2.B,
#        'Alpha': qha2.alpha,
#        'Cv Murnaghan': qha2.cv,
#        'Cp Murnaghan': qha2.cp,
#        }
#df2 = DataFrame(data1, columns= ['Temperature', 'Vmin', 'F', 'Theta',
#                            'B', 'Alpha', 'Cv Murnaghan', 'Cp Murnaghan'])
#export_excel = df.to_excel (r'C:\\Users\\danie\\OneDrive\Documentos\\GitHub\\elasticity\\cp_data.xlsx', index = None, header=True) 
#print (df2)

print('time to run: ', time.time()-a)