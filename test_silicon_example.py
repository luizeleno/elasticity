import numpy as np
from pandas import DataFrame
import matplotlib.pyplot as plt
import elasticity as el
import thermal_expansion as thermal
from scipy import constants
R = constants.gas_constant

#Font size for the graphics
plt.rcParams.update({'font.size': 14})

#Material properties  doi = "10.17188/1190959"
material = "silicon"
structure = "cubic"
density, density_unit = 2.28, "g/cm3"
atomic_mass, atomic_mass_unit = 28.0855, "u"
C11, C12, C44, elastic_constants_unit = 162.07, 63.51, 77.26, "GPa"
chi = 1.

#Experimental data for specific heat at constant pressure
temperatures_exp = np.array([298.15, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 
                        1300, 1400, 1500, 1600, 1687])

cp_exp_J = np.array([20.007,  20.006,  22.258,  23.588,  24.420,  25.050,  25.608,  26.110, 
                    26.569, 26.988,  27.360, 27.707, 28.045,  28.372, 28.674, 28.930])

#Ab initio data
'''
    Calculated by the author.
'''
energies = np.array([-19.19110706, -19.18747301, -19.18192993, -19.19300108, -19.19331307,
                    -19.19219081, -19.18977257])
volumes = np.array([262.146575223, 254.284242613, 246.573276899, 270.169326182, 278.354115558, 
                    286.702563418, 295.224668962])

energies_unit = "Ry"
volumes_unit = "au3"  

#Temperature range
T_min, T_max = 10, 2000
T_amb = 300
T_unit = "K"

#VRH Theory
VRH = el.VRH(structure, C11, C12, C44)
VRH.VRH()
poisson = VRH.PoissonRatio

print(VRH.PoissonRatio, VRH.BulkModulus, VRH.ShearModulus, VRH.YoungModulus)


#EOS fit and plot
eos_murnaghan = thermal.EOS(energies, volumes, method = "Murnaghan")
eos_birch = thermal.EOS(energies, volumes, method = "Birch-Murnaghan")
debye_temp = thermal.Debye(VRH.ShearModulus, VRH.BulkModulus,density, atomic_mass)

print(debye_temp.temp_debye0())

print("Murnaghan fit", eos_murnaghan.fit_eos())
print("Birch fit", eos_birch.fit_eos())

volume_fit = np.linspace(min(volumes), max(volumes), 100)
plt.figure(1)
plt.autoscale(tight=False)
plt.ticklabel_format(useOffset=False)
plt.xlabel("Volume (" + volumes_unit + ")")
plt.ylabel("Energia (" + energies_unit + ")")
plt.plot(volumes, energies, 'ro', label = "Dados experimentais")
plt.plot(volume_fit, eos_murnaghan.murnaghan(volume_fit,*eos_murnaghan.fit_eos()), 
            color = 'g', label = "Murnaghan EOS")
plt.tight_layout()
plt.legend(loc = 9)
plt.savefig('eos_murnaghan_' + material + '.pdf')

plt.figure(2)
plt.autoscale(tight=False)
plt.ticklabel_format(useOffset=False)
plt.xlabel("Volume (" + volumes_unit + ")")
plt.ylabel("Energia (" + energies_unit + ")")
plt.plot(volumes, energies, 'ro', label = "Dados experimentais")
plt.plot(volume_fit, eos_birch.birch_murnaghan(volume_fit,*eos_birch.fit_eos()), 
            color = 'b', label = "Birch-Murnaghan EOS")
plt.tight_layout()
plt.legend(loc = 9)
plt.savefig('eos_birch_' + material + '.pdf')

#Experimental data adjustment
T_exp = []
t_exp_murnaghan = []
t_exp_birch = []
cp_exp = []
cp_simulated_murnaghan = []
cp_simulated_birch = []
for n in range(len(temperatures_exp)):
    T_exp.append(temperatures_exp[n])
    exp_murnaghan = thermal.QuasiHarmonic(temperatures_exp[n],VRH.ShearModulus,VRH.BulkModulus,density, 
                                            atomic_mass, chi, energies, volumes, 
                                            poisson, method="Murnaghan")
    exp_birch = thermal.QuasiHarmonic(temperatures_exp[n],VRH.ShearModulus,VRH.BulkModulus,density, 
                                            atomic_mass, chi, energies, volumes, 
                                            poisson, method="Birch-Murnaghan")
    t_exp_murnaghan.append(temperatures_exp[n]/exp_murnaghan.theta_eq())
    t_exp_birch.append(temperatures_exp[n]/exp_birch.theta_eq())
    cp_exp.append(cp_exp_J[n])
    cp_simulated_murnaghan.append(R*exp_murnaghan.heat_pressure())
    cp_simulated_birch.append(R*exp_birch.heat_pressure())

#Thermal expansion - Murnaghan equation
murnaghan_temperatures = []
murnaghan_volumes = []
murnaghan_energies = []
murnaghan_temp_debye = []
murnaghan_t = []
murnaghan_bulk = []
murnaghan_young = []
murnaghan_cv = []
murnaghan_cp = []
murnaghan_alpha = []
for T in range(T_min, T_max):
    murnaghan_temperatures.append(T)
    exp_murnaghan = thermal.QuasiHarmonic(T,VRH.ShearModulus,VRH.BulkModulus,density, 
                                            atomic_mass, chi, energies, volumes, 
                                            poisson, method="Murnaghan")
    #print(T, exp_murnaghan.V_min())
    murnaghan_volumes.append(exp_murnaghan.V_min())
    murnaghan_energies.append(exp_murnaghan.F(exp_murnaghan.V_min()))
    murnaghan_temp_debye.append(exp_murnaghan.theta_eq())
    murnaghan_t.append(T/exp_murnaghan.theta_eq())
    murnaghan_bulk.append(exp_murnaghan.bulk_modulus())
    murnaghan_cv.append(exp_murnaghan.heat_volume())
    murnaghan_cp.append(exp_murnaghan.heat_pressure())
    murnaghan_young.append(exp_murnaghan.young_modulus(poisson))
    murnaghan_alpha.append(exp_murnaghan.thermal_coef()/10**6)
    if T == 300:
        print("Murnaghan - 300K")
        print("Vmin", exp_murnaghan.V_min())
        print("Fmin", exp_murnaghan.F(exp_murnaghan.V_min()))
        print("Temp Debye", exp_murnaghan.theta_eq())
        print("B", exp_murnaghan.bulk_modulus())
        print("E", exp_murnaghan.young_modulus(poisson))
        print("cv", exp_murnaghan.heat_volume())
        print("cp", exp_murnaghan.heat_pressure())
        print("alpha", exp_murnaghan.thermal_coef())
    elif T==1000:
        print("Murnaghan - 1000K")
        print("Vmin", exp_murnaghan.V_min())
        print("Fmin", exp_murnaghan.F(exp_murnaghan.V_min()))
        print("Temp Debye", exp_murnaghan.theta_eq())
        print("B", exp_murnaghan.bulk_modulus())
        print("E", exp_murnaghan.young_modulus(poisson))
        print("cv", exp_murnaghan.heat_volume())
        print("cp", exp_murnaghan.heat_pressure())
        print("alpha", exp_murnaghan.thermal_coef())

#Thermal expansion - Birch-Murnaghan equation
birch_temperatures = []
birch_volumes = []
birch_energies = []
birch_temp_debye = []
birch_t = []
birch_bulk = []
birch_young = []
birch_cv = []
birch_cp = []
birch_alpha = []
for T in range(T_min, T_max):
    birch_temperatures.append(T)
    exp_birch = thermal.QuasiHarmonic(T,VRH.ShearModulus,VRH.BulkModulus,density, 
                                            atomic_mass, chi, energies, volumes, 
                                            poisson, method="Birch-Murnaghan")
    #print(T, exp_murnaghan.V_min())
    birch_volumes.append(exp_birch.V_min())
    birch_energies.append(exp_birch.F(exp_murnaghan.V_min()))
    birch_temp_debye.append(exp_birch.theta_eq())
    birch_t.append(T/exp_birch.theta_eq())
    birch_bulk.append(exp_birch.bulk_modulus())
    birch_cv.append(exp_birch.heat_volume())
    birch_cp.append(exp_birch.heat_pressure())
    birch_young.append(exp_birch.young_modulus(poisson))
    birch_alpha.append(exp_birch.thermal_coef()/10**6)
    if T == 300:
        print("Birch-Murnaghan - 300K")
        print("Vmin", exp_birch.V_min())
        print("Fmin", exp_birch.F(exp_birch.V_min()))
        print("Temp Debye", exp_birch.theta_eq())
        print("B", exp_birch.bulk_modulus())
        print("E", exp_birch.young_modulus(poisson))
        print("cv", exp_birch.heat_volume())
        print("cp", exp_birch.heat_pressure())
        print("alpha", exp_birch.thermal_coef())
    elif T==1000:
        print("Birch-Murnaghan - 1000K")
        print("Vmin", exp_birch.V_min())
        print("Fmin", exp_birch.F(exp_birch.V_min()))
        print("Temp Debye", exp_birch.theta_eq())
        print("B", exp_birch.bulk_modulus())
        print("E", exp_birch.young_modulus(poisson))
        print("cv", exp_birch.heat_volume())
        print("cp", exp_birch.heat_pressure())
        print("alpha", exp_birch.thermal_coef())

#Plots
plt.figure(3)
plt.autoscale(tight=False)
plt.ticklabel_format(useOffset=False)
plt.xlabel("Temperatura (K)")
plt.ylabel("B (GPa)")
plt.plot(murnaghan_temperatures, murnaghan_bulk, color = 'g', label = "Murnaghan")
plt.plot(birch_temperatures, birch_bulk, color = 'b', label = "Birch-Murnaghan")
plt.tight_layout()
plt.legend(loc = 1)
plt.savefig('bulk_modulus_' + material + '.pdf')

plt.figure(4)
plt.autoscale(tight=False)
plt.ticklabel_format(useOffset=False)
plt.xlabel("Temperatura (K)")
plt.ylabel("E (GPa)")
plt.plot(murnaghan_temperatures, murnaghan_young, color = 'g', label = "Murnaghan")
plt.plot(birch_temperatures, birch_young, color = 'b', label = "Birch-Murnaghan")
plt.tight_layout()
plt.legend(loc = 1)
plt.savefig('young_modulus_' + material + '.pdf')

plt.figure(5)
plt.autoscale(tight=False)
plt.ticklabel_format(useOffset=False)
plt.xlabel("Temperatura (K)")
plt.ylabel(r"Coef. expansão térmica (10$^{-6}$K$^{-1}$)")
plt.plot(murnaghan_temperatures, murnaghan_alpha, color = 'g', label = "Murnaghan")
plt.plot(birch_temperatures, birch_alpha, color = 'b', label = "Birch-Murnaghan")
plt.tight_layout()
plt.legend(loc = 4)
plt.savefig('coef_termico_' + material + '.pdf')

plt.figure(6)
plt.autoscale(tight=False)
plt.ticklabel_format(useOffset=False)
plt.xlabel("Temperatura reduzida")
plt.ylabel("Calor específico (1/R)")
plt.plot(murnaghan_t, murnaghan_cv, color = 'b', label = r"c$_V$")
plt.plot(murnaghan_t, murnaghan_cp, color = 'g', label = r"c$_P$")
#plt.plot(T_exp, cp_exp, color = 'r', label = r"c$_P$ experimental")
plt.tight_layout()
plt.legend(loc = 4)
plt.savefig('specific_heat_murnaghan_' + material + '.pdf')

plt.figure(7)
plt.autoscale(tight=False)
plt.ticklabel_format(useOffset=False)
plt.xlabel("Temperatura reduzida")
plt.ylabel("Calor específico (1/R)")
plt.plot(birch_t, birch_cv, color = 'b', label = r"c$_V$")
plt.plot(birch_t, birch_cp, color = 'g', label = r"c$_P$")
#plt.plot(T_exp, cp_exp, color = 'r', label = r"c$_P$ experimental")
plt.tight_layout()
plt.legend(loc = 4)
plt.savefig('specific_heat_birch_' + material + '.pdf')

plt.figure(8)
plt.autoscale(tight=False)
plt.ticklabel_format(useOffset=False)
plt.xlabel("Temperatura (K)")
plt.ylabel("Temperatura de Debye (K)")
plt.plot(murnaghan_temperatures, murnaghan_temp_debye, color = 'g', label = "Murnaghan")
plt.plot(birch_temperatures, birch_temp_debye, color = 'b', label = "Birch-Murnaghan")
plt.tight_layout()
plt.legend(loc = 1)
plt.savefig('debye_temperature_' + material + '.pdf')

plt.figure(9)
plt.autoscale(tight=False)
plt.ticklabel_format(useOffset=False)
plt.xlabel("Temperatura (K)")
plt.ylabel("Calor específico (J/mol.K)")
plt.plot(T_exp, cp_simulated_murnaghan, color = 'g', label = r"c$_P$ calculado")
plt.plot(T_exp, cp_exp, color = 'r', label = r"c$_P$ experimental")
plt.tight_layout()
plt.legend(loc = 4)
plt.savefig('comparative_specific_heat_murnaghan_' + material + '.pdf')

plt.figure(10)
plt.autoscale(tight=False)
plt.ticklabel_format(useOffset=False)
plt.xlabel("Temperatura (K)")
plt.ylabel("Calor específico [J/mol.K]")
plt.plot(T_exp, cp_simulated_birch, color = 'g', label = r"c$_P$ calculado")
plt.plot(T_exp, cp_exp, color = 'r', label = r"c$_P$ experimental")
plt.tight_layout()
plt.legend(loc = 4)
plt.savefig('comparative_specific_heat_birch_' + material + '.pdf')

#Export data for Excel
data = {'Temperature': T_exp, 
        'Cp experimental': cp_exp,
        'Cp Murnaghan': cp_simulated_murnaghan,
        'Cp Birch': cp_simulated_birch
        }
df = DataFrame(data, columns= ['Temperature', 'Cp experimental', 'Cp Murnaghan', 'Cp Birch'])
#export_excel = df.to_excel (r'C:\Users\danie\OneDrive\Documentos\GitHub\elasticity\cp_data.xlsx', index = None, header=True) 
print (df)