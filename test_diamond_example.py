import numpy as np
from pandas import DataFrame
import matplotlib.pyplot as plt
import elasticity as el
import thermal_expansion as thermal
from scipy import constants
R = constants.gas_constant

#Font size for the graphics
plt.rcParams.update({'font.size': 14})

print(constants.physical_constants["atomic mass constant"][0])
#Material properties  doi = "10.17188/1281384"
material = "diamond"
structure = "cubic"
density, density_unit = 3.5, "g/cm3"
atomic_mass, atomic_mass_unit = 12.0107, "u"
C11, C12, C44, elastic_constants_unit = 1079, 124, 578, "GPa"
chi = 1.

#Experimental data for specific heat at constant pressure
temperatures_exp = np.array([70.16, 75.37, 81.59, 88.65, 96.68, 105.1, 113.04, 125.28, 134.29, 144.1, 153.71, 
                        162.76, 173.33, 181.96, 191.44, 200.94, 211.84, 231.06, 241.09, 252.37, 264.31, 
                        276.61, 287.96, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 550, 600,
                        650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100])

cp_exp_J = np.array([0.0938, 0.124, 0.152, 0.189, 0.248, 0.329, 0.415, 0.578, 0.729, 0.913, 1.118, 1.331,
                    1.616, 1.859, 2.169, 2.491, 2.868, 3.592, 3.843, 4.325, 4.760, 5.263, 5.673, 6.196, 
                    7.030, 7.863, 8.679, 9.475, 10.24, 10.97, 11.67, 12.34, 12.97, 13.57, 14.94, 16.13,
                    17.16, 18.05, 18.83, 19.51, 20.11, 20.65, 21.14, 21.61, 22.06, 22.52])

#Ab initio data
'''
    Calculated by the author.
'''
energies = np.array([-22.76745769, -22.76213240, -22.75431180, -22.77047178, -22.77134749,
                    -22.77024716, -22.76732281])
volumes = np.array([72.2765903507, 70.1081472678, 67.9835155729, 74.4892919379, 76.7466991459, 
                    79.0492590911, 81.39741889])

energies_unit = "Ry"
volumes_unit = "au3"  

#Temperature range
T_min, T_max = 10, 4000
T_amb = 300
T_unit = "K"

#VRH Theory
VRH = el.VRH(structure, C11, C12, C44)
VRH.VRH()
poisson = VRH.PoissonRatio

#EOS fit and plot
eos_murnaghan = thermal.EOS(energies, volumes, method = "Murnaghan")
eos_birch = thermal.EOS(energies, volumes, method = "Birch-Murnaghan")
debye_temp = thermal.Debye(VRH.ShearModulus, VRH.BulkModulus,density, atomic_mass)

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