import numpy as np
import matplotlib.pyplot as plt
import elasticity as el
import thermal_expansion as thermal

#Material properties  doi = "10.17188/1190959"
material = "silicon"
structure = "cubic"
density, density_unit = 2.28, "g/cm3"
atomic_mass, atomic_mass_unit = 28.0855, "u"
C11, C12, C44, elastic_constants_unit = 144, 53, 75, "GPa"
chi = .5

#Experimental data / ab initio
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
T_min, T_max = 10, 1000
T_amb = 300
T_unit = "K"

#VRH Theory
VRH = el.VRH(structure, C11, C12, C44)
VRH.VRH()

#EOS fit and plot
eos_murnaghan = thermal.EOS(energies, volumes, method = "Murnaghan")
eos_birch = thermal.EOS(energies, volumes, method = "Birch-Murnaghan")
debye_temp = thermal.Debye(VRH.ShearModulus, VRH.BulkModulus,density, atomic_mass)
print(debye_temp.theta0)

volume_fit = np.linspace(min(volumes), max(volumes), 100)
plt.figure(1)
plt.xlabel("Volume [" + volumes_unit + "]")
plt.ylabel("Energy [" + energies_unit + "]")
plt.plot(volumes, energies, 'bo', label = "Experimental data")
plt.plot(volume_fit, eos_murnaghan.murnaghan(volume_fit,*eos_murnaghan.fit_eos()), 
            color = 'g', label = "Murnaghan EOS")
plt.plot(volume_fit, eos_birch.birch_murnaghan(volume_fit,*eos_birch.fit_eos()), 
            color = 'r', label = "Birch-Murnaghan EOS")
plt.legend(loc = 1)
plt.show()

#Thermal expansion - Murnaghan equation
murnaghan_temperatures = []
murnaghan_volumes = []
murnaghan_energies = []
murnaghan_temp_debye = []
murnaghan_t = []
murnaghan_bulk = []
murnaghan_cv = []
murnaghan_cp = []
for T in range(T_min, T_max):
    murnaghan_temperatures.append(T)
    exp_murnaghan = thermal.QuasiHarmonic(T,VRH.ShearModulus,VRH.BulkModulus,density, 
                                            atomic_mass, chi, energies, volumes, 
                                            method="Murnaghan")
    #print(T, exp_murnaghan.V_min())
    murnaghan_volumes.append(exp_murnaghan.V_min())
    murnaghan_energies.append(exp_murnaghan.F(exp_murnaghan.V_min()))
    murnaghan_temp_debye.append(exp_murnaghan.theta_eq())
    murnaghan_t.append(T/exp_murnaghan.theta_eq())
    murnaghan_bulk.append(exp_murnaghan.bulk_modulus())
    murnaghan_cv.append(exp_murnaghan.heat_volume())
    murnaghan_cp.append(exp_murnaghan.heat_presure())


#Thermal expansion - Birch-Murnaghan equation
birch_temperatures = []
birch_volumes = []
birch_energies = []
birch_bulk = []
birch_cv = []
birch_cp = []
for T in range(T_min, T_max):
    birch_temperatures.append(T)
    exp_birch = thermal.QuasiHarmonic(T,VRH.ShearModulus,VRH.BulkModulus,density, 
                                            atomic_mass, chi, energies, volumes, 
                                            method="Birch-Murnaghan")
    birch_volumes.append(exp_birch.V_min())
    birch_energies.append(exp_birch.F(exp_birch.V_min()))
    birch_bulk.append(exp_birch.bulk_modulus())
    birch_cv.append(exp_birch.heat_volume())
    birch_cp.append(exp_birch.heat_presure())

#Plots
plt.figure(2)
plt.plot(murnaghan_volumes, murnaghan_energies, color = 'g', label = "Murnaghan")
plt.plot(birch_volumes, birch_energies, color = 'r', label = "Birch-Murnaghan")
plt.show()

plt.figure(3)
plt.plot(murnaghan_t, murnaghan_cv, color = 'b')
plt.plot(murnaghan_t, murnaghan_cp, color = 'g')
plt.show()

plt.figure(4)
plt.plot(murnaghan_volumes, murnaghan_bulk, color = 'g', label = "Murnaghan")
plt.plot(birch_volumes, birch_bulk, color = 'r', label = "Birch-Murnaghan")
plt.show()