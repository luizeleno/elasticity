import numpy as np
import elasticity as el
import thermal_expansion as thermal

energies = np.array([-76.03124099, -76.03638913, -76.04013485, -76.04259556, -76.04385612,
                    -76.04401924, -76.04315491, -76.04135395, -76.03868801, -76.03522839, -76.03103839])
volumes = np.array([65.64719243, 67.74216001, 69.88123056, 72.06486349, 74.2935182, 
                    76.56765409, 78.88773058, 81.25420706, 83.66754295, 86.12819765, 88.63663057])

material_name = "diamond"
structure = "cubic"
density = 3.50
density_unit = "g/cm3"
elastic_constants = (1054, 125, 562)
C11 = 1054
C12 = 125
C44 = 562
elastic_constants_unit = "GPa"
molar_volume = 3.42
molar_volume_unit = "cm3"
T = 300
chi = .5

VRH = el.VRH(structure, C11, C12, C44)

VRH.VRH()

#eos_murnaghan = thermal.EOS(energies, volumes)

#debye = thermal.Debye(VRH.ShearModulus,VRH.BulkModulus,density, molar_volume)

expansion_murnaghan = thermal.QuasiHarmonic(T,VRH.ShearModulus,VRH.BulkModulus,density, 
                                            molar_volume, chi, energies, volumes, 
                                            method="Murnaghan")

expansion_birch = thermal.QuasiHarmonic(T,VRH.ShearModulus,VRH.BulkModulus,density, 
                                        molar_volume, chi, energies, volumes, 
                                        method="Birch-Murnaghan")

print(expansion_murnaghan.V_min())
print(expansion_birch.V_min())

