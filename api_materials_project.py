import pymatgen as mg

type_element = input("What do you want to simulate, a compound (1) or a element (2):" )

print(type_element)
if type_element == "1":
    compound = mg.Composition(input("What is the compound? Example: Fe2O3."))
elif type_element == "2":
    element = mg.Element(input("What is the element? Example: Fe."))
else:
    pass


class APIMaterialProject():

    def __init__(self, density, atomic_mass, n_atoms, V0, energies, volumes, structure):
        self.density = density
        self.atomic_mass = atomic_mass
        self.n_atoms = n_atoms
        self.V0 = V0
        self.energies = energies
        self.volumes = volumes
        self.structure = structure

    def Composition(self, type_, composition):
        if type_ == compound:
            self.composition = mg.Composition(composition)
        elif type_ == element:
            self.composition = mg.Element(composition)
        else:
            pass
        return self.composition

    def IntrinsicProperties(self):
        self.density = self.composition.density
        self.atomic_mass = self.composition.atomic_mass
        self.n_atoms = self.composition.num_atoms
         
        self.structure = self.composition.structure


