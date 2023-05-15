import numpy as np
import numpy.linalg as la


class ElasticityTheory:

    '''
        Generates the compliance matrix, compliance tensor and stiffness matrix based on the given inputs. Each crystal class have a specific number of constants to be given

        Parameters
        ----------
        crystal_class: str
            one of the following crystal classes: isotropic (2), cubic (3), hexagonal (5), trigonal_1 (6), trigonal_2 (7), tetragonal_1 (6), tetragonal_2 (7), orthorhombic (9), monoclinic_1 (13), monoclinic_2 (13), triclinic (21)
        stiffness constants: int
            independent elastic constants, the quantity of constants depends on the crystal class, as indicated between parenthesis above, GPa is the refered unit

        Returns
        ----------
        C: array
            (6x6) numpy array, stiffness matrix based on the given constants and material crystal class
        S: array
            (6x6) numpy array, compliance matrix based on the given constants and material crystal class
        St: array
            (3x3x3x3) numpy array, compliance tensor based on the given constants and material crystal class



    '''

    def __init__(self, crystal_class, *stiff_consts):

        self.crystal_classes = {'isotropic': 2, 'cubic': 3, 'hexagonal': 5,
                                'trigonal_1': 6, 'trigonal_2': 7,
                                'tetragonal_1': 6, 'tetragonal_2': 7,
                                'orthorhombic': 9, 'monoclinic_1': 13,
                                'monoclinic_2': 13, 'triclinic': 21}

        if crystal_class not in self.crystal_classes:
            raise KeyError('Unknown crystal class: %s. Use: %s' % (crystal_class, [s for s in self.crystal_classes]))

        self.crystal_class = crystal_class
        self.nconsts = self.crystal_classes[crystal_class]

        self.C = self.StiffnessMatrix(crystal_class, *stiff_consts)
        self.S = self.ComplianceMatrix(self.C)
        self.St = self.ComplianceTensor(self.S)

    def StiffnessMatrix(self, crystal_class, *stiff_consts):
        '''
        Generates the stiffness matrix based on the given crystal class and independent elastic constants based on voigt notation

        Parameters
        ----------
        crystal_class: str
            one of the following crystal classes: isotropic (2), cubic (3), hexagonal (5), trigonal_1 (6), trigonal_2 (7), tetragonal_1 (6), tetragonal_2 (7), orthorhombic (9), monoclinic_1 (13), monoclinic_2 (13), triclinic (21)
        stiffness constants: int
            independent elastic constants, the quantity of constants depends on the crystal class, as indicated between parenthesis above, GPa is the refered unit

        Returns
        ----------
        C: array
            (6x6) numpy array, stiffness matrix based on the given constants and material crystal class

        '''

        m, n = self.crystal_classes[crystal_class], len(stiff_consts)

        if m != n:
            raise IndexError('Crystal class %s has %d independent constants, but %d were provided' % (crystal_class, m, n))

        if crystal_class == 'isotropic':
            C11, C12 = stiff_consts
            C44 = .5 * (C11 - C12)
            C = np.array([
                [C11, C12, C12,  0.,  0.,  0.],
                [0., C11, C12,  0.,  0.,  0.],
                [0.,  0., C11,  0.,  0.,  0.],
                [0.,  0.,  0., C44,  0.,  0.],
                [0.,  0.,  0.,  0., C44,  0.],
                [0.,  0.,  0.,  0.,  0., C44]
               ]
              )
        elif crystal_class == 'cubic':      # all classes
            C11, C12, C44 = stiff_consts
            C = np.array([
                [C11, C12, C12,  0.,   0.,  0.],
                [0., C11, C12,  0.,   0.,  0.],
                [0.,  0., C11,  0.,   0.,  0.],
                [0.,  0.,  0., C44,   0.,  0.],
                [0.,  0.,  0.,  0.,  C44,  0.],
                [0.,  0.,  0.,  0.,   0., C44]
               ]
              )
        elif crystal_class == 'hexagonal':      # all classes
            C11, C12, C13, C33, C44 = stiff_consts
            C66 = .5 * (C11 - C12)
            C = np.array([
                [C11, C12, C13,  0.,  0.,  0.],
                [0., C11, C13,  0.,  0.,  0.],
                [0.,  0., C33,  0.,  0.,  0.],
                [0.,  0.,  0., C44,  0.,  0.],
                [0.,  0.,  0.,  0., C44,  0.],
                [0.,  0.,  0.,  0.,  0., C66]
               ]
              )
        elif crystal_class == 'tetragonal_1':  # classes 4mm, -42m, 422, 4/mmm
            C11, C12, C13, C33, C44, C66 = stiff_consts
            C = np.array([
                [C11, C12, C13,  0.,  0.,  0.],
                [0., C11, C13,  0.,  0.,  0.],
                [0.,  0., C33,  0.,  0.,  0.],
                [0.,  0.,  0., C44,  0.,  0.],
                [0.,  0.,  0.,  0., C44,  0.],
                [0.,  0.,  0.,  0.,  0., C66]
               ]
              )
        elif crystal_class == 'tetragonal_2':    # classes 4, -4, 4/m
            C11, C12, C13, C33, C44, C66, C16 = stiff_consts
            C = np.array([
                [C11, C12, C13,  0.,  0.,  C16],
                [0., C11, C13,  0.,  0., -C16],
                [0.,  0., C33,  0.,  0.,   0.],
                [0.,  0.,  0., C44,  0.,   0.],
                [0.,  0.,  0.,  0., C44,   0.],
                [0.,  0.,  0.,  0.,  0.,  C66]
               ]
              )
        elif crystal_class == 'trigonal_1':     # classes 32, -3m, 3m
            C11, C12, C13, C14, C33, C44 = stiff_consts
            C66 = .5 * (C11 - C12)
            C = np.array([
                [C11, C12, C13,  C14,  0.,  0.],
                [0., C11, C13, -C14,   0,  0.],
                [0., 0.,  C33,   0.,  0.,  0.],
                [0., 0.,   0.,  C44,  0.,  0.],
                [0., 0.,   0.,   0., C44, C14],
                [0., 0.,   0.,   0.,  0., C66]
               ]
              )
        elif crystal_class == 'trigonal_2':     # classes 3, -3
            C11, C12, C13, C33, C14, C15, C44 = stiff_consts
            C66 = .5 * (C11 - C12)
            C = np.array([
                [C11, C12, C13,  C14,  C15,   0.],
                [0., C11, C13, -C14, -C15,   0.],
                [0.,  0., C33,   0.,   0.,   0.],
                [0.,  0.,  0.,  C44,   0., -C15],
                [0.,  0.,  0.,   0.,  C44,  C14],
                [0.,  0.,  0.,   0.,   0.,  C66]
               ]
              )
        elif crystal_class == 'orthorhombic':       # all classes
            C11, C12, C13, C22, C23, C33, C44, C55, C66 = stiff_consts
            C = np.array([
                [C11, C12, C13,  0.,  0.,  0.],
                [0., C22, C23,  0.,  0.,  0.],
                [0.,  0., C33,  0.,  0.,  0.],
                [0.,  0.,  0., C44,  0.,  0.],
                [0.,  0.,  0.,  0., C55,  0.],
                [0.,  0.,  0.,  0.,  0., C66]
               ]
              )
        elif crystal_class == 'monoclinic_1':       # diad || x2
            C11, C12, C13, C22, C23, C33, C15, C25, C35, C44, C46, C55, C66 = stiff_consts
            C = np.array([
                [C11, C12, C13,  0., C15,  0.],
                [0., C22, C23,  0., C25,  0.],
                [0.,  0., C33,  0., C35,  0.],
                [0.,  0.,  0., C44,  0., C46],
                [0.,  0.,  0.,  0., C55,  0.],
                [0.,  0.,  0.,  0.,  0., C66]
               ]
              )
        elif crystal_class == 'monoclinic_2':       # diad || x3
            C11, C12, C13, C22, C23, C33, C16, C26, C36, C44, C45, C55, C66 = stiff_consts
            C = np.array([
                [C11, C12, C13,  0.,  0., C16],
                [0., C22, C23,  0.,  0., C26],
                [0.,  0., C33,  0.,  0., C36],
                [0.,  0.,  0., C44, C45,  0.],
                [0.,  0.,  0.,  0., C55,  0.],
                [0.,  0.,  0.,  0.,  0., C66]
               ]
              )
        elif crystal_class == 'triclinic':      # all classes
            C11, C12, C13, C14, C15, C16, C22, C23, C24, C25, C26, C33, C34, C35, C36, C44, C45, C46, C55, C56, C66 = stiff_consts
            C = np.array([
                [C11, C12, C13, C14, C15, C16],
                [0., C22, C23, C24, C25, C26],
                [0.,  0., C33, C34, C35, C36],
                [0.,  0.,  0., C44, C45, C46],
                [0.,  0.,  0.,  0., C55, C56],
                [0.,  0.,  0.,  0.,  0., C66]
               ]
              )

        C = C + C.T - np.diag(C.diagonal())     # symmetrize C
        return C

    def ComplianceMatrix(self, C):
        '''
        Generates the compliance matrix based on the given crystal class and independent elastic constants based on voigt notation

        Parameters
        ----------
        crystal_class: str
            one of the following crystal classes: isotropic (2), cubic (3), hexagonal (5), trigonal_1 (6), trigonal_2 (7), tetragonal_1 (6), tetragonal_2 (7), orthorhombic (9), monoclinic_1 (13), monoclinic_2 (13), triclinic (21)
        stiffness constants: int
            independent elastic constants, the quantity of constants depends on the crystal class, as indicated between parenthesis above, GPa is the refered unit

        Returns
        ----------
        S: array
            (6x6) numpy array, compliance matrix based on the given constants and material crystal class
        '''

        return la.inv(C)

    def ComplianceTensor(self, S):
        '''
        Generates the compliance tensor based on the given crystal class and independent elastic constants based on voigt notation

        Parameters
        ----------
        crystal_class: str
            one of the following crystal classes: isotropic (2), cubic (3), hexagonal (5), trigonal_1 (6), trigonal_2 (7), tetragonal_1 (6), tetragonal_2 (7), orthorhombic (9), monoclinic_1 (13), monoclinic_2 (13), triclinic (21)
        stiffness constants: int
            independent elastic constants, the quantity of constants depends on the crystal class, as indicated between parenthesis above, GPa is the refered unit


        Returns
        ----------
        St: array
            (3x3x3x3) numpy array, compliance tensor based on the given constants and material crystal class
        '''


        St = np.zeros((3, 3, 3, 3), dtype=float)

        rel = {'00': 0, '11': 1, '22': 2,
               '12': 3, '21': 3, '02': 4, '20': 4,
               '01': 5, '10': 5}

        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        m_str = str(i) + str(j)
                        n_str = str(k) + str(l)
                        m = rel[m_str]
                        n = rel[n_str]
                        if (m < 3 and n < 3):
                            St[i, j, k, l] = S[m, n]
                        elif (m >= 3 and n >= 3):
                            St[i, j, k, l] = .25 * S[m, n]
                        else:
                            St[i, j, k, l] = .5 * S[m, n]

        return St


class DirectionalProperties(ElasticityTheory):

    '''
        Generates the object that calculate directional properties based on the given input properties

        Parameters
        ----------
        crystal_class: str
            one of the following crystal classes: isotropic (2), cubic (3), hexagonal (5), trigonal_1 (6), trigonal_2 (7), tetragonal_1 (6), tetragonal_2 (7), orthorhombic (9), monoclinic_1 (13), monoclinic_2 (13), triclinic (21)
        stiffness constants: int
            independent elastic constants, the quantity of constants depends on the crystal class, as indicated between parenthesis above, GPa is the refered unit


        Returns
        ----------
        B:  array
            (NxN) np array, Linear compressibility for each direction l
        E: array
            (NxN) np array, Young Modulus for each direction l
        G: array
            (NxNxN) np array | Shear Modulus for l and n directions
        p: array
            (NxNxN) np array | Poisson's Ration for l and n directions
        *N is the number of points specified in the l direction

    '''

    def LinearCompressibility(self, l):
        '''
            Generates the material linear compressibility matrix based on the l direction input

            Parameters
            ----------
            l: array
                Stress direction (NxN) array

            Returns
            ----------
            B: array
                (NxN) np array, Linear compressibility for each direction l
        '''

        B = 0.
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    B += self.St[i, j, k, k] * l[i] * l[j]
        return 1.0 / B

    def YoungModulus(self, l):
        '''
            Generates the material young modulus matrix based on the l direction input

            Parameters
            ----------
            l: array
                Stress direction (NxN) array

            Returns
            ----------
            E: array
                (NxN) np array, Linear compressibility for each direction l
        '''
        E = 0.
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for m in range(3):
                        E += self.St[i, j, k, m] * l[i] * l[j] * l[k] * l[m]
        return 1.0 / E

    def ShearModulus(self, l, n):
        '''
            Generates the material shear modulus matrix based on the l direction input

            Parameters
            ----------
            l: array
                Stress direction (NxN) array
            n: array
                Normal direction (NxNxN) array

            Returns
            ----------
            G: array
                (NxNxN) np array, Linear compressibility for each direction l and normal direction n

        '''

        G = 0.
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for m in range(3):
                        G += self.St[i, j, k, m] * l[i] * n[j] * l[k] * n[m]
        return .25 / G

    def PoissonRatio(self, l, n):
        '''
            Generates the material Poisson Ratio matrix based on the l direction input

            Parameters
            ----------
            l: array
                Stress direction (NxN) array
            n: array
                Normal direction (NxNxN) array

            Returns
            ----------
            PoissonRatio: array
                (NxNxN) np array, Linear compressibility for each direction l and normal direction n

        '''

        p = 0.
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for m in range(3):
                        p += self.St[i, j, k, m] * n[i] * n[j] * l[k] * l[m]
        return -p * self.YoungModulus(l)


class VRH(ElasticityTheory):

    '''
        Voigt-Reuss-Hill approximation to polycristalline elastic properties

        Parameters
        ----------
        crystal_class: str
            one of the following crystal classes: isotropic (2), cubic (3), hexagonal (5), trigonal_1 (6), trigonal_2 (7), tetragonal_1 (6), tetragonal_2 (7), orthorhombic (9), monoclinic_1 (13), monoclinic_2 (13), triclinic (21)
        stiffness constants: int
            independent elastic constants, the quantity of constants depends on the crystal class, as indicated between parenthesis above, GPa is the refered unit


        Returns
        ----------
        BR: int
            Reuss approximation of linear compressibility
        GR: int
            Reuss approximation of shear modulus
        BV: int
            Voigt approximation of linear compressibility
        GV :int
            Voigt approximation of shear modulus
        BH: int
            Hill approximation of linear compressibility
        GH: int
            Hill Reuss approximation of shear modulus
        BulkModulus: int
            VRH approximation of Bulk Modulus using data from the three approximations
        YoungModulus: int
            VRH approximation of Young Modulus using data from the three approximations
        PoissonRatio: int
            VRH approximation of Poisson Ratio using data from the three approximations

    '''

    def Reuss(self):

        sii = self.S[0, 0] + self.S[1, 1] + self.S[2, 2]
        sij = self.S[0, 1] + self.S[0, 2] + self.S[1, 2]
        skl = self.S[3, 3] + self.S[4, 4] + self.S[5, 5]
        self.BR = 1. / (sii + 2. * sij)
        self.GR = 15. / (4. * sii - sij + 3. * skl)

    def Voigt(self):

        cii = self.C[0, 0] + self.C[1, 1] + self.C[2, 2]
        cij = self.C[0, 1] + self.C[0, 2] + self.C[1, 2]
        ckl = self.C[3, 3] + self.C[4, 4] + self.C[5, 5]
        self.BV = (cii + 2. * cij) / 9.
        self.GV = (cii - cij + 3. * ckl) / 15.

    def Hill(self):

        self.BH = .5 * (self.BR + self.BV)
        self.GH = .5 * (self.GR + self.GV)

    def VRH(self):

        self.Reuss()
        self.Voigt()
        self.Hill()

        # Other elastic properties
        B, G = self.BH, self.GH
        self.BulkModulus, self.ShearModulus = B, G
        self.YoungModulus = 9. * B * G / (3. * B + G)
        self.PoissonRatio = (3. * B - 2. * G) / 2. / (3. * B + G)


class DebyeGruneisen(VRH):
    '''
    Class under development
    '''
    pass


class Elasticity(DirectionalProperties, VRH):
    '''
    Class under development
    '''
    pass
