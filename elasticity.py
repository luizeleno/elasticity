import numpy as np
import numpy.linalg as la

class ElasticityTheory:

    '''
        ElasticityTheory class

        input:
            - crystal_class: classified according to Nye (1956).
            One of the following (with no. of independent elastic constants):
                   'isotropic' (2) , 'cubic' (3), 'hexagonal' (5),
                   'trigonal_1' (6), 'trigonal_2' (7), 'tetragonal_1' (6),
                   'tetragonal_2' (7), 'orthorhombic' (9), 'monoclinic_1' (13),
                   'monoclinic_2' (13), 'triclinic' (21)
            - stiffness constants: variable number of constants,
               depending on crystal_class.
               Units: GPa

        output: functions to calculate the Compliance and Stiffness matrices,
           and the Compliance tensor

    '''

    def StiffnessMatrix(self, crystal_class, *stiff_consts):
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
                [0.,  0.,  0.,  0.,  0., C44]])
        elif crystal_class == 'cubic':      # all classes
            C11, C12, C44 = stiff_consts
            C = np.array([
                [C11, C12, C12,  0.,   0.,  0.],
                [0., C11, C12,  0.,   0.,  0.],
                [0.,  0., C11,  0.,   0.,  0.],
                [0.,  0.,  0., C44,   0.,  0.],
                [0.,  0.,  0.,  0.,  C44,  0.],
                [0.,  0.,  0.,  0.,   0., C44]])
        elif crystal_class == 'hexagonal':      # all classes
            C11, C12, C13, C33, C44 = stiff_consts
            C66 = .5 * (C11 - C12)
            C = np.array([
                [C11, C12, C13,  0.,  0.,  0.],
                [0., C11, C13,  0.,  0.,  0.],
                [0.,  0., C33,  0.,  0.,  0.],
                [0.,  0.,  0., C44,  0.,  0.],
                [0.,  0.,  0.,  0., C44,  0.],
                [0.,  0.,  0.,  0.,  0., C66]])
        elif crystal_class == 'tetragonal_1':  # classes 4mm, -42m, 422, 4/mmm
            C11, C12, C13, C33, C44, C66 = stiff_consts
            C = np.array([
                [C11, C12, C13,  0.,  0.,  0.],
                [0., C11, C13,  0.,  0.,  0.],
                [0.,  0., C33,  0.,  0.,  0.],
                [0.,  0.,  0., C44,  0.,  0.],
                [0.,  0.,  0.,  0., C44,  0.],
                [0.,  0.,  0.,  0.,  0., C66]])
        elif crystal_class == 'tetragonal_2':    # classes 4, -4, 4/m
            C11, C12, C13, C16, C33, C44, C66 = stiff_consts
            C = np.array([
                [C11, C12, C13,  0.,  0.,  C16],
                [0., C11, C13,  0.,  0., -C16],
                [0.,  0., C33,  0.,  0.,   0.],
                [0.,  0.,  0., C44,  0.,   0.],
                [0.,  0.,  0.,  0., C44,   0.],
                [0.,  0.,  0.,  0.,  0.,  C66]])
        elif crystal_class == 'trigonal_1':     # classes 32, -3m, 3m
            C11, C12, C13, C14, C33, C44 = stiff_consts
            C66 = .5 * (C11 - C12)
            C = np.array([
                [C11, C12, C13,  C14,  0.,  0.],
                [0., C11, C13, -C14,   0,  0.],
                [0., 0.,  C33,   0.,  0.,  0.],
                [0., 0.,   0.,  C44,  0.,  0.],
                [0., 0.,   0.,   0., C44, C14],
                [0., 0.,   0.,   0.,  0., C66]])
        elif crystal_class == 'trigonal_2':     # classes 3, -3
            C11, C12, C13, C33, C14, C15, C44 = stiff_consts
            C66 = .5 * (C11 - C12)
            C = np.array([
                [C11, C12, C13,  C14,  C15,   0.],
                [0., C11, C13, -C14, -C15,   0.],
                [0.,  0., C33,   0.,   0.,   0.],
                [0.,  0.,  0.,  C44,   0., -C15],
                [0.,  0.,  0.,   0.,  C44,  C14],
                [0.,  0.,  0.,   0.,   0.,  C66]])
        elif crystal_class == 'orthorhombic':       # all classes
            C11, C12, C13, C22, C23, C33, C44, C55, C66 = stiff_consts
            C = np.array([
                [C11, C12, C13,  0.,  0.,  0.],
                [0., C22, C23,  0.,  0.,  0.],
                [0.,  0., C33,  0.,  0.,  0.],
                [0.,  0.,  0., C44,  0.,  0.],
                [0.,  0.,  0.,  0., C55,  0.],
                [0.,  0.,  0.,  0.,  0., C66] ])
        elif crystal_class == 'monoclinic_1':       # diad || x2
            C11, C12, C13, C15, C22, C23, C25, C33, C35, C44, C46, C55, C66 = stiff_consts
            C = np.array([
                [C11, C12, C13,  0., C15,  0.],
                [0., C22, C23,  0., C25,  0.],
                [0.,  0., C33,  0., C35,  0.],
                [0.,  0.,  0., C44,  0., C46],
                [0.,  0.,  0.,  0., C55,  0.],
                [0.,  0.,  0.,  0.,  0., C66]])
        elif crystal_class == 'monoclinic_2':       # diad || x3
            C11, C12, C13, C16, C22, C23, C26, C33, C36, C44, C45, C55, C66 = stiff_consts
            C = np.array([
                [C11, C12, C13,  0.,  0., C16],
                [0., C22, C23,  0.,  0., C26],
                [0.,  0., C33,  0.,  0., C36],
                [0.,  0.,  0., C44, C45,  0.],
                [0.,  0.,  0.,  0., C55,  0.],
                [0.,  0.,  0.,  0.,  0., C66]])
        elif crystal_class == 'triclinic':      # all classes
            C11, C12, C13, C14, C15, C16, C22, C23, C24, C25, C26, C33, C34, C35, C36, C44, C45, C46, C55, C56, C66 = stiff_consts
            C = np.array([
                [C11, C12, C13, C14, C15, C16],
                [0., C22, C23, C24, C25, C26],
                [0.,  0., C33, C34, C35, C36],
                [0.,  0.,  0., C44, C45, C46],
                [0.,  0.,  0.,  0., C55, C56],
                [0.,  0.,  0.,  0.,  0., C66]])

        C = C + C.T - np.diag(C.diagonal())     # symmetrize C
        return C

    def ComplianceMatrix(self, C):

        return la.inv(C)

    def ComplianceTensor(self, S):

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
        Directional properties --- functions to calculate:
            - unidimensional bulk modulus(l)
            - Young's modulus(l)
            - Shear Modulus(l, n)
            #- Poisson's ratio(l, n)

        as a function of strain/stress direction l and normal direction n
    '''

    def BulkModulus(self, l):

        B = 0.
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    B += self.St[i, j, k, k] * l[i] * l[j]
        return 1.0 / B

    def YoungModulus(self, l):

        E = 0.
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for m in range(3):
                        E += self.St[i, j, k, m] * l[i] * l[j] * l[k] * l[m]
        return 1.0 / E

    def ShearModulus(self, l, n):

        G = 0.
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for m in range(3):
                        G += self.St[i, j, k, m] * l[i] * n[j] * l[k] * n[m]
        return .25 / G

    def PoissonRatio(self, l, n):

        p = 0.
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for m in range(3):
                        p += self.St[i, j, k, m] * n[i] * n[j] * l[k] * l[m]
        return -p * self.YoungModulus(l)


class VRH(ElasticityTheory):

    '''
        Voigt-Reuss-Hill approximation
        to polycristalline elastic properties

        The function VHR is the main user interface
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

    pass


class Elasticity(DirectionalProperties, VRH):

    pass

class Crystal(Elasticity):
    def __init__(self):
                self.crystal_classes = {'isotropic': 2, 'cubic': 3, 'hexagonal': 5,
                                        'trigonal_1': 6, 'trigonal_2': 7,
                                        'tetragonal_1': 6, 'tetragonal_2': 7,
                                        'orthorhombic': 9, 'monoclinic_1': 13,
                                        'monoclinic_2': 13, 'triclinic': 21}
                self.crystal_class=""
                self.stiffconsts=[]
                self.archive=''


    def manualinput(self, crystal_class, *stiff_consts):

        if crystal_class not in self.crystal_classes:
            raise KeyError('Unknown crystal class: %s. Use: %s' % (crystal_class, [s for s in self.crystal_classes]))

        self.crystal_class = crystal_class
        self.nconsts = self.crystal_classes[crystal_class]
        self.C = self.StiffnessMatrix(crystal_class, *stiff_consts)
        self.S = self.ComplianceMatrix(self.C)
        self.St = self.ComplianceTensor(self.S)
    def autoinput(self, archive):
        self.archive=archive
        with open(self.archive) as a:
            self.lines=list(a)
            self.c11,self.c12,self.c13,self.c14,self.c15,self.c16=[float(n) for n in list(self.lines[15].split())]
            self.c21,self.c22,self.c23,self.c24,self.c25,self.c26=[float(n) for n in list(self.lines[16].split())]
            self.c31,self.c32,self.c33,self.c34,self.c35,self.c36=[float(n) for n in list(self.lines[17].split())]
            self.c41,self.c42,self.c43,self.c44,self.c45,self.c46=[float(n) for n in list(self.lines[18].split())]
            self.c51,self.c52,self.c53,self.c54,self.c55,self.c56=[float(n) for n in list(self.lines[19].split())]
            self.c61,self.c62,self.c63,self.c64,self.c65,self.c66=[float(n) for n in list(self.lines[20].split())]

            if "Isotropic" in list(self.lines[4].split()):
                self.crystal_class='isotropic'
                self.C = self.StiffnessMatrix(self.crystal_class, self.c11,self.c12)
            if "Cubic" in list(self.lines[4].split()):
                self.crystal_class='cubic'
                self.C = self.StiffnessMatrix(self.crystal_class, self.c11,self.c12,self.c44)
            if "Hexagonal" in list(self.lines[4].split()):
                self.crystal_class='hexagonal'
                self.C = self.StiffnessMatrix(self.crystal_class, self.c11,self.c12,self.c13,self.c33,self.c44)
            if "Trigonal" in list(self.lines[4].split()):
                if 'I' in list(self.lines[4].split()):
                    self.crystal_class='trigonal_1'
                    self.C = self.StiffnessMatrix(self.crystal_class, self.c11,self.c12,self.c13,self.c14,self.c33,self.c44)
                else:
                    self.crystal_class='trigonal_2'
                    self.C = self.StiffnessMatrix(self.crystal_class, self.c11,self.c12,self.c13,self.c14,self.c15,self.c44)
            if "Tetragonal" in list(self.lines[4].split()):
                if 'I' in list(self.lines[4].split()):
                    self.crystal_class='tetragonal _1'
                    self.C = self.StiffnessMatrix(self.crystal_class,self.c11,self.c12,self.c13,self.c33,self.c44,self.c66)
                else:
                    self.crystal_class='tetragonal_2'
                    self.C = self.StiffnessMatrix(self.crystal_class,self.c11,self.c12,self.c13,self.c16,self.c33,self.c44,self.c66)
            if "Monoclinic" in list(self.lines[4].split()):
                if 'I' in list(self.lines[4].split()):
                    self.crystal_class='monoclinic_1'
                    self.C = self.StiffnessMatrix(self.crystal_class, self.c11,self.c12,self.c13,self.c15,self.c22,self.c23,self.c25,self.c33,self.c35,self.c44,self.c55,self.c66)
                else:
                    self.crystal_class='monoclinic_2'
                    self.C = self.StiffnessMatrix(self.crystal_class, self.c11,self.c12,self.c13,self.c16,self.c22,self.c23,self.c26,self.c33,self.c36,self.c44,self.c45,self.c55,self.c66)
            if "Orthorhombic" in list(self.lines[4].split()):
                self.crystal_class='orthorhombic'
                self.C = self.StiffnessMatrix(self.crystal_class, self.c11,self.c12,self.c13,self.c22,self.c23,self.c33,self.c44,self.c55,self.c66)
            if "Triclinic" in list(self.lines[4].split()):
                self.crystal_class='triclinic'
                self.C = self.StiffnessMatrix(self.crystal_class, self.c11,self.c12,self.c13,self.c14,self.c15,self.c16,self.c22,self.c23,self.c24,self.c25,self.c26,

                self.c33,self.c34,self.c35,self.c36,self.c44,self.c45,self.c46,self.c55,self.c56,self.c66)
        self.S = self.ComplianceMatrix(self.C)
        self.St = self.ComplianceTensor(self.S)
