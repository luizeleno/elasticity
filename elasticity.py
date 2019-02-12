import numpy as np
import numpy.linalg as la

'''
    elasticity
    
    A collection of classes and routines to
    help the treatment and presentation
    of single and polycristal elastic properties
    
    Prof. Luiz T. F. Eleno
    Materials Eng. Dept.
    Lorena School of Engineering
    University of Sao Paulo
'''

class ElasticityTheory:
    
    '''
        ElasticityTheory class
        
        input:
            - crystal_class: classified according to Nye (1956). One of the following: 'isotropic', 'cubic', 'hexagonal', 'trigonal_1', 'trigonal_2', 'tetragonal_1', 'tetragonal_2', 'orthorhombic', 'monoclinic_1', 'monoclinic_2', 'triclinic'
            - stiffness constants: variable number of constants, depending on crystal_class.
            
        output: functions to calculate the Compliance and Stiffness matrices, and the Compliance tensor
        
    '''
    
    def __init__(self, crystal_class, *stiff_consts):

        self.crystal_classes = {'isotropic':2, 'cubic':3, 'hexagonal':5, 'trigonal_1':6, 'trigonal_2':7, 'tetragonal_1':6, 'tetragonal_2':7, 'orthorhombic':9, 'monoclinic_1':13, 'monoclinic_2':13, 'triclinic':21}

        if crystal_class not in self.crystal_classes:
            raise KeyError('Unknown crystal class: %s. Use: %s' % (crystal_class, [s for s in self.crystal_classes]))

        self.crystal_class = crystal_class
        self.nconsts = self.crystal_classes[crystal_class]
        
        self.C = self.StiffnessMatrix(crystal_class, *stiff_consts)
        self.S = self.ComplianceMatrix(self.C)
        self.St = self.ComplianceTensor(self.S)
        
    def StiffnessMatrix(self, crystal_class, *stiff_consts):
        
        m, n = self.crystal_classes[crystal_class], len(stiff_consts)
        
        if m != n :
            raise IndexError('Crystal class %s has %d independent constants, but %d were provided' % (crystal_class, m, n))
        
        if crystal_class is 'cubic':
            C11, C12, C44 = stiff_consts
            C=np.array([
                [C11, C12, C12, 0.,   0.,  0.],
                [C12, C11, C12, 0.,   0.,  0.],
                [C12, C12, C11, 0.,   0.,  0.],
                [0.,   0.,  0., C44,  0.,  0.],
                [0.,   0.,  0., 0.,  C44,  0.],
                [0.,   0.,  0., 0.,   0.,  C44]
               ]
              )
        elif crystal_class is 'isotropic':
            C11, C12 = stiff_consts
            C44 = .5 * ( C11 - C12 )
            C=np.array([
                [C11, C12, C12, 0.,   0.,  0.],
                [C12, C11, C12, 0.,   0.,  0.],
                [C12, C12, C11, 0.,   0.,  0.],
                [0.,   0.,  0., C44,  0.,  0.],
                [0.,   0.,  0., 0.,  C44,  0.],
                [0.,   0.,  0., 0.,   0.,  C44]
               ]
              )
        elif crystal_class is 'hexagonal':
            C11, C12, C13, C33, C44 = stiff_consts
            C66 = .5 * ( C11 - C12 )
            C=np.array([
                [C11, C12, C13, 0.,   0.,  0.],
                [C12, C11, C13, 0.,   0.,  0.],
                [C13, C13, C33, 0.,   0.,  0.],
                [0.,   0.,  0., C44,  0.,  0.],
                [0.,   0.,  0., 0.,  C44,  0.],
                [0.,   0.,  0., 0.,   0.,  C66]
               ]
              )
        elif crystal_class is 'tetragonal_1':
            C11, C12, C13, C33, C44, C66 = stiff_consts
            C=np.array([
                [C11, C12, C13, 0.,   0.,  0.],
                [C12, C11, C13, 0.,   0.,  0.],
                [C13, C13, C33, 0.,   0.,  0.],
                [0.,   0.,  0., C44,  0.,  0.],
                [0.,   0.,  0., 0.,  C44,  0.],
                [0.,   0.,  0., 0.,   0.,  C66]
               ]
              )
        elif crystal_class is 'tetragonal_2':
            C11, C12, C13, C33, C44, C66, C16 = stiff_consts
            C=np.array([
                [C11, C12, C13, 0.,  0.,  C16],
                [C12, C11, C13, 0.,  0., -C16],
                [C13, C13, C33, 0.,  0.,  0.],
                [0.,   0.,  0., C44, 0.,  0.],
                [0.,   0.,  0., 0.,  C44, 0.],
                [C16, -C16, 0., 0.,  0.,  C66]
               ]
              )
        
        return C
              
    def ComplianceMatrix(self, C):
        
        return la.inv(C)

    def ComplianceTensor(self, S):

        St = np.zeros((3,3,3,3), dtype=float)

        # 6
        St[0, 0, 0, 0] = S[0, 0]
        St[0, 0, 1, 1] = S[0, 1]
        St[0, 0, 2, 2] = S[0, 2]
        St[0, 0, 1, 2] = .5 * S[0, 3]
        St[0, 0, 0, 2] = .5 * S[0, 4]
        St[0, 0, 0, 1] = .5 * S[0, 5]
        
        # 11
        St[1, 1, 0, 0] = St[0, 0, 1, 1]
        St[2, 2, 0, 0] = St[0, 0, 2, 2]
        St[0, 0, 2, 1] = St[1, 2, 0, 0] = St[2, 1, 0, 0] = St[0, 0, 1, 2]
        St[0, 0, 2, 0] = St[0, 2, 0, 0] = St[2, 0, 0, 0] = St[0, 0, 0, 2]
        St[0, 0, 1, 0] = St[0, 1, 0, 0] = St[1, 0, 0, 0] = St[0, 0, 0, 1]
        
        # 5
        St[1, 1, 1, 1] = S[1, 1]
        St[1, 1, 2, 2] = S[1, 2]
        St[1, 1, 1, 2] = .5 * S[1, 3]
        St[1, 1, 0, 2] = .5 * S[1, 4]
        St[1, 1, 0, 1] = .5 * S[1, 5]
        
        # 10
        St[2, 2, 1, 1] = St[1, 1, 2, 2]
        St[1, 1, 2, 1] = St[1, 2, 1, 1] = St[2, 1, 1, 1] = St[1, 1, 1, 2]
        St[1, 1, 2, 0] = St[0, 2, 1, 1] = St[2, 0, 1, 1] = St[1, 1, 0, 2]
        St[1, 1, 1, 0] = St[0, 1, 1, 1] = St[1, 0, 1, 1] = St[1, 1, 0, 1]
        
        # 4
        St[2, 2, 2, 2] = S[2, 2]
        St[2, 2, 1, 2] = .5 * S[2, 3]
        St[2, 2, 0, 2] = .5 * S[2, 4]
        St[2, 2, 0, 1] = .5 * S[2, 5]
        
        # 9
        St[2, 2, 2, 1] = St[1, 2, 2, 2] = St[2, 1, 2, 2] = St[2, 2, 1, 2]
        St[2, 2, 2, 0] = St[0, 2, 2, 2] = St[2, 0, 2, 2] = St[2, 2, 0, 2]
        St[2, 2, 1, 0] = St[0, 1, 2, 2] = St[1, 0, 2, 2] = St[2, 2, 0, 1]
        
        # 3
        St[1, 2, 1, 2] = .25 * S[3, 3]
        St[1, 2, 0, 2] = .25 * S[3, 4]
        St[1, 2, 0, 1] = .25 * S[3, 5]
        
        # 17
        St[1, 2, 2, 1] = St[2, 1, 1, 2] = St[2, 1, 2, 1] = St[1, 2, 1, 2]
        St[1, 2, 2, 0] = St[2, 1, 0, 2] = St[2, 1, 2, 0] = St[0, 2, 1, 2] = St[0, 2, 2, 1] = St[2, 0, 1, 2] = St[2, 0, 2, 1] = St[1, 2, 0, 2]
        St[1, 2, 1, 0] = St[2, 1, 0, 1] = St[2, 1, 1, 0] = St[0, 1, 1, 2] = St[0, 1, 2, 1] = St[1, 0, 1, 2] = St[1, 0, 2, 1] = St[1, 2, 0, 1]
        
        # 2
        St[0, 2, 0, 2] = .25 * S[4, 4]
        St[0, 2, 0, 1] = .25 * S[4, 5]
        
        # 10
        St[0, 2, 2, 0] = St[2, 0, 0, 2] = St[2, 0, 2, 0] = St[0, 2, 0, 2]
        St[0, 2, 1, 0] = St[2, 0, 0, 1] = St[2, 0, 1, 0] = St[0, 1, 0, 2] = St[0, 1, 2, 0] = St[1, 0, 0, 2] = St[1, 0, 2, 0] = St[0, 2, 0, 1]
        
        # 1
        St[0, 1, 0, 1] = .25 * S[5, 5]
        
        # 3
        St[0, 1, 1, 0] = St[1, 0, 0, 1] = St[1, 0, 1, 0] = St[0, 1, 0, 1]
        
        # total: 81
        
        return St

class DirectionalProperties(ElasticityTheory):
    
    '''
        Directional properties: functions to calculate:
            - unidimensional bulk modulus(l)
            - Young's modulus(l)
            - Shear Modulus(l, n)
            #- Poisson's ratio(l, n)
            
        as a function of crystallographic direction l and n
    '''
        
    def BulkModulus(self, l):
        
        B = 0
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    B += self.St[i, j, k, k] * l[i] * l[j] 
        return 1.0 / B
    
    def YoungModulus(self, l):
        
        E = 0
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for m in range(3):
                        E += self.St[i, j, k, m] * l[k] * l[m] 
        return 1.0 / E

    def ShearModulus(self, l, n):
        
        G = 0
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for m in range(3):
                        G +=self.St[i, j, k, m] * l[i] * n[j] * l[k] * n[m] 
        return .25 / G

    def PoissonRatio(self, l, n):
        
        p = 0
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
        
        sii = S[0,0] + S[1,1] + S[2,2]
        sij = S[0,1] + S[0,2] + S[1,2]
        skl = S[3,3] + S[4,4] + S[5,5]
        self.BR = 1. / ( sii + 2. * sij )
        self.GR = 15. / ( 4. * sii - sij + 3. * skl )
    
    def Voigt(self):
        
        cii = C[0,0] + C[1,1] + C[2,2]
        cij = C[0,1] + C[0,2] + C[1,2]
        ckl = C[3,3] + C[4,4] + C[5,5]
        self.BV = ( cii + 2. * cij ) / 9.
        self.GV = ( cii - cij + 3. * ckl ) / 15.
    
    def VRH(self):
        
        self.Reuss()
        self.Voigt()
        
        B = .5 * (self.BR + self.BV)
        G = .5 * (self.GR + self.GV)
        
        self.BulkModulus, self.ShearModulus = B, G
        self.YoungModulus = 9. * B * G / ( 3. * B + G )
        self.PoissonRatio = ( 3. * B - 2. * G ) / 2. / ( 3. * B + G )

if __name__ == '__main__':
    solido = DirectionalProperties('cubic', 300., 500., 200.)
    print (solido.YoungModulus([1,0,0]))
    #print (solido.ShearModulus([1,0,0], [0,1,0]))
