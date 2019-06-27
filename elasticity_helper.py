import numpy as np
import numpy.linalg as la


class StiffnessConstants:

    '''
        StiffnessConstants class

        input:
            - crystal_class: classified according to Nye (1956).
                One of the following:
                   'isotropic', 'cubic', 'hexagonal',
                   'trigonal_1', 'trigonal_2', 'tetragonal_1',
                   'tetragonal_2', 'orthorhombic', 'monoclinic_1',
                   'monoclinic_2', 'triclinic'
            - compliance constants: variable number of constants,
               depending on crystal_class.
               Units: TPa^{-1}

        output: stiffness constants in GPa

    '''

    def __init__(self, crystal_class, *compl_consts):

        self.crystal_classes = {'isotropic': 2, 'cubic': 3, 'hexagonal': 5,
                                'trigonal_1': 6, 'trigonal_2': 7,
                                'tetragonal_1': 6, 'tetragonal_2': 7,
                                'orthorhombic': 9, 'monoclinic_1': 13,
                                'monoclinic_2': 13, 'triclinic': 21}

        if crystal_class not in self.crystal_classes:
            raise KeyError('Unknown crystal class: %s. Use: %s' % (crystal_class, [s for s in self.crystal_classes]))

        self.crystal_class = crystal_class
        self.nconsts = self.crystal_classes[crystal_class]

        self.S = self.ComplianceMatrix(crystal_class, *compl_consts)
        self.C = self.StiffnessMatrix(self.S)

    def ComplianceMatrix(self, crystal_class, *compl_consts):

        m, n = self.crystal_classes[crystal_class], len(compl_consts)

        if m != n:
            raise IndexError('Crystal class %s has %d independent constants, but %d were provided' % (crystal_class, m, n))

        if crystal_class == 'isotropic':
            S11, S12 = compl_consts
            S44 = .5 * (S11 - S12)
            S = np.array([
                [S11, S12, S12,  0.,  0.,  0.],
                [0., S11, S12,  0.,  0.,  0.],
                [0.,  0., S11,  0.,  0.,  0.],
                [0.,  0.,  0., S44,  0.,  0.],
                [0.,  0.,  0.,  0., S44,  0.],
                [0.,  0.,  0.,  0.,  0., S44]
               ]
              )
        elif crystal_class == 'cubic':      # all classes
            S11, S12, S44 = compl_consts
            S = np.array([
                [S11, S12, S12,  0.,   0.,  0.],
                [0., S11, S12,  0.,   0.,  0.],
                [0.,  0., S11,  0.,   0.,  0.],
                [0.,  0.,  0., S44,   0.,  0.],
                [0.,  0.,  0.,  0.,  S44,  0.],
                [0.,  0.,  0.,  0.,   0., S44]
               ]
              )
        elif crystal_class == 'hexagonal':      # all classes
            S11, S12, S13, S33, S44 = compl_consts
            S66 = .5 * (S11 - S12)
            S = np.array([
                [S11, S12, S13,  0.,  0.,  0.],
                [0., S11, S13,  0.,  0.,  0.],
                [0.,  0., S33,  0.,  0.,  0.],
                [0.,  0.,  0., S44,  0.,  0.],
                [0.,  0.,  0.,  0., S44,  0.],
                [0.,  0.,  0.,  0.,  0., S66]
               ]
              )
        elif crystal_class == 'tetragonal_1':  # classes 4mm, -42m, 422, 4/mmm
            S11, S12, S13, S33, S44, S66 = compl_consts
            S = np.array([
                [S11, S12, S13,  0.,  0.,  0.],
                [0., S11, S13,  0.,  0.,  0.],
                [0.,  0., S33,  0.,  0.,  0.],
                [0.,  0.,  0., S44,  0.,  0.],
                [0.,  0.,  0.,  0., S44,  0.],
                [0.,  0.,  0.,  0.,  0., S66]
               ]
              )
        elif crystal_class == 'tetragonal_2':    # classes 4, -4, 4/m
            S11, S12, S13, S33, S44, S66, S16 = compl_consts
            S = np.array([
                [S11, S12, S13,  0.,  0.,  S16],
                [0., S11, S13,  0.,  0., -S16],
                [0.,  0., S33,  0.,  0.,   0.],
                [0.,  0.,  0., S44,  0.,   0.],
                [0.,  0.,  0.,  0., S44,   0.],
                [0.,  0.,  0.,  0.,  0.,  S66]
               ]
              )
        elif crystal_class == 'trigonal_1':     # classes 32, -3m, 3m
            S11, S12, S13, S14, S33, S44 = compl_consts
            S66 = .5 * (S11 - S12)
            S = np.array([
                [S11, S12, S13,  S14,  0.,  0.],
                [0., S11, S13, -S14,   0,  0.],
                [0., 0.,  S33,   0.,  0.,  0.],
                [0., 0.,   0.,  S44,  0.,  0.],
                [0., 0.,   0.,   0., S44, S14],
                [0., 0.,   0.,   0.,  0., S66]
               ]
              )
        elif crystal_class == 'trigonal_2':     # classes 3, -3
            S11, S12, S13, S33, S14, S15, S44 = compl_consts
            S66 = .5 * (S11 - S12)
            S = np.array([
                [S11, S12, S13,  S14,  S15,   0.],
                [0., S11, S13, -S14, -S15,   0.],
                [0.,  0., S33,   0.,   0.,   0.],
                [0.,  0.,  0.,  S44,   0., -S15],
                [0.,  0.,  0.,   0.,  S44,  S14],
                [0.,  0.,  0.,   0.,   0.,  S66]
               ]
              )
        elif crystal_class == 'orthorhombic':       # all classes
            S11, S12, S13, S22, S23, S33, S44, S55, S66 = compl_consts
            S = np.array([
                [S11, S12, S13,  0.,  0.,  0.],
                [0., S22, S23,  0.,  0.,  0.],
                [0.,  0., S33,  0.,  0.,  0.],
                [0.,  0.,  0., S44,  0.,  0.],
                [0.,  0.,  0.,  0., S55,  0.],
                [0.,  0.,  0.,  0.,  0., S66]
               ]
              )
        elif crystal_class == 'monoclinic_1':       # diad || x2
            S11, S12, S13, S22, S23, S33, S15, S25, S35, S44, S46, S55, S66 = compl_consts
            S = np.array([
                [S11, S12, S13,  0., S15,  0.],
                [0., S22, S23,  0., S25,  0.],
                [0.,  0., S33,  0., S35,  0.],
                [0.,  0.,  0., S44,  0., S46],
                [0.,  0.,  0.,  0., S55,  0.],
                [0.,  0.,  0.,  0.,  0., S66]
               ]
              )
        elif crystal_class == 'monoclinic_2':       # diad || x3
            S11, S12, S13, S22, S23, S33, S16, S26, S36, S44, S45, S55, S66 = compl_consts
            S = np.array([
                [S11, S12, S13,  0.,  0., S16],
                [0., S22, S23,  0.,  0., S26],
                [0.,  0., S33,  0.,  0., S36],
                [0.,  0.,  0., S44, S45,  0.],
                [0.,  0.,  0.,  0., S55,  0.],
                [0.,  0.,  0.,  0.,  0., S66]
               ]
              )
        elif crystal_class == 'triclinic':      # all classes
            S11, S12, S13, S14, S15, S16, S22, S23, S24, S25, S26, S33, S34, S35, S36, S44, S45, S46, S55, S56, S66 = compl_consts
            S = np.array([
                [S11, S12, S13, S14, S15, S16],
                [0., S22, S23, S24, S25, S26],
                [0.,  0., S33, S34, S35, S36],
                [0.,  0.,  0., S44, S45, S46],
                [0.,  0.,  0.,  0., S55, S56],
                [0.,  0.,  0.,  0.,  0., S66]
               ]
              )

        S = S + S.T - np.diag(S.diagonal())     # symmetrize S
        return S

    def StiffnessMatrix(self, S):

        return la.inv(S)
