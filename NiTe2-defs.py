import numpy as np
import elasticity as el

pbesoc = el.Elasticity('trigonal_1', 113.7, 36.6, 27.2, -6.5, 45.7, 11.2)

a = 7.364795908
c = a * 1.379652

# for p in range(1, 6):
#     S = pbesoc.S
#     sigma = [-p, -p, -p, 0., 0., 0.]
#     epsilon = np.dot(S, sigma)
#     a1 = a * (1 + epsilon[0])
#     c1 = c * (1 + epsilon[2]) / a1
#     print(f'celldm(1) = {a1},')
#     print(f'celldm(3) = {c1},')
#
# print()

for t in range(-5, 6):
    S = pbesoc.S
    sigma = [0., 0., t, 0., 0., 0.]
    epsilon = np.dot(S, sigma)
    print(epsilon)
    a1 = a * (1 + epsilon[0])
    c1 = c * (1 + epsilon[2]) / a1
    # print(f'celldm(1) = {a1},')
    # print(f'celldm(3) = {c1},')
