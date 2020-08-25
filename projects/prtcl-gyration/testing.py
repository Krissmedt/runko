import numpy as np

a = np.zeros((4,10,3))

for i in range(0,a.shape[0]):
    for j in range(0,a.shape[1]):
        for k in range(0,a.shape[2]):
            a[i,j,k] += 1
            print(a[i,j,k])


A = np.ravel(a)

print(A)
