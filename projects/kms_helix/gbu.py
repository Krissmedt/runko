import numpy as np

c = 1
gammas = np.linspace(1,10,100,dtype=np.float)

for g in range(0,gammas.shape[0]):
    beta = np.sqrt(1-1./g**2.)
    u = beta*c

    print(beta)
