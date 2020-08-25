import numpy as np
import random as rand
from rvv_functions import *
from rvv_pushers import *

B = np.array([rand.random(),rand.random(),rand.random()])
E = np.array([rand.random(),rand.random(),rand.random()])
x  = np.array([rand.random(),rand.random(),rand.random()])
u_old  = np.array([rand.random(),rand.random(),rand.random()]) - 0.1
c = 1
dt = 0.5

u_new = boris(x,u_old,E,B,dt)
defect = u_new - u_old - (E+np.cross((u_old+u_new)/2,B))*dt
print(defect)
