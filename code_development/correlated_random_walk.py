# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 16:12:56 2021

@author: Romain
"""

# Pure random walk in 1D
from random import seed
from random import random
from matplotlib import pyplot
seed(1)
random_walk = list()
random_walk.append(0) # Starting position
for i in range(1, 1000):
	movement = -1 if random() < 0.5 else 1
	value = random_walk[i-1] + movement
	random_walk.append(value)
pyplot.plot(random_walk)
pyplot.show()


# Random walk with persistence in 1D: Correlated Random Walk
from matplotlib import pyplot
import random
tau = 1.0 # time step
P = 10.0 # persistence: ~how many timesteps between direction changes (from 1 (RW) to +inf (directional))
t_ratio = tau/P
CRW = list()

# Starting conditions
CRW.append(0)
movement0 = 0
movementr = movement0

for i in range(1, 1000):
    rand = random.uniform(0,1)
    movementr = movement0 if rand > t_ratio else -1 if random.uniform(0,1) < 0.5 else 1
    movement0 = movementr      
    value = CRW[i-1]+movementr
    CRW.append(value)
pyplot.plot(CRW)
pyplot.show()



