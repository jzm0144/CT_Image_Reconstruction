from scipy.io import loadmat
import scipy.misc
import os
import sys
import multiprocessing
from pynufft import NUFFT_hsa, NUFFT_cpu
import numpy as np 
import ipdb as ipdb
import matplotlib.pyplot as plt


spoke_range = (np.arange(0, 512) - 256.0 )* np.pi/ 256  # normalized between -pi and pi
M = 512*360
om = np.empty((M,2), dtype = np.float32)


for angle in range(0, 360):
   radian = angle * 2 * np.pi/ 360.0
   spoke_x =  spoke_range * np.cos(radian)
   spoke_y =  spoke_range * np.sin(radian)
   om[512*angle : 512*(angle + 1) ,0] = spoke_x
   om[512*angle : 512*(angle + 1) ,1] = spoke_y


plt.plot(om[:,0], om[:,1],'.')
plt.show()

'''
data = loadmat('project_1_new_data.mat')
kspace = data['kspace']
theta = data['theta']

data = loadmat('ktraj.mat')
ktraj = data['ktraj']

numProjections = kspace.shape[1]
numReadouts = kspace.shape[0]
gridSize = int(numReadouts/2)

print('Number of Projections = ', numProjections)
print('Number of Readout Values = ',numReadouts)

imgGrid = np.zeros((gridSize, gridSize))

myNufft = NUFFT_cpu()
myNufft.plan(om = ktraj,
             Nd = (gridSize, gridSize),
             Kd = (numReadouts, numReadouts),
             Jd = (6,6))

'''
