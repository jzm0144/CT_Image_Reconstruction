from scipy.io import loadmat
import scipy.misc
import os
import sys
import multiprocessing
from pynufft import NUFFT_hsa, NUFFT_cpu
import numpy as np 
import ipdb as ipdb
import matplotlib.pyplot as plt

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
