from scipy.io import loadmat
import scipy.misc
import os
import sys
import multiprocessing
from pynufft import NUFFT_hsa, NUFFT_cpu
import numpy as np 
import ipdb as ipdb
import matplotlib
import matplotlib.pyplot as plt


data = loadmat('project_1_new_data.mat')
kspace = data['kspace']
theta = data['theta']

spoke_range = (np.arange(0, 512) - 256.0 )* np.pi/ 256  # normalized between -pi and pi

M = 512*5120 #M = 512*360

om = np.empty((M,2), dtype = np.float32)

for index in range(0, 5120):
    angle = theta[0, index]
    radian = angle * 2 * np.pi/ 360.0
    spoke_x = spoke_range * np.cos(radian)
    spoke_y = spoke_range * np.sin(radian)

    om[512*index : 512*(index + 1), 0] = spoke_x
    om[512*index : 512*(index + 1), 1] = spoke_y

#plt.plot(om[:,0], om[:,1],'.')
#plt.title("Radial Kspace Trajectory")
#plt.show()


numProjections = kspace.shape[1]
numReadouts = kspace.shape[0]


print('Number of Projections = ', numProjections)
print('Number of Readout Values = ',numReadouts)


myNufft = NUFFT_cpu()
myNufft.plan(om = om,
             Nd = (256,256),
             Kd = (numReadouts, numReadouts),
             Jd = (2,2))

y = kspace.flatten(order='C')
image = myNufft.adjoint(y)
#y = myNufft.forward(image)

#ipdb.set_trace()

plt.subplot(2,2,1)
image0 = myNufft.solve(y, solver='cg',maxiter=50)

#img = image0.real/image0.real.max()
plt.title('Restored image (cg)')
plt.imshow(image0.real,
           cmap=matplotlib.cm.gray,
           norm=matplotlib.colors.Normalize(vmin=0.0, vmax=1))
plt.show()

