import numpy as np
import matplotlib.pyplot as plt
import sys


#######################################################################


data1 = np.loadtxt('model.out')
r1 = data1[:,0]
f1 = data1[:,4]

data2 = np.loadtxt('mesh.out')
r2 = data2[:,0]
f2 = data2[:,4]

plt.plot(f1,r1,'k--',linewidth=1.0)
plt.plot(f2,r2,'r',linewidth=1.0)  



plt.show()
