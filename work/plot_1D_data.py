import numpy as np
import matplotlib.pyplot as plt
import sys


#######################################################################

# set some plotting parameters
plt.rc('xtick', labelsize=18) 
plt.rc('ytick', labelsize=18) 
font = {'size'   : 20}
plt.rc('font', **font)
plt.rcParams['figure.figsize'] = [16, 8]


file = sys.argv[1]
xcs = sys.argv[2]
ycs = sys.argv[3]

data = np.loadtxt(file)

xc = int(xcs)
yc = int(ycs)

x = data[:,xc-1]
y = data[:,yc-1]

plt.semilogx(x,y,'k',linewidth=1.0)
plt.xlabel('non-dimensional wavenumber')
plt.ylabel('gravitational potential response')




plt.show()
