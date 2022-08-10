import numpy as np
import matplotlib.pyplot as plt
import sys


#######################################################################

file = sys.argv[1]
xcs = sys.argv[2]
ycs = sys.argv[3]

data = np.loadtxt(file)

xc = int(xcs)
yc = int(ycs)

x = data[:,xc-1]
y = data[:,yc-1]

plt.plot(x,y,'k',linewidth=1.0)  



plt.show()
