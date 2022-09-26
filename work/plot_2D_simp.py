import numpy as np
import matplotlib.pyplot as plt
import plot_util as pu
import sys


#######################################################################

file1 = sys.argv[1]

pu.plot_2D_data(file1)

plt.gca().invert_yaxis()
plt.show()
