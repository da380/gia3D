import numpy as np
import matplotlib.pyplot as plt
import plot_util as pu
import sys


#######################################################################

file1 = sys.argv[1]
pu.plot_2D_data(file1)




file2 = sys.argv[2]
pu.plot_2D_data(file2)


plt.show()
