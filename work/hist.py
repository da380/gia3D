import numpy as np
import matplotlib.pyplot as plt
import sys


x = np.loadtxt('test.out')
plt.hist(x,100,density = True)
plt.show()
