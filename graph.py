import numpy as np
import matplotlib.pyplot as plt

#data = np.readtxt('random.dat');
data = np.random.random((1000,1))*19+1;

plt.hist(data);
plt.savefig('graph.pdf',format='pdf');
plt.close()
