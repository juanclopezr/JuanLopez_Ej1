import numpy as np
import plt.pyplot as plt

#data = np.readtxt('random.dat');
data = np.random.random((1,1000));

plt.hist(data);
plt.savefig('graph.pdf',format='pdf');
plt.close()
