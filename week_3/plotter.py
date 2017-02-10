import numpy as np
import matplotlib.pyplot as plt

potential = np.genfromtxt('potential.dat')

#xi, yi = np.linspace(0, 255, 256), np.linspace(0, 255, size)

#xi, yi = np.meshgrid(xi, yi)

#dx, dy = np.gradient(-potential)
plt.imshow(potential)
#plt.streamplot(xi, yi, dy, dx, color='black')
#plt.xlim(0, size-1)
#plt.ylim(size-1, 0)
plt.title('Capacitor')
#plt.show()
plt.savefig('V.png')
