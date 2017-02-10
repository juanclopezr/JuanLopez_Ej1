import matplotlib.pyplot as plt
import numpy as np

potential = np.loadtxt('potential.dat')
fieldx = np.loadtxt('fieldx.dat')
fieldy = np.loadtxt('fieldy.dat')

x = np.linspace(-2.5,2.5,256)
y = np.linspace(2.5,-2.5,256)

plt.imshow(potential,extent=[-2.5,2.5,-2.5,2.5])
plt.xlim([-2.5,2.5])
plt.ylim([-2.5,2.5])
plt.xlabel('$x$(cm)')
plt.ylabel('$y$(cm)')
plt.colorbar(label='Voltage(V)')
plt.streamplot(x,y,fieldx,fieldy,color='black')
plt.title('Electric potential and electric field in capacitor')
plt.savefig('placas.pdf')
