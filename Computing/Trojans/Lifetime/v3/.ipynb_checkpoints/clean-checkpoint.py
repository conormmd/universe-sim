import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection


#Params
G = 6.6726e-11; d = 7.7834082e11; J_mass = 1.898130e27; S_mass = 1.98847e30; r_S = d*(J_mass/(J_mass+S_mass)); r_J = d*(S_mass/(J_mass+S_mass)); P = np.sqrt((d**3/(G*(S_mass+J_mass))))*2*np.pi; V_S = (2*np.pi*r_S)/P; V_J = (2*np.pi*r_J)/P

time = np.load('time.npy')*10
dist = np.load('dist.npy')

x = np.zeros(30000)
y=np.zeros(30000)

for i in range(0,30000):
    j = 1000*i
    x[i]=time[j]
    y[i]=dist[j]

plt.plot(x,y)
plt.show()

np.save('x.npy',x)
np.save('y.npy',y)


        
