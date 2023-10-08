import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from functions2 import *

#Params
G = 6.6726e-11; d = 7.7834082e11; J_mass = 1.898130e27; S_mass = 1.98847e30; r_S = d*(J_mass/(J_mass+S_mass)); r_J = d*(S_mass/(J_mass+S_mass)); P = np.sqrt((d**3/(G*(S_mass+J_mass))))*2*np.pi; V_S = (2*np.pi*r_S)/P; V_J = (2*np.pi*r_J)/P

xy = np.load('comoving.npy')
N = len(xy)
dist = np.zeros(N)
time = np.zeros(N)

for i in range(0,N):
    x = xy[i][0]
    y = xy[i][1]
    dist[i]= np.sqrt(x**2 + (y-r_J)**2)/r_J
    time[i]=i*5000/3e7

np.save('time.npy',time)
np.save('dist.npy',dist)

#avg_time, avg_dist = rolling_average(250000,time,dist)

plt.plot(time,dist)
#plt.plot(avg_time,avg_dist)
#plt.xlim(0,1e7)
plt.ylim(0.25,2.2)
plt.show()

