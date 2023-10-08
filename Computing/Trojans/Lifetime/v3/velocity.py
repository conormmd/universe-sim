import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection


#Params
G = 6.6726e-11; d = 7.7834082e11; J_mass = 1.898130e27; S_mass = 1.98847e30; r_S = d*(J_mass/(J_mass+S_mass)); r_J = d*(S_mass/(J_mass+S_mass)); P = np.sqrt((d**3/(G*(S_mass+J_mass))))*2*np.pi; V_S = (2*np.pi*r_S)/P; V_J = (2*np.pi*r_J)/P

time = np.load('time.npy')*10
dist = np.load('dist.npy')
velocity = np.zeros(len(time)-1)
v_t = np.zeros(len(time)-1)
for i in range(1,len(time)-1):
    velocity[i] = (dist[i-1]-dist[i])/(time[i-1]-time[i])
    v_t[i] = (time[i-1]+time[i])/2

np.save('v_x.npy', v_t)
np.save('v_y.npy',velocity)
plt.show()


        
