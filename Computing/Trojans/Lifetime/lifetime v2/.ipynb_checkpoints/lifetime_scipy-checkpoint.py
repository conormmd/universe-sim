import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from scipy.signal import find_peaks
from functions2 import *

#Params
G = 6.6726e-11; d = 7.7834082e11; J_mass = 1.898130e27; S_mass = 1.98847e30; r_S = d*(J_mass/(J_mass+S_mass)); r_J = d*(S_mass/(J_mass+S_mass)); P = np.sqrt((d**3/(G*(S_mass+J_mass))))*2*np.pi; V_S = (2*np.pi*r_S)/P; V_J = (2*np.pi*r_J)/P

time = np.load('time.npy')*10
dist = np.load('dist.npy')
N=len(dist)

for i in range(5000,N-5000):
    #if dist[i-5000]<dist[i-4500]<dist[i-4000]<dist[i-3500]<dist[i-3000]<dist[i-2500]<dist[i-2000]<dist[i-1500]<dist[i-1000]<dist[i-500]<dist[i] and dist[i+5000]<dist[i+4500]<dist[i+4000]<dist[i+3500]<dist[i+3000]<dist[i+2500]<dist[i+2000]<dist[i+1500]<dist[i+1000]<dist[i+500]<dist[i]:
    if dist[i-4000]<dist[i-3500]<dist[i-3000]<dist[i-2500]<dist[i-2000]<dist[i-1500]<dist[i-1000]<dist[i-500]<dist[i] and dist[i+4000]<dist[i+3500]<dist[i+3000]<dist[i+2500]<dist[i+2000]<dist[i+1500]<dist[i+1000]<dist[i+500]<dist[i]:
        plt.scatter(time[i],dist[i],color='blue')
    else:
        continue
plt.show
