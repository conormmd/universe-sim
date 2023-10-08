import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from functions2 import *

#Params
G = 6.6726e-11; d = 7.7834082e11; J_mass = 1.898130e27; S_mass = 1.98847e30; r_S = d*(J_mass/(J_mass+S_mass)); r_J = d*(S_mass/(J_mass+S_mass)); P = np.sqrt((d**3/(G*(S_mass+J_mass))))*2*np.pi; V_S = (2*np.pi*r_S)/P; V_J = (2*np.pi*r_J)/P

time = np.load('time.npy')
dist = np.load('dist.npy')
N = len(dist)

bottom = np.loadtxt('bottom_index.txt')
top = np.loadtxt('top_index.txt')

bot = np.zeros(len(bottom))
t_bot = np.zeros(len(bottom))

for i in range(0,len(bottom)):
    j = int(bottom[i])
    bot[i] = dist[j]
    t_bot[i] = time[j]
up = np.zeros(len(top))
t_up = np.zeros(len(top))
for i in range(0,len(top)):
    j = int(top[i])
    up[i] = dist[j]
    t_up[i] = time[j]
goodup=[]
t_goodup=[]

for i in range(0,len(up)):
    if up[i]>0.75:
        plt.scatter(t_up[i],up[i],color='blue')
    else:
        continue


#plt.xlim(0,1e7)
plt.ylim(0.25,2.2)
plt.show()

