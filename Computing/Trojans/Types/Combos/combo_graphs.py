import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from progressbar import ProgressBar
from functions import *
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix';matplotlib.rcParams['font.family'] = 'STIXGeneral';
plt.rcParams.update({'font.size': 12})
comoving_1 = np.load('comoving_1.npy')
comoving_2 = np.load('comoving_2.npy')
params = np.load('params.npy')
distance_1 = np.load('distance_1.npy')
distance_2 = np.load('distance_2.npy')

ticks = params[0]; a= params[1]
#theta_0 = ast_inits[1]; r_0 = ast_inits[0]
G = 6.6726e-11; d = 7.7834082e11; J_mass = 1.898130e27; S_mass = 1.98847e30; r_S = d*(J_mass/(J_mass+S_mass)); r_J = d*(S_mass/(J_mass+S_mass)); P = np.sqrt((d**3/(G*(S_mass+J_mass))))*2*np.pi; V_S = (2*np.pi*r_S)/P; V_J = (2*np.pi*r_J)/P
M_1 = S_mass;M_2 = J_mass;R = d
L4 = np.array([R*np.sin(np.pi/3),R*np.cos(np.pi/3)])
L5 = np.array([R*np.sin(np.pi/3),R*np.cos(np.pi/3)])

Jup_L5 = np.sqrt(L4[0]**2+(r_J-L4[1])**2)
plot_ticks = np.ndarray([ticks-1])
for i in range(0,ticks-1):
    plot_ticks[i]=i
    
plt.figure(figsize = (8,8))
ax_hs = plt.subplot(6,6,(1,15))
ax_tp = plt.subplot(6,6,(4,18))
ax_ds = plt.subplot(6,6,(19,30))
#plt.scatter(x=0, y=r_J,color='orange', label='Jupiter')
#plt.scatter(x=L4[0], y=L4[1],color='black', label='L4')
ax_tp.plot(comoving_1[:,0], comoving_1[:,1], color='blue')
ax_tp.plot(0,r_J,'o',color='#e36e4b',markersize=5,label='Jupiter');ax_tp.plot(0,0,'o',color='#FDB813',markersize=7,label='Sun')
ax_tp.plot(L5[0],L5[1],'^',color='black',zorder=10,label='L5')
ax_tp.set_xlim(-0.9e12, 0.9e12);ax_tp.set_ylim(-0.9e12, 0.9e12)
ax_tp.set_xticks([]),ax_tp.set_yticks([])
ax_tp.legend()

ax_hs.plot(comoving_2[:,0], comoving_2[:,1], color='red', label='Horseshoe')
ax_hs.plot(0,r_J,'o',color='#e36e4b',markersize=5);ax_hs.plot(0,0,'o',color='#FDB813',markersize=7,label='Sun')
ax_hs.plot(L5[0],L5[1],'^',color='black')
ax_hs.set_xlim(-0.9e12, 0.9e12);ax_hs.set_ylim(-0.9e12, 0.9e12)
ax_hs.set_xticks([]),ax_hs.set_yticks([])


L5_x = np.array([0, ticks])
L5_y = np.array([Jup_L5/r_J, Jup_L5/r_J])
ax_ds.plot(L5_x, L5_y, color='black', label='L5')
ax_ds.plot(plot_ticks, distance_1/r_J, color='blue', label='Tadpole')
ax_ds.plot(plot_ticks, distance_2/r_J, color='red', label='Horseshoe')
ax_ds.set_yticks([0,0.5,1,1.5,2,2.5])
ax_ds.set_yticklabels(['0','','${R}_{J}$','','2${R}_{J}$'])
ax_ds.set_xticks([0,5*P/100000,10*P/100000,15*P/100000,20*P/100000,25*P/100000])
ax_ds.set_xticklabels(['0','5${P}_{J}$','10${P}_{J}$','15${P}_{J}$','20${P}_{J}$','25${P}_{J}$',])
plt.xlim(0,100000)
plt.ylim(0,2.5)
plt.xlabel('Time',fontsize=16)
plt.ylabel('Distance from Jupiter',fontsize=16)
plt.subplots_adjust(hspace=0)
plt.subplots_adjust(wspace=0)
ax_ds.legend()
print(distance_2[25000]/r_J)
plt.show()
#plt.savefig('comparison.png')

    
    
