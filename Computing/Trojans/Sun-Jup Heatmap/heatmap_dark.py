import numpy as np
import matplotlib.pyplot as plt
from functions import *
plt.rc('font', family='serif')
plt.style.use('dark_background')

x = np.load('x.npy')
y = np.load('y.npy')
tracker = np.load('tracker.npy')
G = 6.6726e-11;d = 7.7834082e11;J_mass = 1.898130e27;S_mass = 1.98847e30;r_S = d*(J_mass/(J_mass+S_mass));r_J = d*(S_mass/(J_mass+S_mass));P = np.sqrt((d**3/(G*(S_mass+J_mass))))*2*np.pi;V_S = (2*np.pi*r_S)/P;V_J = (2*np.pi*r_J)/P;theta = 2*np.pi;omega = theta/P 

out = x,y,tracker


plt.figure(figsize = (10,10))

#contourlevelmanual = [-260000000, -259000000, -258000000, -257000000, -256500000, -256100000, -255900000]
contourlevelmanual = [-260000000, -258000000, -256500000, -255950000]
x0 = -0.9e12
x1 = 0.9e12
y0 = -0.9e12
y1 = 0.9e12

ax1 = plt.subplot(4,4,(2,12))

x=out[0]
y=out[1]
Z=out[2]
plt.style.use('dark_background')
heatmap = plt.pcolormesh(x, y, Z, vmin=-265000000, vmax=-255900000, cmap='gray')
plt.contour(x,y,Z, levels=contourlevelmanual, cmap = 'binary')
#plt.colorbar(heatmap)

plt.scatter(0,0,color='#FDB813',s=250)
plt.scatter(0,r_J,color='#e36e4b',s=70)

##########################################################

xx0 = [0,0]
xy0 = [0.835e12,-r_J]

plt.plot(xx0, xy0, linestyle='--', color='deepskyblue')

xL2 = [-0.9e12,0]
yL2 = [0.835e12, 0.835e12]

plt.plot(xL2, yL2, linestyle='--', color='deepskyblue')

xL3 = [-0.9e12,0]
yL3 = [-r_J,-r_J]

plt.plot(xL3, yL3, linestyle='--', color='deepskyblue')

xL1 = [-0.9e12,0]
yL1 = [0.725e12,0.725e12]

plt.plot(xL1, yL1, linestyle='--', color='deepskyblue')

##########################################################

yx0 = [-r_J*np.sin(np.pi/3), r_J*np.sin(np.pi/3)]
yy0 = [r_J*np.cos(np.pi/3), r_J*np.cos(np.pi/3)]

plt.plot(yx0, yy0, linestyle='--', color='red')

xL4 = [-r_J*np.sin(np.pi/3), -r_J*np.sin(np.pi/3)]
yL4 = [-0.9e12, r_J*np.cos(np.pi/3)]

plt.plot(xL4, yL4, linestyle='--', color='red')

xL5 = [r_J*np.sin(np.pi/3), r_J*np.sin(np.pi/3)]
yL5 = [-0.9e12, r_J*np.cos(np.pi/3)]

plt.plot(xL5, yL5, linestyle='--', color='red')

plt.ylim(y0, y1)
plt.xlim(x0, x1)

ax1.set_yticks([])
ax1.set_yticklabels([])
#########################################

ax2 = plt.subplot(4,4,(1,9))

y2 = out[0]
x2 = out[2][:,500]

plt.plot(x2, y2, color='white')

yL2 = [0.835e12, 0.835e12]
xL2 = [-1e10,0]

plt.plot(xL2, yL2, linestyle='--', color='deepskyblue')

yL1 = [0.725e12,0.725e12]
xL1 = [-1e10,0]

plt.plot(xL1, yL1, linestyle='--', color='deepskyblue')

yL3 = [-r_J,-r_J]
xL3 = [-1e10,0]

plt.plot(xL3, yL3, linestyle='--', color='deepskyblue')

plt.xlim(-5e8,-2.5e8)
plt.ylim(y0, y1)

ax2.xaxis.tick_top()
ax2.xaxis.set_label_position('top')

ax2.set_yticks([-0.9e12, -0.6e12, -0.3e12, 0, 0.3e12, 0.6e12,0.9e12,0.725e12,0.835e12,-r_J])
ax2.set_yticklabels(['-0.9', '-0.6','-0.3','0','0.3','0.6','0.9','L1','L2','L3'])

ax2.set_xticks([-5e8,-4e8,-3e8])
ax2.set_xticklabels(['-5','-4','-3'])

plt.ylabel('Y Distance from Barycentre/m$\cdot{10}^{12}$', fontsize=14)
plt.xlabel('$\Phi$/Jkg${}^{-1}\cdot{10}^{-8}$', fontsize=14)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
#########################################

ax3 = plt.subplot(4,4,(14,16))

x3 = out[1]
y3 = out[2][687,:]

plt.plot(x3, y3, color='white')

xL4 = [-r_J*np.sin(np.pi/3), -r_J*np.sin(np.pi/3)]
yL4 = [-1e10,0]

plt.plot(xL4, yL4, linestyle='--', color='red')

xL5 = [r_J*np.sin(np.pi/3), r_J*np.sin(np.pi/3)]
yL5 = [-1e10,0]

plt.plot(xL5, yL5, linestyle='--', color='red')

plt.xlim(x0, x1)
plt.ylim(-4.25e8, -2.5e8)

ax3.set_xticks([-0.9e12, -0.6e12, -0.3e12, 0, 0.3e12, 0.6e12,0.9e12, -r_J*np.sin(np.pi/3), r_J*np.sin(np.pi/3)])
ax3.set_xticklabels(['-0.9', '-0.6','-0.3','0','0.3','0.6','0.9','L4 ',' L5'])

ax3.set_yticks([-4e8,-3.5e8,-3e8,-2.5e8])
ax3.set_yticklabels(['-4','-3.5','-3','-2.5'])

ax3.yaxis.tick_right()
ax3.yaxis.set_label_position('right')

plt.xlabel('X Distance from Barycentre/m$\cdot{10}^{12}$', fontsize=14)
plt.ylabel('$\Phi$/Jkg${}^{-1}\cdot{10}^{-8}$', fontsize=14)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax1.get_yticklabels(), visible=False)

plt.subplots_adjust(hspace=0)
plt.subplots_adjust(wspace=0)

plt.savefig('heatmap_test.png')
