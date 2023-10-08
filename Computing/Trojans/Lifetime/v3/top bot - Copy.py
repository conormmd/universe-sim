import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
#Params
G = 6.6726e-11; d = 7.7834082e11; J_mass = 1.898130e27; S_mass = 1.98847e30; r_S = d*(J_mass/(J_mass+S_mass)); r_J = d*(S_mass/(J_mass+S_mass)); P = np.sqrt((d**3/(G*(S_mass+J_mass))))*2*np.pi; V_S = (2*np.pi*r_S)/P; V_J = (2*np.pi*r_J)/P
rad_hill = (J_mass/(3*S_mass))**(1/3)
x=np.load('x.npy')
y=np.load('y.npy')

v_x = np.load('v_x.npy')
v_y = np.load('v_y.npy')+1
N=len(v_x)
vx = np.zeros(300000)
vy = np.zeros(300000)
plt.rc('font', family='serif')
#plt.style.use('dark_background')
plt.rcParams.update({'font.size': 15})

for i in range(0,300000):
    j = 100*i
    if j == 30000000:
        vx[i]=v_x[-1]
        vy[i]=v_x[-1]
    else:
        vx[i]=v_x[j]
        vy[i]=v_y[j]
#vx = np.abs(vx)
#vy=np.abs(vy)

track_x = [0,200,1650,3750,7885,9450,11222,14465,15343,16716,18763,19450,22718,27786,28587,29473,30973,33460,38060,38747,40189,41549,44873,46603,49490]
track_y = [0.42,0.362,0.352,0.393,0.3898,0.3937,0.3712,0.3739,0.3961,0.4020,0.3925,0.399,0.3973,0.3747,0.411,0.385,0.389,0.378,0.3935,0.41133,0.3717,0.369,0.329064,0.357,0.375]

fig = plt.figure(figsize=(12,6))
ax_bot = plt.subplot(2,4,(5,8))
ax_bot.plot(x,y)
ax_bot.plot([0,50000],[5*rad_hill,5*rad_hill],color='red')
plt.xlim(0,50000)
plt.ylim(0.325,1)
plt.xlabel('Time/${10}^{5}$yr')
plt.ylabel('Jupiter Distance/${R}_{J}$')
ax_bot.set_xticks([0,10000,20000,30000,40000,50000])
ax_bot.set_xticklabels(['0','0.1','0.2','0.3','0.4','0.5'])
ax_bot.set_yticks([0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.50])
ax_bot.set_yticklabels(['','0.35','','0.40','','0.45'])

text = ax_bot.text(0.5,0.2,'5${R}_{H}$',horizontalalignment='left',verticalalignment='top', transform=ax_bot.transAxes,color='red')

ax_top = plt.subplot(2,4,(1,4))
ax_top.plot(vx,vy,color='red')
plt.xlim(0,50000)
plt.ylim(1,1.06)
plt.ylabel('Asteroid Velocity/${V}_{J}$')
ax_top.set_xticks([])
ax_top.set_yticks([1,1.01,1.02,1.03,1.04,1.05,1.06])
ax_top.set_yticklabels(['','','1.02','','1.04','','1.06'])

plt.subplots_adjust(hspace=0)
plt.subplots_adjust(wspace=0)

#plt.savefig('Lifetime_Light.png')

plt.show()
