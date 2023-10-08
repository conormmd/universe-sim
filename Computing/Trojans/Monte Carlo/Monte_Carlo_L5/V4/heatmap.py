import numpy as np
import matplotlib.pyplot as plt
from functions import *
from labellines import labelLine, labelLines

plt.rc('font', family='serif')
plt.style.use('dark_background')
plt.rcParams.update({'font.size': 12})

x = np.load('x.npy')
y = np.load('y.npy')
tracker = np.load('tracker.npy')
G = 6.6726e-11;d = 7.7834082e11;J_mass = 1.898130e27;S_mass = 1.98847e30;r_S = d*(J_mass/(J_mass+S_mass));r_J = d*(S_mass/(J_mass+S_mass));P = np.sqrt((d**3/(G*(S_mass+J_mass))))*2*np.pi;V_S = (2*np.pi*r_S)/P;V_J = (2*np.pi*r_J)/P;theta = 2*np.pi;omega = theta/P 

out = x,y,tracker


plt.figure(figsize = (6.5,13))

#contourlevelmanual = [-260000000, -259000000, -258000000, -257000000, -256500000, -256100000, -255900000]
#contourlevelmanual = [-256500000, -256100000, -255900000]
contourlevelmanual = [-256130000, -255940000]
x0 = -0.9e12
x1 = 0.9e12
y0 = -0.9e12
y1 = 0.9e12

x=out[0]
y=out[1]
Z=out[2]

#plt.rcParams['contour.negative_linestyle'] = 'solid'
#plt.contour(x,y,Z, levels=contourlevelmanual,colors='red')
total=np.load('total.npy')
hit=np.load('hit.npy')
escape=np.load('escape.npy')

for i in range(0,len(total)):
    plt.scatter(total[i][0],total[i][1],color='red',s=3)
for i in range(0,len(hit)):
    plt.scatter(hit[i][0],hit[i][1],color='lime',s=4,zorder=1)
for i in range(0,len(escape)):
    plt.scatter(escape[i][0],escape[i][1],color='white',s=3.5)
plt.scatter(0,0,color='#FDB813',s=250)
plt.scatter(0,r_J,color='#e36e4b',s=70)

plt.scatter(1e18,1e18,color='red',s=10, label='Tadpole')
plt.scatter(1e18,1e18,color='lime',s=10, label='Horseshoe')
plt.scatter(1e18,1e18,color='white',s=10, label='Unstable')

G = 6.6726e-11;d = 7.7834082e11;J_mass = 1.898130e27;S_mass = 1.98847e30;r_S = d*(J_mass/(J_mass+S_mass));r_J = d*(S_mass/(J_mass+S_mass));P = np.sqrt((d**3/(G*(S_mass+J_mass))))*2*np.pi;V_S = (2*np.pi*r_S)/P;V_J = (2*np.pi*r_J)/P;theta = 2*np.pi;omega = theta/P 
hill_rad = r_J*((J_mass/(3*S_mass))**(1/3))
tad_top = np.array([3.82e11,6.75e11])
tad_bot = np.array([6.72e11,-3.88e11])
hs_top = np.array([2.2e11,7.5e11])
segments = 100


x_rad_up=[]
y_rad_up=[]
x_rad_down=[]
y_rad_down=[]
x_hrad = []
y_hrad = []
for i in range(0,segments):
    theta = np.pi/2 - i*np.pi/segments
    x_rad_up.append(r_J*0.95*np.cos(theta))
    y_rad_up.append(r_J*0.95*np.sin(theta))
    x_rad_down.append(r_J*1.05*np.cos(theta))
    y_rad_down.append(r_J*1.05*np.sin(theta))
    x_hrad.append(hill_rad*np.cos(theta))
    y_hrad.append(r_J+hill_rad*np.sin(theta))

#jup_rad_up = plt.plot(x_rad_up,y_rad_up,color='b', label='test2')
#jup_rad_down = plt.plot(x_rad_down,y_rad_down,color='b',label='test')
#labelLines(plt.gca().get_lines(), zorder=2.5)
#jup_hrad = plt.plot(x_hrad, y_hrad,color='r')

#tadpole_top = plt.plot([0,tad_top[0]],[0,tad_top[1]],color='gray',linestyle='--')
#tadpole_bot = plt.plot([0,tad_bot[0]],[0,tad_bot[1]],color='gray',linestyle='--')
#horse_top = plt.plot([0,hs_top[0]],[0,hs_top[1]],color='gray',linestyle='--') 


plt.xticks([0,0.3e12,0.6e12,r_J],['0.0','0.3','0.6','${R}_{J}$'])
#xticklabels(['0.0','0.3','0.6','{R}_{J}'])

plt.yticks([-r_J,-0.6e12,-0.3e12,0,0.3e12,0.6e12,r_J],['${-R}_{J}$','-0.6','-0.3','0.0','0.3','0.6','${R}_{J}$'])
#plt.yticklabels(['{-R}_{J}','-0.6','-0.3','0.0','0.3','0.6','{R}_{J}'])

plt.xlim(0,0.83e12)
plt.ylim(-0.83e12,0.83e12)
plt.ylabel('Y Distance from Barycentre/m$\cdot{10}^{12}$', fontsize=16)
plt.xlabel('X Distance from Barycentre/m$\cdot{10}^{12}$', fontsize=16)
plt.legend(fontsize=13)

plt.savefig('MC.png')
plt.show()
