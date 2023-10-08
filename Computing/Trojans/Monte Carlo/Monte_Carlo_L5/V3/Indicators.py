import numpy as np
import matplotlib.pyplot as plt
from functions import *


G = 6.6726e-11;d = 7.7834082e11;J_mass = 1.898130e27;S_mass = 1.98847e30;r_S = d*(J_mass/(J_mass+S_mass));r_J = d*(S_mass/(J_mass+S_mass));P = np.sqrt((d**3/(G*(S_mass+J_mass))))*2*np.pi;V_S = (2*np.pi*r_S)/P;V_J = (2*np.pi*r_J)/P;theta = 2*np.pi;omega = theta/P 
hill_rad = r_J*((J_mass/(3*S_mass))**(1/3))
plt.figure(figsize = (5,10))
tad_top = np.array([3.82e11,6.75e11]) #1.055808 rad
tad_bot = np.array([6.72e11,-3.88e11])
hs_top = np.array([2.2e11,7.5e11]) #1.2854669 rad
segments = 100
x_rad=[]
y_rad=[]
x_hrad = []
y_hrad = []
for i in range(0,segments):
    theta = np.pi/2 - i*np.pi/segments
    x_rad.append(r_J*np.cos(theta))
    y_rad.append(r_J*np.sin(theta))
    x_hrad.append(hill_rad*np.cos(theta))
    y_hrad.append(r_J+hill_rad*np.sin(theta))

#jup_radius = plt.plot(x_rad,y_rad,color='r')
#jup_hrad = plt.plot(x_hrad, y_hrad,color='r')

#tadpole_top = plt.plot([0,tad_top[0]],[0,tad_top[1]],color='gray',linestyle='--')
#tadpole_bot = plt.plot([0,tad_bot[0]],[0,tad_bot[1]],color='gray',linestyle='--')
#horse_top = plt.plot([0,hs_top[0]],[0,hs_top[1]],color='gray',linestyle='--') 

#plt.xlim(0,1e12)
#plt.ylim(-1e12,1e12)

plt.show()
