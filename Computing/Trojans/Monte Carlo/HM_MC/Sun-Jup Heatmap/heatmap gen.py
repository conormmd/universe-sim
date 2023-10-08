import numpy as np
import matplotlib.pyplot as plt
from functions import *

G = 6.6726e-11;d = 7.7834082e11;J_mass = 1.898130e27;S_mass = 1.98847e30;r_S = d*(J_mass/(J_mass+S_mass));r_J = d*(S_mass/(J_mass+S_mass));P = np.sqrt((d**3/(G*(S_mass+J_mass))))*2*np.pi;V_S = (2*np.pi*r_S)/P;V_J = (2*np.pi*r_J)/P;theta = 2*np.pi;omega = theta/P 
bodies = np.array([
    [[1, S_mass], [-r_S, 0], [0, -V_S]],
    [[2, J_mass], [r_J, 0], [0, V_J]],
    #[[3, 0], [r_r, 0], [0, V_r]],
])

x0 = -0.9e12
x1 = 0.9e12
y0 = -0.9e12
y1 = 0.9e12

radiusx = x1-x0
radiusy = y1-y0
counts = 1000
stepy = radiusy/counts
stepx = radiusx/counts

x,y,tracker= potentialHeatmap(bodies, x0, x1, y0, y1, counts)

np.save('x',x)
np.save('y',y)
np.save('tracker',tracker)



