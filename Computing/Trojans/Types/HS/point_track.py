import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from progressbar import ProgressBar
from functions import *

asteroids = np.load('asteroids.npy')
planets = np.load('planets.npy')
params = np.load('params.npy')
ast_inits = np.load('ast_inits.npy')

ticks = params[0]; a= params[1]
theta_0 = ast_inits[1]; r_0 = ast_inits[0]
G = 6.6726e-11; d = 7.7834082e11; J_mass = 1.898130e27; S_mass = 1.98847e30; r_S = d*(J_mass/(J_mass+S_mass)); r_J = d*(S_mass/(J_mass+S_mass)); P = np.sqrt((d**3/(G*(S_mass+J_mass))))*2*np.pi; V_S = (2*np.pi*r_S)/P; V_J = (2*np.pi*r_J)/P
M_1 = S_mass;M_2 = J_mass;R = d

L4 = np.array([(R/2)*((M_1-M_2)/(M_1+M_2)),(np.sqrt(3)/2)*R])
L5 = np.array([-1*(R/2)*((M_1-M_2)/(M_1+M_2)),(np.sqrt(3)/2)*R])
asteroid = asteroids[0]
asteroid_init = asteroids[0][0]

pos_track = comoving_frame_evolve(np.array([L5[0],L5[1]]), ticks, a, theta_0)

comoving_point = np.ndarray([ticks,2])

for i in range(1, ticks):
    ast_x, ast_y = asteroids[0][i]
    pos_x, pos_y = pos_track[i]
    comoving_point[i][0] = (ast_x-pos_x)
    comoving_point[i][1] = (ast_y-pos_y)

#plt.plot(pos_track[:,0], pos_track[:,1])
#plt.scatter(asteroid[:,0], asteroid[:,1])
#plt.xlim(-0.8e12, 0.8e12)
#plt.ylim(-0.8e12, 0.8e12)
    
plt.plot(comoving_point[:,0], comoving_point[:,1])
plt.show()
    
    
