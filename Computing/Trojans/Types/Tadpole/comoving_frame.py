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
jupiter = planets[1][0]
asteroid = asteroids[0]
asteroid_init = asteroids[0][0]

pos_track_2 = counter_rotate(asteroid, ticks, a)

np.save('comoving', pos_track_2)
np.save('comoving_ticks', ticks)

print(pos_track_2[0])


plt.figure(figsize=(8,8))
#plt.scatter(0,0,color='Red', label='Sun')
plt.scatter(jupiter[0], jupiter[1], color='black', label='Jupiter')
plt.scatter(r_J*np.sin(np.pi/3), r_J*np.cos(np.pi/3), color='red',label='L5')
plt.plot(pos_track_2[:,0], pos_track_2[:,1], color='blue', label='Asteroid Path')
plt.xlim(-0.9e12,0.9e12)
plt.ylim(-0.9e12,0.9e12)
plt.legend()
plt.savefig('comoving.png')
plt.show()
    
    
