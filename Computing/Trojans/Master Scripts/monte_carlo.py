import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from progressbar import ProgressBar
from functions import *
asteroids = np.load('asteroids.npy')
planets = np.load('planets.npy')
params = np.load('params.npy')
ast_inits = np.load('ast_inits.npy')

num_ast = len(asteroids)
ticks = params[0]; a= params[1]

G = 6.6726e-11; d = 7.7834082e11; J_mass = 1.898130e27; S_mass = 1.98847e30; r_S = d*(J_mass/(J_mass+S_mass)); r_J = d*(S_mass/(J_mass+S_mass)); P = np.sqrt((d**3/(G*(S_mass+J_mass))))*2*np.pi; V_S = (2*np.pi*r_S)/P; V_J = (2*np.pi*r_J)/P
M_1 = S_mass;M_2 = J_mass;R = d
L4 = np.array([(R/2)*((M_1-M_2)/(M_1+M_2)),(np.sqrt(3)/2)*R])
L5 = np.array([-1*(R/2)*((M_1-M_2)/(M_1+M_2)),(np.sqrt(3)/2)*R])

jupiter = planets[1][0]

test = counter_rotate_monte_carlo(asteroids,ticks,a)

for i in range(0,len(asteroids)):
    asteroid=asteroids[i]
    plt.scatter(asteroid[0][0],asteroid[0][1],color='black',s=3)
for i in test:
    asteroid=asteroids[i]
    plt.scatter(asteroid[0][0],asteroid[0][1],color='red',s=3)

plt.xlim(0,1e12)
plt.ylim(0,1e12)
plt.legend()
plt.show()
    
