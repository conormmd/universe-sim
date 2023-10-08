import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from progressbar import ProgressBar
from functions import posVector, modVector, gravAcceleration, Grav_Sim, asteroid_generator

asteroids = np.load('asteroids.npy')
planets = np.load('planets.npy')

N_asteroids = len(asteroids)
ticks = len(asteroids[0])
# Good ish 1,11,14,15,16,17,18,25,26,27,28,29,33,34,35
good_asteroids = [37,38,40,42,48,74,92,95]
#37 VERY GOOD
good_init = np.ndarray([len(good_asteroids),2,1])

for i in range(0,len(good_asteroids)):
    good_init[i][0] = asteroids[good_asteroids[i]][0][0]
    good_init[i][1] = asteroids[good_asteroids[i]][0][1]

r_init = np.ndarray([len(good_asteroids),1])

for i in range(0,len(good_asteroids)):
    r_init[i] = np.sqrt(good_init[i][0]**2 + good_init[i][1]**2)

theta_init = np.ndarray([len(good_asteroids),1])

for i in range(0,len(good_asteroids)):
    theta_init[i] = np.arctan(good_init[i][0]/good_init[i][1])

print(np.mean(r_init))
print(np.std(r_init))

print(np.mean(theta_init))
print(np.std(theta_init))
               

