import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from progressbar import ProgressBar
from functions import posVector, modVector, gravAcceleration, Grav_Sim, asteroid_generator

asteroids = np.load('asteroids.npy')
planets = np.load('planets.npy')

N_asteroids = len(asteroids)
ticks = len(asteroids[0])

good_asteroids = [4,7,25,35,49,76,77]

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
               

