import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from progressbar import ProgressBar
from functions import posVector, modVector, gravAcceleration, Grav_Sim, asteroid_generator

asteroids = np.load('asteroids.npy')
planets = np.load('planets.npy')

N_asteroids = len(asteroids)
ticks = len(asteroids[0])

pbar=ProgressBar()
plt.plot(0,0)
plt.plot(0,0)
for j in pbar(range(0,N_asteroids)):
    distance = []
    tick = []
    for i in range(0,ticks-1):
        tick.append(i)
        dist_vector = asteroids[j][i] - planets[1][i]
        dist = modVector(dist_vector)
        distance.append(dist)
    plt.plot(tick, distance, label=str(j))
    plt.ylim(5.5e11,7e11)
    plt.legend()
    plt.show()


