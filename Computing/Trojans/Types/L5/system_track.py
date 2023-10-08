import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from progressbar import ProgressBar
from functions import posVector, modVector, gravAcceleration, Grav_Sim, asteroid_generator

asteroids = np.load('asteroids.npy')
planets = np.load('planets.npy')

N_asteroids = len(asteroids)
ticks = len(asteroids[0])

plt.figure(figsize = (8,8))

pbar=ProgressBar()
plt.plot(planets[0][:,0], planets[0][:,1], label='Sun')
plt.plot(planets[1][:,0], planets[1][:,1], color='black',label='Jupiter')
for i in pbar(range(0,N_asteroids)):
    plt.plot(asteroids[i][:,0], asteroids[i][:,1])


plt.ylabel("Y axis")
plt.xlabel("X axis")

plt.xlim(-0.8e12, 0.8e12)
plt.ylim(-0.8e12, 0.8e12)

plt.legend()

plt.show()
