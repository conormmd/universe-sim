import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from progressbar import ProgressBar
from functions import posVector, modVector, gravAcceleration, Grav_Sim, asteroid_generator

#Params
G = 6.6726e-11; d = 7.7834082e11; J_mass = 1.898130e27; S_mass = 1.98847e30; r_S = d*(J_mass/(J_mass+S_mass)); r_J = d*(S_mass/(J_mass+S_mass)); P = np.sqrt((d**3/(G*(S_mass+J_mass))))*2*np.pi; V_S = (2*np.pi*r_S)/P; V_J = (2*np.pi*r_J)/P

#Planet intials
planet_bodies = np.array([
    [[1, S_mass], [0, -r_S], [-V_S, 0]],
    [[2, J_mass], [0, r_J], [V_J, 0]],
])

#Asteroid gen params
N_asteroids = 10;r0 = 779161721816.5127-800949230.5616066;r1=779161721816.5127+800949230.5616066;theta0 =1.1;theta1=1.2

#Simulation params
ticks = 75000;a = 100000
params = np.array([ticks,a])

#Simulations
bodies, ast_inits = asteroid_generator(N_asteroids, P, r0, r1, theta0, theta1)
asteroids, planets = Grav_Sim(planet_bodies, bodies, ticks, a)

np.save('asteroids', asteroids)
np.save('planets', planets)
np.save('params', params)
np.save('ast_inits', ast_inits)


