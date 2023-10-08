import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from progressbar import ProgressBar
from functions import posVector, modVector, gravAcceleration, Grav_Sim, asteroid_generator

asteroids = np.load('asteroids.npy')
planets = np.load('planets.npy')
params = np.load('params.npy')
ticks=params[0],a=params[1]
G = 6.6726e-11; d = 7.7834082e11; J_mass = 1.898130e27; S_mass = 1.98847e30; r_S = d*(J_mass/(J_mass+S_mass)); r_J = d*(S_mass/(J_mass+S_mass)); P = np.sqrt((d**3/(G*(S_mass+J_mass))))*2*np.pi; V_S = (2*np.pi*r_S)/P; V_J = (2*np.pi*r_J)/P
M_1 = S_mass,M_2 = J_mass,R = d
L4 = np.array([(R/2)*((M_1-M_2)/(M_1+M_2)),(np.sqrt(3)/2)*R])
L5 = np.array([-1*(R/2)*((M_1-M_2)/(M_1+M_2)),(np.sqrt(3)/2)*R])

def comoving_frame_evolve(position, ticks, a):
    pbar=ProgressBar()
    
    G = 6.6726e-11
    P = 374386630
    omega = 2*np.pi/P
    
    position_track = np.ndarray([ticks,2])
    position_track[0] = position

    for t in range(0, ticks-1):
        theta = (t+1)*a*omega
        x_0, y_0 = position_track[t]
        radius = modVector(position)
        x_1 = -1*radius*np.cos(theta)
        y_1 = -1*radius*np.sin(theta)

        position_track[t+1] = x_1, y_1

    return(position_track)

pos_track = comoving_frame_evolve(L4, ticks, a)

