import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def posVector(A,B):
    A_x, A_y = A
    B_x, B_y = B
    x = A_x - B_x
    y = A_y - B_y
    posVector = np.array([x, y])
    return(posVector);
    
def modVector(posVector):
    x, y = posVector
    mod = (x**2 + y**2)**0.5
    return(mod);

def gravAcceleration(mass, posVector, modVector):
    acc = -1*G*mass*posVector/(modVector)**3
    return(acc);

def gravPotential(mass, modVector):
    G = 6.6726e-11
    V = -1*G*mass/modVector
    return(V);
def potentialHeatmap(bodies, x0, x1, y0, y1, counts):
    G = 6.6726e-11;d = 7.7834082e11;J_mass = 1.898130e27;S_mass = 1.98847e30;r_S = d*(J_mass/(J_mass+S_mass));r_J = d*(S_mass/(J_mass+S_mass));P = np.sqrt((d**3/(G*(S_mass+J_mass))))*2*np.pi;V_S = (2*np.pi*r_S)/P;V_J = (2*np.pi*r_J)/P;theta = 2*np.pi;omega = theta/P 
    G = 6.6726e-11
    N_bodies = len(bodies)
    dimension = len(bodies[0][1])
    position = np.ndarray([N_bodies,dimension])
    angle = 2*np.pi
    ang_vel = angle/P
   
    radiusx = x1-x0
    radiusy = y1-y0
    
    stepx = radiusx/counts
    stepy = radiusy/counts
    
    tracker = np.zeros([counts, counts])
    xaxis = np.zeros([counts])
    yaxis = np.zeros([counts])
    for body in range(0,N_bodies):
        position[body] = bodies[body][1]
    
    for x in range(0, counts):
        for y in range(0, counts):
            for body in range(0, N_bodies):
                pos_x = x0 + x*stepx
                pos_y = y0 + y*stepy
                xaxis[x] = pos_x
                yaxis[y] = pos_y
                current_pos = np.array([pos_x,pos_y])
                source_mass = bodies[body][0][1]
                
                pos_source_posVector = posVector(current_pos, position[body])
                pos_source_modVector = modVector(pos_source_posVector)
                gravPot = gravPotential(source_mass, pos_source_modVector)
                
                tracker[x][y] = tracker[x][y] + gravPot
            tracker[x][y] = tracker[x][y] - 0.5*ang_vel**2*modVector(current_pos)**2
    return(xaxis, yaxis, tracker)
