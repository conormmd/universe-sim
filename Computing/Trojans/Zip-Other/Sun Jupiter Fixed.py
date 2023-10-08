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

def Grav_Sim(planet_bodies, bodies, ticks, a):
    
    
    G = 6.6726e-11
    P = 374386630
    omega = 2*np.pi/P
    
    N_planet_bodies = len(planet_bodies)
    dimension = len(planet_bodies[0][1])
    
    planet_position = np.ndarray([N_planet_bodies,ticks,dimension])
    planet_position[0][0] = planet_bodies[0][1]
    planet_position[1][0] = planet_bodies[1][1]
    for t in range(0, ticks-1):
        theta = (t+1)*a*omega
        x_0, y_0 = planet_position[0][t]
        radius = modVector(planet_position[0][0])
        x_1 = -1*radius*np.cos(theta)
        y_1 = -1*radius*np.sin(theta)

        planet_position[0][t+1] = x_1, y_1
        
        x_0, y_0 = planet_position[1][t]
        radius = modVector(planet_position[1][0])
        x_1 = radius*np.cos(theta)
        y_1 = radius*np.sin(theta)

        planet_position[1][t+1] = x_1, y_1
    
    #Establishing parameters
    G = 6.6726e-11
    N_bodies = len(bodies)
    dimension = len(bodies[0][1])
    
    #Establishing tracking array (1->ticks)
    position = np.ndarray([N_bodies,ticks,dimension])
    velocity = np.zeros([N_bodies,dimension])
    
    for body in range(0,N_bodies):
        position[body][0] = bodies[body][1]
        velocity[body] = bodies[body][2]
    
    for i in range(0,ticks-1):
        acceleration = np.zeros([N_bodies,dimension])
        #Calculating total acceleration for each body
        for target in range(0, N_bodies):
            for source in range(0, N_bodies):
                    source_mass = planet_bodies[source][0][1]
                    target_source_posVector = posVector(position[target][i], planet_position[source][i])
                    target_source_modVector = modVector(target_source_posVector)
                    target_accel = gravAcceleration(source_mass, target_source_posVector, target_source_modVector)
                    acceleration[target] = [acceleration[target][0]+target_accel[0],acceleration[target][1]+target_accel[1]]
        
        
        
        
        for target in range(0,N_bodies):
            x_n, y_n = position[target][i]
            Vx_n, Vy_n = velocity[target]
            Ax_n, Ay_n = acceleration[target]
            r_n = np.array([x_n, y_n])
            V_n = np.array([Vx_n, Vy_n])
            A_n = np.array([Ax_n, Ay_n])
            
            z1x = r_n[0] + 0.5*a*V_n[0]
            Vz1x = V_n[0] + 0.5*a*A_n[0]
            z1y = r_n[1] + 0.5*a*V_n[1]
            Vz1y = V_n[1] + 0.5*a*A_n[1]
            
            z1 = np.array([z1x, z1y])
            Vz1 = np.array([Vz1x, Vz1y])
            Az1 = np.array([0, 0])
            
            for source in range(0, N_bodies):
                    source_mass = planet_bodies[source][0][1]
                    target_source_posVector = posVector(z1, planet_position[source][i])
                    target_source_modVector = modVector(target_source_posVector)
                    target_accel = gravAcceleration(source_mass, target_source_posVector, target_source_modVector)
                    Az1 = np.array([Az1[0]+target_accel[0],Az1[1]+target_accel[1]])
            
            z2x = r_n[0] + 0.5*a*Vz1[0]
            Vz2x = V_n[0] + 0.5*a*Az1[0]
            z2y = r_n[1] + 0.5*a*Vz1[1]
            Vz2y = V_n[1] + 0.5*a*Az1[1]
            
            z2 = np.array([z2x, z2y])
            Vz2 = np.array([Vz2x, Vz2y])
            Az2 = np.array([0, 0])
            
            for source in range(0, N_bodies):
                    source_mass = planet_bodies[source][0][1]
                    target_source_posVector = posVector(z2, planet_position[source][i])
                    target_source_modVector = modVector(target_source_posVector)
                    target_accel = gravAcceleration(source_mass, target_source_posVector, target_source_modVector)
                    Az2 = np.array([Az2[0]+target_accel[0],Az2[1]+target_accel[1]]) 
            
            z3x = r_n[0] + a*Vz2[0]
            Vz3x = V_n[0] + a*Az2[0]
            z3y = r_n[1] + a*Vz2[1]
            Vz3y = V_n[1] + a*Az2[1]
            
            z3 = np.array([z3x, z3y])
            Vz3 = np.array([Vz3x, Vz3y])
            Az3 = np.array([0, 0])
            
            for source in range(0, N_bodies):
                    source_mass = planet_bodies[source][0][1]
                    target_source_posVector = posVector(z3, planet_position[source][i])
                    target_source_modVector = modVector(target_source_posVector)
                    target_accel = gravAcceleration(source_mass, target_source_posVector, target_source_modVector)
                    Az3 = np.array([Az3[0]+target_accel[0],Az3[1]+target_accel[1]])
            
            r_n1 = r_n + (a/6)*(V_n + 2*Vz1 + 2*Vz2 + Vz3)
            V_n1 = V_n + (a/6)*(A_n + 2*Az1 + 2*Az2 + Az3)
            
            velocity[target] = V_n1
            position[target][i+1] = r_n1
    
    
    
    
    return(position, planet_position)

G = 6.6726e-11
d = 7.7834082e11

J_mass = 1.898130e27
S_mass = 1.98847e30

r_S = d*(J_mass/(J_mass+S_mass))
r_J = d*(S_mass/(J_mass+S_mass))

P = np.sqrt((d**3/(G*(S_mass+J_mass))))*2*np.pi

V_S = (2*np.pi*r_S)/P
V_J = (2*np.pi*r_J)/P

theta = 2*np.pi
omega = theta/P
##################################
planet_bodies = np.array([
    [[1, S_mass], [-r_S, 0], [0, -V_S]],
    [[2, J_mass], [r_J, 0], [0, V_J]],
])

bodies = np.array([
     [[1, 0], [r_J*np.cos(np.pi/3), r_J*np.sin(np.pi/3)], [-1*V_J*np.sin(np.pi/3), V_J*np.cos(np.pi/3)]],
])

ticks = 374500
a = 1000
#################################
asteroids, planets = Grav_Sim(planet_bodies, bodies, ticks, a)

distance = []
tick = []
for i in range(0,ticks-1):
    tick.append(i)
    dist_vector = asteroids[0][i] - planets[1][i]
    dist = modVector(dist_vector)
    distance.append(dist)

plt.plot(tick, distance)
plt.show()


plt.figure(figsize = (8,8))

plt.plot(planets[0][:,0], planets[0][:,1], label='Sun')
plt.plot(planets[1][:,0], planets[1][:,1], label='Jupiter')
plt.plot(asteroids[0][:,0], asteroids[0][:,1], label='Asteroid')


plt.ylabel("Y axis")
plt.xlabel("X axis")

plt.legend()

plt.show()
