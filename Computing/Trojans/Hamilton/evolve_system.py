import numpy as np
import time
#Params
G = 6.6726e-11; d = 7.7834082e11; J_mass = 1.898130e27; S_mass = 1.98847e30; r_S = d*(J_mass/(J_mass+S_mass)); r_J = d*(S_mass/(J_mass+S_mass)); P = np.sqrt((d**3/(G*(S_mass+J_mass))))*2*np.pi; V_S = (2*np.pi*r_S)/P; V_J = (2*np.pi*r_J)/P
#Asteroid gen params
#N_asteroids = 1;r0=r_J;r1=r_J;theta0=1.1605;theta1=1.1605
N_asteroids = 1;r0=r_J;r1=r_J;theta0=-np.pi/2;theta1=-np.pi/2
#Sim params
ticks = 1000000;a = 5000
##
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
    G = 6.6726e-11
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
        theta = omega*a*(t+1) + np.pi/2
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
    
    #Establishing tracking array (1->ticks)
    position = np.ndarray([N_bodies,ticks,dimension])
    velocity = np.zeros([N_bodies,dimension])
    
    for body in range(0,N_bodies):
        position[body][0] = np.array([bodies[body][0][1],bodies[body][0][2]])
        velocity[body] = np.array([bodies[body][0][3],bodies[body][0][4]])
    
    for i in range(0,ticks-1):
        acceleration = np.zeros([N_bodies,dimension])
        #Calculating total acceleration for each body
        for target in range(0, N_bodies):
            for source in range(0, N_planet_bodies):
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
            
            for source in range(0, N_planet_bodies):
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
            
            for source in range(0, N_planet_bodies):
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
            
            for source in range(0, N_planet_bodies):
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

def asteroid_generator(N_asteroids, P, r0, r1, theta0, theta1):
    
    bodies = np.ndarray([N_asteroids,1,5])
    omega = 2*np.pi/P
    
    for i in range(0, N_asteroids):
        r = np.random.uniform(r0,r1)
        theta = np.random.uniform(theta0, theta1)
        
        v = omega*r
        
        x = r*np.cos(theta)
        y = r*np.sin(theta)
        
        r_vector = np.array([x/r,y/r,0])
        omega_vector = np.array([0,0,1])
        v_vector = np.cross(omega_vector, r_vector)
        
        v_x = np.array(v*v_vector[0])
        v_y = np.array(v*v_vector[1])
        
        bodies[i][0][0] = 0
        
        bodies[i][0][1] = x
        bodies[i][0][2] = y
        
        bodies[i][0][3] = v_x
        bodies[i][0][4] = v_y
        if i == 0:
            ast_inits = np.array([r, theta])
    
    return bodies, ast_inits;
##
#Planet intials
planet_bodies = np.array([
    [[1, S_mass], [0, -r_S], [-V_S, 0]],
    [[2, J_mass], [0, r_J], [V_J, 0]],
])

##
tick = time.time()
#Simulation params
params = np.array([ticks,a])
#Simulations
bodies, ast_inits = asteroid_generator(N_asteroids, P, r0, r1, theta0, theta1)
asteroids, planets = Grav_Sim(planet_bodies, bodies, ticks, a)
##
tock=time.time()
print(tock-tick)
##
ticks = int(ticks/10)
jup = np.ndarray([1,ticks,2])
ast = np.ndarray([1,ticks,2])
for i in range(0,ticks):
	j = 10*i
	jup[0][i]=planets[1][j]
	ast[0][i] = asteroids[0][j]
np.save('/ddn/data/lgxn55/Lifetime/asteroids.npy', ast)
np.save('/ddn/data/lgxn55/Lifetime/planets', jup)
np.save('/ddn/data/lgxn55/Lifetime/params', params)
np.save('/ddn/data/lgxn55/Lifetime/ast_inits', ast_inits)
##


