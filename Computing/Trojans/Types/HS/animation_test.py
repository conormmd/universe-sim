# importing relevant modules & packages
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import rc
G = 6.6726e-11; d = 7.7834082e11; J_mass = 1.898130e27; S_mass = 1.98847e30; r_S = d*(J_mass/(J_mass+S_mass)); r_J = d*(S_mass/(J_mass+S_mass)); P = np.sqrt((d**3/(G*(S_mass+J_mass))))*2*np.pi; V_S = (2*np.pi*r_S)/P; V_J = (2*np.pi*r_J)/P
tick, a = np.load('params.npy')

# set matplotlib figure size
fig = plt.figure(figsize=(8,8))

interval = 1
step = 50
# create subplot figure and axes handlers
ax_comoving = plt.subplot(6,6,(1,15))

plt.title('$Comoving$ $Frame$')
plt.xlim(-1e12,1e12)
plt.ylim(-1e12,1e12)

ax_comoving.set_yticks([])
ax_comoving.set_yticklabels([])
ax_comoving.set_xticks([])
ax_comoving.set_xticklabels([])

comoving = np.load('comoving.npy')
x_comoving = []
y_comoving = []

line_comoving, = ax_comoving.plot(comoving[0][0],comoving[0][1], 
                color='blue', 
                linestyle='-', 
                linewidth=2, 
                markersize=2, label='Asteroid')

sun_comoving = ax_comoving.plot(0,0,'o',color='red',markersize=7,label='Sun')
sun_legend = plt.legend(handles=sun_comoving, loc=[0.003,0.915])
ax_comoving.plot(0,r_J,'o',color='black',markersize=5)

def frameAnimationComoving(i):

    x_comoving.append(comoving[i][0])
    y_comoving.append(comoving[i][1])

    line_comoving.set_xdata(x_comoving)
    line_comoving.set_ydata(y_comoving)

    return line_comoving,

animationComoving = FuncAnimation(fig, # the figure to assign animation too
                            func = frameAnimationComoving, # the frame rendering function
                            frames = np.arange(0,len(comoving),step), # the steps and amount of frames
                            interval = interval) # invertals is the time per frame, in ms
#########################################################
ax_moving = plt.subplot(6,6,(4,18))

plt.title('$Stationary$ $Frame$')
plt.xlim(-1e12, 1e12)
plt.ylim(-1e12, 1e12)

ax_moving.set_yticks([])
ax_moving.set_yticklabels([])
ax_moving.set_xticks([])
ax_moving.set_xticklabels([])

planets = np.load('planets.npy')
asteroids = np.load('asteroids.npy')

sun = planets[0]
x_sun = []
y_sun = []
line_moving_sun, = ax_moving.plot(sun[0][0],sun[0][1], 
                'o',
                color='red',
                markersize=7, label='Sun')
def frameAnimationMoving_Sun(i):
    
    x_sun.append(sun[i][0])
    y_sun.append(sun[i][1])

    line_moving_sun.set_xdata(x_sun)
    line_moving_sun.set_ydata(y_sun)
    
    return line_moving_sun,
animationMoving_Sun = FuncAnimation(fig, # the figure to assign animation too
                            func = frameAnimationMoving_Sun, # the frame rendering function
                            frames = np.arange(0,len(sun),step), # the steps and amount of frames
                            interval = interval) # invertals is the time per frame, in ms
jupiter = planets[1]
x_jupiter = []
y_jupiter = []
line_moving_jupiter, = ax_moving.plot(jupiter[0][0],jupiter[0][1], 
                'o',
                color='black',
                markersize=5, label='Jupiter')
def frameAnimationMoving_Jupiter(i):

    x_jupiter.append(jupiter[i][0])
    y_jupiter.append(jupiter[i][1])

    line_moving_jupiter.set_xdata(x_jupiter[-1])
    line_moving_jupiter.set_ydata(y_jupiter[-1])

    return line_moving_jupiter,

animationMoving_Jupiter = FuncAnimation(fig, # the figure to assign animation too
                            func = frameAnimationMoving_Jupiter, # the frame rendering function
                            frames = np.arange(0,len(sun),step), # the steps and amount of frames
                            interval = interval) # invertals is the time per frame, in ms

asteroid = asteroids[0]
x_asteroid = []
y_asteroid = []
line_moving_asteroid, = ax_moving.plot(comoving[0][0],comoving[0][1], 
                'o',
                color='blue',
                markersize=4, label='Asteroid')
def frameAnimationMoving_Asteroid(i):

    x_asteroid.append(asteroid[i][0])
    y_asteroid.append(asteroid[i][1])

    line_moving_asteroid.set_xdata(x_asteroid[-1])
    line_moving_asteroid.set_ydata(y_asteroid[-1])

    return line_moving_asteroid,
animationMoving_Asteroid = FuncAnimation(fig, # the figure to assign animation too
                            func = frameAnimationMoving_Asteroid, # the frame rendering function
                            frames = np.arange(0,len(sun),step), # the steps and amount of frames
                            interval = interval) # invertals is the time per frame, in ms


hidden_ast = ax_moving.plot(-1e17,-1e17, 'o',color='blue',markersize=4, label='Asteroid')
hidden_jup = ax_moving.plot(-1e17, -1e17, 'o', color='black', markersize=5, label='Jupiter')

ast_legend = plt.legend(handles=hidden_ast, loc=[0.65,0.915])
ax_moving = ax_moving.add_artist(ast_legend)

jup_legend = plt.legend(handles=hidden_jup, loc = [-0.135, 0.915])
###############################################################
ax_distance=plt.subplot(6,6,(19,30))

plt.xlim(0,100000)
plt.ylim(0,1.7e12)
plt.xlabel('$Time$')
plt.ylabel('$Distance$ $from$ $Jupiter$')

ax_distance.set_yticks([])
ax_distance.set_yticklabels([])
ax_distance.set_xticks([])
ax_distance.set_xticklabels([])

distance = np.load('distance.npy')
x_distance = []
y_distance = []

line_distance, = ax_distance.plot(0,distance[0], 
                color='blue', 
                linestyle='-', 
                linewidth=2, 
                markersize=2, label='Comoving Asteroid')

time_text = ax_distance.text(0.91,0.985,'',horizontalalignment='left',verticalalignment='top', transform=ax_distance.transAxes)
time_text.set_text('Time = 10yrs')

def frameAnimationDistance(i):

    if i ==len(comoving):
        plt.clear()
    
    x_distance.append(i)
    y_distance.append(distance[i])
    line_distance.set_xdata(x_distance)
    line_distance.set_ydata(y_distance)

    time_text.set_text('%.3d yrs'%(i*(100000/3.154e7)))
    return line_distance, time_text

animationDistance = FuncAnimation(fig, # the figure to assign animation too
                            func = frameAnimationDistance, # the frame rendering function
                            frames = np.arange(0,len(comoving),step), # the steps and amount of frames
                            interval = interval) # invertals is the time per frame, in ms

plt.subplots_adjust(hspace=0)
plt.subplots_adjust(wspace=0)

plt.show()





