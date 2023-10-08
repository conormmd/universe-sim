import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from progressbar import ProgressBar
from functions import *
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix';matplotlib.rcParams['font.family'] = 'STIXGeneral';
plt.rcParams.update({'font.size': 14})
asteroids = np.load('asteroids.npy');planets = np.load('planets.npy');params = np.load('params.npy');ast_inits = np.load('ast_inits.npy')
num_ast = len(asteroids);ticks = params[0]; a= params[1]
G = 6.6726e-11; d = 7.7834082e11; J_mass = 1.898130e27; S_mass = 1.98847e30; r_S = d*(J_mass/(J_mass+S_mass)); r_J = d*(S_mass/(J_mass+S_mass)); P = np.sqrt((d**3/(G*(S_mass+J_mass))))*2*np.pi; V_S = (2*np.pi*r_S)/P; V_J = (2*np.pi*r_J)/P;M_1 = S_mass;M_2 = J_mass;R = d;L4 = np.array([(R/2)*((M_1-M_2)/(M_1+M_2)),(np.sqrt(3)/2)*R]);L5 = np.array([-1*(R/2)*((M_1-M_2)/(M_1+M_2)),(np.sqrt(3)/2)*R]);jupiter = planets[1][0]
hit,escape = counter_rotate_monte_carlo(asteroids,ticks,a)
plt.figure(figsize = (10,8))
ax_main=plt.subplot(5,6,(1,12))
for i in range(0,len(asteroids)):
    asteroid=asteroids[i]
    plt.scatter(asteroid[0][0],asteroid[0][1],color='black',s=20)
for i in hit:
    asteroid=asteroids[i]
    plt.scatter(asteroid[0][0],asteroid[0][1],color='green',s=20)
for i in escape:
    asteroid=asteroids[i]
    plt.scatter(asteroid[0][0],asteroid[0][1],color='red',s=20)
ax_main.set_yticks([]);ax_main.set_xticks([])
plt.ylim(7.25e11,7.5e11);plt.xlim(2.18e11,2.68e11)
r_J=777598549241.9435
def RadPlot(R):
    N=50
    theta=np.linspace(np.pi/2,0,N)
    x = np.zeros(N);y=np.zeros(N)
    for i in range(0,N):
        x[i]=np.sin(theta[i])*R
        y[i]=np.cos(theta[i])*R
    return(x,y)
def ThPlot(Th):
    N=50;r_J=777598549241.9435
    x=np.linspace(0,r_J*np.cos(Th),N);y=np.linspace(0,r_J*np.sin(Th),N)
    return(x,y)
x,y=ThPlot(1.25);plt.plot(x,y,color='grey',linestyle='dashed');x,y=ThPlot(1.225);plt.plot(x,y,color='grey',linestyle='dashed');x,y=ThPlot(1.275);plt.plot(x,y,color='grey',linestyle='dashed')
x,y=RadPlot(r_J*1.01);plt.plot(x,y,color='blue');x,y=RadPlot(r_J*0.99);plt.plot(x,y,color='blue');x,y=RadPlot(r_J);plt.plot(x,y,linestyle=(0,(5,10)),zorder=0,color='blue')
ax_esc=plt.subplot(5,6,(13,27))
esc_plt = asteroids[1];esc_comoving = counter_rotate(esc_plt,150000,100000);ax_esc.plot(esc_comoving[:,0],esc_comoving[:,1],color='black');ax_esc.plot(0,r_J,'o',color='#e36e4b',markersize=5);ax_esc.plot(0,0,'o',color='#FDB813',markersize=7)
plt.xlim(-1e12,1e12);plt.ylim(-1e12,1e12)
ax_esc.set_yticks([]);ax_esc.set_xticks([])

ax_hit=plt.subplot(5,6,(16,30))
hit_plt = asteroids[2];hit_comoving = counter_rotate(hit_plt,150000,100000);ax_hit.plot(hit_comoving[:,0],hit_comoving[:,1],color='black');ax_hit.plot(0,r_J,'o',color='#e36e4b',markersize=5);ax_hit.plot(0,0,'o',color='#FDB813',markersize=7)
plt.xlim(-1e12,1e12);plt.ylim(-1e12,1e12)
ax_hit.set_yticks([]);ax_hit.set_xticks([])

ax_main.annotate('R$_{\mathrm{J}}$',[2.398e11,7.37e11],rotation=-10,color='blue')
ax_main.annotate('0.99R$_{\mathrm{J}}$',[2.29e11,7.32e11],rotation=-10,color='blue')
ax_main.annotate('1.01R$_{\mathrm{J}}$',[2.56e11,7.43e11],rotation=-10,color='blue')
ax_main.annotate('1.225',[2.595e11,7.3e11],rotation=60,color='gray')
ax_main.annotate('1.25',[2.4e11,7.32e11],rotation=60,color='gray')
ax_main.annotate('1.275',[2.205e11,7.34e11],rotation=60,color='gray')
ax_main.annotate('i',[2.628e11,7.323e11],color='red')
ax_main.annotate('ii',[2.656e11,7.315e11],color='green')
ax_hit.annotate('ii',[-0.92e12,0.9e12],color='green')
ax_esc.annotate('i',[-0.92e12,0.9e12],color='red')
plt.subplots_adjust(hspace=0,wspace=0)
plt.savefig('fringe_MC.png')
plt.show()   
