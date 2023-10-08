import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection


#Params
G = 6.6726e-11; d = 7.7834082e11; J_mass = 1.898130e27; S_mass = 1.98847e30; r_S = d*(J_mass/(J_mass+S_mass)); r_J = d*(S_mass/(J_mass+S_mass)); P = np.sqrt((d**3/(G*(S_mass+J_mass))))*2*np.pi; V_S = (2*np.pi*r_S)/P; V_J = (2*np.pi*r_J)/P

x=np.load('x.npy')
y=np.load('y.npy')
N=len(y)
i_t = []
i_b = []
for i in range(10,N-10):
    if y[i-10]<y[i-5]<y[i-3]<y[i] and y[i+10]<y[i+5]<y[i+3]<y[i]:
        i_t.append(i)
    if y[i-10]>y[i-5]>y[i-3]>y[i] and y[i+10]>y[i+5]>y[i+3]>y[i]:
        i_b.append(i)

x_t = np.zeros(len(i_t))
y_t = np.zeros(len(i_t))
for i in range(0,len(i_t)):
    j = i_t[i]
    x_t[i]=x[j]
    y_t[i]=y[j]

plt.scatter(x_t,y_t)
#plt.xlim(0,50000)
#plt.ylim(0,2.1)
plt.show()
