import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from progressbar import ProgressBar
from functions import *

total=np.load('total.npy')
hit=np.load('hit.npy')
escape=np.load('escape.npy')

for i in range(0,len(total)):
    plt.scatter(total[i][0],total[i][1],color='black',s=3)
for i in range(0,len(hit)):
    plt.scatter(hit[i][0],hit[i][1],color='lime',s=3)
for i in range(0,len(escape)):
    plt.scatter(escape[i][0],escape[i][1],color='red',s=3)

plt.xlim(0,1e12)
plt.ylim(-1e12,1e12)
plt.legend()
plt.show()
    
