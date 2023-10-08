import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
xy = np.load('out.npy')

x=xy[:,0]
y=xy[:,1]

n=len(x)

#xy = xy.reshape(-1, 1, 2)
#segments = np.hstack([xy[:-1], xy[1:]])

#fig, ax = plt.subplots()
#coll = LineCollection(segments, cmap='cool')
#coll.set_array(np.random.random(xy.shape[0]))

#ax.add_collection(coll)
#ax.autoscale_view()

#plt.show()

T=np.linspace(0,1,np.size(x))**2
fig = plt.figure()
ax = fig.add_subplot(111)

# Segment plot and color depending on T
s = 10 # Segment length
for i in range(0,n-s,s):
    ax.plot(x[i:i+s+1],y[i:i+s+1],color=(0.0,0.5,T[i]))
plt.show()

