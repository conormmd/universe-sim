import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from progressbar import ProgressBar
second = 3e7
a=5000
#-pi/2 -> 3e8*a/second
#angle = np.array([-np.pi/2,1.1605,0.08474,0,-1.1397,-1.460987,1.199])
#lifetime = np.array([1.47e8*a/second,1500000*a/second,150000*a/second,15000*a/second,1.86e7*a/second,2.72e7*a/second,8.3e6*a/second])

angle = np.array([-1.5, -1.4, -1.3])
lifetime = np.array([50000,1.8e8*a/second,1.4e8*a/second])
angle = np.array([-np.pi/2,1.1605,0.08474,0,-1.1397,-1.460987,1.199])
lifetime = np.array([1.47e8*a/second,1500000*a/second,150000*a/second,15000*a/second,1.86e7*a/second,2.72e7*a/second,8.3e6*a/second])

print(lifetime)
plt.scatter(angle,lifetime)
plt.xlim(1.571,-1.571)
plt.ylim(0,50000)
plt.show()
