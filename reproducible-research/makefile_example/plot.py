import numpy as np
from matplotlib import pyplot as plt


X = np.load('results/data.npy')
y = np.load('results/fit.npy')
fig, ax = plt.subplots()
ax.plot(y, linewidth=2, c='r')
ax.scatter(np.arange(len(X)), X)
ax.set_title('My awesome plot')
fig.savefig('images/my_plot.png')
