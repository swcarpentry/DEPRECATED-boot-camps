import numpy as np
from sklearn.linear_model import LinearRegression


X = np.load('results/data.npy')
clf = LinearRegression()
clf.fit(np.arange(len(X))[:, np.newaxis], X[:, np.newaxis])
y = clf.predict(np.arange(len(X))[:, np.newaxis])
y.dump('results/fit.npy')
