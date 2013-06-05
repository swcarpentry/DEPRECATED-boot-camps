"""
Creates a noisy dataset
"""

import numpy as np
n = 100
X = 2 * np.arange(100) + np.random.randint(0, 30, (100,))
X.dump('results/data.npy')
