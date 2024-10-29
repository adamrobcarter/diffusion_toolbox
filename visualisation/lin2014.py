import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import common

x = np.array([0.176,0.507,0.842,1.167,1.502,1.846,2.172,2.507,2.842,3.167,3.502,3.837,4.172,4.507,4.842,5.176,5.502,5.837,6.172,6.507,6.842,7.167,7.502,7.842])
y = np.array([0.981,0.972,0.972,0.968,0.97,0.965,0.963,0.961,0.959,0.954,0.947,0.941,0.939,0.939,0.941,0.945,0.945,0.947,0.945,0.936,0.925,0.912,0.908,0.9])

func = lambda t, D : np.exp(-0.16**2 * t * D)
popt, pcov = scipy.optimize.curve_fit(func, x, y)
print(common.format_val_and_unc(popt[0], np.sqrt(pcov[0, 0]), latex=False))
print(func(33.3, *popt))