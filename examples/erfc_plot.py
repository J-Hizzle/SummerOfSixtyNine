# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erfc
# %%
sigma_vals = np.arange(1, 4)
# %%
xvals = np.linspace(-10, 10, 1000)
# %%
for sigma in sigma_vals:
    xinp = xvals/sigma
    yvals = erfc(xinp)
    plt.plot(xvals, yvals, label='{0}'.format(sigma)) 
    if sigma == sigma_vals.max:
        plt.legend()
        plt.show()
# %%
# %%
