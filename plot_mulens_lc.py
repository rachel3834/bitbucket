from os import path
import numpy as np
import matplotlib.pyplot as plt

from pyLIMA import event
from pyLIMA import telescopes
from pyLIMA.models import PSPL_model, FSBL_model, USBL_model, pyLIMA_fancy_parameters
from pyLIMA.fits import DE_fit
from pyLIMA.fits import TRF_fit, DE_fit
from pyLIMA.outputs import pyLIMA_plots

# Configuration
data_dir = './data/'

lc_file = 'ulwdc1_100_W149.txt'
model_file =  'ulwdc1_100_W149_model.txt'

data = np.loadtxt(path.join(data_dir, lc_file))
model = np.loadtxt(path.join(data_dir, model_file))

dt = 2458200.0
ts = data[:,0] - dt
fsize = 18

fig, ax = plt.subplots(figsize=(8,8))

ax.errorbar(ts, data[:,1], yerr=data[:,2], fmt='.', color='b')
ax.plot(model[:,0]-dt, model[:,1], 'k-')

ax.set_xlabel('HJD - '+str(round(dt,0)), fontsize=fsize)
ax.set_ylabel('Mag', fontsize=fsize)
ax.tick_params(axis='x', labelsize=fsize)
ax.tick_params(axis='y', labelsize=fsize)

[xmin, xmax, ymin, ymax] = plt.axis()
plt.gca().invert_yaxis()
plt.grid()

# Create the inset plot
axins = ax.inset_axes([0.5, 0.5, 0.4, 0.4])  # [left, bottom, width, height]
xmin = 320.0
xmax = 340.0
idx1 = np.where(ts >= xmin)[0]
idx2 = np.where(ts <= xmax)[0]
idx = list(set(idx1).intersection(set(idx2)))
axins.errorbar(ts[idx], data[idx,1], yerr=data[idx,2], fmt='.', color='b')
idx1 = np.where(model[:,0]-dt >= xmin)[0]
idx2 = np.where(model[:,0]-dt <= xmax)[0]
idx = list(set(idx1).intersection(set(idx2)))
axins.plot(model[idx,0]-dt, model[idx,1], 'k-')
axins.invert_yaxis()
axins.grid()

plt.tight_layout()
#plt.show()
plt.savefig(path.join(data_dir, lc_file.replace('.txt', '.png')))