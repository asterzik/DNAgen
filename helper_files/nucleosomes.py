from __future__ import division,print_function
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
import sys

if len(sys.argv) != 4:
    raise SyntaxError ('Usage: python nucleosomes.py [datafile] [model_type] [description]')

datafile = sys.argv[1]
pos, en = np.loadtxt(datafile, unpack = True, skiprows = 1)
#pos = pos[:10000]
#en = en[:10000]
spl = UnivariateSpline(pos, en)
spl.set_smoothing_factor(150000)
plt.plot(pos, spl(pos), 'b', lw=3)
#plt.plot(pos, en, 'ro', ms = 0.1)
plt.savefig('results/{}/{}/nuc.png'.format(sys.argv[2],sys.argv[3]))
