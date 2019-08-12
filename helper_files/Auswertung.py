#!/usr/bin/python

from __future__ import absolute_import
from __future__ import print_function, with_statement, division
import numpy as np
from datetime import datetime
import sys
import io
import os
from os import path
sys.path.append(path.dirname( path.abspath(__file__ ) ) )
from moods import motif_match
from utility import gc_content,cpg_oeratio, string_to_array
from plot import plot
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
import matplotlib.ticker as mtick
plt.rcParams.update({'font.size': 18})
#plt.locator_params(axis='y', nbins=6)
#plt.locator_params(axis='x', nbins=3)
if len(sys.argv) !=5:
    raise SyntaxError("Usage: python Auswertung.py [sequence_file] [model_type] [description/'time'] [datafile]")

seq_file = open(sys.argv[1],'r')
model_type = sys.argv[2]
time = sys.argv[3]
if time == 'time':
    time = datetime.now().strftime('%Y-%m-%d_%H:%M:%S')

seq = seq_file.read()
seq_file.close()
print(model_type)
if model_type == 'cnn':
    res = open('results/cnn_result.csv','a')
elif model_type == 'random':
    res = open('results/random_result.csv','a')
elif model_type == 'original':
    res = open('results/original_result.csv','a')
else:
    print('lstm')
    res = open('results/lstm_result.csv','a')

gc = gc_content(seq)
cpg = cpg_oeratio(seq)
hits, tot = motif_match(string_to_array(seq))

pos_org , en_org = np.loadtxt('results/original/s.cerevisiae/s.cerevisiae.fa.nuc', unpack = True, skiprows = 1)
datafile = sys.argv[4]
pos, en = np.loadtxt(datafile, unpack = True, skiprows = 1)
length = len(en) 
assert np.alltrue(pos[:length]==pos_org[:length])
abst = np.sum(abs(en[:length]-en_org[:length]))/length
#pos = pos[:10000]
#en = en[:10000]
spl = UnivariateSpline(pos, en)
spl.set_smoothing_factor(150000)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fig = plt.figure()
plt.gcf().subplots_adjust(bottom=0.13)
plt.gcf().subplots_adjust(left=0.13)
ax = fig.add_subplot(111)
ax.plot(pos,spl(pos), 'black', lw=1)
ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
every_nth = 2
plt.ylabel(r'Nucleosome binding energy')
plt.xlabel('base pair')
#plt.plot(pos, en, 'ro', ms = 0.1)
path = ('results/{}/{}/'.format(sys.argv[2],sys.argv[3]))
if not os.path.exists(path):
    os.mkdir(path)
plt.savefig(path + 'nuc.eps', format='eps', dpi= 1000)
res.write(f'{sys.argv[1]},{np.round(gc,2)},{np.round(cpg,2)},{hits[0]},{hits[1]},{np.round(np.mean(en),1)},{np.round(np.std(en),1)},{np.round(min(en),1)},{np.round(max(en),1)},{abst},')
res.close()
