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
if len(sys.argv) !=4:
    raise SyntaxError("Usage: python Auswertung.py [sequence_file] [model_type] [description/'time']")

seq_file = open(sys.argv[1],'r')
model_type = sys.argv[2]
time = sys.argv[3]
if time == 'time':
    time = datetime.now().strftime('%Y-%m-%d_%H:%M:%S')

seq = seq_file.read()
seq_file.close()

path = './results/{}/{}/'.format(model_type, time)
if not (os.path.isdir(path)):
    os.mkdir(path)
res = open(path + 'results.txt','w')
res.write('Ausertung der Sequenz {} fuer das Netzwerk {} \n \n'.format(sys.argv[1],model_type))

gc = gc_content(seq)
res.write('gc_content: {} \n'.format(gc))

cpg = cpg_oeratio(seq)
res.write('CpG observed to expected ration: {} \n'.format(cpg))

hits, tot = motif_match(string_to_array(seq))
for entry in hits:
    res.write(entry + '\n')
res.write('total: {}'.format(tot))
res.close()
