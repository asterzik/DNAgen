#!/usr/bin/python

#https://github.com/keras-team/keras/blob/master/examples/lstm_text_generation.py

from __future__ import absolute_import
from __future__ import print_function, with_statement, division
import numpy as np
import random
import sys
import io
import os
from os import path
sys.path.append(path.dirname(path.dirname( path.abspath(__file__ ) ) ))
from helper_files.utility import read_fasta_file, string_to_array, get_weights_z

# Input parameters
fasta_input = './../data/chr01.saccharomyces_cerevisiae-JAN-19-2007.fasta'
initial_seq = read_fasta_file(fasta_input)
data = string_to_array(initial_seq)
p = get_weights_z(initial_seq)
for i in range(1,11):
    chrm = np.random.choice(['a', 't', 'g', 'c', 'z'], len(data), p=p)

    c_file = open(f'./../chromosomes/random/random_s.cerevisiae{i}.fa', 'w')
    c_file.write('>chr1_random \n')
    c = ''
    counter = 0
    for i in chrm:
        counter += 1
        c = c + i
        if counter%60 == 0:
           c_file.write(c +'\n')
           c = ''
