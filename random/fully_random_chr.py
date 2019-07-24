#!/usr/bin/python


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
for i in range(10):
    chrm = np.random.choice(['a', 't', 'g', 'c'], len(data))
    c_file = open(f'../chromosomes/random/fully_random_s.cerevisiae{i}.fa', 'w')
    c_file.write('>chr1_random \n')
    c = ''
    counter = 0
    for i in chrm:
        counter += 1
        c = c + i
        if counter%60 == 0:
            c_file.write(c +'\n')
            c = ''
