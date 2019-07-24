
from __future__ import division
from __future__ import print_function
import os,fnmatch
import numpy as np

from Bio import SeqIO
import re

from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder

# difference sequence (http://mathonline.wikidot.com/difference-sequences)
# markov permutation
# random permutation
# sequence alignment (Edit Distance, Levenshtein Distance, Alignment Distance)
#        Alignment Matrix (Needleman-Wunsch Algorithm) Traceback
#        Similarity
#        Smith-Waterman Algorithm
#        Substitution/Similarity Matrices
#        Heuristic Multiple Alignment: Progressive Alignment
#        Guide tree

#
# Function to read in FASTA file
def read_fasta_file(input):

    seq = ''
    for seq_record in SeqIO.parse(input, "fasta"):
        seq = seq + str(seq_record.seq[:])
    
    return seq

#
# function to convert a DNA sequence string to a numpy array
# converts to lower case, changes any non 'acgt' characters to 'z'
def string_to_array(my_string):
    my_string = my_string.lower()
    my_string = re.sub('[^acgt]', 'z', my_string)
    my_array = np.array(list(my_string))
    return my_array

#

# function to encode a DNA sequence string as an ordinal vector
# returns a numpy vector with a=0.25, c=0.50, g=0.75, t=1.00, z=0.00
def ordinal_encoder(my_array):

    # create a label encoder with 'acgtz' alphabet
    from sklearn.preprocessing import LabelEncoder
    label_encoder = LabelEncoder()
    label_encoder.fit(np.array(['a','c','g','t','z']))

    integer_encoded = label_encoder.transform(my_array)
    float_encoded = integer_encoded.astype(float)
    float_encoded[float_encoded == 0] = 0.25 # A
    float_encoded[float_encoded == 1] = 0.50 # C
    float_encoded[float_encoded == 2] = 0.75 # G
    float_encoded[float_encoded == 3] = 1.00 # T
    float_encoded[float_encoded == 4] = 0.00 # anything else, z
    return float_encoded


#
# function to one-hot encode a DNA sequence string
# non 'acgt' bases (n) are 0000
# returns a L x 4 numpy array

def one_hot_encoder(my_array):

    # create a label encoder with 'acgtz' alphabet
    label_encoder = LabelEncoder()
    label_encoder.fit(np.array(['a','c','g','t','z']))

    integer_encoded = label_encoder.transform(my_array)
    onehot_encoder = OneHotEncoder(sparse=False, dtype=int, categories=[range(6)])
    integer_encoded = integer_encoded.reshape(len(integer_encoded), 1)
    onehot_encoded = onehot_encoder.fit_transform(integer_encoded)
    onehot_encoded = np.delete(onehot_encoded, -1, 1)
    
    return onehot_encoded


def ordinal_one_hot_encoder(my_array):
    '''Not being used so far'''

    label_encoder = LabelEncoder()
    label_encoder.fit(np.array([0,1,2,3,4]))

    integer_encoded = label_encoder.transform(my_array)
    onehot_encoder = OneHotEncoder(sparse=False, dtype=int, categories=[range(6)])
    integer_encoded = integer_encoded.reshape(len(integer_encoded), 1)
    onehot_encoded = onehot_encoder.fit_transform(integer_encoded)
    onehot_encoded = np.delete(onehot_encoded, -1, 1)
    
    return onehot_encoded


def mutate(dna):
    '''Takes a string of bases and replaces one of the bases by another random base. Returns a string'''
    dna_list = list(dna)
    mutation_site = random.randint(0, len(dna_list) - 1)
    dna_list[mutation_site] = random.choice(list('ATCG'))
    return ''.join(dna_list)


def get_base_counts(dna):
    '''Takes a string of bases and returns a dictionary containing the count for every base'''
    counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
    dna = dna.upper()
    for base in dna:
        if base in {'A','T','G','C'}:
            counts[base] += 1
    
    return counts
def get_base_counts_z(dna):
    '''Takes a string of bases and returns a dictionary containing the count for every base'''
    counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'Z':0}
    dna = dna.upper()
    for base in dna:
        if base in {'A','T','G','C','Z'}:
            counts[base] += 1
    
    return counts

def get_weights(dna):
    dic = get_base_counts(dna)
    a = dic['A']
    t = dic['T']
    g = dic['G']
    c = dic['C']
    
    su = a+t+g+c
        
    return [a/su,t/su,g/su,c/su]

def get_weights_z(dna):
    dic = get_base_counts_z(dna)
    a = dic['A']
    t = dic['T']
    g = dic['G']
    c = dic['C']
    z = dic['Z']
    
    su = a+t+g+c+z
        
    return [a/su,t/su,g/su,c/su,z/su]
def gc_content(seq):
    '''Takes a string of bases as input and prints count for every base as well as gc content. No outoput'''

    counts = get_base_counts(seq)
    a = counts['A']
    t = counts['T']
    g = counts['G']
    c = counts['C']
    
    return (g+c)*100.0/(a+t+g+c)

def gc_content_data(data):
    '''Takes an array of bases and returns the gc_content'''
    
    counts = {'a': 0, 't': 0, 'g': 0, 'c': 0, 'z': 0}
    
    for base in data:
        counts[base] += 1


    a = counts['a']
    t = counts['t']
    g = counts['g']
    c = counts['c']
    
    gc_content=(g+c)*100.0/(a+t+g+c)
   
    #print ("gc_content= %f" %(gc_content))

    return gc_content
def gc_count(data):
    '''Takes an array, returns g+c'''
    counts = {'a': 0, 't': 0, 'g': 0, 'c': 0, 'z': 0}
    
    for base in data:
        counts[base] += 1


    a = counts['a']
    t = counts['t']
    g = counts['g']
    c = counts['c']
    
    return g + c
 
def cpg_oeratio(dna):
    
    """ Takes DNA string and returns cpg content in %"""
    data = string_to_array(dna)
    observed = 0
    for i in range(len(data)-1):
        if data[i]=='c' and data[i+1]=='g':
            observed += 1
    dic = get_base_counts(dna)
    a = dic['A']
    t = dic['T']
    g = dic['G']
    c = dic['C']
    expected = c*g/len(data)
    return observed/expected
