from __future__ import absolute_import
from __future__ import print_function, with_statement, division
from copy import deepcopy
import numpy as np

def pre(data, model, seq_length = 500, salt_pepper = 1, sp_amount = 0.05):
    ''' Takes the original sequence as np.array, name of the model and
    optionally seq_length, salt_pepper and sp_amount.
    salt_pepper defines how many augmented copies will be made.
    sp_amount is the percentage of 'mutations' that is going to be made'''
    
    #Salt and Pepper Augmentation
    sequences = []
    next_base = [] 
    for i in range(salt_pepper):
        data1 = deepcopy(data)
        mask = np.random.choice([True,False],len(data),p=[sp_amount,1-sp_amount])
        vals = np.random.choice(['a','c','g','t'],np.sum(mask))
        np.place(data1,mask,vals)
    
    #Creating inputs (sequences) and targets (next_base) for the network
        for i in range(0, len(data1) - seq_length):
            sequences.append(data1[i: i + seq_length])
            next_base.append(data1[i + seq_length])

    return sequences, next_base


