#!/usr/bin/python

#https://github.com/keras-team/keras/blob/master/examples/lstm_text_generation.py

from __future__ import absolute_import
from __future__ import print_function, with_statement, division
import keras
from keras.models import load_model
import numpy as np
import random
import sys
import io
from datetime import datetime
import tensorflow as tf
from keras.callbacks import LambdaCallback
import os
from os import path
sys.path.append(path.dirname( path.abspath(__file__ ) ) )
#from helper_files.moods import motif_match
from helper_files.utility import read_fasta_file, string_to_array, gc_content, gc_content_data, get_weights
from helper_files.plot import plot
#from preprocessing import pre

# Input parameters
fasta_input = 'data/chr01.saccharomyces_cerevisiae-JAN-19-2007.fasta'
initial_seq = read_fasta_file(fasta_input)
data = string_to_array(initial_seq)
#data = data[:100]
if len(sys.argv) != 6:
    raise SyntaxError('Usage: python generate_stateful.py [name of model] [number of episodes to train] [Saved_model or None] [description] [seq_length]')

model_type = sys.argv[1]
num_epochs = int(sys.argv[2])
saved_model = sys.argv[3]
time = sys.argv[4]
if time == 'time':
    time = datetime.now().strftime('%Y-%m-%d_%H:%M:%S')
seq_length = int(sys.argv[5])
if model_type == 'vae_lstm':
    seq_length = len(data)



# Hyper-parameters

args =  {
  "hidden_size": 1,
  "seq_length": seq_length,
  "learning_rate": 1e-1,
  "vocab_size": 4,
  "num_epochs": num_epochs,
  "batch_size": 1,
  "kernel_size": 147,
  "num_filters": 5,
  "patience": 100,
  "min_delta": 0,
  "interm_dim": 10000,
  "latent_dim": 1000
}

#hidden_size   = args["hidden_size"]
learning_rate = args["learning_rate"]
vocab_size    = args["vocab_size"]
batch_size    = args["batch_size"]
patience      = args["patience"]
min_delta     = args["min_delta"]

# Kommandozeilennavigation der Modelle:
print('Using {} model'.format(model_type))

if model_type == 'lstm':
    from models.lstm_stateful import Model
    logdir = "logs/lstm/"

else:
    raise ValueError('Unknown model type')
if saved_model == 'None':
    model = Model(args)
else:
    model = Model('./saved_models/{}/{}'.format(model_type,saved_model), load = True)
    print('Using saved model:{}'.format(saved_model))
logdir = logdir + str(time)
if not (os.path.isdir(logdir)):
    os.mkdir(logdir)

if num_epochs >0:
    print("Learning for", num_epochs, " epochs")
    sequences = []
    next_base = []
    for i in range(0, len(data) - seq_length):
        sequences.append(data[i: i + seq_length])
        next_base.append(data[i + seq_length])
    num_val = int(0.1*len(sequences))

    base_indices = dict((b, i) for i, b in enumerate(np.array(['a','c','g','t'])))#enumerate(np.array(['a','c','g','t','z'])))
    indices_base = dict((i, b) for i, b in enumerate(np.array(['a','c','g','t'])))# enumerate(np.array(['a','c','g','t','z'])))

    inputs = np.zeros((len(sequences), seq_length, vocab_size), dtype=np.bool)
    targets = np.zeros((len(sequences), vocab_size), dtype=np.bool)
    for i, sequence in enumerate(sequences):
        for t, char in enumerate(sequence):
            inputs[i, t, base_indices[char]] = 1
        if model_type != 'vae_lstm':
            targets[i, base_indices[next_base[i]]] = 1
    
    val_inputs = inputs[num_val:]
    val_targets = targets[num_val:]
    inputs = inputs[:num_val]
    targets = targets[:num_val]


# Callbacks

    savedir = 'saved_models/'+str(model_type)+'/'+str(time)
    if not (os.path.isdir(savedir)):
        os.mkdir(savedir)

    loss = []
    valloss = []
    acc = []
    valacc = []
    log = open(logdir + '/history.csv', 'w')
    log.write('epoch, loss, acc, val_loss, val_acc \n')
    log.close()
    for i in range(num_epochs):
        log = open(logdir + '/history.csv', 'a')
        model.model.reset_states()
        mean_loss = []
        mean_acc = []
        mean_valloss =[]
        mean_valacc = []
        for batch in range(len(inputs)):
            c_loss, c_acc = model.model.train_on_batch(np.reshape(inputs[batch],(1,1,4)), np.reshape(targets[batch],(1,4)))
            mean_loss.append(c_loss)
            mean_acc.append(c_acc)
        model.model.reset_states()
        for j in range(len(val_inputs)):
            v_loss, v_acc = model.model.test_on_batch(np.reshape(val_inputs[j],(1,1,4)), np.reshape(val_targets[j],(1,4)))
            mean_valloss.append(v_loss)
            mean_valacc.append(v_acc)
        model.model.reset_states()
        valloss = np.mean(mean_valloss)
        log.write('{},{},{},{},{}\n'.format(i,np.mean(mean_loss),np.mean(mean_acc),valloss,np.mean(mean_valacc)))
        log.close()
        model.model.save(savedir + '/checkpoint-{}-{}.h5'.format(i, valloss))
        print(loss, valloss, acc, valacc)




dnafile = open('chromosomes/{}/{}.fa'.format(model_type, time), 'w')
base_indices = dict((b, i) for i, b in enumerate(np.array(['a','c','g','t'])))
indices_base = dict((i, b) for i, b in enumerate(np.array(['a','c','g','t'])))
dnafile.write('>seq1\n')

#Create random start sequence and generate new nucleotides.
#Subsequently append sequence and cut first base of.
#Use new sequence for next prediction step.
seq = np.random.choice(['a', 't', 'g', 'c'], seq_length, p=get_weights(initial_seq))
size = 0
print('len(data) = {}'.format(len(data)))
while size < len(data):
    for i in range(len(seq)):
        x_pred = np.zeros((1, seq_length, vocab_size))
        for t, base in enumerate(seq):
            x_pred[0, t, base_indices[base]] = 1.
    size += 1
    preds = model.model.predict(x_pred, verbose=0)[0] #this is just because it is a nested list [[0,1,...]]
    next_index = model.sample(preds, 1)
    next_base = indices_base[next_index]
    # get first base to the rear and shift everything else to the left, order is right here, double checked.
    seq = np.roll(seq, seq_length - 1)
    seq[seq_length - 1] = next_base
    dnafile.write(str(next_base))
    if size % 60 == 0:
        dnafile.write('\n')
dnafile.close()
print("finish generation")
