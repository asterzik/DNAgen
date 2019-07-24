#!/usr/bin/python

#Some of the preprocessing and final generation was done according to https://github.com/keras-team/keras/blob/master/examples/lstm_text_generation.py
#make it compatible with python2
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
from helper_files.utility import read_fasta_file, string_to_array, gc_content, gc_content_data, get_weights
from helper_files.plot import plot

# Input parameters
fasta_input = 'data/chr01.saccharomyces_cerevisiae-JAN-19-2007.fasta'
initial_seq = read_fasta_file(fasta_input)
data = string_to_array(initial_seq)
data = data[:2200]
if len(sys.argv) != 6:
    raise SyntaxError('Usage: python generate.py [name of model] [number of episodes to train] [Saved_model or None] [description] [seq_length]')

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
  "seq_length": seq_length,
  "vocab_size": 4,
  "num_epochs": num_epochs,
  "batch_size": 50,
  "patience": 50,
  "min_delta": 0,
}

vocab_size    = args["vocab_size"]
batch_size    = args["batch_size"]
patience      = args["patience"]
min_delta     = args["min_delta"]

print('Using {} model'.format(model_type))

if model_type == 'cnn':
    from models.cnn import Model
    logdir = "logs/cnn/"
elif model_type == 'lstm':
    from models.lstm import Model
    logdir = "logs/lstm/"
else:
    raise ValueError('Unknown model type')

if saved_model == 'None':
    model = Model(args)
else:
    model = Model('saved_models/{}/{}'.format(model_type,saved_model), load = True)
    print('Using saved model:{}'.format(saved_model))

#Set path for tensorboard
logdir = logdir + str(time)
if not (os.path.isdir(logdir)):
    os.mkdir(logdir)

if num_epochs >0:
    print("Learning for", num_epochs, " epochs")

    # Cut data into overlapping sequences
    sequences = []
    next_base = []
    for i in range(0, len(data) - seq_length):
        sequences.append(data[i: i + seq_length])
        next_base.append(data[i + seq_length])
    
    #one-hot enconding
    base_indices = dict((b, i) for i, b in enumerate(np.array(['a','c','g','t'])))
    indices_base = dict((i, b) for i, b in enumerate(np.array(['a','c','g','t'])))

    inputs = np.zeros((len(sequences), seq_length, vocab_size), dtype=np.bool)
    targets = np.zeros((len(sequences), vocab_size), dtype=np.bool)
    for i, sequence in enumerate(sequences):
        for t, char in enumerate(sequence):
            inputs[i, t, base_indices[char]] = 1
            targets[i, base_indices[next_base[i]]] = 1



# Callbacks

    savedir = 'saved_models/'+str(model_type)+'/'+str(time)
    if not (os.path.isdir(savedir)):
        os.mkdir(savedir)
    save_callback = keras.callbacks.ModelCheckpoint(filepath = savedir + '/checkpoint-{epoch:02d}-{val_loss:.2f}.h5', mode='auto', period=1)
    tensorboard_callback = keras.callbacks.TensorBoard(log_dir=logdir)
    #monitor validation loss to prevent over fitting, optional 
    early_stopping = keras.callbacks.EarlyStopping(monitor='val_loss', min_delta=min_delta, patience=patience)

# Training
print('start training')
model.model.fit(inputs, targets,validation_split = 0.2,
          batch_size = batch_size,
          epochs=num_epochs, verbose = 0,
          callbacks=[save_callback, tensorboard_callback,early_stopping])


# Generation of Chromosome
print('start generation')
dnafile = open('chromosomes/{}/{}.fa'.format(model_type, time), 'w')
base_indices = dict((b, i) for i, b in enumerate(np.array(['a','c','g','t'])))
indices_base = dict((i, b) for i, b in enumerate(np.array(['a','c','g','t'])))
dnafile.write('>seq1\n')
#Create random start sequence and generate new nucleotides.
#Subsequently append sequence and cut first base of.
#Use new sequence for next prediction step.
seq = np.random.choice(['a', 't', 'g', 'c'], seq_length, p=get_weights(initial_seq))
size = 0
while size < len(data):
    for i in range(len(seq)):
        x_pred = np.zeros((1, seq_length, vocab_size))
        for t, base in enumerate(seq):
            x_pred[0, t, base_indices[base]] = 1.
    size += 1
    preds = model.model.predict(x_pred, verbose=0)[0] #this is just because it is a nested list [[0,1,...]]
    next_index = model.sample(preds)
    next_base = indices_base[next_index]
    # get first base to the rear and shift everything else to the left, order is right here, double checked.
    seq = np.roll(seq, seq_length - 1)
    seq[seq_length - 1] = next_base
    dnafile.write(str(next_base))
    if size % 60 == 0:
        dnafile.write('\n')
dnafile.close()
print("finish generation")
