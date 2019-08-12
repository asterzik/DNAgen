# DNAgen

This repository contains some of the code that I used during my bachelor thesis.
It can be used for predicting DNA sequences of chromosomal size using neural networks. I included two example LSTMs, one CNN and one VAE in this repository. They can be adjusted easily.

generate.py generates a chromosome based on S. cerevisiae's chromosome 1.

generate_augmented.py augments the dataset by salt and pepper augmentation.
Both of these can be used for LSTM and CNN.

generate_stateful generates a chromosome with a stateful LSTM network.

The generate_incl_vae.py file could use a VAE for prediction as well, but it does not work for DNA of chromosomal size.
Some of the code was inspired by 
https://github.com/keras-team/keras/blob/master/examples/lstm_text_generation.py.

I also included a testpipeline for testing the resulting chromosomes on their quality.
This makes use of several other software:

MOODS: https://github.com/jhkorhonen/MOODS

EMBOSS newcpgreport: http://emboss.sourceforge.net/apps/cvs/emboss/apps/newcpgreport.html

NucEnerGen: http://nucleosome.rutgers.edu/nucenergen/


The network training was done with Keras on python 2 and testing on python 3.
