#https://github.com/twairball/keras_lstm_vae/blob/master/lstm_vae/vae.py

from __future__ import absolute_import
from __future__ import print_function, with_statement, division
import numpy as np
from keras.layers import Dense, Input, Lambda, LSTM, RepeatVector, TimeDistributed
from keras.activations import softmax
from keras.optimizers import Adam
from keras.models import Model as K_Model
from keras.models import load_model
from keras import objectives
from keras import backend as K


class Model():
    def __init__(self, args, load = np.bool_(False)):
        if not load:
            self.args = args

            seq_length    = args["seq_length"]
            learning_rate = args["learning_rate"]
            vocab_size    = args["vocab_size"]
            num_epochs    = args["num_epochs"]
            batch_size    = args["batch_size"]
            interm_dim    = 1000
            latent_dim    = 100

            # Reparameterization Trick
            def sampling(args): 
                """Reparameterization trick by sampling from an isotropic unit Gaussian.
             Arguments
                args (tensor): mean and log of variance of Q(z|X)
             Returns
                z (tensor): sampled latent vector
            """
                z_mean, z_log_var = args
                batch = K.shape(z_mean)[0] #TODO is this really the batchsize?
                print('batchsize sampling: %b', batch)
                epsilon = K.random_normal(shape=(batch, latent_dim))
                return z_mean + K.exp(0.5 * z_log_var) * epsilon



            shape = (seq_length, vocab_size,)
            inp = Input(shape)
            h = LSTM(interm_dim)(inp)

            z_mean = Dense(latent_dim, name='z_mean')(h)
            z_log_var = Dense(latent_dim, name='z_log_var')(h)
            z = Lambda(sampling, output_shape=(latent_dim,), name='z')([z_mean, z_log_var])

            decoder_h = LSTM(interm_dim, return_sequences=True)
            decoder_mean = LSTM(vocab_size, return_sequences=True, activation = softmax)

            h_decoded = RepeatVector(seq_length)(z)
            h_decoded = decoder_h(h_decoded)

            h_mean = decoder_mean(h_decoded)
            output = TimeDistributed(Dense(vocab_size, activation='softmax'))(h_mean)
            
            #end-to-end autoencoder
            self.model = K_Model(inp, output)
	    # generator, from latent space to reconstructed inputs
            decoder_input = Input((latent_dim,))
            
            _h_decoded = RepeatVector(seq_length)(decoder_input)
            _h_decoded = decoder_h(_h_decoded)
            
            _h_mean = decoder_mean(_h_decoded)
            _output = TimeDistributed(Dense(vocab_size, activation='softmax'))(_h_mean)
            self.generator = K_Model(decoder_input, _output)

            def vae_loss(inp,output):
                xent_loss = objectives.mse(inp, output)
                kl_loss = - 0.5 * K.mean(1 + z_log_var - K.square(z_mean) - K.exp(z_log_var))
                return xent_loss + kl_loss

            self.model.compile(optimizer='rmsprop', loss=vae_loss)
            #self.model.summary()
            #This is only to avoid: AttributeError: 'Model' object has no attribute
            #'stateful_metric_names'. Compilaton here doesn't change the prediction
            # (see https://github.com/keras-team/keras/issues/9394)
            self.generator.compile('sgd','mse')
        else:
            self.model = load_model(args)



