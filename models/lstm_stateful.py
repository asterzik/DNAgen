import numpy as np
from keras.models import Sequential
from keras.layers import Dense
from keras.layers import LSTM
from keras.models import load_model
from keras.optimizers import RMSprop
import random

class Model():
    def __init__(self, args, load = False):
        if not load:
            self.args = args
        
            hidden_size   = args["hidden_size"]
            seq_length    = args["seq_length"]
            learning_rate = args["learning_rate"]
            vocab_size    = args["vocab_size"]
            num_epochs    = args["num_epochs"]

            self.model = Sequential()
            self.model.add(LSTM(400, input_shape=(seq_length,vocab_size),batch_size = 1, stateful = True))
            self.model.add(Dense(vocab_size, activation='softmax'))

            self.optimizer = RMSprop(lr=0.001)
            self.model.compile(loss='categorical_crossentropy', optimizer=self.optimizer, metrics=['acc'])
            print(self.model.summary())
        else:
            self.model = load_model(args)

    def sample(self, preds, temperature=1.0):
        # helper function to sample an index from a probability array
        preds = np.asarray(preds).astype('float64')
#        preds = np.log(preds) / temperature
 #       exp_preds = np.exp(preds)
        preds = preds / np.sum(preds)
        probas = np.random.multinomial(1, preds, 1)
        return np.argmax(probas)
