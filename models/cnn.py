import numpy as np
from keras.models import Sequential
from keras.layers import Dense
from keras.layers import Conv1D
from keras.layers import MaxPooling1D
from keras.layers import Flatten
from keras.layers import LeakyReLU
from keras.layers import Dropout
from keras.models import load_model
from keras.optimizers import Adam
import random
from keras import backend as K
import tensorflow as tf

K.tensorflow_backend._get_available_gpus()

class Model():
    def __init__(self, args, load = False):
        if not load:
            self.args = args
        
            seq_length    = args["seq_length"]
            vocab_size    = args["vocab_size"]

            self.model = Sequential()
            self.model.add(Conv1D(filters = 64 , kernel_size = 3, \
                   activation='relu', padding='same', input_shape=(seq_length,vocab_size)))
            self.model.add(Conv1D(filters = 64 , kernel_size = 3, \
                    activation='relu', padding='same'))
            #self.model.add(LeakyReLU(0.1))
            self.model.add(MaxPooling1D())
            self.model.add(Conv1D(128, 3, activation='relu'))
            self.model.add(Conv1D(128, 3, activation='relu'))
            #self.model.add(LeakyReLU(0.1))
            self.model.add(Dropout(0.5))
            self.model.add(MaxPooling1D())
            #self.model.add(LeakyReLU(0.1))
            self.model.add(Flatten())
            self.model.add(Dense(64, activation='relu'))
            self.model.add(Dense(4, activation='softmax'))
            self.optimizer = Adam(lr=0.001)
            self.model.compile(loss='categorical_crossentropy', optimizer=self.optimizer, metrics=['acc'])
            print(self.model.summary())
        else:
            self.model = load_model(args)

    def sample(self, preds):
        # helper function to sample an index from a probability array
        preds = np.asarray(preds).astype('float64')
        preds = preds / np.sum(preds)
        probas = np.random.multinomial(1, preds, 1)
        return np.argmax(probas)
