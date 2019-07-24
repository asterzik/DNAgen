import numpy as np
import matplotlib
matplotlib.use('Agg')
from datetime import datetime
import matplotlib.pyplot as plt

def plot(model_type, gc_content, CTCF_count_original, CTCF_count_generated, path = './gen_dna/figs/'):
    '''Plots the specified quantities with matplotlib.pyplot.
    It's possible to define a custom path'''
    
    time = datetime.now().strftime('%Y-%m-%d_%H:%M:%S')


    gc = np.asarray(gc_content)
    co = np.asarray(CTCF_count_original)
    cg = np.asarray(CTCF_count_generated)

    ##################
    #
    # Save to txt files
    #
    ##################
    myfile = open('./gen_dna/figs/{}/{}.txt'.format(model_type, time), 'w')
    myfile.write('gc_content: \n {} \n \n count_original: \n {} \n \n count_generated: \n {}'.format(gc, co, cg))
    myfile.close()

    ##################
    #
    # Create pictures
    #
    ##################
    plt.plot(np.arange(len(gc)),gc, 'bo')
    plt.title('gc content generated')
    plt.ylabel('gc_content')
    plt.xlabel('epochs')
    plt.savefig(path +'gc' + str(time) +'.png')
    plt.close()

    plt.plot(np.arange(len(cg)),cg, 'mo', label = 'generated')
    plt.plot(np.arange(len(co)),co, 'go', label = 'original')
    plt.legend()
    plt.title('CTCF count')
    plt.ylabel('CTCF count')
    plt.xlabel('epochs')
    plt.savefig(path +'count_CTCF' +str(time)+'.png')

