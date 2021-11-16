# IO defs for loading pathway data.
# Zilin Song, 23 AUG 2021
# 

import numpy

def load_ds(sysname, pathname,):
    '''Load various normalized/rescaled datasets and the corresponding labels from 
    sysname:    toho_amp, toho_cex;
    pathname:   r1ae, r2ae;
    '''
    # directory settings.
    prefx    = '../2.datasets/4.conclude_ds/fin_ds'
    xname    = 'x'
    xpohanme = 'x_onehot'
    yname    = 'y'

    xdir     = f'{prefx}/{sysname}.{pathname}.{xname}.npy'
    xpohdir  = f'{prefx}/{sysname}.{pathname}.{xpohanme}.npy'
    ydir     = f'{prefx}/{sysname}.{pathname}.{yname}.npy'

    x    = numpy.load(xdir)
    xpoh = numpy.load(xpohdir)
    y    = numpy.expand_dims(numpy.load(ydir), 1)

    return x, xpoh, y
