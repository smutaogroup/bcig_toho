# IO defs for loading pathway data.
# Zilin Song, 29 AUG 2021
# 

import numpy
from tensorflow import keras

def load_ds(sysname, pathname):
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

def load_model(sysname, pathname):
    '''Load the trained machine learning models
    '''
    tag = 0 if sysname == 'toho_amp' and pathname == 'r1ae' else \
          1 if sysname == 'toho_cex' and pathname == 'r1ae' else \
          2 if sysname == 'toho_amp' and pathname == 'r2ae' else \
          3 if sysname == 'toho_cex' and pathname == 'r2ae' else \
          4 if sysname == 'both'     and pathname == 'r1ae' else \
          5 if sysname == 'both'     and pathname == 'r2ae' else \
          None

    moddir = f'../3.dwnn/cpu_mod{tag}/model.h5'
    mod = keras.models.load_model(moddir)

    return mod

def load_pred_barriers(sysname, pathname):

    tag = 0 if sysname == 'toho_amp' and pathname == 'r1ae' else \
          1 if sysname == 'toho_cex' and pathname == 'r1ae' else \
          2 if sysname == 'toho_amp' and pathname == 'r2ae' else \
          3 if sysname == 'toho_cex' and pathname == 'r2ae' else \
          4 if sysname == 'both'     and pathname == 'r1ae' else \
          5 if sysname == 'both'     and pathname == 'r2ae' else \
          None

    ypreddir = f'../3.dwnn/cpu_mod{tag}/y_pred.npy'
    ypred = numpy.load(ypreddir)

    ypred_paths = numpy.reshape(ypred, (-1, 50))
    ypred_barrier = []
    for pid in range(ypred_paths.shape[0]):
        barrier = -999.

        for repid in range(ypred_paths.shape[1]):
            delta_ener = ypred_paths[pid][repid] - ypred_paths[pid][0]
            if delta_ener >= barrier:
                barrier = delta_ener
        
        ypred_barrier.append(barrier)
    
    return numpy.asarray(ypred_barrier)

def load_grad(sysname, pathname, fgname):
    '''Load the computed gradient.
    '''
    graddir = f'./fg_grad/{sysname}.{pathname}.{fgname}.npy'
    grad = numpy.load(graddir)

    return grad

def load_intergrad(sysname, pathname, fgname):
    '''Load the path-wise cumsum of replica-wise intergrated gradient.
    '''
    graddir = f'./fg_intergrad/{sysname}.{pathname}.{fgname}.npy'
    grad = numpy.load(graddir)

    return grad
