# IO defs for loading pathway data.
# Zilin Song, 21 AUG 2021
# 

import numpy

def load_ds(sysname, pathname,):
    '''Load various normalized/rescaled datasets and the corresponding labels from 
    sysname:    toho_amp, toho_cex;
    pathname:   r1ae, r2ae;
    '''
    # directory settings.
    prefx    = '/users/zilins/scratch/2.proj_toho2lig_acy/2.datasets/5.train_test_split/fin_ds'

    xtrain_name    = 'x_train'
    xpohtrain_anme = 'x_onehot_train'
    ytrain_name    = 'y_train'

    xtest_name    = 'x_test'
    xpohtest_anme = 'x_onehot_test'
    ytest_name    = 'y_test'

    xtrain_dir    = f'{prefx}/{sysname}.{pathname}.{xtrain_name}.npy'
    xpohtrain_dir = f'{prefx}/{sysname}.{pathname}.{xpohtrain_anme}.npy'
    ytrain_dir    = f'{prefx}/{sysname}.{pathname}.{ytrain_name}.npy'

    xtest_dir     = f'{prefx}/{sysname}.{pathname}.{xtest_name}.npy'
    xpohtest_dir  = f'{prefx}/{sysname}.{pathname}.{xpohtest_anme}.npy'
    ytest_dir     = f'{prefx}/{sysname}.{pathname}.{ytest_name}.npy'

    xtrain    = numpy.load(xtrain_dir)
    xpohtrain = numpy.load(xpohtrain_dir)
    ytrain    = numpy.load(ytrain_dir)

    xtest    = numpy.load(xtest_dir)
    xpohtest = numpy.load(xpohtest_dir)
    ytest    = numpy.load(ytest_dir)

    return xtrain, xpohtrain, ytrain, xtest, xpohtest, ytest