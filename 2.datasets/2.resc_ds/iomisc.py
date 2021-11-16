# IO defs for loading pathway data.
# Zilin Song, 21 AUG 2021
# 

import numpy

def load_rel_x(sysname, pathname, ):
    '''Load various relative datasets and the corresponding labels from 
        ../1.merge_ds
    sysname:    toho_amp, toho_cex;
    pathname:   r1ae, r2ae;
    '''
    _basedir = '../1.merge_ds/conf_ds'
    _datname = 'rel_x'
    _lblname = 'rel_xlbl'

    _datdir = f'{_basedir}/{sysname}.{pathname}.{_datname}.npy'
    _lbldir = f'{_basedir}/{sysname}.{pathname}.{_lblname}.npy'

    dat = numpy.load(_datdir)
    lbl = numpy.load(_lbldir)

    return dat, lbl

def load_rel_y(sysname, pathname, ):
    '''Load the energies from
    ../1.merge_ds/conf_ds
    '''
    _basedir = '../1.merge_ds/ener_ds'
    _datname = 'rel_y'
    _lblname = 'rel_ylbl'

    _datdir = f'{_basedir}/{sysname}.{pathname}.{_datname}.npy'
    _lbldir = f'{_basedir}/{sysname}.{pathname}.{_lblname}.npy'

    dat = numpy.load(_datdir)
    lbl = numpy.load(_lbldir)

    return dat, lbl
