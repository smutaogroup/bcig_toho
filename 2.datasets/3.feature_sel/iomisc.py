# IO defs for loading pathway data.
# Zilin Song, 22 AUG 2021
# 

import numpy

def load_x(sysname, pathname, whichds='raw', ):
    '''Load various normalized/rescaled datasets and the corresponding labels from 
    sysname:    toho_amp, toho_cex;
    pathname:   r1ae, r2ae;
    whichds:    norm, resc, raw
    '''
    # directory settings.
    # _datprefx = '../2.norm_ds/conf_ds_{}'.format(whichds)
    
    _datprefx = '../1.merge_ds/conf_ds' if whichds ==  'raw' else \
                '../2.resc_ds/conf_ds'  if whichds == 'resc' else \
                None
    _lblrepfx = '../1.merge_ds/conf_ds'
    _datname = 'rel_x' if whichds == 'raw' else 'x'
    _lblname = 'rel_xlbl'

    _datdir = '{0}/{1}.{2}.{3}.npy'.format(_datprefx, sysname, pathname, _datname, )
    _lbldir = '{0}/{1}.{2}.{3}.npy'.format(_lblrepfx, sysname, pathname, _lblname, )

    if sysname == 'both':
        return numpy.load(_datdir)
    else:
        return numpy.load(_datdir), numpy.load(_lbldir)

def load_y(sysname, pathname, ):
    '''Load the energies from
    ../1.merge_ds/conf_ds
    '''
    _basedir = '../1.merge_ds/ener_ds'
    _datname = 'rel_y'
    _lblname = 'rel_ylbl'

    _datdir = '{0}/{1}.{2}.{3}.npy'.format(_basedir, sysname, pathname, _datname, )
    _lbldir = '{0}/{1}.{2}.{3}.npy'.format(_basedir, sysname, pathname, _lblname, )

    dat = numpy.load(_datdir)
    lbl = numpy.load(_lbldir)

    return dat, lbl
