# IO defs for loading pathway data.
# Zilin Song, 21 AUG 2021
# 

import numpy

def load_raw_ds(sysname, pathname, dsname, ):
    '''Load various raw datasets and the corresponding labels from 
        ../0.mk_ds
    sysname:    toho_amp, toho_cex;
    pathname:   r1ae, r2ae;
    dsname:     rxc, chembonds, hbonds, ener.
    '''
    _basedir = '../0.mk_ds'
    _dsdir   = 'rawds_{0}'.format(dsname)
    _datname = '{0}_distmat'.format(dsname)   if dsname in ['rxc', 'chembonds', 'hbonds', ] else \
               '{0}mat'.format(dsname)        if dsname in ['ener', ]                       else \
                None

    _lblname = '{0}_distlabel'.format(dsname) if dsname in ['rxc', 'chembonds', 'hbonds', ] else \
               '{0}label'.format(dsname)      if dsname in ['ener', ]                       else \
                None

    if _datname == None or _lblname == None: raise ValueError('Unexpected dsname: {0}'.format(dsname))

    _datdir = '{0}/{1}/{2}.{3}.{4}.npy'.format(_basedir, _dsdir, sysname, pathname, _datname, )
    _lbldir = '{0}/{1}/{2}.{3}.{4}.npy'.format(_basedir, _dsdir, sysname, pathname, _lblname, )

    dat = numpy.load(_datdir)
    lbl = numpy.load(_lbldir)

    return dat, lbl
