# Merge features and make relative values.
# Zilin Song, 21 AUG 2021
# 

import numpy, iomisc

def make_relative_ds(ds, ref_offset, nrep=50):
    '''Process the ds into relative dist data.
    ds_elabels      -> labels of replica energies.
    ref_offset = 47 -> last replica (AE) as reference.
    '''
    relative_ds = []
    for rid in range(ds.shape[0]):
        if rid % nrep == 0:
            refdistrow = ds[ rid + ref_offset ].copy()    # use Last replica (AE) as the reference.
        # refdistrow = ds[ref_offset].copy()              # use one replica as the global reference.
        distrow = numpy.round_(ds[rid] - refdistrow, 8, )
        relative_ds.append(distrow)

    return numpy.asarray(relative_ds)

def merge_ds(sysname, pathname, ):
    '''Merge the x and y datasets. 
    Also make them relative data to the AE state.
    '''
    # outputs
    y, ylbl = iomisc.load_raw_ds(sysname, pathname, 'ener', )
    ref_offset = 49  # last replica is ref 
    rel_y = make_relative_ds(y, ref_offset)

    numpy.save('./ener_ds/{0}.{1}.rel_y.npy'.format(sysname, pathname, ), rel_y)
    numpy.save('./ener_ds/{0}.{1}.rel_ylbl.npy'.format(sysname, pathname, ), ylbl)
    
    # features
    x_rc, xlbl_rc = iomisc.load_raw_ds(sysname, pathname, 'rxc', )
    x_cb, xlbl_cb = iomisc.load_raw_ds(sysname, pathname, 'chembonds', )
    x_hb, xlbl_hb = iomisc.load_raw_ds(sysname, pathname, 'hbonds', )
    # merge features
    x     = numpy.concatenate((   x_rc,    x_cb,    x_hb, ), axis=1, )
    xlbl  = numpy.concatenate((xlbl_rc, xlbl_cb, xlbl_hb, ), axis=0, )
    rel_x = make_relative_ds(x, ref_offset)

    numpy.save('./conf_ds/{0}.{1}.rel_x.npy'.format(sysname, pathname, ), rel_x)
    numpy.save('./conf_ds/{0}.{1}.rel_xlbl.npy'.format(sysname, pathname, ), xlbl)
    numpy.save('./raw_ds/{0}.{1}.raw_x.npy'.format(sysname,pathname,), x)
    numpy.save('./raw_ds/{0}.{1}.raw_xlbl.npy'.format(sysname, pathname, ), xlbl)

for s in ['toho_amp', 'toho_cex']:
    for p in ['r1ae', 'r2ae', ]:
        merge_ds(s, p, )
