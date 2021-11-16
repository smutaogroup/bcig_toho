# Normalize the merged features.
# Zilin Song, 21 AUG 2021
# 

import numpy, iomisc

def find_range(x, abs_vals=False, ):
    '''Find the min and max values of the dataset.
    '''
    if not abs_vals:
        return numpy.amin(x, axis=0, ), numpy.amax(x, axis=0, )
    else:
        return numpy.amax(numpy.abs(x), axis=0, )

def rescale(x, x_min, x_max, resc_x_min=-1., resc_x_max=1): 
    '''rescale the feautres. 
    '''
    x_rescale = (x - x_min) / (x_max - x_min) * (resc_x_max - resc_x_min) + resc_x_min
    if numpy.amax(x_rescale) > resc_x_max or numpy.amin(x_rescale) < resc_x_min: 
        raise ValueError('Rescaled values out of range.\nCheck feature alignment.')
    else:
        return x_rescale

def resc_ds(sysname, pathname, nrep=50):
    '''Produce rescaled ds.
    '''
    # xall_rel = numpy.concatenate( (iomisc.load_rel_x('toho_amp', pathname, )[0], iomisc.load_rel_x('toho_cex', pathname, )[0]), axis=0)

    # load ds, x_rel = ds that was made relative to acylenzymes states.  
    x_rel = iomisc.load_rel_x(sysname, pathname, )[0]
    npath = 200
    
    # range of the variables for each feature dim.
    # x_abs_max = find_range(xall_rel, abs_vals=True)
    # x_max = find_range(x_rel, abs_vals=True)
    x_resc = []
    onehot = []
    for p in range(npath):
        xpath = x_rel[p*nrep:(p+1)*nrep]
        # x_min, x_max = xpath[0], xpath[-1]
        # x_min, x_max = find_range(xpath)
        x_max = find_range(xpath, abs_vals=True)
        xpath_resc = rescale(xpath, -x_max, x_max)
        x_resc.extend(xpath_resc)

        # pid labels.
        onehot_path = numpy.zeros( (npath) )
        onehot_path[p] = 1
        onehot.extend([onehot_path for i in range(nrep)])

    x_resc=numpy.asarray(x_resc)
    print(x_resc.shape)
    numpy.save('./conf_ds/{0}.{1}.x.npy'.format(sysname, pathname, ), x_resc, )
    numpy.save('./conf_ds/{0}.{1}.x_onehot.npy'.format(sysname, pathname, ), onehot, )

    # onehot-encoding for both ds
    onehot=[]
    for p in range(npath*2):
        # pid labels.
        onehot_path = numpy.zeros( (npath*2) )
        onehot_path[p] = 1
        onehot.extend([onehot_path for i in range(nrep)])
    
    numpy.save('./conf_ds/both.{0}.x_onehot.npy'.format(pathname), onehot)

for s in ['toho_amp', 'toho_cex']:
    for p in ['r1ae', 'r2ae']:
        resc_ds(s, p, )
