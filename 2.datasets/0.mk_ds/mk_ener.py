# IO options for loading pathway data.
# Zilin Song, 20 AUG 2021
# 

import iomisc, numpy

def extract_ener(sysname, pathname, theoryname='b3lyp', ):
    '''Extract all 198 pathways for one mechanism.
    '''
    # directories.
    _fw_edir, _bw_edir, _fw_bdir, _bw_bdir = iomisc.enerdirs(sysname, pathname, theoryname)

    # pathway energies loading...
    rep_enermat = []
    rep_enerlabel = []

    for _edir, direction in [(_fw_edir, 'fw', ), (_bw_edir, 'bw'), ]:
        lines = open(_edir, 'r', ).readlines()
        for line in lines:
            words = line.split()
            
            if words[0][0:2] == 'f:':
                rep_enermat.append(float(words[2]))
                rep_enerlabel.append('{0}.{1}.{2}.path{3}.rep{4}'.format(
                        sysname, pathname, direction, str(words[0].split(':')[1]), str(words[1].split(':')[1]), 
                    )
                )

    # pathway barriers loading...
    path_barrier_mat = []

    for _bdir, direction in [(_fw_bdir, 'fw', ), (_bw_bdir, 'bw'), ]:

        lines = open(_bdir, 'r', ).readlines()
        for line in lines:
            words = line.split()
            path_barrier_mat.append(float(words[1]))

    return rep_enermat, rep_enerlabel, path_barrier_mat

def process_ener(sysname, pathname, ):
    '''Produce ener datasets.
    '''
    # energies 
    emat, elabel, barriermat = extract_ener(sysname, pathname, )

    # save ds: energies. Note that enerpathlabel should be exactly the same with distpathlabel. 
    numpy.save('./rawds_ener/{0}.{1}.enermat.npy'.format(sysname, pathname.replace('-', '')), emat, )
    numpy.save('./rawds_ener/{0}.{1}.enerlabel.npy'.format(sysname, pathname.replace('-', '')), elabel, )
    numpy.save('./rawds_ener/{0}.{1}.barriermat.npy'.format(sysname, pathname.replace('-', '')), barriermat, )
    
if __name__ == "__main__":
    for s in ['toho_amp', 'toho_cex', ]:
        for p in ['r1ae', 'r2ae', ]:
            process_ener(s, p, )
