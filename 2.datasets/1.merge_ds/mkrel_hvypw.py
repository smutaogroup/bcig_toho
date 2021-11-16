# Merge features and make relative values.
# Zilin Song, 21 AUG 2021
# 
 
import numpy
ampr1 = numpy.load('../0.mk_ds/rawds_hvypw/toho_amp.r1ae.hvypw_distmat.npy')
ampr2 = numpy.load('../0.mk_ds/rawds_hvypw/toho_amp.r2ae.hvypw_distmat.npy')
cexr1 = numpy.load('../0.mk_ds/rawds_hvypw/toho_cex.r1ae.hvypw_distmat.npy')
cexr2 = numpy.load('../0.mk_ds/rawds_hvypw/toho_cex.r2ae.hvypw_distmat.npy')

pw_all = numpy.concatenate( (ampr1, cexr1, ampr2, cexr2), axis=0, )

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

rel_pw_all = make_relative_ds(pw_all, 49)

print(rel_pw_all.shape)

numpy.save('./hvypw_rel/hvypw_rel_ds.npy', rel_pw_all)