# group individual replicas into sequences of 3-replicas
# Zilin Song, 22 AUG 2021
# 

from re import X
import numpy, iomisc

nrep   = 50
npath  = 200
seqlen = 3
nseq_per_path = nrep - seqlen + 1

for sysname in ['toho_amp', 'toho_cex']:
    for pathname in ['r1ae', 'r2ae']:
        x, xlbl, pid = iomisc.load_x(sysname, pathname)
        y, ylbl      = iomisc.load_y(sysname, pathname)
        
        with open(f'./fin_ds/{sysname}.{pathname}.xlbl.dat', 'w') as logout:
            for i in range(xlbl.shape[0]):
                logout.write(f'{i:>4}\t{xlbl[i]}\n')

        numpy.save(f'./fin_ds/{sysname}.{pathname}.x.npy', numpy.asarray(x))
        numpy.save(f'./fin_ds/{sysname}.{pathname}.xlbl.npy', numpy.asarray(xlbl, dtype=str))
        numpy.save(f'./fin_ds/{sysname}.{pathname}.x_onehot.npy', numpy.asarray(pid))
        numpy.save(f'./fin_ds/{sysname}.{pathname}.y.npy', numpy.asarray(y))
        numpy.save(f'./fin_ds/{sysname}.{pathname}.ylbl.npy', numpy.asarray(ylbl, dtype=str))

# sysname == 'both'
for pathname in ['r1ae', 'r2ae']:
    x, pid = iomisc.load_x('both', pathname)
    
    y_amp, ylbl_amp = iomisc.load_y('toho_amp', pathname)
    y_cex, ylbl_cex = iomisc.load_y('toho_cex', pathname)

    y = numpy.concatenate((y_amp, y_cex), axis=0)
    ylbl = numpy.concatenate((ylbl_amp, ylbl_cex), axis=0)

    numpy.save(f'./fin_ds/both.{pathname}.x.npy', numpy.asarray(x))
    numpy.save(f'./fin_ds/both.{pathname}.x_onehot.npy', numpy.asarray(pid))
    numpy.save(f'./fin_ds/both.{pathname}.y.npy', numpy.asarray(y))
    numpy.save(f'./fin_ds/both.{pathname}.ylbl.npy', numpy.asarray(ylbl, dtype=str))
