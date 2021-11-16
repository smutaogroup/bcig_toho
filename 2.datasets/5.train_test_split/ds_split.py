# group individual replicas into sequences of 3-replicas
# Zilin Song, 22 AUG 2021
# 

from re import X
import numpy, iomisc

nrep   = 50
npath  = 200
seqlen = 3
nseq_per_path = nrep - seqlen + 1

test_seed=[]
for p in range(npath):
    seed = numpy.asarray(numpy.random.choice(50, 5, replace=False))
    test_seed.append(seed+p*50)
    print(seed+p*50)
    
    # numpy.asarray(numpy.random.choice(10000, size=1000, replace=False, )) # random seed.
test_seed = numpy.asarray(test_seed).flatten()
print(test_seed)

for sysname in ['toho_amp', 'toho_cex']:
    for pathname in ['r1ae', 'r2ae']:
        x, xlbl, pid = iomisc.load_x(sysname, pathname)
        y, ylbl      = iomisc.load_y(sysname, pathname)
        
        with open(f'./fin_ds/{sysname}.{pathname}.xlbl.dat', 'w') as logout:
            for i in range(xlbl.shape[0]):
                logout.write(f'{i:>4}\t{xlbl[i]}\n')

        numpy.save(f'./fin_ds/{sysname}.{pathname}.x_test.npy', numpy.asarray(          x[test_seed]))
        numpy.save(f'./fin_ds/{sysname}.{pathname}.y_test.npy', numpy.asarray(          y[test_seed]))
        numpy.save(f'./fin_ds/{sysname}.{pathname}.x_onehot_test.npy', numpy.asarray( pid[test_seed]))
        numpy.save(f'./fin_ds/{sysname}.{pathname}.x_train.npy', numpy.asarray(       numpy.delete(  x, test_seed, axis=0)))
        numpy.save(f'./fin_ds/{sysname}.{pathname}.y_train.npy', numpy.asarray(       numpy.delete(  y, test_seed, axis=0)))
        numpy.save(f'./fin_ds/{sysname}.{pathname}.x_onehot_train.npy', numpy.asarray(numpy.delete(pid, test_seed, axis=0)))
        numpy.save(f'./fin_ds/{sysname}.{pathname}.xlbl.npy', numpy.asarray(     xlbl, dtype=str))
        numpy.save(f'./fin_ds/{sysname}.{pathname}.ylbl.npy', numpy.asarray(     ylbl, dtype=str))

# sysname == 'both'
for pathname in ['r1ae', 'r2ae']:
    x, pid = iomisc.load_x('both', pathname)
    
    y_amp, ylbl_amp = iomisc.load_y('toho_amp', pathname)
    y_cex, ylbl_cex = iomisc.load_y('toho_cex', pathname)

    y = numpy.concatenate((y_amp, y_cex), axis=0)
    ylbl = numpy.concatenate((ylbl_amp, ylbl_cex), axis=0)

    numpy.save(f'./fin_ds/both.{pathname}.x_test.npy', numpy.asarray(          x[test_seed]))
    numpy.save(f'./fin_ds/both.{pathname}.y_test.npy', numpy.asarray(          y[test_seed]))
    numpy.save(f'./fin_ds/both.{pathname}.x_onehot_test.npy', numpy.asarray( pid[test_seed]))
    numpy.save(f'./fin_ds/both.{pathname}.x_train.npy', numpy.asarray(       numpy.delete(  x, test_seed, axis=0)))
    numpy.save(f'./fin_ds/both.{pathname}.y_train.npy', numpy.asarray(       numpy.delete(  y, test_seed, axis=0)))
    numpy.save(f'./fin_ds/both.{pathname}.x_onehot_train.npy', numpy.asarray(numpy.delete(pid, test_seed, axis=0)))
    numpy.save(f'./fin_ds/both.{pathname}.ylbl.npy', numpy.asarray(ylbl, dtype=str))
