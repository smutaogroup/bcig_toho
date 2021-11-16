# Integrate the feature gradients path-wise and exponentially average their sum.
# Zilin Song, 29 Aug 2021
# 

import iomisc, numpy

'''Summation of the Boltzmann factor-weighted fg contribution.
'''

for pathname in ['r1ae', 'r2ae', ]:
    for sysname in ['toho_amp', 'toho_cex', 'both']:
        fgroups = ['fg{}'.format(str(i)) for i in [1, 2, 3, 4,    6, 7, ] ] if pathname == 'r1ae' \
             else ['fg{}'.format(str(i)) for i in [1, 2,       5, 6, 7, ] ]

        for fg in fgroups:
                
            fg_grad_per_path = iomisc.load_grad(sysname, pathname, fg)      # shape (200,48) for toho_amp/cex
            # also sum the last acyl-enzyme replica.
            fg_grad_path_integrated = numpy.sum(numpy.abs(numpy.cumsum(fg_grad_per_path, axis=1)), axis=1) + numpy.abs(numpy.sum(fg_grad_per_path, axis=1))   # shape (200)    for toho_amp/cex
            numpy.save(f'./fg_intergrad/{sysname}.{pathname}.{fg}.npy', fg_grad_path_integrated)
