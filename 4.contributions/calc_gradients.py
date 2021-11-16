# Compute the gradients at each replica.
# Zilin Song, 29 AUG 2021
# 

import iomisc, numpy, multiprocessing

def get_path(ds, pathid, nrep):
    '''Return the pathway selected via pathid.
    '''
    return ds[pathid*nrep:(pathid+1)*nrep].copy()

def get_fgmask(fgname, pathname, xrow):
    '''Return feature indices w.r.t. names of feature group.
        fg1, S70-LIGAND bond forming;
        fg2, LIGAND bond breaking;
      # fg3, S70-Wat proton transfer;
      # fg4, Wat-E166 proton transfer;
      * fg5, S70-K73 proton transfer;
        fg6, K73-S130 proton transfer;
        fg7, S130-LIGAND proton transfer.
        
    #: only in r1ae pathways;
    *: only in r2ae pathways.
    '''
    # definitions of feature groups.
    #  idx defined corresponding to 
    #       /users/zilins/scratch/2.proj_toho2lig_acy/2.datasets/4.conclude_ds/fin_ds/*dat
    group_dict = {
        'fg1': [10,                 ], 
        'fg2': [ 9,                 ], 
        'fg3': [ 0,  1,             ], 
        'fg4': [ 2,  3,             ], 
        'fg6': [ 4,  5, 12,         ], 
        'fg7': [ 6,  7,  8, 11, 13, ], 
    } if pathname == 'r1ae' else {
        'fg1': [ 8,                 ], 
        'fg2': [ 7,                 ], 
        'fg5': [ 0,  1,             ], 
        'fg6': [ 2,  3,             ], 
        'fg7': [ 4,  5,  6,  9, 10, ], 
    } if pathname == 'r2ae' else None

    # feature indices in fgroup
    indices = group_dict[fgname]

    # make mask: selected are 1, unselected are 0.
    _mask = numpy.zeros(xrow.shape[0])
    _mask[indices] = 1

    # also show the unmasked labels
    _masked_labels = xrow[_mask == 1]
    return _mask, _masked_labels

def get_gradient_along_path(x_path, xpoh_path, model, fg_mask, pert=0.01):
    '''Compute the gradient per path of n replicas.
    '''
    # Note that the end-point replicas are discarded: range(50) - [1,48]
    # The gradients at end-points (R/TS) are considered zero.
    xm = x_path[0: -2]
    x  = x_path[1: -1]
    xp = x_path[2:   ]
    xpoh = xpoh_path[1: -1]

    # Compute perturbations on target features belong to the feature group.
    # Note that other features are retained with a perturbation of zero.
    d_xp_masked = pert * (xp-x) * fg_mask       #  forward pert
    d_xm_masked = pert * (x-xm) * fg_mask       # backward pert

    # Perturbed x.
    x_dxm, x_dxp = x-d_xm_masked, x+d_xp_masked

    # Predicted perturbed x.
    f_x_dxm = numpy.squeeze( model.predict( (x_dxm, xpoh) ), axis=1)
    f_x     = numpy.squeeze( model.predict( (x,     xpoh) ), axis=1)
    f_x_dxp = numpy.squeeze( model.predict( (x_dxp, xpoh) ), axis=1)

    # Partial gradients with numpy.
    _grad = []
    _grad = f_x_dxp - f_x_dxm
    
    return _grad / 2. / pert

def compute_gradients(sysname, pathname, nrep=50):
    '''Compute the path-wise numerical gradients and store it as numpy.array().
    '''
    npath = 400 if sysname == 'both' else 200
    x, xpoh, y = iomisc.load_ds(sysname, pathname)
    model = iomisc.load_model(sysname, pathname)

    # definition of feature groups.
    fgroups = ['fg{}'.format(str(i)) for i in [1, 2, 3, 4,    6, 7, ] ] if pathname == 'r1ae' \
         else ['fg{}'.format(str(i)) for i in [1, 2,       5, 6, 7, ] ]

    # iterate over each feature group.
    for fg in fgroups:
        fg_mask = get_fgmask(fg, pathname, x[0])[0]

        fg_gradients = []       # The feature group gradients per path of nrep replicas. 
        
        # iteratively compute the fg gradients per path.
        for pathid in range(npath):
            x_path, xpoh_path = get_path(x, pathid, nrep), get_path(xpoh, pathid, nrep)

            fg_grad = get_gradient_along_path(x_path, xpoh_path, model, fg_mask)
            fg_gradients.append(fg_grad)
            # verbose.
            print(f'Done: {sysname:>10} {pathname} {fg} path{str(pathid):3} : shape({len(fg_gradients)}, {len(fg_gradients[0])})', flush=True)

        numpy.save(f'./fg_grad/{sysname}.{pathname}.{fg}.npy', fg_gradients)

for s in ['toho_amp', 'toho_cex', 'both']:
    for p in ['r1ae', 'r2ae']:
        p = multiprocessing.Process(target=compute_gradients, args=(s, p))
        p.start()
