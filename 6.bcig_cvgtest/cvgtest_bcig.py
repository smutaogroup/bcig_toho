# Convergence test of the BCIG contributions.
# Zilin Song, Jan 3 2022
# 

import iomisc, numpy, multiprocessing, sys
from sklearn import mixture

def get_path(ds, pathid, nrep):
    '''Return the pathway selected via pathid.'''
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
    f_x_dxp = numpy.squeeze( model.predict( (x_dxp, xpoh) ), axis=1)
    # f_x     = numpy.squeeze( model.predict( (x,     xpoh) ), axis=1)

    # Partial gradients with numpy.
    _grad = f_x_dxp - f_x_dxm
    
    return _grad / 2. / pert

def boltzmann_prob(arr, normalize=True):
    '''Return the Boltzmann probabilities of each array element.'''
    Kb = 0.001987204258 # kcal/mol/K, Boltzmann Const equals R - gas constant.
    T  = 310            # K         , Temperature
    bolzm_prob = numpy.exp( (-arr) / Kb / T )

    if normalize: # return the normalized Boltzmann prob.
        return bolzm_prob / numpy.sum(bolzm_prob)
    else: 
        return bolzm_prob

def gaussian_prob(arr, normalize=True):
    '''Return the Gaussion probabilities of each array element.'''
    mean = numpy.mean(arr)
    std  = numpy.std(arr)

    gauss_prob = 1 / (std * numpy.sqrt(numpy.pi*2)) * numpy.exp( -0.5*numpy.square( (arr-mean)/std ) )
    
    if normalize: # return the normalized (discrete) Gaussian probs
        return gauss_prob / numpy.sum(gauss_prob)
    else: 
        return gauss_prob

def gaussian_mixture_prob(arr, normalize=True):
    """Density estimation using Gaussian Mixture Models."""
    arr = numpy.reshape(arr, (-1, 1))
    gmm = mixture.GaussianMixture(n_components=2)
    gmm.fit(arr)
    gmm_prob = numpy.exp(gmm.score_samples(arr))
    if normalize: # return the normalized GMM prob.
        return gmm_prob / numpy.sum(gmm_prob)
    else: 
        return gmm_prob

def compute_bcig(sysname, pathname, fg, ipath, nrep=50):
    '''Compute the path-wise numerical gradients of npaths selected with numpy.rand.'''
    x, xpoh, y = iomisc.load_ds(sysname, pathname)
    model = iomisc.load_model(sysname, pathname)
    fg_mask = get_fgmask(fg, pathname, x[0])[0]
    fg_gradients = []

    for pathid in ipath:
        x_path, xpoh_path = get_path(x, pathid, nrep), get_path(xpoh, pathid, nrep)
        fg_grad = get_gradient_along_path(x_path, xpoh_path, model, fg_mask)
        fg_gradients.append(fg_grad)
        # verbose.
        # print(f'Done: {sysname:>10} {pathname} {fg} path{str(pathid):3} :' \
        #     'shape({len(fg_gradients)}, {len(fg_gradients[0])})', flush=True)

    cummulative_integrads = numpy.sum(numpy.abs(numpy.cumsum(fg_gradients, axis=1)), axis=1) \
                          + numpy.abs(numpy.sum(fg_gradients, axis=1))
    barriers = iomisc.load_pred_barriers(sysname, pathname)[ipath]
    # gauss_prob = gaussian_prob(barriers)
    gauss_prob = gaussian_mixture_prob(barriers)
    bolzm_prob = boltzmann_prob(barriers)

    bcig = numpy.sum(cummulative_integrads * gauss_prob * bolzm_prob)
    return bcig

def bootstrapping_bcig(sysname, pathname, npath, ntrial=10):
    """Perform bootstrapping on bcig convergence tests."""
    fgroups = ['fg{}'.format(str(i)) for i in [1, 2, 3, 4,    6, 7, ] ] if pathname == 'r1ae' \
         else ['fg{}'.format(str(i)) for i in [1, 2,       6, 7, 5, ] ] # note that reordered 6, 7, 5, for plotting.

    bcig = []

    for fg in fgroups:
        fg_bcig = []

        for _ in range(ntrial):
            ipath = numpy.random.choice(200, size=npath, replace=False)
            bcig_ipath = compute_bcig(sysname, pathname, fg, ipath) * npath # scale by the number of pathways so that the magnitude are compatible.
            # verbose.
            print(f'Done: {sysname:>10} {pathname} {fg} npath_{npath} trail_{_} {bcig_ipath}', flush=True)

            fg_bcig.append(bcig_ipath)

        bcig.append(fg_bcig)
    
    numpy.save(f'./bcigs/{sysname}.{pathname}.bootstrap_{npath}.npy', numpy.asarray(bcig))

if __name__ == '__main__':
    ds_idx = 3 # int(sys.argv[1])       # which dataset.
    sysnames  = ['toho_amp', 'toho_cex', 'toho_amp', 'toho_cex', ]
    pathnames = ['r1ae',     'r1ae',     'r2ae',     'r2ae',     ]
    sysname   = sysnames[ds_idx]
    pathname  = pathnames[ds_idx]

    for npath in [10*i+100 for i in range(0, 11)]:
        p = multiprocessing.Process(target=bootstrapping_bcig, args=(sysname, pathname, npath))
        p.start()
