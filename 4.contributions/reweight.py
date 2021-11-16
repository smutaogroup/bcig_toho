# Integrate the feature gradients path-wise and exponentially average their sum.
# Zilin Song, 30 Aug 2021
# 

import iomisc, numpy

'''Summation of the Boltzmann factor-weighted fg contribution.
'''
def gaussian_prob(arr, normalize=True):
    '''Return the Gaussion probabilities of each array element.
    '''
    # Note that in a complete sampled Gaussian, median = mean. 
    # Due to the under-sampling, the average of mean and median values should NOT be taken:
    #     It biases the real mean towards the pseudo-median.
    mean = numpy.mean(arr)
    std  = numpy.std(arr)

    gauss_prob = 1 / (std * numpy.sqrt(numpy.pi*2)) * numpy.exp( -0.5*numpy.square( (arr-mean)/std ) )
    
    if normalize: # return the normalized (discrete) Gaussian probs
        return gauss_prob / numpy.sum(gauss_prob)
    else: 
        return gauss_prob

def boltzmann_prob(arr, normalize=True):
    '''Return the Boltzmann probabilities of each array element.
    '''
    Kb = 0.001987204258 # kcal/mol/K, Boltzmann Const equals R - gas constant.
    T  = 310            # K         , Temperature

    # take prob relative to the lowest barrier
    # raw barriers give very small prob, leading numerical errors.
    b_min = numpy.min(arr)
    bolzm_prob = numpy.exp( (b_min-arr) / Kb / T )

    if normalize: # return the normalized Boltzmann prob.
        return bolzm_prob / numpy.sum(bolzm_prob)
    else: 
        return bolzm_prob

def merged_reweighting_probabiligy(arr):
    '''Merge the boltzmann and gaussian to one expresion
    '''
    mean = numpy.mean(arr)
    std  = numpy.std(arr)
    b_min = numpy.min(arr)
    Kb = 0.001987204258 # kcal/mol/K, Boltzmann Const. equals R - gas constant.
    T  = 310            # K         , Temperature
    prob =  1 / (std * numpy.sqrt(numpy.pi*2)) * numpy.exp( (b_min-arr) / Kb / T - numpy.square( (arr-mean)/std ) / 2. )
    return prob / numpy.sum(gaussian_prob(arr, normalize=False)) / numpy.sum(boltzmann_prob(arr, normalize=False))

def reweight_contribution(barriers, cumsum_igs):
    '''reweight the contribution of each pathway using Gaussian probability,
    and then calculate the Boltzmann weighted contribution.
    '''
    gauss_prob = gaussian_prob(barriers)
    bolzm_prob = boltzmann_prob(barriers)
    merged_prob = merged_reweighting_probabiligy(barriers)

    # Normalization factor = numpy.sum(bolzm_prob) * numpy.sum(gauss_prob), 
    # included in the weights so as to avoid numerical problems.
    
    # Followings to ensure correctness.
    weighted_cumsum_igs = numpy.sum(cumsum_igs * gauss_prob * bolzm_prob)  # step-wise 
    # weighted_cumsum_igs = numpy.sum(cumsum_igs * merged_prob)  # from equation 3 
    return weighted_cumsum_igs

logout = open('contribution.log', 'w')

for pathname in ['r1ae', 'r2ae', ]:
    for sysname in ['toho_amp', 'toho_cex', 'both']:
        fgroups = ['fg{}'.format(str(i)) for i in [1, 2, 3, 4,    6, 7, ] ] if pathname == 'r1ae' \
             else ['fg{}'.format(str(i)) for i in [1, 2,       5, 6, 7, ] ]
        
        logout.write(f'\n{sysname}: {pathname}\n\n')

        logstr1 = ''
        logstr2 = ''
        for fg in fgroups:
            cumsum_igs = iomisc.load_intergrad(sysname, pathname, fg)
            barriers   = iomisc.load_pred_barriers(sysname, pathname)

            if sysname != 'both':   # toho_amp or toho_cex, individually
                reweighted_cumsum_igs = reweight_contribution(barriers, cumsum_igs)
                logstr1 +=f'\t{fg}\t\t{reweighted_cumsum_igs}\n'

            else:                   # both. need to first split to toho_amp, toho_cex subsets.
                amp_reweighted_cumsum_igs = reweight_contribution(barriers[:200], cumsum_igs[:200])
                logstr1 += f'\tamp_{fg}\t\t{amp_reweighted_cumsum_igs}\n'
                cex_reweighted_cumsum_igs = reweight_contribution(barriers[200:], cumsum_igs[200:])
                logstr2 += f'\tcex_{fg}\t\t{cex_reweighted_cumsum_igs}\n'
        
        logout.write(logstr1 +'\n'+ logstr2)


            
