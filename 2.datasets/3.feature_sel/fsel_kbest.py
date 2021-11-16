# Normalize the merged features.
# Zilin Song, 22 AUG 2021
# 

import numpy, iomisc
from numpy.lib.function_base import average
from numpy.core.fromnumeric import var
from sklearn.feature_selection import SelectKBest, mutual_info_regression

def kbest_mi_selection(x, y, nfeatures, ):
    # feature selection.
    selector = SelectKBest(score_func=mutual_info_regression, k=nfeatures, )
    selector.fit(x, y, )
    selected_idx = selector.get_support(indices=True)
    
    return selected_idx

def aver_mi_per_path(x, y, npath=200, nrep=50):
    '''compute the averaged MI per path'''
    mi = []
    for p in range(npath):
        xpath  = x[p*nrep: p*nrep+nrep]
        ypath  = y[p*nrep: p*nrep+nrep]
        mipath = mutual_info_regression(xpath, ypath)
        mi.append(mipath)

    mi = numpy.asarray(mi)
    print(mi.shape)
    aver_mi = numpy.average(mi, axis=0)
    return aver_mi

def aver_var_per_path(x, thresh, npath=200, nrep=50):
    variance = []
    for p in range(npath):
        xpath  = x[p*nrep: p*nrep+nrep]
        varpath = numpy.var(xpath, axis=0)
        variance.append(varpath)

    variance = numpy.asarray(variance)
    print(variance.shape)
    aver_variance = numpy.average(variance, axis=0)
    sel_idx = numpy.argwhere(aver_variance>thresh)

    # print(sel_idx)
    return sel_idx

def select_features(pathname, ):
    '''Feature selection on relative ds. 
    '''
    # load raw data.
    amp_x, amp_xlbl,  = iomisc.load_x('toho_amp', pathname, 'raw')
    cex_x, cex_xlbl,  = iomisc.load_x('toho_cex', pathname, 'raw')
    amp_y             = iomisc.load_y('toho_amp', pathname, )[0]
    cex_y             = iomisc.load_y('toho_cex', pathname, )[0]

    # feature selection: remove low variances with raw data.
    lv_amp_idx = aver_var_per_path(amp_x, 0.03)
    lv_cex_idx = aver_var_per_path(cex_x, 0.03)
    lv_idx = numpy.unique(numpy.concatenate( (lv_amp_idx, lv_cex_idx) ) ) # merge and concat retained features.
    # lv_idx = numpy.intersect1d(lv_amp_idx, lv_cex_idx)

    # apply low variance selection.
    amp_x, amp_xlbl,  = iomisc.load_x('toho_amp', pathname, 'raw')
    cex_x, cex_xlbl,  = iomisc.load_x('toho_cex', pathname, 'raw')
    lv_amp_x    = numpy.take(  amp_x, lv_idx, axis=1)
    lv_cex_x    = numpy.take(  cex_x, lv_idx, axis=1)
    lv_amp_xlbl = numpy.take(amp_xlbl, lv_idx)
    lv_cex_xlbl = numpy.take(cex_xlbl, lv_idx)
    print('{}'.format(lv_cex_xlbl.shape))

    # feature selection: based on mutual informaiton.
    nfeatures = 15 if pathname == 'r1ae' else 11
    mi_amp_idx = kbest_mi_selection(lv_amp_x, amp_y, nfeatures)
    mi_cex_idx = kbest_mi_selection(lv_cex_x, cex_y, nfeatures)
    # mi_idx = numpy.unique(numpy.concatenate( (mi_amp_idx, mi_cex_idx, ) ) )
    mi_idx = numpy.intersect1d(mi_amp_idx, mi_cex_idx)
    if pathname == 'r1ae':mi_idx = numpy.sort(numpy.append(mi_idx, 8)) # feature 8 is manually appended: Ser130 HG1 - beta-lactam O10B/O12. H atoms are always selected.
    
    # apply selections (lv_idx, mi_idx) to the unselected ds. 
    amp_x, amp_xlbl,  = iomisc.load_x('toho_amp', pathname, 'resc')
    cex_x, cex_xlbl,  = iomisc.load_x('toho_cex', pathname, 'resc')
    
    # apply lv_idx
    lv_amp_x    = numpy.take(   amp_x, lv_idx, axis=1)
    lv_cex_x    = numpy.take(   cex_x, lv_idx, axis=1) 
    lv_amp_xlbl = numpy.take(amp_xlbl, lv_idx)
    lv_cex_xlbl = numpy.take(cex_xlbl, lv_idx)
    # apply mi_idx
    mi_amp_x    = numpy.take(lv_amp_x,    mi_idx, axis=1)
    mi_cex_x    = numpy.take(lv_cex_x,    mi_idx, axis=1)
    mi_amp_xlbl = numpy.take(lv_amp_xlbl, mi_idx)
    mi_cex_xlbl = numpy.take(lv_cex_xlbl, mi_idx)

    print('{}'.format(mi_amp_xlbl))
    print('{}'.format(mi_cex_xlbl))
    print('{}\n'.format(mi_cex_xlbl.shape))

    # save dat.
    numpy.save('./conf_ds/{0}.{1}.x.npy'.format(   'toho_amp', pathname, ), mi_amp_x   )
    numpy.save('./conf_ds/{0}.{1}.x.npy'.format(   'toho_cex', pathname, ), mi_cex_x   )
    numpy.save('./conf_ds/{0}.{1}.x.npy'.format(       'both', pathname, ), numpy.concatenate( (mi_amp_x, mi_cex_x ), axis=0) )
    numpy.save('./conf_ds/{0}.{1}.xlbl.npy'.format('toho_amp', pathname, ), mi_amp_xlbl)
    numpy.save('./conf_ds/{0}.{1}.xlbl.npy'.format('toho_cex', pathname, ), mi_cex_xlbl)

if __name__ == '__main__':
    select_features('r1ae')
    select_features('r2ae')
