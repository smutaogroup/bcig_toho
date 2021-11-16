# IO options for loading pathway data.
# Zilin Song, 20 AUG 2021
# 

import iomisc, numpy, multiprocessing, dist_compute

def extract_rxc(mda_universe, mda_universe_label, pathname, nrep=50, ):
    '''Extract all reaction coordinates in one RPM coordinate of nrep replicas..
    '''
    rep_distlabel       = []      # Labels
    path_rxc_distmat    = []
    path_rxc_labelmat   = []

    for repid in range(1, nrep+1):
        rep_distlabel.append('{0}.rep{1}'.format(mda_universe_label, str(repid)))
        path_rxc_distrow, path_rxc_labelrow = dist_compute.dist_rx(mda_universe, repid, pathname, )
        path_rxc_distmat.append(path_rxc_distrow)
        path_rxc_labelmat.append(path_rxc_labelrow)
    
    return path_rxc_distmat, path_rxc_labelmat

def check_label(labelmat, ):
	'''Check for consistent labeling.
	'''
	for i in range(len(labelmat)):
		for j in range(len(labelmat[i])):
			if labelmat[0][j] != labelmat[i][j]:
				print('Inequal dist label detected: {0} vs {1} @ f{2}r{3}\n'.format(labelmat[0][j], labelmat[i][j], str(i), str(j), ))
				exit()

def process_conf(sysname, pathname, ):
    '''Process all conformations.
    '''
    rxc_distmat     = []
    rxc_distlabel   = []

    for pid in [i for i in range(1, 101, )]:
        u = iomisc.load_path(sysname, pathname, 'fw', pid, )
        u_label = '{0}.{1}.{2}.path{3}'.format(sysname, pathname, 'fw', str(pid), ) 
        print(u_label, flush=True, )

        path_distmat, path_distlabel = extract_rxc(u, u_label, pathname, )
        rxc_distmat += path_distmat
        rxc_distlabel += path_distlabel
    
    for pid in [i for i in range(1, 101, )]:
        u = iomisc.load_path(sysname, pathname, 'bw', pid, )
        u_label = '{0}.{1}.{2}.path{3}'.format(sysname, pathname, 'bw', str(pid), ) 
        print(u_label, flush=True, )

        path_distmat, path_distlabel = extract_rxc(u, u_label, pathname, )
        rxc_distmat += path_distmat
        rxc_distlabel += path_distlabel

    rxc_distmat = numpy.asarray(rxc_distmat)
    
    print('Finished dist extraction: rxc_distmat.shape = {0}.\n Checking labels...'.format(rxc_distmat.shape))
    check_label(rxc_distlabel)
    print('Labels_checked.\nDONE.')

    numpy.save('./rawds_rxc/{0}.{1}.rxc_distmat.npy'.format(sysname, pathname.replace('-', '')), rxc_distmat, )
    numpy.save('./rawds_rxc/{0}.{1}.rxc_distlabel.npy'.format(sysname, pathname.replace('-', '')), rxc_distlabel[0], )

def main():

    jobs = []

    p0 = multiprocessing.Process(target=process_conf, args=('toho_amp', 'r1ae'))
    p1 = multiprocessing.Process(target=process_conf, args=('toho_amp', 'r2ae'))
    p2 = multiprocessing.Process(target=process_conf, args=('toho_cex', 'r1ae'))
    p3 = multiprocessing.Process(target=process_conf, args=('toho_cex', 'r2ae'))

    jobs.append(p0)
    jobs.append(p1)
    jobs.append(p2)
    jobs.append(p3)

    for p in jobs:
        p.start()

if __name__ == "__main__":
	main()
