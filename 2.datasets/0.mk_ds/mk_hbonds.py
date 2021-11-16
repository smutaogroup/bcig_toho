# IO options for loading pathway data.
# Zilin Song, 20 AUG 2021
# 

import iomisc, numpy, multiprocessing, dist_compute

def launcher_find_hbonds():
    '''Simply a parallel launcher for the nested def.
    '''

    def find_hbonds(sysname, pathname, direction, nrep=50, ):
        '''Detect all possible Hbonds, output as a log .npy .
        Format of the log: donor:hydrogen or hydrogen:acceptor. 
        NOTE that all labels in reaction coordinates are of the format: hydrogen:heavy
        '''
        hbonds_labels = []
        
        for pid in [i for i in range(1, 101, )]:
            u = iomisc.load_path(sysname, pathname, direction, pid, )
            u_label = '{0}.{1}.{2}.path{3}'.format(sysname, pathname, direction, str(pid), ) 
            print(u_label, flush=True, )

            for repid in range(1, nrep+1):
                rep_labels = dist_compute.hbond_detect(u, repid, )
                hbonds_labels.append(rep_labels)

        numpy.save('./rawds_hbonds_detectlabels/detect.{0}.{1}.{2}.hbonds_labels.npy'.format(sysname, pathname, direction, ), 
                    numpy.asarray(hbonds_labels, dtype=object),  # suppress warnings on numpy arrays of inhomogenous shape.
                )

    jobs = []

    p0 = multiprocessing.Process(target=find_hbonds, args=('toho_amp', 'r1ae', 'fw'))
    p1 = multiprocessing.Process(target=find_hbonds, args=('toho_amp', 'r2ae', 'fw'))
    p2 = multiprocessing.Process(target=find_hbonds, args=('toho_cex', 'r1ae', 'fw'))
    p3 = multiprocessing.Process(target=find_hbonds, args=('toho_cex', 'r2ae', 'fw'))
    p4 = multiprocessing.Process(target=find_hbonds, args=('toho_amp', 'r1ae', 'bw'))
    p5 = multiprocessing.Process(target=find_hbonds, args=('toho_amp', 'r2ae', 'bw'))
    p6 = multiprocessing.Process(target=find_hbonds, args=('toho_cex', 'r1ae', 'bw'))
    p7 = multiprocessing.Process(target=find_hbonds, args=('toho_cex', 'r2ae', 'bw'))

    jobs.append(p0)
    jobs.append(p1)
    jobs.append(p2)
    jobs.append(p3)
    jobs.append(p4)
    jobs.append(p5)
    jobs.append(p6)
    jobs.append(p7)

    for p in jobs:
        p.start()

def get_hbond_labels():
    '''Extract all hydrogen bonds labels.
    This label is then used for replica-wise selection.
    '''
    def find_unique(sysname):
        '''Find unique labels in H-bonds
        '''
        hblabels = []
        for p in ['r1ae', 'r2ae', ]:

            for d in ['fw', 'bw', ]:
                lbl_list = numpy.load(
                    './rawds_hbonds_detectlabels/detect.{0}.{1}.{2}.hbonds_labels.npy'.format(
                        sysname, p, d, ), 
                    allow_pickle=True,).tolist()

                for pathlbls in lbl_list:   # for each replca in pathway.

                    for lbl in pathlbls:    # for each selected label in(from) replica
                        revlbl = '{0}:{1}'.format(lbl.split(':')[1], lbl.split(':')[0], )

                        if not (lbl in hblabels or revlbl in hblabels): # filter out identical atom pairs
                            hblabels.append(lbl)
        return hblabels

    amp_labels = find_unique('toho_amp')
    cex_labels = find_unique('toho_cex')

    uni_amplbl, uni_cexlbl = dist_compute.unify_hbond_labels(amp_labels, cex_labels)

    print(len(uni_amplbl))
    print(len(uni_cexlbl))

    for i in range(len(uni_amplbl)):
        print('{}___{}'.format(uni_amplbl[i], uni_cexlbl[i]), flush=True)

    return uni_amplbl, uni_cexlbl

def extract_hbonds(mda_universe, mda_universe_label, bond_labels, nrep=50, ):
    '''Extract all hydrogen donor/acceptor distances according to atom pair 
    specified in bond_labels
    '''
    rep_distlabel = []

    # hydrogen bonds.
    path_bonds_distmat  = []
    path_bonds_labelmat = []
    for repid in range(1, nrep+1):
        rep_distlabel.append('{0}.rep{1}'.format(mda_universe_label, str(repid)))
        path_bonds_distrow, path_bonds_labelrow = dist_compute.dist_hbonds(mda_universe, repid, bond_labels, )
        path_bonds_distmat.append(path_bonds_distrow)
        path_bonds_labelmat.append(path_bonds_labelrow)
    
    return path_bonds_distmat, path_bonds_labelmat

def check_label(labelmat, ):
	'''Check for consistent labeling.
	'''
	for i in range(len(labelmat)):
		for j in range(len(labelmat[i])):
			if labelmat[0][j] != labelmat[i][j]:
				print('Inequal dist label detected: {0} vs {1} @ f{2}r{3}\n'.format(labelmat[0][j], labelmat[i][j], str(i), str(j), ))
				exit()

def process_hbonds(sysname, pathname, bondlabels, ):
    '''Process all Hbond conformations.
    '''
    bonds_distmat   = []
    bonds_distlabel = []

    for pid in [i for i in range(1, 101, )]:
        u = iomisc.load_path(sysname, pathname, 'fw', pid, )
        u_label = '{0}.{1}.{2}.path{3}'.format(sysname, pathname, 'fw', str(pid), )
        print(u_label, flush=True, )

        path_distmat, path_distlabel = extract_hbonds(u, u_label, bondlabels, )
        bonds_distmat += path_distmat
        bonds_distlabel += path_distlabel

    for pid in [i for i in range(1, 101, )]:
        u = iomisc.load_path(sysname, pathname, 'bw', pid, )
        u_label = '{0}.{1}.{2}.path{3}'.format(sysname, pathname, 'bw', str(pid), ) 
        print(u_label, flush=True, )

        path_distmat, path_distlabel = extract_hbonds(u, u_label, bondlabels, )
        bonds_distmat += path_distmat
        bonds_distlabel += path_distlabel

    bonds_distmat = numpy.asarray(bonds_distmat)
    
    print('Finished dist extraction: hbonds_distmat.shape = {0}.\n Checking labels...'.format(bonds_distmat.shape))
    check_label(bonds_distlabel)
    print('Labels_checked.\nDONE.')

    numpy.save('./rawds_hbonds/{0}.{1}.hbonds_distmat.npy'.format(sysname, pathname), bonds_distmat, )
    numpy.save('./rawds_hbonds/{0}.{1}.hbonds_distlabel.npy'.format(sysname, pathname), bonds_distlabel[0], )

if __name__ == "__main__":
    #launcher_find_hbonds()    # initial run,
    #exit()

    amp_lbl, cex_lbl = get_hbond_labels()

    jobs = []

    p0 = multiprocessing.Process(target=process_hbonds, args=('toho_amp', 'r1ae', amp_lbl, ))
    p1 = multiprocessing.Process(target=process_hbonds, args=('toho_amp', 'r2ae', amp_lbl, ))
    p2 = multiprocessing.Process(target=process_hbonds, args=('toho_cex', 'r1ae', cex_lbl, ))
    p3 = multiprocessing.Process(target=process_hbonds, args=('toho_cex', 'r2ae', cex_lbl, ))

    jobs.append(p0)
    jobs.append(p1)
    jobs.append(p2)
    jobs.append(p3)

    for p in jobs:
        p.start()
