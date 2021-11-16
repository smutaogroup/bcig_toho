# IO options for loading pathway data.
# Zilin Song, 20 AUG 2021
# 

import iomisc, numpy, multiprocessing, dist_compute

def extract_conf(mda_universe, mda_universe_label, nrep=50, ):
    '''Extract all bonds (chemical/hydrogen) in one RPM coordinate of nrep replicas..
    '''
    rep_distlabel       = []      # Labels

    # chemical bonds.
    path_bonds_distmat  = []
    path_bonds_labelmat = []
    for repid in range(1, nrep+1):
        print(mda_universe_label+'.rep'+str(repid), flush=True, )

        rep_distlabel.append('{0}.rep{1}'.format(mda_universe_label, str(repid)))
        path_bonds_distrow, path_bonds_labelrow = dist_compute.dist_bonds(mda_universe, repid, )
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

def process_chembonds(sysname, pathname, ):
    '''Process all conformations.
    '''
    bonds_distmat     = []
    bonds_distlabel   = []

    for pid in [i for i in range(1, 101, )]:
        u = iomisc.load_path(sysname, pathname, 'fw', pid, )
        u_label = '{0}.{1}.{2}.path{3}'.format(sysname, pathname, 'fw', str(pid), ) 
        print(u_label, flush=True, )

        path_distmat, path_distlabel = extract_conf(u, u_label, )
        bonds_distmat += path_distmat
        bonds_distlabel += path_distlabel
    
    for pid in [i for i in range(1, 101, )]:
        u = iomisc.load_path(sysname, pathname, 'bw', pid, )
        u_label = '{0}.{1}.{2}.path{3}'.format(sysname, pathname, 'bw', str(pid), ) 
        print(u_label, flush=True, )

        path_distmat, path_distlabel = extract_conf(u, u_label, )
        bonds_distmat += path_distmat
        bonds_distlabel += path_distlabel

    bonds_distmat = numpy.asarray(bonds_distmat)
    
    print('Finished dist extraction: chembonds_distmat.shape = {0}.\n Checking labels...'.format(bonds_distmat.shape))
    check_label(bonds_distlabel)
    print('Labels_checked.\nDONE.')

    numpy.save('./rawds_chembonds/{0}.{1}.chembonds_distmat.npy'.format(sysname, pathname), bonds_distmat, )
    numpy.save('./rawds_chembonds/{0}.{1}.chembonds_distlabel.npy'.format(sysname, pathname), bonds_distlabel[0], )

def main():

    jobs = []

    p0 = multiprocessing.Process(target=process_chembonds, args=('toho_amp', 'r1ae'))
    p1 = multiprocessing.Process(target=process_chembonds, args=('toho_amp', 'r2ae'))
    p2 = multiprocessing.Process(target=process_chembonds, args=('toho_cex', 'r1ae'))
    p3 = multiprocessing.Process(target=process_chembonds, args=('toho_cex', 'r2ae'))

    jobs.append(p0)
    jobs.append(p1)
    jobs.append(p2)
    jobs.append(p3)

    for p in jobs:
        p.start()

if __name__ == "__main__":
    main()
