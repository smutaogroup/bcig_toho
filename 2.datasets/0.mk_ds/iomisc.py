# IO defs for loading pathway data.
# Zilin Song, 20 AUG 2021
# 

import MDAnalysis as mda

def basedir():
    '''An ugly way to put a global var.
    '''
    return '/users/zilins/scratch/2.proj_toho2lig_acy/1.sampling'

def enerdirs(sysname, pathname, theoryname, ):
    '''Build the directories to the energies.
    
    sysname: 	toho_amp, toho_cex;
    pathname:	r1ae, r2ae;
    theoryname: dftb3, b3lyp.

    One pathname will return two directories, e.g.:
    /users/zilins/scratch/2.proj_toho2lig_acy/1.sampling/ 0 .toho_amp. r1/9.paths.sp/         0.dftb3   /path.ene
    /users/zilins/scratch/2.proj_toho2lig_acy/1.sampling/ 2 .toho_amp. ae/9a.r1.paths.sp/     0.dftb3   /path.ene

    Will return 2 str: fwdir & bwdir
    '''
    # sysid:    number before sysname;
    if sysname == 'toho_amp':
        _bw_sysid = 2
    elif sysname == 'toho_cex':
        _bw_sysid = 5
    else:
        print('Invalid sysname: ' + sysname)

    if pathname == 'r1ae':
        _fw_sysid = _bw_sysid - 2
        _fw_pathname = 'r1/9.paths.sp'
        _bw_pathname = 'ae/9a.r1.paths.sp'
    elif pathname == 'r2ae':
        _fw_sysid = _bw_sysid - 1
        _fw_pathname = 'r2/9.paths.sp'
        _bw_pathname = 'ae/9b.r2.paths.sp'
    else:
        print('Invalid pathname: ' + pathname)
    
    if theoryname == 'dftb3':
        _theoryname = '0.dftb3'
    elif theoryname == 'b3lyp':
        _theoryname = '1.b3lyp'
    else:
        print('Invalid theoryname: ' + theoryname)
    
    _bdir = basedir()
    _fw_edir = '{0}/{1}.{2}.{3}/{4}/pathraw.ene'.format(_bdir, str(_fw_sysid), sysname, _fw_pathname, _theoryname )
    _bw_edir = '{0}/{1}.{2}.{3}/{4}/pathraw.ene'.format(_bdir, str(_bw_sysid), sysname, _bw_pathname, _theoryname )

    _fw_barrierdir = '{0}/{1}.{2}.{3}/{4}/pathbarrier.ene'.format(_bdir, str(_fw_sysid), sysname, _fw_pathname, _theoryname )
    _bw_barrierdir = '{0}/{1}.{2}.{3}/{4}/pathbarrier.ene'.format(_bdir, str(_bw_sysid), sysname, _bw_pathname, _theoryname )

    return _fw_edir, _bw_edir, _fw_barrierdir, _bw_barrierdir

def load_path(sysname, pathname, whichdirection, pathid, ):
    '''Build the directories to the conformations; 
    Load and return a MDAnalysis.Universe instance.
    
    sysname: 		toho_amp, toho_cex;
    pathname:		r1ae, r2ae;
    whichdirection:	fw, bw;
    pathid:     	a number from 1-100, 

    Conformation directories:
    /users/zilins/scratch/2.proj_toho2lig_acy/1.sampling/3.toho_cex. r1/8.    dftb.paths  /path_opt/toho_cex.path_f52.cor
    /users/zilins/scratch/2.proj_toho2lig_acy/1.sampling/4.toho_cex .r2/8.    dftb.paths  /path_opt/toho_cex.path_f52.psf
    /users/zilins/scratch/2.proj_toho2lig_acy/1.sampling/5.toho_cex. ae/8b.r2.dftb.paths  /path_opt/toho_cex.path_f65.cor

    Will return a MDAnalysis.Universe object of the path.
    '''

    # sysid:    number before sysname;
    if sysname == 'toho_amp':
        _sysid = 0
    elif sysname == 'toho_cex':
        _sysid = 3
    else:
        print('Invalid sysname: ' + sysname)

    if   pathname == 'r1ae' and whichdirection == 'fw':
        _sysid += 0
        _pathname = 'r1/8.dftb.paths'
    elif pathname == 'r2ae' and whichdirection == 'fw':
        _sysid += 1
        _pathname = 'r2/8.dftb.paths'
    elif pathname == 'r1ae' and whichdirection == 'bw':
        _sysid += 2
        _pathname = 'ae/8a.r1.dftb.paths'
    elif pathname == 'r2ae' and whichdirection == 'bw':
        _sysid += 2
        _pathname = 'ae/8b.r2.dftb.paths'
    else:
        print('Invalid pathname: ' + sysname + '\nInvalid whichdirection ' + whichdirection)
    
    _bdir = basedir()

    _psfdir = '{0}/{1}.{2}.{3}/path_opt/{2}.path_f{4}.psf'.format(_bdir, _sysid, sysname, _pathname, str(pathid))
    _cordir = '{0}/{1}.{2}.{3}/path_opt/{2}.path_f{4}.cor'.format(_bdir, _sysid, sysname, _pathname, str(pathid))

    return mda.Universe(_psfdir, _cordir, topology_format='PSF', format='CRD', )
