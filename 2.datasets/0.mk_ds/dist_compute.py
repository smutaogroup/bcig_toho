# Defs for extracting features.
# Zilin Song, 25 AUG 2021
# 

from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

def interatomic_dist(atom_i, atom_j, ):
    '''Compute the interatomic distance between 2 atoms.
    '''
    return	round(
        ((atom_i.position[0] - atom_j.position[0])**2 + \
         (atom_i.position[1] - atom_j.position[1])**2 + \
         (atom_i.position[2] - atom_j.position[2])**2   \
        ) ** 0.5, 8)

def dist_rx(mda_universe, repid, pathname, ):
    '''Compute and return all reaction coordinates.
    These selections have to be hard-coded, one way or the other.

    dx:  sel_atom_i    sel_atom_j
    d0:  Ser70-HG1  =  Ser70-OG
    d1:  Ser70-HG1  = Wat433-OH2
    d2:  Wat433-H2  = Wat433-OH2
    d3:  Wat433-H2  = Glu165-OE1
    d3a: Wat433-H2  = Glu165-OE2
    d4:  Lys73-HZ2  =  Ser70-OG
    d5:  Lys73-HZ2  =  Lys73-NZ
    d6:  Lys73-HZ1  =  Lys73-NZ
    d7:  Lys73-HZ1  = Ser130-OG
    d8:  Ser130-HG1 = Ser130-OG
    d9:  Ser130-HG1 = Lig285-N
    d10: Ser70-OG   = Lig285-C  -> this entry is included in the chembonds selection;
    d11: Lig285-C   = Lig285-N  -> this entry is included in the chembonds selection.
    d12: Ser130-HG1 = Lig285-CO2-O
    '''
    _selbase = 'segid Q{0} and'.format(str(repid))
    atom_i_sel = [
        ' (resid  69 and resname  SER and name HG1)',	# d0
        ' (resid  69 and resname  SER and name HG1)',	# d1
        ' (resid 433 and resname TIP3 and name  H2)',	# d2
        ' (resid 433 and resname TIP3 and name  H2)',	# d3
        ' (resid 433 and resname TIP3 and name  H2)',	# d3a
        ' (resid  72 and resname  LYS and name HZ2)',	# d4
        ' (resid  72 and resname  LYS and name HZ2)',	# d5
        ' (resid  72 and resname  LYS and name HZ1)',	# d6
        ' (resid  72 and resname  LYS and name HZ1)',	# d7
        ' (resid 129 and resname  SER and name HG1)',	# d8
        ' (resid 129 and resname  SER and name HG1)',	# d9
        # missing SER130-HG1				# d12
    ]
    atom_j_sel = [
        ' (resid  69 and resname  SER and name  OG)',	# d0
        ' (resid 433 and resname TIP3 and name OH2)',	# d1
        ' (resid 433 and resname TIP3 and name OH2)',	# d2
        ' (resid 165 and resname  GLU and name OE1)',	# d3
        ' (resid 165 and resname  GLU and name OE2)',	# d3a
        ' (resid  69 and resname  SER and name  OG)',	# d4
        ' (resid  72 and resname  LYS and name  NZ)',	# d5
        ' (resid  72 and resname  LYS and name  NZ)',	# d6
        ' (resid 129 and resname  SER and name  OG)',	# d7
        ' (resid 129 and resname  SER and name  OG)',	# d8
        # missing Lig285-N				# d9
        # missing Lig285-CO2-O				# d12
    ]
    
    # Determine ligand name. 
    lig_sel= 'segid Q{0} and resid 285 and (resname {1} or resname {2})'.format(str(repid), 'AMP', 'CEX')
    lig = mda_universe.select_atoms(lig_sel)
    if lig.n_residues != 1 or (not lig.residues[0].resname in ['AMP', 'CEX', ]):
        print('LIG selected {0} residues or selected residue name {1}: check dist_rx().'.format(str(lig.n_residues), str(lig.residues[0])))
        exit()
    
    # Append missing ligand rx atoms.
    if lig.residues[0].resname == 'AMP':
        atom_i_sel.append(' (resid 129 and resname SER and name HG1)')
        atom_j_sel.append(' (resid 285 and resname {0} and name  N4)'.format(lig.residues[0].resname))
        atom_j_sel.append(' (resid 285 and resname {0} and name O12)'.format(lig.residues[0].resname))
    elif lig.residues[0].resname == 'CEX':
        atom_i_sel.append(' (resid 129 and resname SER and name HG1)')
        atom_j_sel.append(' (resid 285 and resname {0} and name  N5)'.format(lig.residues[0].resname))
        atom_j_sel.append(' (resid 285 and resname {0} and name O10B)'.format(lig.residues[0].resname))
    else:
        print('Wrong residue[0] name at resid 285: ' + str(lig.residue[0].resname))
        exit()

    _rxcsel = [0, 1, 2, 3, 4, 7, 8, 9, 10, 11, ] if pathname == 'r1ae' else [5, 6, 7, 8, 9, 10, 11, ]
    _distmat = []
    _distlbl = []
    for i in _rxcsel:
        group0 = mda_universe.select_atoms(_selbase + atom_i_sel[i])
        group1 = mda_universe.select_atoms(_selbase + atom_j_sel[i])

        if group0.n_atoms == 1 and group1.n_atoms == 1:
            atom_i = group0.atoms[0]
            atom_j = group1.atoms[0]
            d = interatomic_dist(atom_i, atom_j, )
            l =	'{0}.{1}.{2}:{3}.{4}.{5}'.format(
                atom_i.residue.resname, atom_i.residue.resid, atom_i.name,
                atom_j.residue.resname, atom_j.residue.resid, atom_j.name,
            )
            _distmat.append(d)
            _distlbl.append(l)
        else: 
            raise ValueError('Pairwise selection retrived more than 1 atoms in each group: Check dist_compute.dist_rx()\n')
            
    return _distmat, _distlbl

def qmhvy_selection(mda_universe, repid, ):
    '''Atom selection for heavy atoms in the QM region.
    return (selected_AtomGroup, ligand_name, )
    '''
    qmhvy_sel = 'segid Q{0} and ('  \
                '    (resid  69 and (name  CB or name  OG))'  \
                ' or (resid  72 and (name  CB or name  CG or name  CD or name CE or name NZ))'  \
                ' or (resid 129 and (name  CB or name  OG))'  \
                ' or (resid 165 and (name  CB or name  CG or name  CD or name OE1 or name OE2))'  \
                ' or (resid 169 and (name  CB or name  CG or name OD1 or name ND2))'  \
                ' or (resid 233 and (name  CB or name  CG or name  CD or name  CE or name NZ))'  \
                ' or (resid 234 and (name  CB or name OG1 or name CG2))'  \
                ' or (resid 236 and (name  CB or name  OG))'  \
                ' or (resid 433 and (name OH2 ))'  \
                ')'.format(str(repid))
    qmhvy = mda_universe.select_atoms(qmhvy_sel)

    # Should be 10 residues with 53 atoms.
    if qmhvy.n_residues != 9 or qmhvy.n_atoms != 29:
        print([i for i in qmhvy.atoms])
        raise ValueError(
            'QMHVY selected {0}/{1} residues/atoms, should be 9/29: check qmhvy_selection().'.format(
                qmhvy.n_residues, qmhvy.n_atoms,)
        )

    # select heavy atoms on ligands, excluding the fuzed ring structure.
    ## select ligand for determine name: AMP or CEX.
    lig_selbase = 'segid Q{0} and resid 285 and (resname {1} or resname {2})'.format(
        str(repid), 'AMP', 'CEX')
    lig = mda_universe.select_atoms(lig_selbase)

    ## sanity check for selected ligands
    if lig.n_residues != 1 or (not lig.residues[0].resname in ['AMP', 'CEX', ]):
        raise ValueError(
            'LIG selected {0} residues or wrong residue name {1}: check qmhvy_selection().'.format(
                str(lig.n_residues), str(lig.residues[0]))
        )

    # ordered ligand selection.
    lig_ordered = None

    # NOTE that atoms on ligands are aligned for selection by order, 
    #           otherwise selected atoms are ranked by incremental indices.
    # Also NOTE that 3 carbon atoms on the fuzed rings are discarded due to alignment reasons.
    if lig.residues[0].resname == 'AMP':
        lig_ordered	= mda_universe.select_atoms(lig_selbase + ' and name  N4')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name  C5')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name  C6')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name  C7')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name  O8')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name N14')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name C15')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name O16')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name C17')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name N34')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name C18')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name C19')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name C20')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name C21')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name C22')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name C23')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name  S1')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name  C3')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name C11')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name O12')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name O13')
    elif lig.residues[0].resname == 'CEX':
        lig_ordered	= mda_universe.select_atoms(lig_selbase + ' and name  N5')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name  C6')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name  C7')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name  C8')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name  O8')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name N11')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name C12')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name O12')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name C13')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name N20')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name C14')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name C15')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name C16')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name C17')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name C18')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name C19')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name  S1')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name  C4')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name C10')  \
                    + mda_universe.select_atoms(lig_selbase + ' and name O10B') \
                    + mda_universe.select_atoms(lig_selbase + ' and name O10A')
    else:
        raise ValueError(
            'LIG selected {0} residues or wrong residue name {1}: check qmhvy_selection().'.format(
                str(lig.n_residues), str(lig.residues[0]))
        )

    qmhvy += lig_ordered

    if qmhvy.n_residues == 10 and qmhvy.n_atoms == 50: 
        return qmhvy, lig.residues[0].resname
    else: 
        print([i for i in qmhvy.atoms])
        raise ValueError(
            'QMHVY selected {0}/{1} residues/atoms, should be 10/51: check qmhvy_selection().'.format(
                qmhvy.n_residues, qmhvy.n_atoms, )
            )

def dist_hvypw(mda_universe, repid, ):
    '''Compute the pair wise distances b/w heavy atoms.
    '''
    qmhvy_atoms, ligname = qmhvy_selection(mda_universe, repid, )
    _distmat = []
    _distlbl = []
    
    for i in range(qmhvy_atoms.n_atoms-1):
        for j in range(i+1, qmhvy_atoms.n_atoms):
            # all possible (selected atom pairs)
            atom_i = qmhvy_atoms[i]
            atom_j = qmhvy_atoms[j]
            d = interatomic_dist(atom_i, atom_j)
            l ='{0}.{1}.{2}:{3}.{4}.{5}'.format(
                atom_i.residue.resname, atom_i.residue.resid, atom_i.name,
                atom_j.residue.resname, atom_j.residue.resid, atom_j.name,
            )
            _distmat.append(d)
            _distlbl.append(l) 
    return _distmat, _distlbl

def dist_bonds(mda_universe, repid, ):
    '''Compute and return all chemical bonding distances.
    NOTE: the S70-OG -- LIG-C acyl bond is appended manually.
    NOTE: C9/C10 in AMP and C3/C9 in CEX are _always_ unselected.
    '''
    qmhvy_atoms, ligname = qmhvy_selection(mda_universe, repid, )
    all_bonds = qmhvy_atoms.bonds

    _distmat = []
    _distlbl = []

    # Determine atoms. Note that this aligns the order of the bond lengths,
    # according to the order in qmhvy_atoms selection.
    # NOTE this is an i/o bottle neck, 
    #      but I don't care since I am using this only once.
    for i in range(qmhvy_atoms.n_atoms-1):
        for j in range(i+1, qmhvy_atoms.n_atoms):
            # all possible (selected) atom pairs
            atom_i = qmhvy_atoms[i]
            atom_j = qmhvy_atoms[j]

            for bond in all_bonds:                      # all bond list between selected atoms.
                if atom_i in bond and atom_j in bond:   # detect if bonded.
                    d = interatomic_dist(atom_i, atom_j)
                    l = '{0}.{1}.{2}:{3}.{4}.{5}'.format(
                        atom_i.residue.resname, atom_i.residue.resid, atom_i.name,
                        atom_j.residue.resname, atom_j.residue.resid, atom_j.name,
                    )
                    _distmat.append(d)
                    _distlbl.append(l)
                    
    # bonds from reaction.
    atom_ser70og = mda_universe.select_atoms(
        'segid Q{0} and (resid  69 and name  OG)'.format(str(repid))
    )
    atom_ligc = mda_universe.select_atoms(
        'segid Q{0} and resid 285 and name C7'.format(str(repid))
    ) if ligname == 'AMP' else mda_universe.select_atoms(
        'segid Q{0} and resid 285 and name C8'.format(str(repid))
    )
    
    ## sanity check for selected bondings
    if atom_ser70og.n_atoms != 1 or (not atom_ser70og.residues[0].resname in ['SER', ]):
        raise ValueError(
            'LIG selected {0} residues or wrong residue name {1}: check dist_bonds().'.format(
                str(atom_ser70og.n_atoms), str(atom_ser70og.residues[0]))
        )
    elif atom_ligc.n_atoms != 1 or (not atom_ligc.residues[0].resname in ['AMP', 'CEX', ]):
        raise ValueError(
            'LIG selected {0} residues or wrong residue name {1}: check dist_bonds().'.format(
                str(atom_ligc.n_atoms), str(atom_ligc.residues[0]))
        )
    
    _distmat.append(interatomic_dist(atom_ser70og.atoms[0], atom_ligc.atoms[0], ))
    _distlbl.append('{0}.{1}.{2}:{3}.{4}.{5}'.format(
                        atom_ser70og.atoms[0].residue.resname, 
                        atom_ser70og.atoms[0].residue.resid, 
                        atom_ser70og.atoms[0].name,
                        atom_ligc.atoms[0].residue.resname, 
                        atom_ligc.atoms[0].residue.resid, 
                        atom_ligc.atoms[0].name,
                        )
                    )

    return _distmat, _distlbl

def hbond_detect(mda_universe, repid, ):
    '''Detects all hbonds on that replica, defined as Donor-Acceptor distances.
    save the labels of hbonds as an inhomogeneous numpy.array()
    '''
    # atoms to be selected in residues.
    qmres_sel = '(segid Q{0} and ('  \
                '    (resid  69 and (name  CB or name  OG))'  \
                ' or (resid  72 and (name  CB or name  CG or name  CD or name CE or name NZ))'  \
                ' or (resid 129 and (name  CB or name  OG))'  \
                ' or (resid 165 and (name  CB or name  CG or name  CD or name OE1 or name OE2))'  \
                ' or (resid 169 and (name  CB or name  CG or name OD1 or name ND2))'  \
                ' or (resid 233 and (name  CB or name  CG or name  CD or name  CE or name NZ))'  \
                ' or (resid 234 and (name  CB or name OG1 or name CG2))'  \
                ' or (resid 236 and (name  CB or name  OG))'  \
                ' or (resid 433 and (name OH2 ))'  \
                ') )'.format(str(repid))

    # see if intended selection
    qmres = mda_universe.select_atoms(qmres_sel)
    if qmres.n_residues != 9 or qmres.n_atoms != 29:
        print([i for i in qmres.atoms])
        raise ValueError(
            'qmres selected {0}/{1} residues/atoms, should be 9/29: check hbond_detect().'.format(
                qmres.n_residues, qmres.n_atoms,)
        )

    # atoms to be selected in ligands.
    # NOTE that the selection is not aligned (due to the limitations in HBA).
    qmlig_sel = ' or (segid Q{0} and resid 285 and (resname AMP or resname CEX) and (not name H*) )'.format(
        str(repid), 
    )
    
    # merge selection string
    allhvy_sel = qmres_sel + qmlig_sel

    # test if intended selection.
    ## NOTE difference carbons are retained: amp-C9/C10, cex-C4-C9
    allhvy = mda_universe.select_atoms(allhvy_sel)    
    if allhvy.n_residues != 10 or allhvy.n_atoms != 53: 
        print([i for i in allhvy.atoms])
        raise ValueError(
            'allhvy selected {0}/{1} residues/atoms, should be 10/51: check hbond_detect().'.format(
                allhvy.n_residues, allhvy.n_atoms, )
            )

    # hydrogens to be selected.
    allh_sel = 'segid Q{0} and name H* and ('  \
                '    (resid  69 and resname SER)' \
                ' or (resid  72 and resname LYS)' \
                ' or (resid 129 and resname SER)' \
                ' or (resid 165 and resname GLU)' \
                ' or (resid 169 and resname ASN)' \
                ' or (resid 233 and resname LYS)' \
                ' or (resid 234 and resname THR)' \
                ' or (resid 236 and resname SER)' \
                ' or (resid 433 and resname TIP3)' \
                ')'.format(str(repid))

    # Hydrogen bonding analysis.
    hba = HBA(universe=mda_universe, 
                 donors_sel=allhvy_sel, 
                 hydrogens_sel=allh_sel, 
                 acceptors_sel=allhvy_sel,  
                )
    hba.run()
    results = hba.hbonds

    # Labels of Hbonds.
    ## NOTE format  Donor:Acceptor
    labels = []
    for hb in results:
        donor_atom = mda_universe.atoms[int(hb[1])]
        accep_atom = mda_universe.atoms[int(hb[3])]

        da_l = '{0}.{1}.{2}:{3}.{4}.{5}'.format(
                donor_atom.residue.resname, donor_atom.residue.resid, donor_atom.name,
                accep_atom.residue.resname, accep_atom.residue.resid, accep_atom.name,
            )
        ad_l = '{0}.{1}.{2}:{3}.{4}.{5}'.format(
                accep_atom.residue.resname, accep_atom.residue.resid, accep_atom.name,
                donor_atom.residue.resname, donor_atom.residue.resid, donor_atom.name,
            )

        if not (da_l in labels or ad_l in labels):      # filter repeated entries.
            labels.append(da_l)
        
    return labels

def unify_hbond_labels(amp_labels, cex_labels):
    '''Here return unified labels for replica-wise feature extraction. 
    NOTE ligand atoms always at l1.
    This is mostly for feature alignment, e.g.:

    uni_amp_labels___uni_cex_labels

    SER.69.OG:TIP3.433.OH2___SER.69.OG:TIP3.433.OH2
    LYS.72.NZ:SER.129.OG___LYS.72.NZ:SER.129.OG
    LYS.72.NZ:SER.69.OG___LYS.72.NZ:SER.69.OG
    SER.129.OG:AMP.285.O12___SER.129.OG:CEX.285.O10B
    ASN.169.ND2:GLU.165.OE1___ASN.169.ND2:GLU.165.OE1
    LYS.233.NZ:AMP.285.O12___LYS.233.NZ:CEX.285.O10B
    THR.234.OG1:AMP.285.O13___THR.234.OG1:CEX.285.O10A
    TIP3.433.OH2:ASN.169.OD1___TIP3.433.OH2:ASN.169.OD1
    TIP3.433.OH2:GLU.165.OE1___TIP3.433.OH2:GLU.165.OE1
    SER.129.OG:AMP.285.N4___SER.129.OG:CEX.285.N5
    THR.234.OG1:AMP.285.O12___THR.234.OG1:CEX.285.O10B
    ASN.169.ND2:AMP.285.N34___ASN.169.ND2:CEX.285.N20
    LYS.233.NZ:SER.129.OG___LYS.233.NZ:SER.129.OG
    LYS.72.NZ:GLU.165.OE2___LYS.72.NZ:GLU.165.OE2
    LYS.72.CE:GLU.165.OE2___LYS.72.CE:GLU.165.OE2
    SER.236.OG:AMP.285.O13___SER.236.OG:CEX.285.O10A
    TIP3.433.OH2:GLU.165.OE2___TIP3.433.OH2:GLU.165.OE2
    LYS.233.NZ:THR.234.OG1___LYS.233.NZ:THR.234.OG1
    SER.236.OG:AMP.285.N34___SER.236.OG:CEX.285.N20
    TIP3.433.OH2:GLU.165.CD___TIP3.433.OH2:GLU.165.CD
    GLU.165.OE1:AMP.285.O16___GLU.165.OE1:CEX.285.O12
    SER.129.OG:AMP.285.O13___SER.129.OG:CEX.285.O10A
    ASN.169.ND2:GLU.165.OE2___ASN.169.ND2:GLU.165.OE2
    TIP3.433.OH2:AMP.285.O16___TIP3.433.OH2:CEX.285.O12
    ASN.169.ND2:TIP3.433.OH2___ASN.169.ND2:TIP3.433.OH2
    LYS.233.NZ:AMP.285.O13___LYS.233.NZ:CEX.285.O10A
    SER.236.OG:AMP.285.O12___SER.236.OG:CEX.285.O10B
    '''
    uni_amp_labels = []
    uni_cex_labels = []

    # corresponding amp atoms to cex atoms by order in both lists. 
    # These lists has to be hard-coded.
    at_amp = ['N4', 'C5', 'C6', 'C7', 'O8', 'N14', 'C15', 'O16', 'C17', 'N34', 'C18', 'C19', 'C20', 'C21', 'C22', 'C23', 'S1', 'C2', 'C3', 'C11',  'O12',  'O13', ]
    at_cex = ['N5', 'C6', 'C7', 'C8', 'O8', 'N11', 'C12', 'O12', 'C13', 'N20', 'C14', 'C15', 'C16', 'C17', 'C18', 'C19', 'S1', 'C2', 'C4', 'C10', 'O10B', 'O10A', ]

    # First loop through amp_labels
    for al in amp_labels:
        l = al.split(':')
        l0, l1 = l[0], l[1]

        if l1.split('.')[0] != 'AMP':   # irrelavent to ligand, add directly.
            uni_amp_labels.append(al)
            uni_cex_labels.append(al)
        
        else:                           # ligand relavent
            # get at_amp indices
            at_amp_i = [ i for i in range(len(at_amp)) if at_amp[i] == l1.split('.')[2] ]

            cl = '{0}:{1}'.format(
                l0, 
                'CEX.285.{0}'.format(at_cex[at_amp_i[0]])
            )
            uni_amp_labels.append(al)
            uni_cex_labels.append(cl)

    # Then loop through cex_labels
    for cl in cex_labels:
        l = cl.split(':')
        l0, l1 = l[0], l[1]

        if l1.split('.')[0] != 'CEX':   # irrelavent to ligand, check if already selected.
            rev_cl = '{0}:{1}'.format(l1, l0)

            if not (cl in uni_cex_labels or rev_cl in uni_cex_labels):  # not selected, add.
                uni_amp_labels.append(cl)
                uni_cex_labels.append(cl)
        
        else:                           # relavent to ligand, check if already selected.
            if not cl in uni_cex_labels:
                at_cex_i = [ i for i in range(len(at_cex)) if at_cex[i] == l1.split('.')[2] ]
                al = '{0}:{1}'.format(
                    l0, 
                    'AMP.285.{0}'.format(at_amp[at_cex_i[0]])
                )
                
                uni_amp_labels.append(al)
                uni_cex_labels.append(cl)

    return uni_amp_labels, uni_cex_labels

def dist_hbonds(mda_universe, repid, da_labels, ):
    '''Compute and return all inter-donor/acceptor pairwise distances specified in da_labels
    '''
    _distmat = []
    _distlbl = []

    for pairlabel in da_labels:
        atomlabels = pairlabel.split(':')
        l0, l1 = atomlabels[0], atomlabels[1]

        group0 = mda_universe.select_atoms('segid Q{0} and resid {1} and resname {2} and name {3}'.format(
            str(repid), l0.split('.')[1], l0.split('.')[0], l0.split('.')[2], )
        )

        group1 = mda_universe.select_atoms('segid Q{0} and resid {1} and resname {2} and name {3}'.format(
            str(repid), l1.split('.')[1], l1.split('.')[0], l1.split('.')[2], )
        )

        if group0.n_atoms == 1 and group1.n_atoms == 1:

            atom_i = group0.atoms[0]
            atom_j = group1.atoms[0]
            d = interatomic_dist(atom_i, atom_j, )

            _distmat.append(d)
            _distlbl.append(pairlabel)

        else: 
            raise ValueError('Pairwise selection retrived more than 1 atoms in each group: Check dist_compute.dist_hbonds()\n')

    return _distmat, _distlbl
