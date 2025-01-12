#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
This script takes in the vasp poscar input file  and converts it into the QE scf input format
'''


import glob
import os
import sys
import numpy as np
import math
import shutil
from parameters import    atom_info,Pseudodir,Molecule,ibrav,ecutwfc,ecutrho,conv_thr,diagonalization,mixing_beta,electron_maxstep,k_points_unitcell,van_der_Waals,vdw_corr                                                             




list_names = []
for filename in sorted(glob.glob('Raman-POSCAR.*')):
    name1 = filename.split('.vasp')
    name2= name1[0].split('Raman-POSCAR.')[1]
    with open(filename) as f:
        content = f.readlines()
        entry1= content[2].split()
        entry2= content[3].split()
        entry3= content[4].split()
        entry4=content[5].split()
        entry5=content[6].split()
        atyp_uniq=entry4
        ntyp=len(atyp_uniq)
        
        Atoms=(entry4[0]+' ')*int(entry5[0]) # Creating a list of all Atoms (not necessary unique) in the syste
        for jj in range(len(entry4)-1):
            Atoms=Atoms+(entry4[jj+1]+' ')*int(entry5[jj+1])
        Atoms=Atoms.split()
        nat=len(Atoms)
    with open('scf.'+name2+'.'+'in','w') as scf:
        scf.write('&CONTROL' + '\n')
        scf.write('  calculation  = "scf",' + '\n')
        scf.write('  restart_mode  = "from_scratch",' + '\n')
        scf.write('  verbosity  = "high",' + '\n')
        scf.write('  pseudo_dir   = "' + str(Pseudodir) + '",' + '\n')
        scf.write('  outdir  = "./",' + '\n')
        scf.write('  prefix       = "' + str(Molecule) + '",' + '\n')
        scf.write('  tstress  = .true.,' + '\n')
        scf.write('  tprnfor  = .true.,' + '\n')
        scf.write('  nstep  = 1000,' + '\n')
        scf.write('/' + '\n')
        scf.write('&SYSTEM' + '\n')
        scf.write('  ibrav = 0,' + '\n')     
        scf.write('  nat       = ' + str(nat) + ',' + '\n')
        scf.write('  ntyp       = ' + str(ntyp) + ',' + '\n')
        scf.write('  ecutwfc   = ' + str(ecutwfc) + ',' + '\n')
        scf.write('  ecutrho   = ' + str(ecutrho) + ',' + '\n')
        if van_der_Waals==True:
          scf.write('  vdw_corr   = ' + str(vdw_corr) + ',' + '\n')  
        scf.write('/' + '\n')
        scf.write('&ELECTRONS' + '\n')
        scf.write('  conv_thr    = ' + str(conv_thr) + ',' + '\n')
        scf.write('  diagonalization       = "' + str(diagonalization) + '",' + '\n')
        scf.write('  mixing_beta    = ' + str(mixing_beta) + ',' + '\n')
        scf.write('  electron_maxstep    = ' + str(electron_maxstep) + ',' + '\n')
        scf.write('/' + '\n')
        scf.write('ATOMIC_SPECIES' + '\n')
        potential_info = ' '
        for at in atyp_uniq:
          atomic_mass = atom_info[at][0]
          potential_info = atom_info[at][1]
          scf.write('  ' + str(at) + '  ' + str(atomic_mass) + '  ' + str(potential_info) + '\n')
        scf.write('K_POINTS (automatic)' + '\n')
        scf.write(str(k_points_unitcell) + '\n')        
        scf.write('CELL_PARAMETERS (angstrom)' + '\n')  # As phonopy uses AU (bohr) unit for length when interfaced with QE, the unit of  CELL_PARAMETERS is set to bohr
        scf.write(entry1[0]+' '+entry1[1]+' '+entry1[2]+ '\n')
        scf.write(entry2[0]+' '+entry2[1]+' '+entry2[2]+ '\n')
        scf.write(entry3[0]+' '+entry3[1]+' '+entry3[2]+ '\n')
    
        scf.write('ATOMIC_POSITIONS (crystal)' + '\n')
        kk=0
        for Z in Atoms:  # Getting the crystal atomic coordinates from the poscar and write them into scf.in 
            entry6=content[kk+8].split()
            scf.write(Z+' '+entry6[0]+' '+entry6[1]+' '+entry6[2]+ '\n')
            kk=kk+1
    
        scf.close()
        













