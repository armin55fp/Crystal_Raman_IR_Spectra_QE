# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 10:49:27 2022

@author: armin.taheri
This code creates the scf.in input files,the header.in file, and supercells for calculation of the IR/Raman spectra using phononpy+QE. Also, it creates three folders called part1, part2, and part3 and divided all the supercell files between them and submits the DFT calculations to TerraHPC. It also creates a BORN_cal folder and submits the dielectric and Born charge calculation using pw,x and ph.x modules in Quantum Espresso. All DFT setting will be read from parameter.py file. The only external file that this code needs is a file called  "relaxed_atomic_position.txt" which contains the relaxed atomic position (x,y,z) of the ...

Before running the code you need to do: >> conda activate phonopy
"""

import glob
import os
import sys
import numpy as np
import math
import shutil
from parameters import Pseudodir,atom_info,Molecule,ibrav,nat,ntyp,ecutwfc,ecutrho,conv_thr,diagonalization,mixing_beta,electron_maxstep,atyp_uniq,k_points_unitcell,k_points_supercell,tr2_ph,N_proc,nmix_ph,niter_ph,lattice1,lattice2,lattice3,Na,Nb,Nc,van_der_Waals,vdw_corr,Type_of_atomic_position



cwd = os.getcwd() # Saving the current working directory address in cwd 


# Creating the "scf.in" file based on the settings provided before.

with open('scf.in', 'w') as scf:  
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
    scf.write('CELL_PARAMETERS (angstrom)' + '\n')
    scf.write(lattice1[0]+' '+lattice1[1]+' '+lattice1[2]+ '\n')
    scf.write(lattice2[0]+' '+lattice2[1]+' '+lattice2[2]+ '\n')
    scf.write(lattice3[0]+' '+lattice3[1]+' '+lattice3[2]+ '\n')
    scf.write('ATOMIC_POSITIONS (crystal)' + '\n')
    scf.write('\n')
    all_atoms = []
    coord_x = []
    coord_y = []
    coord_z = []
    with open("relaxed_atomic_position.txt") as f:

        for line in f:
            linestrip = line.strip()
            entry = linestrip.split()
            if Type_of_atomic_position==1: # relaxed_atomic_position.txt with ()
              all_atoms.append(entry[1])
              coord_x.append(entry[6])
              coord_y.append(entry[7])
              coord_z.append(entry[8])
            elif Type_of_atomic_position==2:  # relaxed_atomic_position.txt without ()
              all_atoms.append(entry[0])
              coord_x.append(entry[1])
              coord_y.append(entry[2])
              coord_z.append(entry[3])               
    jj=0       
    for position in coord_x:
        coord_x[jj]=float(coord_x[jj])
        coord_y[jj]=float(coord_y[jj])
        coord_z[jj]=float(coord_z[jj])
        jj=jj+1
        
                    
    i=0    
    for a in all_atoms:
         
         scf.write(all_atoms[i]+ '\t' +str(coord_x[i]) + '\t' + str(coord_y[i]) + '\t' + str(coord_z[i]) + '\n')
         i+=1

    scf.close()



# Creating and submitting the BORN effective charge calculations

os.system('mkdir BORN_cal')
os.system('scp scf.in BORN_cal')
with open('ph.in', 'w') as ph:  # Creating the "ph.in" file based on the settings provided before for the Born effective charge calcualtion
    ph.write('           ' +'phonons in '+Molecule+ '\n')
    ph.write('  &inputph' + '\n')
    ph.write('  outdir  = "./",' + '\n')
    ph.write('  prefix       = "' + str(Molecule) + '",' + '\n')
    ph.write('  tr2_ph    = ' + str(tr2_ph) + ',' + '\n')
    ph.write('  epsil  = .true.,' + '\n')
    ph.write('  fildyn = "dyn",' + '\n')
    ph.write('  trans  = .false.,' + '\n')
    ph.write('  nmix_ph    = ' + str(nmix_ph) + ',' + '\n')
    ph.write('  niter_ph    = ' + str(niter_ph) + ',' + '\n')
    ph.write('/' + '\n')
    ph.write('0 0 0'+'\n')
ph.close()
os.system('scp ph.in BORN_cal')
with open('scfph_sub', 'w') as scfph:# Creating the submission batch in TerraHPC for solvent calculations on 16 cores
        scfph.write('#!/bin/bash' + '\n')
        scfph.write('#SBATCH -N 1' + '\n')
        scfph.write('#SBATCH --tasks-per-node='+str(N_proc)+'\n')
        scfph.write('#SBATCH --mem-per-cpu=7G' + '\n')
        scfph.write('#SBATCH --job-name='+ Molecule+'_BORN' + '\n')
        scfph.write('\n')
        scfph.write('\n')
        scfph.write('module purge'+'\n')
        scfph.write('module load ohpc'+'\n')
        scfph.write('module load quantum-espresso-7.0-gcc-8.3.0-6pkcl6c'+'\n')
        scfph.write('\n')
        scfph.write('\n')
        scfph.write('\n')        
        scfph.write('mpirun'+'  '+'pw.x'+'  '+ '< scf.in > scf.out'+'\n')
        scfph.write('mpirun'+'  '+'ph.x'+'  '+ '< ph.in > ph.out'+'\n')
        scfph.write('\n')
scfph.close()
os.system('scp scfph_sub BORN_cal')
os.system('rm scfph_sub')
os.chdir('BORN_cal')
#os.system('sbatch scfph_sub')
os.chdir(cwd)
print('scf and ph calculations for getting the born effective charges have been successfully submitted')


## Creating the supercell files using phonopy


os.system("phonopy --qe -d --dim="+'"'+str(Na)+' '+str(Nb)+' '+str(Nc)+'"'+ " -c scf.in -v")


## Creating the header.in file

Nsup=nat*Na*Nb*Nc # calculation of the number of atoms in supercell
with open('header.in', 'w') as header:
    header.write('&CONTROL' + '\n')
    header.write('  calculation  = "scf",' + '\n')
    header.write('  pseudo_dir   = "' + str(Pseudodir) + '",' + '\n')
    header.write('  outdir  = "./",' + '\n')
    header.write('  tstress  = .true.,' + '\n')
    header.write('  tprnfor  = .true.,' + '\n')
    header.write('  disk_io  = ".none.",' + '\n')
    header.write('/' + '\n')
    header.write('&SYSTEM' + '\n')
    header.write('  ibrav = 0,' + '\n') 
    header.write('  nat       = ' +  str(Nsup) + ',' + '\n')
    header.write('  ntyp       = ' + str(ntyp) + ',' + '\n')
    header.write('  ecutwfc   = ' +  str(ecutwfc) + ',' + '\n')
    header.write('  ecutrho   = ' +  str(ecutrho) + ',' + '\n')
    if van_der_Waals==True:
      scf.write('  vdw_corr   = ' + str(vdw_corr) + ',' + '\n')  
    header.write('/' + '\n')
    header.write('&ELECTRONS' + '\n')
    header.write('  conv_thr    = ' + str(conv_thr) + ',' + '\n')
    header.write('  diagonalization       = "' + str(diagonalization) + '",' + '\n')
    header.write('  mixing_beta    = ' + str(mixing_beta) + ',' + '\n')
    header.write('  electron_maxstep    = ' + str(electron_maxstep) + ',' + '\n')
    header.write('/' + '\n')
    header.write('K_POINTS (automatic)' + '\n')
    header.write(str(k_points_supercell) + '\n')
    
    header.close()



## Creating three folders part1, part2, and part3 and dividing the supercell files between them.

os.system('ls')

N1=int(input("Please enter the number of displaced supercells created by phonopy: "))  # Number of the created supercells by phonopy. double check if this number is actually the same as the number of the created supercells.

command="for i in {001.."+str(N1)+"};do cat header.in supercell-$i.in >| "+Molecule+"-$i.in; done"

os.system(command)    

os.system("mkdir part1 part2 part3")

N2=math.floor(N1/3)

command2="for N in {001.."+str(N2)+"}; do scp "+Molecule+"-$N.in part1; done"        

os.system(command2)

N3=N2+1

command3="for N in {"+str(N3).zfill(3)+".."+str(N2*2)+"}; do scp "+Molecule+"-$N.in part2; done"

os.system(command3)

N4=2*N2+1

command4="for N in {"+str(N4).zfill(3)+".."+str(N1)+"}; do scp "+Molecule+"-$N.in part3; done"   

os.system(command4)

#------------------------------------------------------

# Creating the submission batch script for part1


with open('scf1_sub', 'w') as scf1:
        scf1.write('#!/bin/bash' + '\n')
        scf1.write('#SBATCH -N 1' + '\n')
        scf1.write('#SBATCH --tasks-per-node='+str(N_proc)+'\n')
        scf1.write('#SBATCH --mem-per-cpu=7G' + '\n')
        scf1.write('#SBATCH --job-name='+ Molecule+'p1' + '\n')
        scf1.write('\n')
        scf1.write('\n')
        scf1.write('module purge'+'\n')
        scf1.write('module load ohpc'+'\n')
        scf1.write('module load quantum-espresso-7.0-gcc-8.3.0-6pkcl6c'+'\n')
        scf1.write('\n')
        scf1.write('\n')
        scf1.write('\n')
        scf1.write('for N in {'+'001'+'..'+str(N2).zfill(3)+'}'+'\n')
        scf1.write('do'+'\n')
        scf1.write('mpirun'+'  '+ 'pw.x'+'  '+ '< '+Molecule+'-$N.in'+''+'>'+''+Molecule+'-$N.out'  + '\n')
        scf1.write('done'+'\n')
        scf1.write('\n')
scf1.close()

#---------------------------------------------------------

# Creating the submission batch script for part2

with open('scf2_sub', 'w') as scf2:
        scf2.write('#!/bin/bash' + '\n')
        scf2.write('#SBATCH -N 1' + '\n')
        scf2.write('#SBATCH --tasks-per-node='+str(N_proc)+'\n')
        scf2.write('#SBATCH --mem-per-cpu=7G' + '\n')
        scf2.write('#SBATCH --job-name='+ Molecule+'p2' + '\n')
        scf2.write('\n')
        scf2.write('\n')
        scf2.write('module purge'+'\n')
        scf2.write('module load ohpc'+'\n')
        scf2.write('module load quantum-espresso-7.0-gcc-8.3.0-6pkcl6c'+'\n')
        scf2.write('\n')
        scf2.write('\n')
        scf2.write('\n')
        scf2.write('for N in {'+str(N3).zfill(3)+'..'+str(N2*2).zfill(3)+'}'+'\n')
        scf2.write('do'+'\n')
        scf2.write('mpirun'+'  '+ 'pw.x'+'  '+ '< '+Molecule+'-$N.in'+''+'>'+''+Molecule+'-$N.out'  + '\n')
        scf2.write('done'+'\n')
        scf2.write('\n')
scf2.close()

#----------------------------------------------------------

# Creating the submission batch script for part3

with open('scf3_sub', 'w') as scf3:
        scf3.write('#!/bin/bash' + '\n')
        scf3.write('#SBATCH -N 1' + '\n')
        scf3.write('#SBATCH --tasks-per-node='+str(N_proc)+'\n')
        scf3.write('#SBATCH --mem-per-cpu=7G' + '\n')
        scf3.write('#SBATCH --job-name='+ Molecule+'p3' + '\n')
        scf3.write('\n')
        scf3.write('\n')
        scf3.write('module purge'+'\n')
        scf3.write('module load ohpc'+'\n')
        scf3.write('module load quantum-espresso-7.0-gcc-8.3.0-6pkcl6c'+'\n')
        scf3.write('\n')
        scf3.write('\n')
        scf3.write('\n')
        scf3.write('for N in {'+str(N4).zfill(3)+'..'+str(N1).zfill(3)+'}'+'\n')
        scf3.write('do'+'\n')
        scf3.write('mpirun'+'  '+ 'pw.x'+'  '+ '< '+Molecule+'-$N.in'+''+'>'+''+Molecule+'-$N.out'  + '\n')
        scf3.write('done'+'\n')
        scf3.write('\n')
        scf3.write('sleep 2400'+'\n')
        scf3.write('cd ..'+'\n')
        #scf3.write('source /mnt/cephfs/home/armin5/ARMIN_env/bin/activate'+ '\n')
        #scf3.write('export PYTHONPATH=${PYTHONPATH}:/mnt/cephfs/home/armin5/ARMIN_env/bin/Phonopy-Spectroscopy-master'+ '\n')
        #scf3.write('export PATH=${PATH}:/mnt/cephfs/home/armin5/ARMIN_env/bin/Phonopy-Spectroscopy-master/Scripts'+ '\n')
        #scf3.write('python3 IR_Raman_Step2.py'+ '\n')
        #scf3.write('\n')
        
scf3.close()



#Submiiting all three supercell simulations to terraHPC

os.system('scp scf1_sub part1')
os.system('scp scf2_sub part2')
os.system('scp scf3_sub part3')
os.system('rm scf1_sub scf2_sub scf3_sub ph.in')
os.chdir('part1')
#os.system('sbatch scf1_sub')
os.chdir(cwd)
os.chdir('part2')
#os.system('sbatch scf2_sub')
os.chdir(cwd)
os.chdir('part3')
#os.system('sbatch scf3_sub')
os.chdir(cwd)
#print('All the supercell calculations have been successfully submitted')
