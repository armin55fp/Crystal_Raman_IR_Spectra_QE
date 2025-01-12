# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 10:49:27 2022

@author: armin.taheri

This code gathers all the output of all the supercell DFT calculations, creates the BORN file from the ph.out file, creates a Result folder and  calculates the IR Spectra in a frequency  range from Nu_min to Nu_max. 
Make sure to update the Molecule name ("Molecule"), number of atoms ("natoms") and also the minimum and maximum frequency in the plot (Nu_min and Nu_max).

Also, the following commands must be perfomrd IN THE LINUX SHELL before running the code:

command 1:   conda activate phonopy
command 2:   export PYTHONPATH=${PYTHONPATH}:/mnt/cephfs/home/armin5/ARMIN_env/bin/Phonopy-Spectroscopy-master
command 3:   export PATH=${PATH}:/mnt/cephfs/home/armin5/ARMIN_env/bin/Phonopy-Spectroscopy-master/Scripts

space group websites:

http://gernot-katzers-spice-pages.com/character_tables/D3.html
https://cryst.ehu.es/cgi-bin/cryst/programs/nph-table?from=raman
http://gernot-katzers-spice-pages.com/character_tables/index.html

"""
import glob
import os
import sys
import numpy as np
import math
import shutil
from parameters import Molecule,nat,tr2_ph,N_proc,nmix_ph,niter_ph,Na,Nb,Nc



cwd = os.getcwd()
os.chdir(cwd)
band_no=555555555555  # Just a random number to get into the while loop for getting Raman-active modes! 
NN1=0 # will take care of number of generated displaced supercells for the dielectric calculations by the spectroscopy code

N1=int(input("Please enter the number of displaced supercells created by phonopy in the previous steps: "))

 

         

# Creting the force set file  and the mesh.yaml file based on thesupercell DFT outputs

os.system('conda activate phonopy')

os.system('phonopy --qe -f '+Molecule+'-{001..'+str(N1).zfill(3)+'}.out --nac')

# creating the setting.conf file 
with open('settings.conf', 'w') as settings:  
    settings.write('DIM ='+str(Na)+' '+str(Nb)+' '+str(Nc)+ '\n')
    settings.write('FC_SYMMETRY = .TRUE.' + '\n')
    settings.write('FORCE_CONSTANTS = WRITE' + '\n')
    settings.write('HDF5 = .TRUE.' + '\n')
    settings.write('NAC = .TRUE.' + '\n')
    settings.close()

os.system('phonopy  --qe -c scf.in settings.conf')

os.system('phonopy --qe  --dim='+'"'+str(Na)+' '+str(Nb)+' '+str(Nc)+'"' +' --readfc --hdf5 --fc-symmetry --mesh="1 1 1" --eigenvectors -c scf.in')
os.system(' phonopy --dim='+'"'+str(Na)+' '+str(Nb)+' '+str(Nc)+'"'+' --readfc --hdf5 --fc-symmetry --irreps="0 0 0" ')

os.system('cat irreps.yaml')


whichmode=int(input("Do you want to consider all modes for dielectric calculations (1), or wish to only specify Raman active modes(2)? Enter either (1) or (2):  " ))

if whichmode==1:
 bands=""
 for x in range(1,3*nat+1):
  bands=bands+' '+str(x)
 NN1=3*nat
elif whichmode==2:
 bands=""
 while band_no!=0:
   band_no=int(input("Please enter each Raman active mode number one by one followed by Enter. Enter 0 after the last one:  "))
   if band_no!=0:
    bands=bands+' '+str(band_no)
    NN1=NN1+1  # getting the number of generated displaced supercells for dielectric calculation 
                                                                                                 
print("Raman active band numbers are: "+bands+','+"Total number of Raman active bands is: "+str(NN1)) 



os.system('phonopy-raman -d --bands'+'=" '+bands+' " ')
os.system('python3 poscar2qe.py')


os.system('mkdir dielectric_cal')
os.system('scp scf.*.*.in dielectric_cal')
os.chdir('dielectric_cal')
os.system('mkdir part1 part2 part3')


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

os.system('scp ph.in part1')
os.system('scp ph.in part2')
os.system('scp ph.in part3')



N2=math.floor(NN1/3)
bandlist=bands.split()
bandlist1=bandlist[0:N2]
bandlist2=bandlist[N2:2*N2]
bandlist3=bandlist[2*N2:]


for item in bandlist1:
  command1="scp scf."+item.zfill(4)+".001.in part1"
  command2="scp scf."+item.zfill(4)+".002.in part1"
  os.system(command1)
  os.system(command2)
  
for item in bandlist2:
  command3="scp scf."+item.zfill(4)+".001.in part2"
  command4="scp scf."+item.zfill(4)+".002.in part2"
  os.system(command3)
  os.system(command4)
  
for item in bandlist3:
  command5="scp scf."+item.zfill(4)+".001.in part3"
  command6="scp scf."+item.zfill(4)+".002.in part3"
  os.system(command5)
  os.system(command6)

commandd1=""
commandd2=""
commandd3=""

for item in bandlist1:
  commandd1=commandd1+item.zfill(4)+","
for item in bandlist2:
  commandd2=commandd2+item.zfill(4)+","
for item in bandlist3:
  commandd3=commandd3+item.zfill(4)+","
commandd1=commandd1.rstrip(commandd1[-1]) # Getting rid of the las "," in the command  
commandd2=commandd2.rstrip(commandd2[-1])  
commandd3=commandd3.rstrip(commandd3[-1])  

## Creating the submission batch script for part1

with open('scf1_sub', 'w') as scf1:
        scf1.write('#!/bin/bash' + '\n')
        scf1.write('#SBATCH -N 1' + '\n')
        scf1.write('#SBATCH --tasks-per-node='+str(N_proc)+'\n')
        scf1.write('#SBATCH --mem-per-cpu=7G' + '\n')
        scf1.write('#SBATCH --job-name='+ Molecule+'d1' + '\n')
        scf1.write('\n')
        scf1.write('\n')
        scf1.write('module purge'+'\n')
        scf1.write('module load ohpc'+'\n')
        scf1.write('module load quantum-espresso-7.0-gcc-8.3.0-6pkcl6c'+'\n')
        scf1.write('\n')
        scf1.write('\n')
        scf1.write('\n')
        scf1.write('for N in {'+commandd1+'}'+'\n')
        scf1.write('do'+'\n')
        scf1.write('for i in {001,002}'+'\n')
        scf1.write('do'+'\n')
        scf1.write('mpirun'+'  '+ 'pw.x'+'  '+ '< '+'scf.$N.$i.in'+''+'>'+' '+'scf.$N.$i.out'  + '\n')
        scf1.write('mpirun'+'  '+ 'ph.x'+'  '+ '< '+'ph.in'+''+'>'+' '+'ph.$N.$i.out'  + '\n')
        scf1.write('done'+'\n')
        scf1.write('done'+'\n')
        scf1.write('\n')
scf1.close()

## Creating the submission batch script for part2

with open('scf2_sub', 'w') as scf2:
        scf2.write('#!/bin/bash' + '\n')
        scf2.write('#SBATCH -N 1' + '\n')
        scf2.write('#SBATCH --tasks-per-node='+str(N_proc)+'\n')
        scf2.write('#SBATCH --mem-per-cpu=7G' + '\n')
        scf2.write('#SBATCH --job-name='+ Molecule+'d2' + '\n')
        scf2.write('\n')
        scf2.write('\n')
        scf2.write('module purge'+'\n')
        scf2.write('module load ohpc'+'\n')
        scf2.write('module load quantum-espresso-7.0-gcc-8.3.0-6pkcl6c'+'\n')
        scf2.write('\n')
        scf2.write('\n')
        scf2.write('\n')
        scf2.write('for N in {'+commandd2+'}'+'\n')
        scf2.write('do'+'\n')
        scf2.write('for i in {001,002}'+'\n')
        scf2.write('do'+'\n')
        scf2.write('mpirun'+'  '+ 'pw.x'+'  '+ '< '+'scf.$N.$i.in'+''+'>'+' '+'scf.$N.$i.out'  + '\n')
        scf2.write('mpirun'+'  '+ 'ph.x'+'  '+ '< '+'ph.in'+''+'>'+' '+'ph.$N.$i.out'  + '\n')
        scf2.write('done'+'\n')
        scf2.write('done'+'\n')
        scf2.write('\n')
scf2.close()

## Creating the submission batch script for part3

with open('scf3_sub', 'w') as scf3:
        scf3.write('#!/bin/bash' + '\n')
        scf3.write('#SBATCH -N 1' + '\n')
        scf3.write('#SBATCH --tasks-per-node='+str(N_proc)+'\n')
        scf3.write('#SBATCH --mem-per-cpu=7G' + '\n')
        scf3.write('#SBATCH --job-name='+ Molecule+'d3' + '\n')
        scf3.write('\n')
        scf3.write('\n')
        scf3.write('module purge'+'\n')
        scf3.write('module load ohpc'+'\n')
        scf3.write('module load quantum-espresso-7.0-gcc-8.3.0-6pkcl6c'+'\n')
        scf3.write('\n')
        scf3.write('\n')
        scf3.write('\n')
        scf3.write('for N in {'+commandd3+'}'+'\n')
        scf3.write('do'+'\n')
        scf3.write('for i in {001,002}'+'\n')
        scf3.write('do'+'\n')
        scf3.write('mpirun'+'  '+ 'pw.x'+'  '+ '< '+'scf.$N.$i.in'+''+'>'+' '+'scf.$N.$i.out'  + '\n')
        scf3.write('mpirun'+'  '+ 'ph.x'+'  '+ '< '+'ph.in'+''+'>'+' '+'ph.$N.$i.out'  + '\n')
        scf3.write('done'+'\n')
        scf3.write('done'+'\n')
        scf3.write('\n')
 #       scf3.write('sleep 2400'+'\n')
 #       scf3.write('cd ..'+'\n')
 #       scf3.write('cd ..'+'\n')
 #       scf3.write('source /mnt/cephfs/home/armin5/ARMIN_env/bin/activate'+ '\n')
 #       scf3.write('export PYTHONPATH=${PYTHONPATH}:/mnt/cephfs/home/armin5/ARMIN_env/bin/Phonopy-Spectroscopy-master'+ '\n')
 #       scf3.write('export PATH=${PATH}:/mnt/cephfs/home/armin5/ARMIN_env/bin/Phonopy-Spectroscopy-master/Scripts'+ '\n')
 #       scf3.write('python3 IR_Raman_Step3.py'+ '\n')
 #       scf3.write('\n')
scf3.close()


os.system('scp scf1_sub part1')
os.system('scp scf2_sub part2')
os.system('scp scf3_sub part3')
os.system('rm scf1_sub scf2_sub scf3_sub')

#os.chdir('part1')
##os.system('sbatch scf1_sub')
#os.chdir(cwd)
#os.chdir('dielectric_cal')
#os.chdir('part2')
##os.system('sbatch scf2_sub')
#os.chdir(cwd)
#os.chdir('dielectric_cal')
#os.chdir('part3')
##os.system('sbatch scf3_sub')
#os.chdir(cwd)
