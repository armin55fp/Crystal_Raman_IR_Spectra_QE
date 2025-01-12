# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 10:49:27 2022

@author: armin.taheri

This code gathers all the output of all the supercell DFT calculations, creates the BORN file from the ph.out file, creates a Result folder and  calculates the IR Spectra in a frequency  range from Nu_min to Nu_max. 
Make sure to update the Molecule name ("Molecule"), number of atoms ("natoms") and also the minimum and maximum frequency in the plot (Nu_min and Nu_max).

Also, the following commands must be perfomrd IN THE LINUX SHELL before running the code:

command 1:   source /mnt/cephfs/home/armin5/ARMIN_env/bin/activate
command 2:   export PYTHOATH=${PYTHONPATH}:/mnt/cephfs/home/armin5/ARMIN_env/bin/Phonopy-Spectroscopy-master
command 3:   export PATH=${PATH}:/mnt/cephfs/home/armin5/ARMIN_env/bin/Phonopy-Spectroscopy-master/Scripts

"""
import glob
import os
import sys
import numpy as np
import math
import shutil
from parameters import Molecule,nat


nonequivalent=[1,2] # A list to keep track the number of non-equivalent atoms in the BORN file

N1=5555555555555555  # just a random number to get into the while loop!

mode=int(input("Do you want to consider all atoms of the unitcell in the BORN file (1), or wish to only consider non-equivalent atoms (2)? Enter either (1) or (2):  " ))


if mode==2:
 while N1!=0:
   N1=int(input("Please enter the number associated with non-equivalent atoms one by one and then press enter. enter 0 after the last one:  "))
   if N1!=0:
    nonequivalent.append(N1+2)





cwd = os.getcwd()
os.system('mkdir Result')
ResAdd=str(cwd)+'/Result'
os.chdir('part1')
os.system('scp'+' '+ Molecule+'-*.out'+' '+str(cwd))
os.chdir(cwd)
os.chdir('part2')
os.system('scp'+' '+ Molecule+'-*.out'+' '+str(cwd))
os.chdir(cwd)
os.chdir('part3')
os.system('scp'+' '+ Molecule+'-*.out'+' '+str(cwd))
os.chdir(cwd)


#############################################################################################
os.chdir('BORN_cal')
shutil.copy('ph.out',ResAdd)
os.chdir(cwd)
os.chdir('Result')

# Extracting  the dielectric constant and Born charges from ph.out

dielectricx = []
dielectricy = []
dielectricz = []
Born_x = []
Born_y = []
Born_z = []

#  The dielectric constant


k=0
with open('ph.out','r') as f:
    lines=f.readlines()
    for line in lines:
        linestrip = line.strip()
        entry = line.split()
        if 'Dielectric constant in cartesian axis' in line:
            N=k
        k+=1

j=0
with open('ph.out','r') as f:
    lines=f.readlines()
    for line in lines:
        linestrip = line.strip()
        entry = linestrip.split()
        if N+2 <= j <= N+4:

            dielectricx.append(entry[1])
            dielectricy.append(entry[2])
            dielectricz.append(entry[3])
        j+=1


#The  born effective charges

kk=0
with open('ph.out','r') as f:
    lines=f.readlines()
    for line in lines:
        linestrip = line.strip()
        entry = line.split()
        if 'Effective charges (d Force / dE) in cartesian axis without acoustic sum rule applied (asr)' in line:
            NN=kk
        kk+=1

jj=0
with open('ph.out','r') as f:
    lines=f.readlines()
    for line in lines:
        linestrip = line.strip()
        entry = linestrip.split()
        if NN+3 <= jj <= NN+(4*nat+1):
            if (jj-NN)%(4) != 2:
                #print(entry)

                Born_x.append(entry[2])
                Born_y.append(entry[3])
                Born_z.append(entry[4])
        jj+=1
#----------------------------------------------
# Wrting the dielectric constant and the Born eggective charges in the BORN file


p=0
# dielectic constant
with open('BORN', 'w') as born:
     born.write('2'+'\n')
     born.write(dielectricx[0]+'\t'+dielectricy[0]+'\t'+dielectricz[0]+'\t'+dielectricx[1]+'\t'+dielectricy[1]+'\t'+dielectricz[1]+'\t'+dielectricx[2]+'\t'+dielectricy[2]+'\t'+dielectricz[2]+'\n')
# Born charges
     for p in range(nat):
         born.write(Born_x[3*p]+' '+Born_y[3*p]+' '+Born_z[3*p]+' '+Born_x[3*p+1]+' '+Born_y[3*p+1]+' '+Born_z[3*p+1]+' '+Born_x[3*p+2]+' '+Born_y[3*p+2]+' '+Born_z[3*p+2]+'\n')

born.close()


#-------------------------------------------------
os.system('scp'+' '+ 'BORN'+' '+str(cwd))
os.chdir(cwd)

if mode==2:
 todelete=''
 for L in range(1,nat+3):
   if L not in nonequivalent:
     todelete=todelete+str(L)+'d'+';'
 command="sed -i.bak  -e"+" "+"'"+todelete+"'"+' '+"BORN"  

 os.system(command)
    
    
                                              
 
 
  










