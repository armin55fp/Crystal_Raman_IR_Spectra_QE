# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 10:49:27 2022

@author: armin.taheri

This code uses Phonopy-Spectroscopy package to create final Raman and IR spectra results.
command 1:   source /mnt/cephfs/home/armin5/ARMIN_env/bin/activate
command 2:   export PYTHONPATH=${PYTHONPATH}:/mnt/cephfs/home/armin5/ARMIN_env/bin/Phonopy-Spectroscopy-master
command 3:   export PATH=${PATH}:/mnt/cephfs/home/armin5/ARMIN_env/bin/Phonopy-Spectroscopy-master/Scripts

"""
import glob
import os
import sys
import numpy as np
import math
import shutil

w_min=float(input("Please enter the lowest frequency in the plot in cm^-1:  "))
w_max=float(input("Please enter the miximum frequency in the plot in cn^-1:  "))


cwd = os.getcwd()
os.chdir(cwd)
os.chdir('dielectric_cal')
os.chdir('part1')
os.system('scp ph.*.*.out'+' '+str(cwd))
os.chdir(cwd)
os.chdir('dielectric_cal')
os.chdir('part2')
os.system('scp ph.*.*.out'+' '+str(cwd))
os.chdir(cwd)
os.chdir('dielectric_cal')
os.chdir('part3')
os.system('scp ph.*.*.out'+' '+str(cwd))
os.chdir(cwd)


os.system('source /mnt/cephfs/home/armin5/ARMIN_env/bin/activate')
os.system('phonopy-ir --ir_reps --spectrum_range='+'"'+str(w_min)+' '+str(w_max)+'"' )
os.system('python3 qe2outcar.py')
os.system('phonopy-raman -r OUTCAR.*')
os.system('phonopy-raman -p --ir_reps --spectrum_range='+'"'+str(w_min)+' '+str(w_max)+'"')
print ('done!')
