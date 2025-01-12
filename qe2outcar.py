#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
This script takes in the output from QE and converts into the VASP OUTCAR format
'''

import glob
import os
import sys
import numpy as np
import math


atom_info = {'C': [12.010, 4.000], 'H': [1.0078, 1.000], 'O': [15.999, 6.000],
             'S': [32.065, 6.000], 'N': [14.0067, 5.000], 'Na': [22.9898, 1.000],
             'Ca': [40.078, 2.000], 'Al': [26.981539, 3.000], 'Si': [28.0855, 4.000],
             'Mg': [24.305, 2.000], 'P': [30.98376, 5.000], 'K': [39.0983, 1.000],
             'Se': [78.96, 6.000],'Ti':[47.867,4.000]}  # is ZVAL of Ti actually 4?  

list_names = []
for filename in sorted(glob.glob('ph.*.*.out')):
    name1 = filename.split('.out')[0]
    name2=name1.split('ph')
    list_names.append(name2)

    content = open(filename,"r")
    lines = content.readlines()
    out=open("OUTCAR"+name2[1], "w")
    sys.stdout = out

    atom_list = []
    num_atoms = -1
# Copying elements and masses

    for index,line in enumerate(lines):
        if 'site n.' in line:
            j=index+1
            while 'Computing dynamical' not in lines[j]:
                num_atoms = num_atoms + 1
                temp = lines[j]
                temp = temp.strip()
                atom_list.append(temp[6:8])
                j=j+1
            break
    atom_list = list(filter(None,atom_list))
    atom_list = [x.strip(' ') for x in atom_list]
    res = {}

    for i in atom_list:
        res[i] = atom_list.count(i)
    for key in res:
        out.write("   POMASS =   " + str(atom_info[key][0])  + ";" + "ZVAL   =    " + str(atom_info[key][1]) + "    mass and valenz" + '\n')
    out.write("   number of dos      NEDOS =    301   number of ions     NIONS =      " + str(num_atoms) + '\n')
    print("   ions per type =               ", end="")
    for key in res:
        print(str(res[key])+"   ", sep=' ', end='', flush=True)
    print('\n')
 
    for index,line in enumerate(lines):
        count1 = 0
        count2 = 0
        if 'Effective charges' in line:
            out.write(' BORN EFFECTIVE CHARGES (in e, cummulative output)' + '\n')
            out.write(' -------------------------------------------------' + '\n')
            for k in range(index+2, index+(4*num_atoms)+2):
                if k==index+2:
                    count1+=1
                    count2 =0
                    out.write(" ion    " + str(count1) + '\n')
                elif (k-index-2)%4==0:
                    count1+=1
                    count2 =0
                    out.write(" ion    " + str(count1) + '\n')
                else:
                    count2+=1
                    temp = lines[k]
                    temp = temp.strip()
                    temp = temp[8:]
                    temp = temp.rstrip(temp[-1])
                    out.write("    " + str(count2) + "    " + temp + '\n')
            break
    out.write('\n')
    out.write('\n')
    out.write("--------------------------------------------------------------------------------------------------------" + '\n')
    out.write('\n')
  
#Copying Dielectric constants
    for index,line in enumerate(lines):
        if 'Dielectric constant in cartesian axis' in line:
            out.write(" MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in DFT)" + '\n')
            out.write("------------------------------------------------------" + '\n')
            for i in range(index+2,index+5):
                temp = lines[i]
                temp = temp.strip()
                temp = temp.rstrip(temp[-1])
                temp = temp.lstrip(temp[0])
                out.write(temp+'\n')
                if i==index+4:
                    out.write("------------------------------------------------------" + '\n')
            break 

