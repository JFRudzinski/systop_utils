#!/usr/bin/env python

import os

py_interp = '/Users/joedinski/miniconda3/envs/mosdef37/bin/python3.7'
path = '/Users/joedinski/work/Projects/FAIRmat/LAMMPS_DEV/MD_Overview/KA/1_xyz_file/'
sys_nm = 'KA_1_xyz_file'
top_fnm = path+'log.lammps'
conf_fnm = path+'pos_vel.xyz'
top_format = 'DATA'
conf_format = 'XYZ'

os.system(py_interp+' '+'../../Get_sys-top_MDAnalysis-to-MBUILD_w-Overview_general.py'+' '+sys_nm+' '+top_fnm+' '+conf_fnm+' '+top_format+' '+conf_format)