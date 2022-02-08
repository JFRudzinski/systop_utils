#!/usr/bin/env python

import os

py_interp = '/Users/joedinski/miniconda3/envs/mosdef37/bin/python3.7'
path = '/Users/joedinski/work/Projects/FAIRmat/nomad/dependencies/parsers/lammps/tests/data/hexane_cyclohexane/'
sys_nm = 'hexane_cyclohexane'
top_fnm = path+'data.'+sys_nm
conf_fnm = path+sys_nm+'_nvt_20ps.lammpstrj'
top_format = 'DATA'
conf_format = 'LAMMPSDUMP'

os.system(py_interp+' '+'../Get_sys-top_MDAnalysis-to-MBUILD_w-Overview_general.py'+' '+sys_nm+' '+top_fnm+' '+conf_fnm+' '+top_format+' '+conf_format)