#!/usr/bin/env python

import os

py_interp = '/Users/joedinski/miniconda3/envs/mosdef37/bin/python3.7'
path = '/Users/joedinski/work/Projects/FAIRmat/nomad/dependencies/parsers/lammps/tests/data/methane_xyz/'
sys_nm = 'methane_xyz'
top_fnm = path+'data.64xmethane_from_restart'
conf_fnm = path+'64xmethane-nvt.xyz'
top_format = 'DATA'
conf_format = 'XYZ'

os.system(py_interp+' '+'../Get_sys-top_MDAnalysis-to-MBUILD_w-Overview_general.py'+' '+sys_nm+' '+top_fnm+' '+conf_fnm+' '+top_format+' '+conf_format)