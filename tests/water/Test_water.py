#!/usr/bin/env python

import os

py_interp = '/Users/joedinski/miniconda3/envs/mosdef37/bin/python3.7'

path = '/Users/joedinski/work/Projects/FAIRmat/nomad/dependencies/parsers/gromacs/tests/data/water/'
sys_nm = 'water'
top_fnm = path+'topol.tpr'
conf_fnm = path+'conf.gro'
top_format = 'None'
conf_format = 'None'

os.system(py_interp+' '+'../Get_sys-top_MDAnalysis-to-MBUILD_w-Overview_general.py'+' '+sys_nm+' '+top_fnm+' '+conf_fnm+' '+top_format+' '+conf_format)