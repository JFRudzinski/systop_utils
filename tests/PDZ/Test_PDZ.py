#!/usr/bin/env python

import os

py_interp = '/Users/joedinski/miniconda3/envs/mosdef37/bin/python3.7'

path = '/Users/joedinski/work/Projects/FAIRmat/GRO_DEV/MD_Overview/PDZ/PDZ2_SAP-97/no_ligand/'
sys_nm = '2awx'
top_fnm = path+sys_nm+'.tpr'
conf_fnm = path+sys_nm+'.NPT-out.gro'
top_format = 'None'
conf_format = 'None'

os.system(py_interp+' '+'../Get_sys-top_MDAnalysis-to-MBUILD_w-Overview_general.py'+' '+sys_nm+' '+top_fnm+' '+conf_fnm+' '+top_format+' '+conf_format)