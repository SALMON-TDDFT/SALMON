#!/usr/bin/env python
#
#   Copyright 2017 ARTED developers
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#

Description = """
ARTED to SALMON input file converter;
This script helps to translates the ordinally ARTED input files to
SALMON-TDDFT (v.1.0.0) input file formats.
"""

Usage = """%prog -t sc|ms < ARTED_file > SALMON_file"""

import sys
from optparse import OptionParser
from collections import deque, defaultdict
from math import pi

param_list_sc = [
    'entrance_option', 'Time_shutdown', 'backup', 'entrance_iter', 'SYSname', 
    'directory', 'functional', 'cval', 'propagator', 'ps_format', 'PSmask_option', 
    'alpha_mask', 'gamma_mask', 'eta_mask', 'aL', 'ax', 'ay', 'az', 'Sym', 
    'crystal_structure', 'Nd', 'NLx', 'NLy', 'NLz', 'NQx', 'NQy', 'NQz', 
    'file_kw', 'NEwald', 'aEwald', 'KbTev', 'NB', 'Nelec', 'FSset_option', 
    'Ncg', 'Nmemory_MB', 'alpha_MB', 'NFSset_start', 'NFSset_every', 'Nscf', 
    'ext_field', 'Longi_Trans', 'MD_option', 'AD_RHO', 'Nt', 'dt', 'dAc', 
    'Nomega', 'domega', 'AE_shape', 'IWcm2_1', 'tpulsefs_1', 'omegaev_1', 
    'phi_CEP_1', 'Epdir_1_x', 'Epdir_1_y', 'Epdir_1_z', 'IWcm2_2', 
    'tpulsefs_2', 'omegaev_2', 'phi_CEP_2', 'Epdir_2_x', 'Epdir_2_y', 
    'Epdir_2_z', 'T1_T2fs', 'NI', 'NE', 
]

param_list_ms = [
    'entrance_option', 'Time_shutdown', 'backup', 'entrance_iter', 'SYSname', 
    'directory', 'functional', 'cval', 'propagator', 'ps_format', 'PSmask_option', 
    'alpha_mask', 'gamma_mask', 'eta_mask', 'aL', 'ax', 'ay', 'az', 'Sym', 
    'crystal_structure', 'Nd', 'NLx', 'NLy', 'NLz', 'NQx', 'NQy', 'NQz', 
    'FDTDdim', 'TwoD_shape', 'NX_m', 'NY_m', 'HX_m', 'HY_m', 'NKsplit', 
    'NXYsplit', 'NXvacL_m', 'NXvacR_m', 'NEwald', 'aEwald', 'KbTev', 'NB', 
    'Nelec', 'FSset_option', 'Ncg', 'Nmemory_MB', 'alpha_MB', 'NFSset_start', 
    'NFSset_every', 'Nscf', 'MD_option', 'AD_RHO', 'Nt', 'dt', 'dAc', 
    'Nomega', 'domega', 'AE_shape', 'IWcm2_1', 'tpulsefs_1', 'omegaev_1', 
    'phi_CEP_1', 'Epdir_1_x', 'Epdir_1_y', 'Epdir_1_z', 'IWcm2_2', 
    'tpulsefs_2', 'omegaev_2', 'phi_CEP_2', 'Epdir_2_x', 'Epdir_2_y', 
    'Epdir_2_z', 'T1_T2fs', 'NI', 'NE', 
]

au_fs = +2.418884e-02
au_ev = +2.721139e+01

element_list = [
    'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 
    'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 
    'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 
    'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 
    'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 
    'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 
    'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 
]

ps_suffix_tbl = {"ky": "_rps.dat", "abinit": ".fhi", 'fhi': ".cpi"}

aeshape_tbl = {"asin2cos": "Acos2", "esin2cos": "Ecos2"}

restart_tbl = {"new": "new", "reentrance": "restart"}

def str_f(lit):
    return lit.strip().replace("'", "").replace('"', '')

def float_f(lit):
    return float(lit.replace("d", "e").replace("D", "E"))

def int_f(lit):
    return int(round(float_f(lit)))

def repr_f(item):
    if type(item) is str:
        return "'%s'" % str_f(item)
    elif type(item) is int:
        return "%d" % item
    elif type(item) is float:
        return ("%e" % item).replace("e", "d")
    elif type(item) is list:
        return "%s" % ",".join([repr_f(x) for x in item])



parser = OptionParser(description=Description, usage=Usage)
parser.add_option("-t", "--type", metavar="type_calc", dest="type_calc",
                  default="", type=str, help="calculation type: 'sc' or 'ms'")
opts, args = parser.parse_args()

if opts.type_calc == "sc":
    param_q = deque(param_list_sc)
elif opts.type_calc == "ms":
    param_q = deque(param_list_ms)
else:
    parser.print_help(sys.stderr)
    sys.exit(-1)

buff = deque([])
for line in sys.stdin:
    buff += line.split("!")[0].replace(",", " ").split()


##########################################################################
ARTED = {}
while param_q:
    param_name = param_q.popleft()
    if (param_name == "file_kw"
            and 0 < int_f(ARTED["NQx"])
            and 0 < int_f(ARTED["NQy"])
            and 0 < int_f(ARTED["NQz"])):
        continue
    else:
        ARTED[param_name]=buff.popleft()

NE = int_f(ARTED["NE"])
NI = int_f(ARTED["NI"])
PS = ps_suffix_tbl[str_f(ARTED["ps_format"]).lower()]

izatom_list = []
pseudo_file_list = [] 
for i in range(NE):
    iz = int_f(buff.popleft())
    izatom_list.append(iz)
    pseudo_file_list.append(element_list[iz-1] + PS)
    
lloc_ps_list = []
for i in range(NE):
    lref = int_f(buff.popleft())
    lloc_ps_list.append(lref)

position_list = []
for i in range(NI):
    ia = int_f(buff.popleft())
    rx = float_f(buff.popleft())
    ry = float_f(buff.popleft())
    rz = float_f(buff.popleft())
    ii = int_f(buff.popleft())
    position_list.append(
        (element_list[izatom_list[ii-1]-1], rx, ry, rz, ii)
    )
    
##########################################################################
SALMON = defaultdict(dict)

SALMON['calculation']['calc_mode'] = 'GS_RT'
SALMON['calculation']['use_ehrenfest_md'] = str_f(ARTED['MD_option']).lower()
#SALMON['calculation']['use_force'] 
#SALMON['calculation']['use_geometry_opt'] 

SALMON['control']['restart_option'] = restart_tbl[
    str_f(ARTED['entrance_option']).lower()
]
SALMON['control']['backup_frequency'] = int_f(ARTED['backup'])
SALMON['control']['time_shutdown'] = int_f(ARTED['Time_shutdown'])
SALMON['control']['sysname'] = str_f(ARTED['SYSname'])
#SALMON['control']['dump_filename']

SALMON['units']['unit_system'] = "a.u."

SALMON['system']['iperiodic'] = 3
#SALMON['system']['ispin']

aL = float_f(ARTED['aL'])
ax = float_f(ARTED['ax'])
ay = float_f(ARTED['ay'])
az = float_f(ARTED['az'])

SALMON['system']['al'] = [aL*ax, aL*ay, aL*az]
SALMON['system']['isym'] = int_f(ARTED['Sym'])
SALMON['system']['crystal_structure'] = str_f(ARTED['crystal_structure'])
SALMON['system']['nstate'] = int_f(ARTED['NB'])
SALMON['system']['nelec'] = int_f(ARTED['Nelec'])
SALMON['system']['nelem'] = NE
SALMON['system']['natom'] = NI
temperature = float_f(ARTED['KbTev']) / au_ev
if 0.0 <= temperature:
    SALMON['system']['temperature'] = temperature
else:
    SALMON['system']['temperature'] = -1.0

SALMON['pseudo']['pseudo_file'] = pseudo_file_list
#SALMON['pseudo']['Lmax_ps']
SALMON['pseudo']['Lloc_ps'] = lloc_ps_list
SALMON['pseudo']['iZatom'] = izatom_list
SALMON['pseudo']['psmask_option'] = str_f(ARTED['PSmask_option'])
SALMON['pseudo']['alpha_mask'] = float_f(ARTED['alpha_mask'])
SALMON['pseudo']['gamma_mask'] = float_f(ARTED['gamma_mask'])
SALMON['pseudo']['eta_mask'] = float_f(ARTED['eta_mask'])

SALMON['functional']['xc'] = str_f(ARTED['functional'])
SALMON['functional']['cval'] = float_f(ARTED['cval'])

#SALMON['rgrid']['dl']
SALMON['rgrid']['num_rgrid'] = [
    int_f(ARTED["NLx"]), int_f(ARTED["NLy"]), int_f(ARTED["NLz"])
]

SALMON['kgrid']['num_kgrid'] = [
    int_f(ARTED["NQx"]), int_f(ARTED["NQy"]), int_f(ARTED["NQz"])
]
if "file_kw" in ARTED:
    SALMON['kgrid']['file_kw'] = str_f(ARTED['file_kw'])

SALMON['tgrid']['nt'] = int_f(ARTED['Nt'])
SALMON['tgrid']['dt'] = float_f(ARTED['dt'])

SALMON['propagation']['propagator'] = str_f(ARTED['propagator'])

SALMON['scf']['ncg'] = int_f(ARTED['Ncg'])
SALMON['scf']['nmemory_mb'] = int_f(ARTED['Nmemory_MB'])
SALMON['scf']['alpha_mb'] = float_f(ARTED['alpha_MB'])
SALMON['scf']['fsset_option'] = str_f(ARTED['FSset_option'])
SALMON['scf']['nfsset_start'] = int_f(ARTED['NFSset_start'])
SALMON['scf']['nfsset_every'] = int_f(ARTED['NFSset_every'])
SALMON['scf']['nscf'] = int_f(ARTED['Nscf'])

if opts.type_calc == "ms" or str_f(ARTED["ext_field"]).lower() == "lf":
    aeshape = aeshape_tbl[str_f(ARTED['AE_shape']).lower()]
    tw1 = float_f(ARTED['tpulsefs_1']) / au_fs
    tw2 = float_f(ARTED['tpulsefs_2']) / au_fs
    omega1 = float_f(ARTED['omegaev_1']) / au_ev
    omega2 = float_f(ARTED['omegaev_2']) / au_ev
    old_cep1 = float_f(ARTED['phi_CEP_1'])
    old_cep2 = float_f(ARTED['phi_CEP_2'])
    new_cep1 = (old_cep1 + omega1 * tw1 / (4.0 * pi)) % 1
    new_cep2 = (old_cep2 + omega2 * tw2 / (4.0 * pi)) % 1
    SALMON['emfield']['ae_shape1'] = aeshape
    SALMON['emfield']['ae_shape2'] = aeshape
    SALMON['emfield']['pulse_tw1'] = tw1
    SALMON['emfield']['pulse_tw2'] = tw2
    SALMON['emfield']['omega1'] = omega1
    SALMON['emfield']['omega2'] = omega2
    SALMON['emfield']['phi_cep1'] = new_cep1
    SALMON['emfield']['phi_cep2'] = new_cep2
    #SALMON['emfield']['amplitude1']
    #SALMON['emfield']['amplitude2']
    SALMON['emfield']['rlaser_int_wcm2_1'] = float_f(ARTED['IWcm2_1'])
    SALMON['emfield']['rlaser_int_wcm2_2'] = float_f(ARTED['IWcm2_2'])
    SALMON['emfield']['t1_t2'] = float_f(ARTED['T1_T2fs']) / au_fs

SALMON['emfield']['epdir_re1'] = [
    float_f(ARTED['Epdir_1_x']),
    float_f(ARTED['Epdir_1_y']),
    float_f(ARTED['Epdir_1_z']),
]
SALMON['emfield']['epdir_re2'] = [
    float_f(ARTED['Epdir_2_x']),
    float_f(ARTED['Epdir_2_y']),
    float_f(ARTED['Epdir_2_z']),
]
#SALMON['emfield']['epdir_im1'] 
#SALMON['emfield']['epdir_im2'] 

SALMON['analysis']['projection_option'] = str_f(ARTED['AD_RHO']).lower()
#SALMON['analysis']['out_dos']
#SALMON['analysis']['out_dos_start']
#SALMON['analysis']['out_dos_end']
#SALMON['analysis']['iout_dos_nenergy']
#SALMON['analysis']['out_dos_smearing']
#SALMON['analysis']['out_dos_method'] 
#SALMON['analysis']['out_dos_fshift'] 
#SALMON['analysis']['out_dns'] 
#SALMON['analysis']['out_dns_rt'] 
#SALMON['analysis']['out_dns_rt_step'] 
#SALMON['analysis']['format3d'] 
#SALMON['analysis']['numfiles_out_3d'] 

SALMON['ewald']['newald'] = int_f(ARTED['NEwald'])
SALMON['ewald']['aewald'] = float_f(ARTED['aEwald'])

if opts.type_calc == "ms":
    SALMON['calculation']['use_ms_maxwell'] = 'y'
    SALMON['multiscale']['fdtddim'] = str_f(ARTED['FDTDdim'])
    SALMON['multiscale']['twod_shape'] = str_f(ARTED['TwoD_shape'])
    SALMON['multiscale']['nx_m'] = int_f(ARTED['NX_m'])
    SALMON['multiscale']['ny_m'] = int_f(ARTED['NY_m'])
    #SALMON['multiscale']['nz_m']
    SALMON['multiscale']['hx_m'] = float_f(ARTED['HX_m'])
    SALMON['multiscale']['hy_m'] = float_f(ARTED['HY_m'])
    #SALMON['multiscale']['hz_m']
    SALMON['multiscale']['nksplit'] = int_f(ARTED['NKsplit'])
    SALMON['multiscale']['nxysplit'] = int_f(ARTED['NXYsplit'])
    SALMON['multiscale']['nxvacl_m'] = int_f(ARTED['NXvacL_m'])
    SALMON['multiscale']['nxvacr_m'] = int_f(ARTED['NXvacR_m'])
elif opts.type_calc == "sc":
    SALMON['calculation']['use_ms_maxwell'] = 'n'
    SALMON['emfield']['trans_longi'] = str_f(ARTED['Longi_Trans']).lower()
    if str_f(ARTED["ext_field"]).lower() == "lr":
        SALMON['emfield']['ae_shape1'] = 'impulse'
        SALMON['emfield']['e_impulse'] = float_f(ARTED['dAc'])
        SALMON['analysis']['nenergy'] = int_f(ARTED['Nomega'])
        SALMON['analysis']['de'] = float_f(ARTED['domega'])

##########################################################################

for group in sorted(SALMON.keys()):
    sys.stdout.write("&%s\n" % group)
    for key in sorted(SALMON[group].keys()):
        sys.stdout.write("\t%s=%s\n" % (key, repr_f(SALMON[group][key])))
    sys.stdout.write("/\n")

sys.stdout.write("&atomic_red_coor\n")
for item in position_list:
    sys.stdout.write("\t'%s' %.5f %.5f %.5f %d\n" % tuple(item))
sys.stdout.write("/\n")
