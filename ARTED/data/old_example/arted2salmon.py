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
SALMON-TDDFT (preliminary version) input file formats.
"""

Usage = """%prog -t sc|ms < arted_inputfile > salmon_inputfile"""

    
import sys
from optparse import OptionParser
from collections import deque

input_list_sc = [
    "entrance_option", "Time_shutdown", "backup_frequency","entrance_iter", "SYSname",
    "directory", "functional", "cval", "propagator","ps_format", "PSmask_option",
    "alpha_mask", "gamma_mask", "eta_mask", "aL", "ax", "ay", "az", "Sym",
    "crystal_structure", "Nd", "NLx", "NLy", "NLz", "NKx", "NKy", "NKz",
    "file_kw", "NEwald", "aEwald", "KbTev", "NB", "Nelec", "FSset_option",
    "Ncg", "Nmemory_MB", "alpha_MB", "NFSset_start", "NFSset_every", "Nscf",
    "ext_field", "Longi_Trans", "MD_option", "AD_RHO", "Nt", "dt", "dAc",
    "Nomega", "domega", "AE_shape", "IWcm2_1", "tpulsefs_1", "omegaev_1",
    "phi_CEP_1", "Epdir_1", "IWcm2_2", "tpulsefs_2", "omegaev_2", "phi_CEP_2",
    "Epdir_2", "T1_T2fs", "NI", "NE",
]

input_list_ms = [
    "entrance_option", "Time_shutdown", "backup_frequency", "entrance_iter", "SYSname",
    "directory", "functional", "cval", "propagator", "ps_format", "PSmask_option",
    "alpha_mask", "gamma_mask", "eta_mask", "aL", "ax", "ay", "az", "Sym",
    "crystal_structure", "Nd", "NLx", "NLy", "NLz", "NKx", "NKy", "NKz",
    "file_kw", "FDTDdim", "TwoD_shape", "NX_m", "NY_m", "HX_m", "HY_m",
    "NKsplit", "NXYsplit", "NXvacL_m", "NXvacR_m", "NEwald", "aEwald",
    "KbTev", "NB", "Nelec", "FSset_option", "Ncg", "Nmemory_MB", "alpha_MB",
    "NFSset_start", "NFSset_every", "Nscf", "MD_option", "AD_RHO", "Nt", "dt",
    "dAc", "Nomega", "domega", "AE_shape", "IWcm2_1", "tpulsefs_1",
    "omegaev_1", "phi_CEP_1", "Epdir_1", "IWcm2_2", "tpulsefs_2", "omegaev_2",
    "phi_CEP_2", "Epdir_2", "T1_T2fs", "NI", "NE",
]

output_list = [
    ["group_function", [
        "cfunction",
    ]],
    ["control", [
        "entrance_option", "Time_shutdown", "backup_frequency", "entrance_iter",
        "SYSname", "directory",
    ]],
    ["system", [
        "functional", "cval", "aL", "ax", "ay", "az", "Sym",
        "crystal_structure", "NB", "Nelec", "ext_field", "MD_option",
        "AD_RHO", "NE", "NI",
    ]],
    ["propagation", [
        "propagator",
    ]],
    ["rgrid", [
        "Nd", "NLx", "NLy", "NLz",
    ]],
    ["kgrid", [
        "NKx", "NKy", "NKz", "file_kw",
    ]],
    ["tstep", [
        "Nt", "dt",
    ]],
    ["pseudo", [
        "ps_format", "PSmask_option", "alpha_mask", "gamma_mask", "eta_mask",
    ]],
    ["electrons", [
        "NEwald", "aEwald", "KbTev", "Ncg", "Nmemory_MB", "alpha_MB",
        "FSset_option", "NFSset_start", "NFSset_every", "Nscf",
    ]],
    ["incident", [
        "Longi_Trans", "dAc", "AE_shape", "IWcm2_1", "tpulsefs_1",
        "omegaev_1", "phi_CEP_1", "Epdir_1", "IWcm2_2", "tpulsefs_2",
        "omegaev_2", "phi_CEP_2", "Epdir_2", "T1_T2fs",
    ]],
    ["response", [
        "Nomega", "domega",
    ]],
    ["multiscale", [
        "FDTDdim", "TwoD_shape", "NX_m", "NY_m", "HX_m", "HY_m", "NKsplit",
        "NXYsplit", "NXvacL_m", "NXvacR_m",
    ]],
]

parser = OptionParser(description=Description, usage=Usage)
parser.add_option("-t", "--type", metavar="MODE", dest="MODE",
                  default="", type=str, help="calculation type: 'sc' or 'ms'")
opts, args = parser.parse_args()

if opts.MODE.lower() == "sc":
    input_q = deque(input_list_sc)
    data = {"cfunction": "'singlecell'"}
elif opts.MODE.lower() == "ms":
    input_q = deque(input_list_ms)
    data = {"cfunction": "'multiscale'"}
else:
    parser.print_help(sys.stderr)
    sys.exit(-1)

buff = deque([])
for line in sys.stdin:
    temp = line.split("!")[0]
    temp = temp.replace(",", " ")
    buff += temp.split()

while input_q:
    param = input_q.popleft()
    if (param == "file_kw") and (int(data["NKx"]) > 0):
        continue
    elif param in ["Epdir_1", "Epdir_2"]:
        data[param] = ",".join([buff.popleft() for i in xrange(3)])
    else:
        data[param] = buff.popleft()

for group, param_list in output_list:
    print("&%s" % group)
    for param in param_list:
        if param in data:
            print("\t%s=%s" % (param, data[param]))
    print("\t/")

NE = int(data["NE"])
NI = int(data["NI"])

Zatom = []
for i in range(NE):
    Zatom.append(buff.popleft())

Lref = []
for i in range(NE):
    Lref.append(buff.popleft())

print("&atomic_spiecies")
for i in range(NE):
    print("\t%d\t%s\t%s" % (i + 1, Zatom[i], Lref[i]))
print("\t/")

print("&atomic_positions")
for i in range(NI):
    print("\t".join([""] + [buff.popleft() for j in range(5)]))
print("\t/")
