#! /usr/bin/env python
#
#   Copyright 2017 SALMON developers
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

#
# Autotools (configure) like script with Python.
#
from optparse import OptionParser, OptionGroup
import os

SOURCE_DIR = os.path.dirname(__file__)

def on_or_off(v) :
  if v:
    return 'on'
  else:
    return 'off'

def debug_or_release(v) :
  if v:
    return 'Debug'
  else:
    return 'Release'

def add_option(dic, name, var) :
  if var is not None:
    dic[name] = on_or_off(var)

def add_env(dic, name, var) :
  if var is not None:
    dic[name] = var

usage  = "usage: %prog [options]"
parser = OptionParser(usage)

parser.add_option('-n', '--dry-run', action='store_true', default=False, dest='dry_run', help='don\'t actually run.')
parser.add_option('-v', '--verbose', action='store_true', default=False, dest='verbose', help='show verbose messages.')

group = OptionGroup(parser, 'Build target')
group.add_option('-a', '--arch',   action='store', default=None, dest='arch',   help='cross compile mode. ARCH format should be <COMPILER>-<SYSTEM>')
group.add_option('--prefix',       action='store', default=None, dest='prefix', help='package install prefix.')
parser.add_option_group(group)

group = OptionGroup(parser, 'Library options')
group.add_option('--enable-mpi',        action='store_true',  dest='mpi')
group.add_option('--disable-mpi',       action='store_false', dest='mpi',       help='enable/disable MPI parallelization.')
group.add_option('--enable-scalapack',  action='store_true',  dest='scalapack')
group.add_option('--disable-scalapack', action='store_false', dest='scalapack', help='disable/disable computations with ScaLAPACK library.')
parser.add_option_group(group)

group = OptionGroup(parser, 'Optimization options')
group.add_option('--enable-reduce-for-manycore',  action='store_true',  dest='reduce_manycore')
group.add_option('--disable-reduce-for-manycore', action='store_false',  dest='reduce_manycore', help='enable/disable reduction code optimization for many-core processor.')
parser.add_option_group(group)

group = OptionGroup(parser, 'Optimization options for stencil computations')
group.add_option('--old-stencil',           action='store_false', dest='stencil_optimized', help='use old implementations.')
group.add_option('--explicit-vec',          action='store_true',  dest='explicit_vec',      help='enable explicit vectorization. it requires --simd-set option to be set.')
group.add_option('--compiler-vec',          action='store_false', dest='explicit_vec',      help='entrust optimization to a compiler.')
group.add_option('--simd-set',              action='store',       dest='simd',              help='specifies SIMD instruction set. (e.g. AVX, AVX_512, HPC_ACE2...)')
group.add_option('--enable-swp',            action='store_true',  dest='swp')
group.add_option('--disable-swp',           action='store_false', dest='swp',               help='enable/disable software prefetch in the explicit vec.')
group.add_option('--enable-array-padding',  action='store_true',  dest='padding')
group.add_option('--disable-array-padding', action='store_false', dest='padding',           help='enable/disable array padding for the cache utilization.')
group.add_option('--enable-domain-pow2',    action='store_true',  dest='domain_two')
group.add_option('--disable-domain-pow2',   action='store_false', dest='domain_two',        help='enable/disable whether the optimization assumes that 3-D domain size is power of two.')
group.add_option('--enable-loop-blocking',  action='store_true',  dest='loop_blocking')
group.add_option('--disable-loop-blocking', action='store_false', dest='loop_blocking',     help='enable/disable loop blocking.')
parser.add_option_group(group)

group = OptionGroup(parser, 'Debug options')
group.add_option('-d', '--debug', action='store_true', default=False, dest='debug', help='enable debug build.')
group.add_option('--nvtx',        action='store_true',                dest='nvtx',  help='use NVIDIA Tools Extention Library.')
group.add_option('--hpsi_test',   action='store_true',                dest='hpsi_test',  help='use joint hpsi subroutine (test).')
parser.add_option_group(group)

(options, args) = parser.parse_args()

### check options
dict = {}
if options.arch is not None:
  dict['CMAKE_TOOLCHAIN_FILE']     = options.arch.lower()
if options.prefix is not None:
  dict['CMAKE_INSTALL_PREFIX']     = options.prefix
dict['CMAKE_BUILD_TYPE']           = debug_or_release(options.debug)
dict['CMAKE_VERBOSE_MAKEFILE']     = on_or_off(options.verbose)

add_option(dict, 'USE_MPI',             options.mpi)
add_option(dict, 'USE_SCALAPACK',       options.scalapack)

add_option(dict, 'OPT_STENCIL',         options.stencil_optimized)
add_option(dict, 'DOMAIN_IS_POW2',      options.domain_two)
add_option(dict, 'ARRAY_PADDING',       options.padding)
add_option(dict, 'EXPLICIT_VEC',        options.explicit_vec)
add_option(dict, 'LOOP_BLOCKING',       options.loop_blocking)
add_option(dict, 'SW_PREFETCH',         options.swp)
add_option(dict, 'REDUCE_FOR_MANYCORE', options.reduce_manycore)

add_option(dict, 'USE_NVTX',            options.nvtx)
add_option(dict, 'HPSI_TEST',           options.hpsi_test)
if options.simd is not None:
  dict['SIMD_SET'] = options.simd.upper()

define = ''
for k,v in dict.items():
  define = '{0} -D {1}={2}'.format(define, k, v)

env = ''
for var in args:
  (k, v) = var.split('=', 1)
  env = '{0} {1}="{2}"'.format(env, k, v)

### configuration
comm = '{2} cmake {0} {1}'.format(define, SOURCE_DIR, env)
print('    $ %s' % comm)
os.system(comm)
