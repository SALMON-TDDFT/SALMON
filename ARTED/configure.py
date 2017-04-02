#! /usr/bin/env python
#
#   Copyright 2016 ARTED developers
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
parser.add_option_group(group)

group = OptionGroup(parser, 'Optimization options')
group.add_option('--old-stencil',         action='store_false', dest='stencil_optimized', help='use old implementation of the stencil computation code.')
group.add_option('--explicit-vec',        action='store_true',  dest='explicit_vec',      help='enable explicit vectorization in the stencil computation with C-lang.')
group.add_option('--compiler-vec',        action='store_false', dest='explicit_vec',      help='defer to compiler vectorization in the stencil computation with Fortran90.')
group.add_option('--simd-set',            action='store',       dest='simd',              help='specifies SIMD instruction set. (e.g. AVX, HPC_ACE2...)')
group.add_option('--enable-swp',          action='store_true',  dest='swp',               help='enable software prefetch in the explicit vec.')
group.add_option('--disable-swp',         action='store_false', dest='swp',               help='disable software prefetch in the explicit vec.')
group.add_option('--array-padding',       action='store_true',  dest='padding',           help='array padding applied to the stencil computation.')
group.add_option('--domain-pow2',         action='store_true',  dest='domain_two',        help='3-D domain size is power of two.')
group.add_option('--loop-blocking',       action='store_true',  dest='loop_blocking',     help='loop blocking applied to the stencil computation.')
group.add_option('--opt-current',         action='store_true',  dest='current_optimized', help='enable the current routine optimization in RT.')
group.add_option('--reduce-manycore',     action='store_true',  dest='reduce_manycore',   help='enable reduction code optimization for many-core processor.')
parser.add_option_group(group)

group = OptionGroup(parser, 'Debug options')
group.add_option('-d', '--debug', action='store_true', default=False, dest='debug', help='enable debug build.')
group.add_option('--papi',        action='store_true',                dest='papi',  help='use PAPI profiling (SC only).')
group.add_option('--nvtx',        action='store_true',                dest='nvtx',  help='use NVIDIA Tools Extention Library.')
parser.add_option_group(group)

(options, args) = parser.parse_args()

### check options
dict = {}
if options.arch is not None:
  dict['CMAKE_TOOLCHAIN_FILE']     = options.arch.lower()
dict['CMAKE_BUILD_TYPE']           = debug_or_release(options.debug)
dict['CMAKE_VERBOSE_MAKEFILE']     = on_or_off(options.verbose)

add_option(dict, 'USE_PAPI',                   options.papi)
add_option(dict, 'USE_NVTX',                   options.nvtx)
add_option(dict, 'OPT_STENCIL',                options.stencil_optimized)
add_option(dict, 'OPT_CURRENT',                options.current_optimized)
add_option(dict, 'DOMAIN_IS_POW2',             options.domain_two)
add_option(dict, 'ENABLE_ARRAY_PADDING',       options.padding)
add_option(dict, 'ENABLE_EXPLICIT_VEC',        options.explicit_vec)
add_option(dict, 'ENABLE_LOOP_BLOCKING',       options.loop_blocking)
add_option(dict, 'ENABLE_SWPREFETCH',          options.swp)
add_option(dict, 'ENABLE_REDUCE_FOR_MANYCORE', options.reduce_manycore)
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
print '    $', comm
os.system(comm)
