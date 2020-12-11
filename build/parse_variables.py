#!/usr/bin/env python3


# Copyright (C) 2020 Martin Lueders 
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA.
#



import sys
import os.path
import glob
import json
import getopt
from variables import Variables


def usage():
    """print usage information."""

    print('Usage: parse_variables [options]')
    print('')
    print('Opdions:')
    print('  -h, --help              Prints this help and exits.')
    print('  -v, --versbose          Print names of parsed files.')
    print('  -s, --srcdir=DIR        Name of the top source directory')
    print('  -b, --builddir=DIR      Name of the top build directory')
    print('  -d, --definitions=NAME  Name of the JSON definitions file')
    print('  --enable-varinfo, --disable-varinfo: enable/disable the output of "varinfo" and "varinfo_orig". Default: enabled' )
    print('  --enable-headers, --disable-headers: enable/disable the output of the header files. Default: enabled.')
    print('  --enable-variables, --disable-variables: enable/disable the output of the "variables" file. Default: enabled.')
    print('  --enable-json, --disable-json: enable/disable the generation of the JSON definitions file. Default: enabled')
    print('  --enable-untested, --disable-untested: enable/disable the generation of the "untested.txt" file. Default: disabled')


# set up defaults:

top_srcdir = "../"
top_builddir = "./"
variable_defs_name = ''

write_varinfo = True
write_headers = True
write_variables_file = True
write_json = True
write_untested = False

verbose = False

# Parse command line options:

# change to resemble mk_varinfo.pl behaviour:
# -s top_srcdir
# -b top_builddir

try:
    options, args = getopt.getopt(sys.argv[1:], "s:b:d:vh", ['srcdir:','builddir:','definitions:','verbose','help',
        'enable-varinfo', 'enable-headers', 'enable-variables', 'enable-json', 'enable-untested',
        'disable-varinfo','disable-headers','disable-variables','disable-json', 'disable-untested'
        ])
except  getopt.GetoptError:
    print('Error !')
    usage()
    sys.exit(-1)

for (opt, arg) in options:
    if opt in ['-h', '--help']:
        usage()
        sys.exit(0)
    if opt in ['-s', '--top_srcdir']:
        top_srcdir = arg
    if opt in ['-b', '--top_builddir']:
        top_builddir = arg
    elif opt in ['-d','--definitions']:
        variable_defs_name = arg
    elif opt in ['-v','--verbose']:
        verbose = True
    elif opt in ['--enable-varinfo']:
        write_varinfo = True
    elif opt in ['--enable-headers']:
        write_headers = True 
    elif opt in ['--enable-variables']:
        write_variables = True
    elif opt in ['--enable-json']:
        write_json = True
    elif opt in ['--enable-untested']:
        write_untested = True
    elif opt in ['--disable-varinfo']:
        write_varinfo = False
    elif opt in ['--disable-headers']:
        write_headers = False
    elif opt in ['--disable-variables']:
        write_variables = False
    elif opt in ['--disable-json']:
        write_json = False
    elif opt in ['--disable-untested']:
        write_untested = False
        


testdir  = top_srcdir+"/testsuite/"
srcdir   = top_srcdir+"/src/"
sharedir = top_builddir+"/share/"
includedir = top_builddir+'/src/include/'

if variable_defs_name == '':
    variable_defs_name = sharedir+'/varinfo.json'



variables = Variables(sources=srcdir, tests=testdir)
print('Parsing sources complete.')

if write_json:
    variables.export_json(variable_defs_name)

if write_varinfo:
    print('Generating varinfo files.')
    variables.write_varinfo(filename=sharedir+'varinfo_orig', filterHTML=False)
    variables.write_varinfo(filename=sharedir+'varinfo', filterHTML=True)

if write_headers:
    print('Generating header files.')
    variables.write_defaults_header(filename=includedir+'defaults.h')
    variables.write_options_header(filename=includedir+'options.h')

if write_variables_file:
    print('Generating variables file.')
    variables.write_variables_file(filename=sharedir+'variables')

if write_untested:
    untested_file = 'untested.txt'
else:
    untested_file = None


print('\n')
print('Number of variables:                                         ', variables.length())
print('Number of variables: not explicitely referenced in test files', variables.number_of_untested(without_default=False, filename=untested_file))
print('Number of variables: not referenced in test files            ', variables.number_of_untested())