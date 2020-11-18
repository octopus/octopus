#!/usr/bin/env python3

import sys
import os.path
import glob
import json
import getopt
from variables import Variables



# set up defaults:

top_srcdir = "../"
top_builddir = "../"
variable_defs_name = ''

write_varinfo = True
write_headers = True
write_variables_file = True
verbose = False

# Parse command line options:

# change to resemble mk_varinfo.pl behaviour:
# -s top_srcdir
# -b top_builddir

try:
    options, args = getopt.getopt(sys.argv[1:], "s:b:d:v", ['srcdir:','builddir:','definitions:','verbose',
        'enable-varinfo', 'enable-heders', 'enable-variables', 'enable-json',
        'disable-varinfo','disable-heders','disable-variables','disable-json',
        ])
except  getopt.GetoptError:
    sys.exit(-1)

for (opt, arg) in options:
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
    elif opt in ['--disable-varinfo']:
        write_varinfo = False
    elif opt in ['--disable-headers']:
        write_headers = False
    elif opt in ['--disable-variables']:
        write_variables = False
    elif opt in ['--disable-json']:
        write_json = False
        


testdir  = top_srcdir+"testsuite/"
srcdir   = top_srcdir+"src/"
sharedir = top_builddir+"share/"
includedir = top_builddir+'src/include/'

if variable_defs_name == '':
    variable_defs_name = sharedir+'varinfo.json'



variables = Variables(sources=srcdir, tests=testdir)

variables.export_json(variable_defs_name)

if write_varinfo:
    variables.write_varinfo(filename=sharedir+'varinfo_ORIG', filterHTML=False)
    variables.write_varinfo(filename=sharedir+'varinfo', filterHTML=True)

if write_headers:
    variables.write_defaults_header(filename=includedir+'defaults.h')
    variables.write_options_header(filename=includedir+'options.h')

if write_variables_file:
    variables.write_variables_file(filename=sharedir+'variables')



print('\n')
print('Number of variables:                                         ', variables.length())
print('Number of variables: not explicitely referenced in test files', variables.number_of_untested(without_default=False, filename='untested.txt'))
print('Number of variables: not referenced in test files            ', variables.number_of_untested())


# Finally dump the data into the JSON file.

