#!/usr/bin/env python3

import sys
import os.path
import glob
import json
import getopt
from variables import Variables



# set up directories and file names:

top_srcdir = "../"
top_builddir = "../"
verbose = False
variable_defs_name = ''

# Parse command line options:

# change to resemble mk_varinfo.pl behaviour:
# -s top_srcdir
# -b top_builddir

try:
    options, args = getopt.getopt(sys.argv[1:], "s:b:d:v", ['srcdir','builddir','definitions','verbose' ])
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
        


testdir  = top_srcdir+"testsuite/"
srcdir   = top_srcdir+"src/"
sharedir = top_builddir+"share/"
includedir = top_builddir+'src/include/'

if variable_defs_name == '':
    variable_defs_name = sharedir+'varinfo.json'



variables = Variables(sources=srcdir, tests=testdir)

variables.export_json(variable_defs_name)

variables.print_varinfo(filename=sharedir+'new_varinfo_ORIG', filterHTML=False)
variables.print_varinfo(filename=sharedir+'new_varinfo', filterHTML=True)

variables.print_defaults_header(filename='defaults.h')
variables.print_options_header(filename='options.h')

variables.print_variables(filename='variables')



print('\n')
print('Number of variables:                                         ', variables.length())
print('Number of variables: not explicitely referenced in test files', variables.number_of_untested(without_default=False, filename='untested.txt'))
print('Number of variables: not referenced in test files            ', variables.number_of_untested())


# Finally dump the data into the JSON file.

