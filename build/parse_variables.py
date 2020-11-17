#!/usr/bin/env python3

import sys
import os.path
import glob
import json
import getopt



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

if variable_defs_name == '':
    variable_defs_name = sharedir+'varinfo.json'

varinfo_json = open(variable_defs_name,'w')

variables = dict()
variable_list = list()

# Parse source files:

srcfiles = glob.glob(srcdir+"*/*.F90")

parse_calls = []

for f in srcfiles:
    if verbose:
        print("Parsing "+f)

    with open(f,'r') as source:

        # first parse calls in the source:

        line = source.readline()
        line_number = 1
        while line != '':
            if 'call parse_variable(' in line:
                parse_line= line.strip()
                if parse_line[-1] is '&':
                    line = source.readline()
                    line_number += 1
                    parse_line = parse_line[:-1] + line.strip()
                parse_calls.append( (parse_line, line_number) )
            line = source.readline()
            line_number += 1

        source.seek(0)

        # parse variable descriptions in the source:

        for line in source:
            if "!%" in line:
                words = line.split()

                if "!%variable" in line.lower():
                    parsing_mode = 0
                    var_name = words[1]
                    var_options = []
                    var_default = None

                if "!%type" in line.lower():
                    parsing_mode = 0
                    if len(words)>1:
                        var_type = words[1]
                    else:
                        var_type = ""

                if "!%default" in line.lower():
                    parsing_mode = 0
                    if len(words)>1:
                        var_default = words[1:]
                    else:
                        var_default=""

                if "!%section" in line.lower():
                    parsing_mode = 0
                    if len(words)>1:
                        var_section = words[1]
                    else:
                        var_section = ""

                if "!%option" in line.lower():
                    if len(words)>1:
                        tmp = dict()
                        tmp['Name'] = words[1]
                        parsing_mode = 1
                        if len(words) > 2:
                            tmp['Value'] = words[2]
                        var_options.append(tmp)

                if "!%description" in line.lower():
                    parsing_mode = 2
                    if(len(words)>1):
                        var_description = [' '.join(words[1:])]
                    else:
                        var_description = []

                if "!% " in line.lower():
                    if parsing_mode is 1:
                        var_options[-1]['Description'] = line.replace('!%','').strip()
                    if parsing_mode is 2:
                        var_description.append(line.replace('!%','').strip())

                if "!%end" in line.lower():
                    parsing_mode = 0
                    variables[var_name.lower()] = {
                        'Name': var_name,
                        'Type':var_type, 
                        'Default':var_default, 
                        'Section':var_section, 
                        'Options':var_options, 
                        'Description':var_description, 
                        'Testfiles':[],
                        'Sourcefile':f.replace(srcdir,'')
                        }

    source.close()

# associate calls with the variable description:
#
#    (we seperated this to relax the requirement that the description is before the call)

for var in variables.keys():

    for call in parse_calls:
        if var in call[0]:
            variables[var]['CallLine'] = call[0]
            variables[var]['LineNumber'] = call[1]
        


files = glob.glob(testdir+"*/*.inp")

# Parse test input files:
for f in files:
    if verbose:
        print('Parsing '+f+':')
    input = open(f,'r')
    for line in input:
        words = line.split()
        if len(words)>0:
            var =  words[0].strip().lower()
            if var in variables.keys():
                variables[var]['Testfiles'].append(f.replace(testdir,''))
    input.close()


untested = 0

for v in variables.keys():
    if len( variables[v]['Testfiles'] ) is 0:
        untested += 1

print('\n')
print('Number of variables:                                         ', len(variables))
print('Number of variables: not explicitely referenced in test files', untested)


# Finally dump the data into the JSON file.

print(json.dumps(variables, sort_keys=True, indent=2), file=varinfo_json)
varinfo_json.close()

