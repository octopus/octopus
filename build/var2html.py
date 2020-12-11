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


variables_top_dir = 'Variables'



try:
    options, args = getopt.getopt(sys.argv[1:], "i:o:d:")
except  getopt.GetoptError:
    usage()
    sys.exit(-1)

input_file_name = '-'

for (opt, arg) in options:
    if opt == '-i':
        input_file_name = arg
    elif opt == '-o':
        output_file_name = arg
    elif opt in ['-d', 'definitions']:
        variable_defs_name = arg
        
if input_file_name is '-':
    if len(args)>0:
        input_file_name = args[0]


# Read Variable definitions from JSON:

variables = Variables()
variables.import_json(variable_defs_name)


sections = dict()

for var in variables.keys():

    variable = variables[var]

    full_section = variable['Section'].split('::')
    section = full_section[0]
    if len(full_section) > 1:
        subsection = full_section[1]
    else:
        subsection = ''

    if section not in sections.keys():
        sections[section] = dict()

    if subsection in sections[section].keys():
        sections[section][subsection].add(var)
    else:
        sections[section][subsection] = {var}



if not os.path.isdir(variables_top_dir):
    os.mkdir(variables_top_dir)

top_index = open(variables_top_dir+'/_index.md','w')

print("""
---
Title: Variables Overview
---

{{% children depth=5 %}}

""", file=top_index)

file = sys.stdout


for section in sorted(sections):


    section_dir = variables_top_dir+'/'+section.replace(' ','_')
    if not os.path.isdir(section_dir):
        os.mkdir(section_dir)
    
    section_index = open(section_dir+'/_index.md', 'w')

    print('---', file=section_index)
    print('Title: '+section,file=section_index)
    print('---\n', file=section_index)
    print('{{% children depth=5 %}}', file=section_index)


    for sub in sections[section]:
        if sub:

            subsection_dir = section_dir+'/'+sub.replace(' ','_')
            if not os.path.isdir(subsection_dir):
                os.mkdir(subsection_dir)

            subsection_index = open(subsection_dir+'/_index.md', 'w')

            print('---', file=subsection_index)
            print('Title: '+sub,file=subsection_index)
            print('---\n', file=subsection_index)
            print('{{% children depth=5 %}}', file=subsection_index)

        else:

            subsection_dir = section_dir


        for var in sections[section][sub]:

            variables_file = open(subsection_dir+'/'+ var+'.md','w')
            print('---', file=variables_file)
            print('Title: '+variables[var]['Name'],file=variables_file)
            if sub:
                print('tags: ["Variables", "'+section+'", "'+sub+'"]', file=variables_file)
            else:
                print('tags: ["Variables", "'+section+'"]', file=variables_file)
            print('---\n', file=variables_file)


            variables.write_variable_md(var, file=variables_file)
            print('<hr>\n', file=variables_file)

            print('{{%expand "Source information"%}}',file=variables_file)
            sourcefile = variables[var]['Sourcefile']
            if 'LineNumber' in variables[var]:
                linenumber = ': '+str(variables[var]['LineNumber'])
            else:
                linenumber = ''
            print('<a href="https://gitlab.com/octopus-code/octopus/-/blob/10.3/src/'+sourcefile+'">'+sourcefile+'</a> '+linenumber  ,file=variables_file)
            if 'CallLine' in variables[var]:
                print('```Fortran',file=variables_file)
                print(' '+variables[var]['CallLine'],file=variables_file)
                print('```',file=variables_file)

            print('{{%/expand%}}\n',file=variables_file)


            if 'Testfiles' in variables[var]:
                if len(variables[var]['Testfiles']) > 0:
                    print('{{%expand "Featured in Testfiles"%}}',file=variables_file)
                    print('<ul>',file=variables_file)
                    for testfile in variables[var]['Testfiles']:
                        print('<li><a href="https://gitlab.com/octopus-code/octopus/-/blob/10.3/testsuite/'+testfile+'">'+testfile+'</a></li>',file=variables_file)
                    print('</ul>',file=variables_file)
                    print('{{%/expand%}}',file=variables_file)
            variables_file.close()

        if sub:
            subsection_index.close()

    section_index.close()

top_index.close()

