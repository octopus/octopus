#!/usr/bin/env python3

import sys
import getopt
import json
import re


def annotate(str, mode='OctopusWiki'):
    """
    Return annotated string for variable name.

    Arguments:
    str:  string containing a word of the input file
    mode: annotation mode. Currently implemented:
    - 'OctopusWiki': '{{variable|key|section}}'
    - 'Test': 'Function(key, section)'
    """

    key = str.strip()

    if mode == 'OctopusWiki':
        section = variables[key.lower()]['Section'].split('::')[0]
        return '{{variable|'+key+'|'+section+'}}'

    # Here we can add other formats, e.g. for MkDocs

    if mode == 'Test':
        section = variables[key.lower()]['Section']
        return 'Function('+key+','+section+')'    



implemented_modes = ['OctopusWiki', 'Test']

# Some default names:
 
input_file_name = '-'
output_file_name = '-'
variable_defs_name = '../share/varinfo.json'
mode = 'OctopusWiki'


# Parse command line options:

try:
    options, args = getopt.getopt(sys.argv[1:], "i:o:v:m:")
except  getopt.GetoptError:
    print(' annotate_input -i <inputfile> -o <outputfile> -v <variables-definition>')
    sys.exit(-1)

for (opt, arg) in options:
    if opt == '-i':
        input_file_name = arg
    elif opt == '-o':
        output_file_name = arg
    elif opt in ['-d', 'definitions']:
        variable_defs_name = arg
    elif opt == '-m':
        if arg in implemented_modes:
            mode = arg
        else:
            print('Unsopported mode: \''+arg+'\'. Try one of: '+', '.join(implemented_modes))
            sys.exit(-1)
        
if input_file_name is '-':
    if len(args)>0:
        input_file_name = args[0]


# Read Variable definitions from JSON:

varinfo_json = open(variable_defs_name,'r')
variables = json.load(varinfo_json)
varinfo_json.close()

# Open input and output files:

if input_file_name is '-':
    input = sys.stdin
else:
    try:
        input = open(input_file_name,'r')
    except FileNotFoundError:
        print(input_file_name+' not found.')
        sys.exit(-1)

if output_file_name is '-':
    output = sys.stdout
else:
    output = open(output_file_name, 'w')


# Process input file:


for line in input:
    # Extract indentation:
    lead = ''
    if len(line.strip())>1:
        lead = line[0:re.search('\S', line).start()]
        
    words = line.split()
    for i in range(0,len(words)):
        # Remove prefix '%' and save for later:
        if words[i].strip()[0] =='%':
            prefix = '%'
            words[i] = words[i].replace('%','')
        else:
            prefix=''
        # Annotate, if word is defined variable:
        if words[i].lower() in variables.keys():
            words[i] = annotate(words[i], mode=mode)
        # Re-attach prefix:
        words[i] = prefix+words[i]
    # Joins words and indent by one space (for tutorial formatting) 
    new_line = ' ' + lead + ' '.join(words)
    print(new_line, file=output)


input.close()
output.close()
