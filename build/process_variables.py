import sys
import os.path
import json
import re 

varinfo_json = open('../share/varinfo.json','r')
variables = json.load(varinfo_json)
varinfo_json.close()

# Define some functions, we will use later:

def cleanhtml(raw_html):
    """This function removes HTML tags from the string."""
    cleanr = re.compile('<.*?>')
    cleantext = re.sub(cleanr, '', raw_html)
    return cleantext

def leavehtml(raw_html):
    """This is a dummy function which leaves the HTML tags in the string."""
    return raw_html

def extract_arg(string):
    """Extract the number x from a 'bit(x)' string."""
    try:
        return int(string.replace('bit(','').replace(')',''))
    except TypeError:
        print(string+' was not of the form bit(X)')

def bit(x):
    """Calculate 2^x if x is an integer."""
    if type(x) is int:
        return 2**x
    else: 
        print('bit argument must be integer')
        raise TypeError

def is_number(s):
    """Returns True is string is a number. """
    return s.replace('.','',1).isdigit()

def print_varinfo(variables, file=sys.stdout, filter=leavehtml):
    """Print the varinfo files with optional filtering of the HTML tags.

       Arguments:
       - file: file object into which the variable info will be written
       - filter: a function object: the default leavehtml is a dummy, which does not change the string.
                 in order to filter the HTML tags, pass file=cleanhtml.
    """

    for var in variables.keys():
    
        variable = variables[var]
        print('Variable '+variable['Name'], file=file)
        print('Type '+variable['Type'], file=file)
        if( variable['Default']):
            print('Default '+' '.join(variable['Default']), file=file)
        print('Section '+variable['Section'], file=file)
        print('Description', file=file)
        for desc in variable['Description']:
            print(' '+ filter(desc), file=file)
        for option in variable['Options']:
            value = option.get('Value')
            print('Option '+option['Name'] + ' ', value, file=file )
            # print( type( option.get('Description')))
            if 'Description' in option:
                print(' ', filter(option.get('Description')) , file=file)
        print('END\n', file=file)
    

def print_defaults_header(variables, file):
    """Generate the defaults.h file."""

    for var in sorted(variables.keys()):
        current = variables[var]
        if 'Default' in current:
            if current['Default'] and current['Type'] in ('integer', 'real','complex'):
                value = current['Default'][0]
                for option in current['Options']:
                    if (value.lower() in option['Name'].lower()) and ('Value' in option):
                        value = option['Value']
                if is_number(value):
                    value = value+'_8'
                print('#define DEFAULT__'+var.upper()+' ('+value+')', file=file)


def print_options_header(variables, file):
    for var in sorted(variables.keys()):
        current = variables[var]
        for option in sorted(current['Options'],key=lambda opt: opt['Name']):
            if 'Value' in option:
                value = option['Value']
                if 'bit(' in value:
                    value = str(bit(extract_arg(value)))
                print('#define OPTION__'+var.upper()+'__'+option['Name'].upper()+' ('+value+'_8)', file=file)


def print_variables(variables, file):
    """Generate the variables file."""
    
    option_values = dict()
    for var in variables.keys():
        for option in variables[var]['Options']:
            if 'Value' in option:
                value = option['Value']
                if 'bit(' in value:
                    value = str(bit(extract_arg(value)))
                option_values[option['Name']] = value
    
    for option in sorted(option_values.keys()):
        print(option+' = '+option_values[option], file=file)
    


varinfo = open('../share/new_varinfo','w')
varinfo_ORIG = open('../share/new_varinfo_ORIG','w')

file_variables = open('../share/new_variables','w')
file_defaults_header = open('defaults.h','w')
file_options_header = open('options.h','w')


print_varinfo(variables, file=varinfo_ORIG, filter=leavehtml)
print_varinfo(variables, file=varinfo, filter=cleanhtml)

print_variables(variables, file_variables)
print_defaults_header(variables, file_defaults_header)
print_options_header(variables, file_options_header)


varinfo.close()
varinfo_ORIG.close()

file_variables.close()
file_defaults_header.close()
file_options_header.close()




