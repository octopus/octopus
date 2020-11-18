#
#

import sys
import glob
import json
import re


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


class Variables:
    """
    The Variables class is a wrapper around the dictionary, which contains the full info about all Octopus variables.

    """

    def __init__(self, verbose=False, sources=None, tests=None, json=None):
        """
        Constructor of the Variables class.

        Without arguments, an empty dictionary is created. The optional argument verbose controls whether the names of parsed files is echoed to stdout.
        
        When the optional arguments sources or tests are given, the corresponding files are parsed. 
        The arguments should be the directories, in which the sources or tests can be found.
        Alternatively, the dictionary can be imported from a JSON file (no validation is performed).
        """

        if json and sources:
            print('Variables cannot be initialized from source files and a json file at the same time.')
            raise RuntimeError

        self.variables = dict()
        self.verbose = verbose
    
        if sources:
            self.parse_sources(sources)

        if tests:
            self.parse_testfiles(tests)

        if json:
            self.import_json(variable_defs_filename=json)


    def __getitem__(self, key):
        """
        Implement the [] operator for the vatriable class.
        """

        try:
            return self.variables[key]
        except KeyError:
            print('Warning: Key '+key+' not found.', file=stderr)
            return None


    def keys(self):
        """
        Return the list of keys.
        """

        return self.variables.keys()


    def length(self):
        """
        Return the number of elements in the dictionary.
        """

        return len(self.variables)


    def parse_sources(self, srcdir):
        """
        Parse all file in srcdir/*/*.F90 and fill in dictionary.
        """
    
        srcfiles = glob.glob(srcdir+"*/*.F90")
    
        for f in srcfiles:
            if self.verbose:
                print("Parsing "+f)
        
            with open(f,'r') as source:

                parse_calls = []
                keys_in_file = []

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
                            keys_in_file.append(var_name.lower())
                            self.variables[var_name.lower()] = {
                                'Name': var_name,
                                'Type':var_type, 
                                'Default':var_default, 
                                'Section':var_section, 
                                'Options':var_options, 
                                'Description':var_description, 
                                'Testfiles':[],
                                'Sourcefile':f.replace(srcdir,'')
                                }

            # associate calls with the variable description:
            #
            #    (we seperated this to relax the requirement that the description is before the call)

            for key in keys_in_file:
            
                for call in parse_calls:
                    if key in call[0].lower():
                        self.variables[key]['CallLine'] = call[0]
                        self.variables[key]['LineNumber'] = call[1]
                
        

            source.close()
        
    
    def parse_testfiles(self, testdir):

        files = glob.glob(testdir+"*/*.inp")

        for f in files:
            if self.verbose:
                print('Parsing '+f+':')
            input = open(f,'r')
            for line in input:
                words = line.split()
                if len(words)>0:
                    key =  words[0].strip().lower()
                    if key in self.variables.keys():
                        self.variables[key]['Testfiles'].append(f.replace(testdir,''))
            input.close()


    def export_json(self, variable_defs_filename):

        varinfo_json = open(variable_defs_filename,'w')
        print(json.dumps(self.variables, sort_keys=True, indent=2), file=varinfo_json)
        varinfo_json.close()


    def import_json(self, variable_defs_filename):

        varinfo_json = open(variable_defs_filename,'w')
        self.variables = json.load(varinfo_json)
        varinfo_json.close()


    def write_varinfo(self, filename='-', filterHTML=True):
        """Write the varinfo files with optional filtering of the HTML tags.
    
           Arguments:
           - file: file object into which the variable info will be written
           - filter: a function object: the default leavehtml is a dummy, which does not change the string.
                     in order to filter the HTML tags, pass file=cleanhtml.
        """

        if filterHTML:
            filter = cleanhtml
        else:
            filter = leavehtml

        if filename=='-':
            file = sys.stdout
        else:
            file = open(filename, 'w')

        for key in self.variables.keys():
        
            variable = self.variables[key]
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
    
        if file != sys.stdout:
            file.close()


    def write_defaults_header(self, filename):
        """Generate the defaults.h file."""
    
        if filename == '-':
            file = sys.stdout
        else:
            file = open(filename, 'w')

        for key in sorted(self.variables.keys()):
            current = self.variables[key]
            if 'Default' in current:
                if current['Default'] and current['Type'] in ('integer', 'real','complex'):
                    value = current['Default'][0]
                    for option in current['Options']:
                        if (value.lower() in option['Name'].lower()) and ('Value' in option):
                            value = option['Value']
                    if is_number(value):
                        value = value+'_8'
                    print('#define DEFAULT__'+key.upper()+' ('+value+')', file=file)

        if file != sys.stdout:
            file.close()    


    def write_options_header(self, filename):
        """Generate the options.h file."""

        file = open(filename, 'w')

        for key in sorted(self.variables.keys()):
            current = self.variables[key]
            for option in sorted(current['Options'],key=lambda opt: opt['Name']):
                if 'Value' in option:
                    value = option['Value']
                    if 'bit(' in value:
                        value = str(bit(extract_arg(value)))
                    print('#define OPTION__'+key.upper()+'__'+option['Name'].upper()+' ('+value+'_8)', file=file)

        file.close()


    def write_variables_file(self, filename):
        """Generate the variables file."""
        
        file = open(filename,'w')
    
        option_values = dict()
        for key in self.variables.keys():
            for option in self.variables[key]['Options']:
                if 'Value' in option:
                    value = option['Value']
                    if 'bit(' in value:
                        value = str(bit(extract_arg(value)))
                    option_values[option['Name']] = value
        
        for option in sorted(option_values.keys()):
            print(option+' = '+option_values[option], file=file)
    
        file.close()



    def number_of_untested(self, without_default=True, filename=None):

        untested = 0

        if filename:
            if filename=='-':
                file=sys.stdout
            else:
                file=open(filename, 'w')


        for key in self.variables.keys():
            current = self.variables[key] 
            if len( current['Testfiles'] ) is 0:
                if without_default:
                    if not current['Default']:
                        untested += 1
                        if filename:
                            print(current['Name'],file=file)
                else:
                    untested += 1
                    if filename:
                        print(current['Name'],file=file)

        if filename:
            if filename != '-':
                file.close()

        return untested