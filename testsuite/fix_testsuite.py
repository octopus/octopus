#!/usr/bin/env python

##########################################################################
# BerkeleyGW, Copyright (c) 2012, The Regents of the University of
# California, through Lawrence Berkeley National Laboratory (subject to
# receipt of any required approvals from the U.S. Dept. of Energy).
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
# 
# (1) Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
# 
# (2) Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
# 
# (3) Neither the name of the University of California, Lawrence
# Berkeley National Laboratory, U.S. Dept. of Energy nor the names of
# its contributors may be used to endorse or promote products derived
# from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# You are under no obligation whatsoever to provide any bug fixes,
# patches, or upgrades to the features, functionality or performance of
# the source code ("Enhancements") to anyone; however, if you choose to
# make your Enhancements available either publicly, or directly to
# Lawrence Berkeley National Laboratory, without imposing a separate
# written license agreement for such Enhancements, then you hereby grant
# the following license: a  non-exclusive, royalty-free perpetual
# license to install, use, modify, prepare derivative works, incorporate
# into other computer software, distribute, and sublicense such
# enhancements or derivative works thereof, in binary and source code
# form.
##########################################################################

# Felipe Homrich da Jornada, September 2011
# imported from BerkeleyGW r4576, r5630

# Python script for updating a test file given its result.
# The script searches for all *.test files recursively, and updates those that
# were affected by the result.
# It can also help one develop a new test file: just create the test skeleton
# and the fields you would like to match, and enter zero for the matched
# values. The script will update these numbers accordingly.

# USAGE: go to the testsuite directory and type ./fix_testsuite.py test_result
#   where test_result contains the results from the test file.

# WARNING: this script updates all test results. Don't blindly accept all
#   changes! You might have meant to change only one component (A), but you
#   might have introduced a bug in component B! So, take a look at all numbers
#   before you "fix" the test result!

import sys
import os
import fnmatch
import re

if len(sys.argv)!=2:
	sys.exit('Usage: %s test_result'%(sys.argv[0]))

f_result = open(sys.argv[1])
f_in  = None
f_out = None

#Return the next line containing a "match" string from the test file.
#All lines before the match are copied to the .test.new file
def get_next_match(match_name):
	global f_in, f_out

	while 1:
		line = f_in.readline()
		if not line: break
                match = re.search('match(\s)*;(\s)*%s'%re.escape(match_name), line)
                if match:
			return line
		f_out.write(line)
	raise Exception('Read past end of file, could not find match!')

#Update the current .test and .test.new files. Finishes copying the rest of the
#.test file into .test.new, and closes/reopens new files.
def refresh_files(new_in=None, new_out=None):
	global f_in, f_out, n_cor, n_fine

	#Finishes copying the .test => .test.new
	if not f_in is None:
		f_out.write(f_in.read())
		print('     corrections: %d/%d'%(n_cor, n_cor+n_fine))

	if not(f_in is None):
		f_in.close()
	if not(f_out is None):
		f_out.close()
	if not ((new_in is None) or (new_out is None)):
		f_in  = open(new_in)
		f_out = open(new_out, 'w')


#Mimics the "find" utility
def recursive_glob(rootdir='.', pattern='*'):
	return [os.path.join(rootdir, filename)
		for rootdir, dirnames, filenames in os.walk(rootdir)
		for filename in filenames
		if fnmatch.fnmatch(filename, pattern)]

#Get test name given a test file
re_obj = re.compile(r'Test\s*:\s*(.*)')
def get_test_name(fname):
	with open(fname) as f:
		str = f.read()
		try:
			test_name = re_obj.search(str).group(1)
		except:
			test_name = ''
		return '***** %s *****'%(test_name.strip())

#Find all tests and gets all titles
test_files = recursive_glob('.', '*.test')
test_names = [get_test_name(fname) for fname in test_files]
print('Found %d test files:'%(len(test_files)))
for test_file in test_files:
	print('      %s'%(test_file))
print('')

ignored_tests=[]

ntot_fine=0
ntot_cor=0
n_fine=0
n_cor=0
test_exists=False
for line in f_result.readlines():
	if line[0:18] == 'Using test file  :':
		# There is a new test! Look for the file...
		cur_file = line[19:-1].strip()
		test_exists = False
		if (os.path.isfile(cur_file)):
			for tname,tfile in zip(test_names, test_files):
				if os.path.samefile(cur_file, tfile):
					refresh_files(tfile, tfile+".new")
					print('>> Entering test: %s'%tname)
					print('          output: %s'%tfile+".new")
					n_fine=0
					n_cor=0
					test_exists=True
					break
		if not test_exists:
			ignored_tests.append(cur_file)
	
	if not test_exists: continue

        match = re.search('Match\s(.*):', line)
	if match:
		match_name = match.group(1)
	if line[:20] == '   Calculated value ':
		calc_val = line[21:]
	if line[:20] == '   Reference value  ':
		ref_val = line[21:]

	if line[:10] == ' Execution':
		continue

	if '[   OK   ]' in line:
		n_fine += 1; ntot_fine += 1

	elif '[  FAIL  ]' in line:
		n_cor += 1; ntot_cor += 1
		buf = get_next_match(match_name)
		idx = buf.rfind(';')
		rest = buf[idx+1:].rstrip()
		spaces = ' '*rest.count(' ')
		buf = '%s;%s%s\n'%(buf[:idx], spaces, calc_val.strip())
		f_out.write(buf)

refresh_files()
f_result.close()

print('\nSummary:')
print('  %d values were fine'%(ntot_fine))
print('  %d values were corrected'%(ntot_cor))
print('  %d test file(s) was/were ignored:'%(len(ignored_tests)))
for ignored in ignored_tests:
	print('      %s'%(ignored))
print('')

#Workaround to get raw_input and input
try: input = raw_input
except: pass

def get_answer(question, default=False):
	ans=input(question).lower().strip()
	try:
		return ans[0]=='y'
	except:
		return default

if not get_answer('Would you like to see the diff now? (Y/n) ', True):
	sys.exit()

for tfile in test_files:
	if (os.path.isfile(tfile+'.new')):
                os.system('diff -U 0 %s %s.new'%(tfile,tfile))
                print('')

print('If you are happy with the diff, would you like me to update the files now?')
if not get_answer('I\'ll not commit anything, so you can still "svn diff" to see the changes. (Y/n) ', True):
	sys.exit()

for tfile in test_files:
	if (os.path.isfile(tfile+'.new')):
		os.system('mv %s.new %s'%(tfile,tfile))

print('\nAll done! Enjoy the time that I saved you!\n')
