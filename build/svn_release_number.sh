#!/bin/sh
# 
# $Id$

# this is meant to be called from octopus/src/Makefile.am
if [ -x "$(which svn)" ]; then
	svn info | grep Revision | awk -F: '{print $2}' | tr -d [:space:]
else
	find . -type f -exec grep --exclude="./ChangeLog" --binary-files=without-match '$Id:' \{\} \; | tr -d \#\!\* | awk '{print $3,"["$2,$4"]"}' | grep -v 'qw(' | sort | tail -1
fi
