#!/bin/sh
# 
# $Id$

# this is meant to be called from octopus/src/Makefile.am
cd `dirname "$0"`/..
if [ -x "$(which svn)" ] && svn info > /dev/null 2>&1 ; then
	svn info | grep Revision | awk -F: '{print $2}' | tr -d [:space:]
else
	find . -type f ! -name ChangeLog ! -name \*.svn\* \
	    ! -name \*.o ! -name \*.a ! -name \*.so \
	    -exec grep '$Id:' \{\} \; | \
	    tr -d \#\!\* | awk '{print $3,"["$2,$4"]"}' | \
	    grep -v 'qw(' | sort | tail -1
fi
