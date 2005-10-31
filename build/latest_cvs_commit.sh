#!/bin/sh
# 
# $Id$

# this is meant to be called from octopus/src/Makefile.am
cd ../ && find . -type f -exec grep --exclude="./ChangeLog" --binary-files=without-match 'Exp \$' \{\} \; | tr -d \#\! | awk '{print $4" "$5" "$2" "$3}' | sort | tail -1
