#!/bin/bash

## Copyright (C) 2002-2006 J. Alberdi-Rodriguez
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.
##
## $Id$

#save the current directory
cd `dirname $0`
dox_dir=$PWD
cd ../..
oct_base_dir=$PWD

# go to the source directory
cd $oct_base_dir

echo "Configuring"
build/mk_varinfo.pl -s . -b .
# to produce options.h

# set all -DHAVE options
sed s'|#undef HAVE_|#define HAVE_|' config.h.in | grep '^#define' > src/include/config_F90.h
# make dummy headers
touch src/include/mpi.h
touch src/include/fcs_fconfig.h

cd $oct_base_dir/src

# FIXME: set FCCPP and version in Doxyfile

echo "Creating documentation"

# Doxyfile is for version 1.8.6
cp $dox_dir/Doxyfile .
cp $dox_dir/octopus.png .
doxygen

# clean up and move results
mv doxygen_doc $dox_dir/
rm Doxyfile octopus.png

cd ..

# undo modifications to this directory
cd src/include
rm config_F90.h
make config_F90.h
rm mpi.h fcs_fconfig.h
cd ../..

echo "Doxygen documentation created in $dox_dir/doxygen_doc/html"
