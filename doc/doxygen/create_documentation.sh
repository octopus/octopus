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
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
## 02111-1307, USA.
##
## $Id$

#save the current directory
dox_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $dox_dir
cd ../..
oct_base_dir=$PWD

# go to the source directory
cd $oct_base_dir/src

# iterate all folders
for folder in $(ls)
do
    if [ -d $folder ]; then
	echo "Creating documentation for $folder"
	cd $oct_base_dir/src/$folder
	cp $dox_dir/Doxyfile .
	#change project name
	sed -i "s/Octopus/Octopus_$folder/" Doxyfile
	#copy to the new location
	sed -i "s|doxygen_doc|$dox_dir/$folder|" Doxyfile
	#call doxygen
	doxygen >/dev/null 2>/dev/null
	
        # go back to the source directory, to iterate next folder
	cd $oct_base_dir/src
    fi
done
cd ..
echo "Doxygen documentation created in doc/doxygen/index.html"
