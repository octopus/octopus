#!/usr/bin/env bash
#
# Copyright (C) 2016 D. Strubbe
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
# $Id$

# Used by buildbot to determine if result of configure, encoded in config.h,
# has changed from the reference result.

BUILDER="$1"
BRANCH="$2"

if [ "x$BRANCH" == "x" ]; then
    BRANCH=trunk
fi
    
echo ""
echo "Checking config.h against reference:"

# show commands in terminal output
set -x

# download ref_file with svn
svn export --force http://www.tddft.org/svn/octopus/buildbot/config.h/$BRANCH/$BUILDER

# these two fields will generally be different, and that is fine
diff -I '^#define BUILD_TIME' -I '^#define LATEST_SVN' config.h $BUILDER
