## Process this file with automake to produce Makefile.in

## Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

bin_SCRIPTS = \
	annotate_input.py

dta1dir = $(pkgdatadir)
dist_dta1_DATA = variables.py

EXTRA_DIST =                   \
	mk_varinfo.pl          \
	parse_variables.py     \
	variables.py           \
	annotate_input.py      \
	mk_functionals_list.pl \
	mk_kfunctionals_list.pl \
	preprocess.pl          \
	var2html.pl            \
	git_commit_hash.sh     \
	true_compiler.sh       \
	octopus.spec.in

DISTCLEANFILES = \
	octopus.spec
