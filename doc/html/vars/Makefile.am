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
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
## 02111-1307, USA.
##
## $Id: Makefile.am 2440 2006-09-29 13:44:03Z marques $

docdir = $(pkgdatadir)/doc/vars

dist_doc_DATA = \
	alpha.js \
	sections.js \
	vars_index.html \
	Execution.html \
	Calculation_Modes.html \
	Hamiltonian.html \
	Linear_Response.html \
	Math.html \
	Mesh.html \
	Output.html \
	SCF.html \
	States.html \
	System.html \
	Time-Dependent.html \
	Open_Boundaries.html \
	Utilities.html

CLEANFILES = $(dist_doc_DATA)

alpha.js: RUN_SCRIPT
sections.js: RUN_SCRIPT
vars_index.html: RUN_SCRIPT
Execution.html: RUN_SCRIPT
Calculation_Modes.html: RUN_SCRIPT
Hamiltonian.html: RUN_SCRIPT
Linear_Response.html: RUN_SCRIPT
Math.html: RUN_SCRIPT
Mesh.html: RUN_SCRIPT
Output.html: RUN_SCRIPT
SCF.html: RUN_SCRIPT
States.html: RUN_SCRIPT
System.html: RUN_SCRIPT
Time-Dependent.html: RUN_SCRIPT
Open_Boundaries.html: RUN_SCRIPT
Utilities.html: RUN_SCRIPT

RUN_SCRIPT:
	$(top_srcdir)/build/var2html.pl -b $(top_builddir)
