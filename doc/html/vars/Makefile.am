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

dist_doc_DATA =     \
	alpha.js    \
	sections.js \
	vars_index.html \
	Generalities.html \
	Geometry_Optimization.html \
	Generator_Coordinates.html \
	Hamiltonian.html \
	Linear_Response.html \
	Math.html \
	Mesh.html \
	Optimal_Control.html \
	Output.html \
	SCF.html \
	States.html \
	System.html \
	Time_Dependent.html \
	Transport.html \
	Unoccupied_States.html \
	Utilities.html

CLEANFILES = $(dist_doc_DATA)

alpha.js: RUN_SCRIPT
sections.js: RUN_SCRIPT
vars_index.html: RUN_SCRIPT
Generalities.html: RUN_SCRIPT
Geometry_Optimization.html: RUN_SCRIPT
Generator_Coordinates.html: RUN_SCRIPT
Hamiltonian.html: RUN_SCRIPT
Linear_Response.html: RUN_SCRIPT
Math.html: RUN_SCRIPT
Mesh.html: RUN_SCRIPT
Optimal_Control.html: RUN_SCRIPT
Output.html: RUN_SCRIPT
SCF.html: RUN_SCRIPT
States.html: RUN_SCRIPT
System.html: RUN_SCRIPT
Time_Dependent.html: RUN_SCRIPT
Transport.html: RUN_SCRIPT
Unoccupied_States.html: RUN_SCRIPT
Utilities.html: RUN_SCRIPT

RUN_SCRIPT:
	$(top_srcdir)/build/var2html.pl -b $(top_builddir)
