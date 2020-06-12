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

docdir = $(pkgdatadir)/doc/vars

dist_doc_DATA = \
	alpha.js \
	sections.js \
	vars_index.html \
	Atomic_Orbitals.html \
	ClassicalParticles.html \
	Execution.html \
	Calculation_Modes.html \
	Hamiltonian.html \
	Linear_Response.html \
	Math.html \
	MaxwellStates.html \
	Mesh.html \
	Output.html \
	SCF.html \
	States.html \
	System.html \
	Time-Dependent.html \
	Utilities.html

CLEANFILES = $(dist_doc_DATA)

alpha.js: RUN_SCRIPT
sections.js: RUN_SCRIPT
vars_index.html: RUN_SCRIPT
Atomic_Orbitals.html: RUN_SCRIPT
Execution.html: RUN_SCRIPT
Calculation_Modes.html: RUN_SCRIPT
Hamiltonian.html: RUN_SCRIPT
Linear_Response.html: RUN_SCRIPT
Math.html: RUN_SCRIPT
Mesh.html: RUN_SCRIPT
Output.html: RUN_SCRIPT
SCF.html: RUN_SCRIPT
States.html: RUN_SCRIPT
MaxwellStates.html: RUN_SCRIPT
System.html: RUN_SCRIPT
Time-Dependent.html: RUN_SCRIPT
ClassicalParticles.html: RUN_SCRIPT
Utilities.html: RUN_SCRIPT

RUN_SCRIPT:
	$(top_srcdir)/build/var2html.pl -b $(top_builddir)
