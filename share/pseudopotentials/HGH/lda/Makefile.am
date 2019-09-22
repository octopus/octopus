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

hghdir = $(pkgdatadir)/pseudopotentials/HGH/lda

dist_hgh_DATA = \
	Ag.hgh \
	Al.hgh \
	Ar.hgh \
	As.hgh \
	At.hgh \
	Au.hgh \
	Ba.hgh \
	Be.hgh \
	B.hgh \
	Bi.hgh \
	Br.hgh \
	Ca.hgh \
	Cd.hgh \
	C.hgh \
	C.hgh.spacing \
	Cl.hgh \
	Co.hgh \
	Cr.hgh \
	Cs.hgh \
	Cu.hgh \
	D.hgh  \
	Fe.hgh \
	F.hgh \
	Ga.hgh \
	Ge.hgh \
	He.hgh \
	He.hgh.spacing \
	Hg.hgh \
	H.hgh \
	I.hgh \
	In.hgh \
	Ir.hgh \
	K.hgh \
	Kr.hgh \
	Li.hgh \
	Mg.hgh \
	Mn.hgh \
	Mo.hgh \
	Na.hgh \
	Nb.hgh \
	Ne.hgh \
	N.hgh \
	N.hgh.spacing \
	Ni.hgh \
	O.hgh \
	Os.hgh \
	Pb.hgh \
	Pd.hgh \
	P.hgh \
	Po.hgh \
	Pt.hgh \
	Rb.hgh \
	Re.hgh \
	Rh.hgh \
	Rn.hgh \
	Ru.hgh \
	Sb.hgh \
	Sc.hgh \
	Se.hgh \
	S.hgh \
	Si.hgh \
	Sn.hgh \
	Sr.hgh \
	T.hgh  \
	Ta.hgh \
	Tc.hgh \
	Te.hgh \
	Ti.hgh \
	Tl.hgh \
	V.hgh \
	W.hgh \
	Xe.hgh \
	Y.hgh \
	Zn.hgh \
	Zr.hgh


# In order to get `make check' working with a `VPATH' build (e. g. to use
# the `distcheck' target) the pseudopotential files have to be accessible
# in the build tree.
check:
	for f in $(dist_hgh_DATA); do \
	  if [ ! -e $$f -o $$f -ot $(srcdir)/$$f ]; then \
	    cp $(srcdir)/$$f $$f; \
	  fi \
	done

check-short: check
check-long: check

clean-local:
	@if [ -n "$(top_builddir)" -a "$(top_builddir)" != "$(top_srcdir)" ]; then \
	  rm -f $(dist_hgh_DATA); \
	fi

.PHONY: check check-short check-long
