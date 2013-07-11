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
## $Id$

hghdir = $(pkgdatadir)/PP/HGH

dist_hgh_DATA = \
	configurations \
	Ag.hgh \
	Ag_sc.hgh \
	Al.hgh \
	Ar.hgh \
	As.hgh \
	At.hgh \
	Au.hgh \
	Au_sc.hgh \
	Ba.hgh \
	Ba_sc.hgh \
	Be.hgh \
	Be_sc.hgh \
	B.hgh \
	Bi.hgh \
	Br.hgh \
	Ca.hgh \
	Ca_sc.hgh \
	Cd.hgh \
	Cd_sc.hgh \
	Ce_sc.hgh \
	C.hgh \
	Cl.hgh \
	Co.hgh \
	Co_sc.hgh \
	Cr.hgh \
	Cr_sc.hgh \
	Cs.hgh \
	Cs_sc.hgh \
	Cu.hgh \
	Cu_sc.hgh \
	Dy_sc.hgh \
	Er_sc.hgh \
	Eu_sc.hgh \
	Fe.hgh \
	Fe_sc.hgh \
	F.hgh \
	Ga.hgh \
	Ga_sc.hgh \
	Gd_sc.hgh \
	Ge.hgh \
	He.hgh \
	Hf_sc.hgh \
	Hg.hgh \
	Hg_sc.hgh \
	H.hgh \
	Ho_sc.hgh \
	I.hgh \
	In.hgh \
	In_sc.hgh \
	Ir.hgh \
	Ir_sc.hgh \
	K.hgh \
	Kr.hgh \
	K_sc.hgh \
	La_sc.hgh \
	Li.hgh \
	Li_sc.hgh \
	Lu_sc.hgh \
	Mg.hgh \
	Mg_sc.hgh \
	Mn.hgh \
	Mn_sc.hgh \
	Mo.hgh \
	Mo_sc.hgh \
	Na.hgh \
	Na_sc.hgh \
	Nb.hgh \
	Nb_sc.hgh \
	Nd_sc.hgh \
	Ne.hgh \
	N.hgh \
	Ni.hgh \
	Ni_sc.hgh \
	O.hgh \
	Os.hgh \
	Os_sc.hgh \
	Pb.hgh \
	Pd.hgh \
	Pd_sc.hgh \
	P.hgh \
	Pm_sc.hgh \
	Po.hgh \
	Pr_sc.hgh \
	Pt.hgh \
	Pt_sc.hgh \
	Rb.hgh \
	Rb_sc.hgh \
	Re.hgh \
	Re_sc.hgh \
	Rh.hgh \
	Rh_sc.hgh \
	Rn.hgh \
	Ru.hgh \
	Ru_sc.hgh \
	Sb.hgh \
	Sc.hgh \
	Sc_sc.hgh \
	Se.hgh \
	S.hgh \
	Si.hgh \
	Sm_sc.hgh \
	Sn.hgh \
	Sr.hgh \
	Sr_sc.hgh \
	Ta.hgh \
	Ta_sc.hgh \
	Tb_sc.hgh \
	Tc.hgh \
	Tc_sc.hgh \
	Te.hgh \
	Ti.hgh \
	Ti_sc.hgh \
	Tl.hgh \
	Tl_sc.hgh \
	Tm_sc.hgh \
	V.hgh \
	V_sc.hgh \
	W.hgh \
	W_sc.hgh \
	Xe.hgh \
	Yb_sc.hgh \
	Y.hgh \
	Y_sc.hgh \
	Zn.hgh \
	Zn_sc.hgh \
	Zr.hgh \
	Zr_sc.hgh


# In order to get `make check' working with a `VPATH' build (e. g. to use
# the `distcheck' target) the pseudopotential files have to be accessible
# in the build tree.
check:
	for f in $(dist_hgh_DATA); do \
	  if [ ! -e $$f -o $$f -ot $(srcdir)/$$f ]; then \
	    cp $(srcdir)/$$f $$f; \
	  fi \
	done

check-full: check

clean-local:
	@if [ -n "$(top_builddir)" -a "$(top_builddir)" != "$(top_srcdir)" ]; then \
	  rm -f $(dist_hgh_DATA); \
	fi

.PHONY: check check-full
