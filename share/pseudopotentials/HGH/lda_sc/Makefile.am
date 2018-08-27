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

hghdir = $(pkgdatadir)/pseudopotentials/HGH/lda_sc

dist_hgh_DATA = \
	Ag_sc.hgh \
	Au_sc.hgh \
	Ba_sc.hgh \
	Be_sc.hgh \
	Ca_sc.hgh \
	Cd_sc.hgh \
	Ce_sc.hgh \
	Co_sc.hgh \
	Cr_sc.hgh \
	Cs_sc.hgh \
	Cu_sc.hgh \
	Dy_sc.hgh \
	Er_sc.hgh \
	Eu_sc.hgh \
	Fe_sc.hgh \
	Ga_sc.hgh \
	Gd_sc.hgh \
	Hf_sc.hgh \
	Hg_sc.hgh \
	Ho_sc.hgh \
	In_sc.hgh \
	Ir_sc.hgh \
	K_sc.hgh \
	La_sc.hgh \
	Li_sc.hgh \
	Lu_sc.hgh \
	Mg_sc.hgh \
	Mn_sc.hgh \
	Mo_sc.hgh \
	Na_sc.hgh \
	Nb_sc.hgh \
	Nd_sc.hgh \
	Ni_sc.hgh \
	Os_sc.hgh \
	Pd_sc.hgh \
	Pm_sc.hgh \
	Pr_sc.hgh \
	Pt_sc.hgh \
	Rb_sc.hgh \
	Re_sc.hgh \
	Rh_sc.hgh \
	Ru_sc.hgh \
	Sc_sc.hgh \
	Sm_sc.hgh \
	Sr_sc.hgh \
	Ta_sc.hgh \
	Tb_sc.hgh \
	Tc_sc.hgh \
	Ti_sc.hgh \
	Tl_sc.hgh \
	Tm_sc.hgh \
	V_sc.hgh \
	W_sc.hgh \
	Yb_sc.hgh \
	Y_sc.hgh \
	Zn_sc.hgh \
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

check-short: check
check-long: check

clean-local:
	@if [ -n "$(top_builddir)" -a "$(top_builddir)" != "$(top_srcdir)" ]; then \
	  rm -f $(dist_hgh_DATA); \
	fi

.PHONY: check check-short check-long
