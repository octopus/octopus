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
## $Id: Makefile.am 2995 2007-06-13 17:49:22Z xavier $

include $(srcdir)/../common-rules.make

EXTRA_PROGRAMS = oct-cross-section oct-rsf oct-harmonic-spectrum \
                 oct-broad oct-vibrational oct-xyz-anim oct-center-geom oct-help \
                 oct-oscillator-strength

bin_PROGRAMS = @UTILITY_PROGRAMS@

oct_cross_section_SOURCES = cross_section.F90
oct_cross_section_LDADD   = @LIBS_LAPACK@ @LIBS_BLAS@ $(core_LIBS) $(external_LIBS)

oct_oscillator_strength_SOURCES = oscillator_strength.F90
oct_oscillator_strength_LDADD   = @LIBS_LAPACK@ @LIBS_BLAS@ $(core_LIBS) $(external_LIBS)

oct_rsf_SOURCES = rotational_strength.F90
oct_rsf_LDADD   = @LIBS_LAPACK@ @LIBS_BLAS@ $(core_LIBS) $(external_LIBS)

oct_harmonic_spectrum_SOURCES = harmonic_spectrum.F90
oct_harmonic_spectrum_LDADD   = @LIBS_LAPACK@ @LIBS_BLAS@ $(core_LIBS) $(external_LIBS)

oct_xyz_anim_SOURCES = xyzanim.F90
oct_xyz_anim_LDADD   = $(all_LIBS)

oct_broad_SOURCES = broad.F90
oct_broad_LDADD =  $(core_LIBS) $(external_LIBS)

oct_vibrational_SOURCES = phonon_spectrum.F90
oct_vibrational_LDADD =  $(core_LIBS) $(external_LIBS)

oct_center_geom_SOURCES = centergeom.F90
oct_center_geom_LDADD = $(all_LIBS)

oct_help_SOURCES = help.F90
oct_help_LDADD =  $(core_LIBS)

depend:
	@makedepf90 $(core_SRCS) $(oct_cross_section_SOURCES) $(oct_rsf_SOURCES) $(oct_harmonic_spectrum_SOURCES) $(oct_xyz_anim_SOURCES) $(oct_broad_SOURCES) $(oct_vibrational_SOURCES) $(oct_center_geom_SOURCES) $(oct_atoms_magnet_SOURCES) $(oct_help_SOURCES)
