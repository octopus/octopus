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

EXEEXT =

EXTRA_PROGRAMS = oct-propagation_spectrum oct-rotatory_strength oct-harmonic-spectrum \
                 oct-casida_spectrum oct-vibrational oct-xyz-anim oct-center-geom oct-help \
                 oct-oscillator-strength \
		 oct-dielectric-function \
		 oct-liquid
##		 oct-photoelectron_spectrum

bin_PROGRAMS = @UTILITY_PROGRAMS@

oct_propagation_spectrum_SOURCES = propagation_spectrum.F90
oct_propagation_spectrum_LDADD   = $(all_LIBS)

oct_oscillator_strength_SOURCES = oscillator_strength.F90
oct_oscillator_strength_LDADD   = $(all_LIBS)

oct_dielectric_function_SOURCES = dielectric_function.F90
oct_dielectric_function_LDADD   = $(all_LIBS)

oct_rotatory_strength_SOURCES = rotatory_strength.F90
oct_rotatory_strength_LDADD   = $(all_LIBS)

oct_harmonic_spectrum_SOURCES = harmonic_spectrum.F90
oct_harmonic_spectrum_LDADD   = $(all_LIBS)

##oct_photoelectron_spectrum_SOURCES = photoelectron_spectrum.F90
##oct_photoelectron_spectrum_LDADD   = @LIBS_LAPACK@ @LIBS_BLAS@ $(core_LIBS) $(external_LIBS)

oct_xyz_anim_SOURCES = xyzanim.F90
oct_xyz_anim_LDADD   = $(all_LIBS)

oct_casida_spectrum_SOURCES = casida_spectrum.F90
oct_casida_spectrum_LDADD =  $(core_LIBS) $(external_LIBS)

oct_vibrational_SOURCES = vibrational.F90
oct_vibrational_LDADD =  $(core_LIBS) $(external_LIBS)

oct_center_geom_SOURCES = centergeom.F90
oct_center_geom_LDADD = $(all_LIBS)

oct_liquid_SOURCES = liquid.F90
oct_liquid_LDADD = $(all_LIBS)

oct_help_SOURCES = help.F90
oct_help_LDADD =  $(core_LIBS)

depend:
	@makedepf90 $(core_SRCS) $(oct_dielectric_function_SOURCES) $(oct_propagation_spectrum_SOURCES) $(oct_rotatory_strength_SOURCES) $(oct_harmonic_spectrum_SOURCES) $(oct_xyz_anim_SOURCES) $(oct_casida_spectrum_SOURCES) $(oct_vibrational_SOURCES) $(oct_center_geom_SOURCES) $(oct_atoms_magnet_SOURCES) $(oct_help_SOURCES)
## $(oct_photoelectron_spectrum_SOURCES)

propagation_spectrum.o : propagation_spectrum.F90 
rotatory_strength.o : rotatory_strength.F90 
harmonic_spectrum.o : harmonic_spectrum.F90 
xyzanim.o : xyzanim.F90 
casida_spectrum.o : casida_spectrum.F90 
phonon_spectrum.o : phonon_spectrum.F90 
centergeom.o : centergeom.F90 
help.o : help.F90 
dielectric_function.o : dielectric_function.F90
liquid.o : liquid.F90
##photoelectron_spectrum.o : photoelectron_spectrum.F90
