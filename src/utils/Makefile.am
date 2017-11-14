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

include $(srcdir)/../common-rules.make

EXEEXT =

EXTRA_PROGRAMS = oct-propagation_spectrum \
	         oct-harmonic-spectrum \
                 oct-casida_spectrum      \
                 oct-vibrational_spectrum \
                 oct-infrared_spectrum    \
                 oct-xyz-anim \
                 oct-center-geom \
                 oct-help \
                 oct-oscillator-strength \
                 oct-dielectric-function \
                 oct-photoelectron_spectrum \
                 oct-convert \
                 oct-local_multipoles \
                 oct-conductivity \
                 oct-floquet_observables \
                 oct-wannier90_interface

bin_PROGRAMS = @UTILITY_PROGRAMS@

oct_propagation_spectrum_SOURCES = propagation_spectrum.F90
oct_propagation_spectrum_LDADD   = $(top_builddir)/src/td/libtd.a \
                                   $(top_builddir)/src/hamiltonian/libhamiltonian.a \
                                   $(top_builddir)/src/states/libstates.a \
                                   $(top_builddir)/src/frozen/libfrozen.a \
                                   $(top_builddir)/src/grid/libgrid.a \
                                   $(top_builddir)/src/ions/libions.a \
                                   $(top_builddir)/src/species/libspecies.a \
                                   $(top_builddir)/src/math/libmath.a \
                                   $(top_builddir)/src/basic/libbasic.a \
                                   $(other_LIBS)

oct_oscillator_strength_SOURCES = oscillator_strength.F90
oct_oscillator_strength_LDADD    = $(top_builddir)/src/td/libtd.a \
                                   $(top_builddir)/src/hamiltonian/libhamiltonian.a \
                                   $(top_builddir)/src/states/libstates.a \
                                   $(top_builddir)/src/frozen/libfrozen.a \
                                   $(top_builddir)/src/grid/libgrid.a \
                                   $(top_builddir)/src/ions/libions.a \
                                   $(top_builddir)/src/species/libspecies.a \
                                   $(top_builddir)/src/math/libmath.a \
                                   $(top_builddir)/src/basic/libbasic.a \
                                   $(other_LIBS)

oct_dielectric_function_SOURCES = dielectric_function.F90
oct_dielectric_function_LDADD    = $(top_builddir)/src/td/libtd.a \
                                   $(top_builddir)/src/hamiltonian/libhamiltonian.a \
                                   $(top_builddir)/src/states/libstates.a \
                                   $(top_builddir)/src/frozen/libfrozen.a \
                                   $(top_builddir)/src/grid/libgrid.a \
                                   $(top_builddir)/src/ions/libions.a \
                                   $(top_builddir)/src/species/libspecies.a \
                                   $(top_builddir)/src/math/libmath.a \
                                   $(top_builddir)/src/basic/libbasic.a \
                                   $(other_LIBS)

oct_harmonic_spectrum_SOURCES = harmonic_spectrum.F90
oct_harmonic_spectrum_LDADD      = $(top_builddir)/src/td/libtd.a \
                                   $(top_builddir)/src/hamiltonian/libhamiltonian.a \
                                   $(top_builddir)/src/states/libstates.a \
                                   $(top_builddir)/src/frozen/libfrozen.a \
                                   $(top_builddir)/src/grid/libgrid.a \
                                   $(top_builddir)/src/ions/libions.a \
                                   $(top_builddir)/src/species/libspecies.a \
                                   $(top_builddir)/src/math/libmath.a \
                                   $(top_builddir)/src/basic/libbasic.a \
                                   $(other_LIBS)

oct_photoelectron_spectrum_SOURCES = photoelectron_spectrum.F90
oct_photoelectron_spectrum_LDADD = $(top_builddir)/src/td/libtd.a \
                                   $(top_builddir)/src/system/libsystem.a \
                                   $(top_builddir)/src/hamiltonian/libhamiltonian.a \
                                   $(top_builddir)/src/states/libstates.a \
                                   $(top_builddir)/src/frozen/libfrozen.a \
                                   $(top_builddir)/src/poisson/libpoisson.a \
                                   $(top_builddir)/src/grid/libgrid.a \
                                   $(top_builddir)/src/ions/libions.a \
                                   $(top_builddir)/src/species/libspecies.a \
                                   $(top_builddir)/src/math/libmath.a \
                                   $(top_builddir)/src/basic/libbasic.a \
                                   $(other_LIBS)

oct_xyz_anim_SOURCES = xyzanim.F90
oct_xyz_anim_LDADD               = $(top_builddir)/src/ions/libions.a \
                                   $(top_builddir)/src/species/libspecies.a \
                                   $(top_builddir)/src/math/libmath.a \
                                   $(top_builddir)/src/basic/libbasic.a \
                                   $(other_LIBS)

oct_casida_spectrum_SOURCES = casida_spectrum.F90
oct_casida_spectrum_LDADD        = $(top_builddir)/src/ions/libions.a \
                                   $(top_builddir)/src/species/libspecies.a \
                                   $(top_builddir)/src/math/libmath.a \
                                   $(top_builddir)/src/basic/libbasic.a \
                                   $(other_LIBS)

oct_vibrational_spectrum_SOURCES = vibrational.F90
oct_vibrational_spectrum_LDADD   = $(top_builddir)/src/td/libtd.a \
                                   $(top_builddir)/src/hamiltonian/libhamiltonian.a \
                                   $(top_builddir)/src/states/libstates.a \
                                   $(top_builddir)/src/frozen/libfrozen.a \
                                   $(top_builddir)/src/grid/libgrid.a \
                                   $(top_builddir)/src/ions/libions.a \
                                   $(top_builddir)/src/species/libspecies.a \
                                   $(top_builddir)/src/math/libmath.a \
                                   $(top_builddir)/src/basic/libbasic.a \
                                   $(other_LIBS)

oct_infrared_spectrum_SOURCES = infrared.F90
oct_infrared_spectrum_LDADD      = $(top_builddir)/src/grid/libgrid.a \
                                   $(top_builddir)/src/ions/libions.a \
                                   $(top_builddir)/src/species/libspecies.a \
                                   $(top_builddir)/src/math/libmath.a \
                                   $(top_builddir)/src/basic/libbasic.a \
                                   $(other_LIBS)

oct_center_geom_SOURCES = centergeom.F90
oct_center_geom_LDADD            = $(top_builddir)/src/ions/libions.a \
                                   $(top_builddir)/src/species/libspecies.a \
                                   $(top_builddir)/src/math/libmath.a \
                                   $(top_builddir)/src/basic/libbasic.a \
                                   $(other_LIBS)

oct_convert_SOURCES = convert.F90
oct_convert_LDADD = $(top_builddir)/src/td/libtd.a \
                                   $(top_builddir)/src/system/libsystem.a \
                                   $(top_builddir)/src/hamiltonian/libhamiltonian.a \
                                   $(top_builddir)/src/states/libstates.a \
                                   $(top_builddir)/src/frozen/libfrozen.a \
                                   $(top_builddir)/src/poisson/libpoisson.a \
                                   $(top_builddir)/src/grid/libgrid.a \
                                   $(top_builddir)/src/ions/libions.a \
                                   $(top_builddir)/src/species/libspecies.a \
                                   $(top_builddir)/src/math/libmath.a \
                                   $(top_builddir)/src/basic/libbasic.a \
                                   $(other_LIBS)

oct_local_multipoles_SOURCES = local_multipoles.F90 \
                                local_write.F90 

oct_local_multipoles_LDADD = $(top_builddir)/src/td/libtd.a \
                                   $(top_builddir)/src/system/libsystem.a \
                                   $(top_builddir)/src/hamiltonian/libhamiltonian.a \
                                   $(top_builddir)/src/states/libstates.a \
                                   $(top_builddir)/src/frozen/libfrozen.a \
                                   $(top_builddir)/src/poisson/libpoisson.a \
                                   $(top_builddir)/src/grid/libgrid.a \
                                   $(top_builddir)/src/ions/libions.a \
                                   $(top_builddir)/src/species/libspecies.a \
                                   $(top_builddir)/src/math/libmath.a \
                                   $(top_builddir)/src/basic/libbasic.a \
                                   $(other_LIBS)

oct_conductivity_SOURCES = conductivity.F90
oct_conductivity_LDADD  = $(top_builddir)/src/td/libtd.a \
                                   $(top_builddir)/src/hamiltonian/libhamiltonian.a \
                                   $(top_builddir)/src/states/libstates.a \
                                   $(top_builddir)/src/frozen/libfrozen.a \
                                   $(top_builddir)/src/grid/libgrid.a \
                                   $(top_builddir)/src/ions/libions.a \
                                   $(top_builddir)/src/species/libspecies.a \
                                   $(top_builddir)/src/math/libmath.a \
                                   $(top_builddir)/src/basic/libbasic.a \
                                   $(other_LIBS)


oct_help_SOURCES = help.F90
oct_help_LDADD                   = $(top_builddir)/src/basic/libbasic.a \
                                   $(other_LIBS)

oct_floquet_observables_SOURCES = floquet_observables.F90 
oct_floquet_observables_LDADD = $(top_builddir)/src/td/libtd.a \
                                   $(top_builddir)/src/scf/libscf.a \
                                   $(top_builddir)/src/system/libsystem.a \
                                   $(top_builddir)/src/hamiltonian/libhamiltonian.a \
                                   $(top_builddir)/src/states/libstates.a \
                                   $(top_builddir)/src/frozen/libfrozen.a \
                                   $(top_builddir)/src/poisson/libpoisson.a \
                                   $(top_builddir)/src/grid/libgrid.a \
                                   $(top_builddir)/src/ions/libions.a \
                                   $(top_builddir)/src/species/libspecies.a \
                                   $(top_builddir)/src/math/libmath.a \
                                   $(top_builddir)/src/basic/libbasic.a \
                                   $(other_LIBS)

oct_wannier90_interface_SOURCES = wannier90_interface.F90
oct_wannier90_interface_LDADD = $(top_builddir)/src/td/libtd.a \
                                   $(top_builddir)/src/scf/libscf.a \
                                   $(top_builddir)/src/system/libsystem.a \
                                   $(top_builddir)/src/hamiltonian/libhamiltonian.a \
                                   $(top_builddir)/src/states/libstates.a \
                                   $(top_builddir)/src/frozen/libfrozen.a \
                                   $(top_builddir)/src/poisson/libpoisson.a \
                                   $(top_builddir)/src/grid/libgrid.a \
                                   $(top_builddir)/src/ions/libions.a \
                                   $(top_builddir)/src/species/libspecies.a \
                                   $(top_builddir)/src/math/libmath.a \
                                   $(top_builddir)/src/basic/libbasic.a \
                                   $(other_LIBS)


# Include files - keep in alphabetical order
noinst_HEADERS =                         \
	wannier90_interface_defs_from_pwscf.F90 


# debug files from MacOS gcc with -g
clean-local :
	rm -rf *.dSYM


# This is to generate the dependencies. It makes use of the makedepf90
# program of Erik Edelmann
# (http://www.helsinki.fi/~eedelman/makedepf90.html). It should be run
# from time to time (just typing make depend within the src
# directory), substituting the dependencies below by the resulting new
# dependencies. (For normal compilation of octopus, there is
# absolutely no need of installing makedepf90).
depend:
	@makedepf90 $(core_SRCS) $(oct_dielectric_function_SOURCES) $(oct_propagation_spectrum_SOURCES) $(oct_harmonic_spectrum_SOURCES) $(oct_xyz_anim_SOURCES) $(oct_casida_spectrum_SOURCES) $(oct_vibrational_spectrum_SOURCES) $(oct_infrared_spectrum_SOURCES) $(oct_center_geom_SOURCES) $(oct_help_SOURCES) $(oct_photoelectron_spectrum_SOURCES) $(oct_convert_SOURCES) $(oct_local_multipoles_SOURCES) $(oct_floquet_observables_SOURCES)


# DO NOT EDIT DIRECTLY --- use make depend to generate the dependencies and copy-paste the result here

dielectric_function.o : dielectric_function.F90 
propagation_spectrum.o : propagation_spectrum.F90 
harmonic_spectrum.o : harmonic_spectrum.F90 
xyzanim.o : xyzanim.F90 
casida_spectrum.o : casida_spectrum.F90 
vibrational.o : vibrational.F90 
infrared.o : infrared.F90 
centergeom.o : centergeom.F90 
help.o : help.F90 
photoelectron_spectrum.o : photoelectron_spectrum.F90 
convert.o : convert.F90 
local_multipoles.o : local_multipoles.F90 local_write.o 
local_write.o : local_write.F90 
floquet_observables.o: floquet_observables.F90
wannier90_interface.o: wannier90_interface.F90
