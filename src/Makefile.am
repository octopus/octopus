bin_PROGRAMS = octopus

octopus_SOURCES = global.F90 liboct.F90 io.F90 pulpo.F90 units.F90 math.F90 \
	spline.F90 fft.F90 ps_xc.F90 ps_kb.F90 ps.F90 specie.F90 atom.F90 \
	mesh.F90 hartree.F90 states.F90 xc.F90 systm.F90 lasers.F90 hamiltonian.F90 \
	eigen.F90 mix.F90 lcao.F90 scf.F90 unocc.F90 static_pol.F90 \
	td.F90 run.F90 main.F90

EXTRA_octopus_SOURCES = \
	mesh1D.F90 mesh1D_create.F90 mesh1D_inc.F90 ps1D.F90 hartree1D.F90 \
	mesh3D.F90 mesh3D_create.F90 mesh3D_inc.F90 ps3D.F90 hartree3D.F90 \
	xc_pot.F90 xc_LDA.F90 xc_GGA.F90 xc_KLI.F90 xc_KLI_x.F90 xc_KLI_SIC.F90 \
	hamiltonian_inc.F90 h_external_pot.F90 h_forces.F90 \
	eigen_cg1.F90 eigen_cg2.F90 \
	states_inc.F90 \
	td_init.F90 td_rti.F90 \
	mix_inc.F90

octopus_LDADD = $(LIB_FFTW) $(LIB_GSL) $(LIB_LAPACK) ../liboct/liboct.a

CLEANFILES = *~ *.bak
