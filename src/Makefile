###############################################################################
#                                                                             #
#                         Makefile for TDDFT                                  #
#                                                                             #
###############################################################################
PROG_NAME_3D = octopus.x
PROG_NAME_1D = octopus_1d.x

#############################################################################
# modules                                                                   # 
#############################################################################
SRC_ALL = fdf.F90 global.F90 pulpo.F90 units.F90 math.F90 spline.F90 fft.F90

SRC_1D = ps_xc.F90 ps_kb.F90 ps1D.F90 specie.F90 atom.F90 mesh1D.F90 \
	hartree1D.F90 states.F90 xc.F90 systm.F90 hamiltonian.F90 \
	lcao.F90 eigen.F90 mix.F90 scf.F90 static_pol.F90\
	lasers.F90 td.F90 run.F90 main.F90

SRC_3D = ps_xc.F90 ps_kb.F90 ps.F90 specie.F90 atom.F90 mesh3D.F90 \
	hartree.F90 states.F90 xc.F90 systm.F90 hamiltonian.F90 \
	lcao.F90 eigen.F90 mix.F90 scf.F90 static_pol.F90 \
	lasers.F90 td.F90 run.F90 main.F90

OBJS_1D = $(SRC_ALL:.F90=.o) $(SRC_1D:.F90=.o) 
OBJS_3D = $(SRC_ALL:.F90=.o) $(SRC_3D:.F90=.o) 

SPECTRUM_SRC = global.F90 fdf.F90 spectrum.F90 calc_spectrum.F90
SPECTRUM_OBJ = $(SPECTRUM_SRC:.F90=.o)

PLAN_SRC = global.F90 fdf.F90 psm_consts.F90 geom.F90 mesh.F90 \
	math.F90 fermi.F90 charge_plane.F90
PLAN_OBJ = $(PLAN_SRC:.F90=.o)

###############################################################################
# Rules
###############################################################################
%.o:%.F90
	$(CPP) $(CPPFLAGS) $< > $*.f90
	$(FC) $(FFLAGS) -o $@ -c $*.f90
	@rm -f $*.f90 

###############################################################################
# Makefile options........                                                    #
###############################################################################
all: what $(OBJS_3D)
	$(FC) $(LDFLAGS) $(OBJS_3D) \
	  $(TDDFTLIBS) -o $(EXEC_DIR)/$(PROG_NAME_3D)

1d: what  $(OBJS_1D)
	$(FC) $(LDFLAGS) $(OBJS_1D) \
	  $(TDDFTLIBS) -o $(EXEC_DIR)/$(PROG_NAME_1D)

what:
	@echo
	@echo " TDDFT compilation for $(ARCH) architecture......"
	@echo

spectrum: $(SPECTRUM_OBJ)
	$(FC) $(LDFLAGS) $(SPECTRUM_OBJ) \
	  $(UTIL_LIBS) -o $(EXEC_DIR)/spectrum.x

plan: $(PLAN_OBJ)
	$(FC) $(LDFLAGS) $(PLAN_OBJ) \
	  $(UTIL_LIBS) -o $(EXEC_DIR)/plan.x


###############################################################################
# Cleaning, taring, etc                                                       #
###############################################################################
touchfast:
	touch $(ALL_SRCS)

clean:
	@echo "==> Cleaning objects"
	rm -f *.o *.mod *~
	rm -f *.p.* *.op

###############################################################################
# Individual dependencies. Please update consistently...                      #
###############################################################################
mesh3D.o: mesh3D_create.F90 mesh3D_inc.F90
xc.o: xc_pot.F90 xc_LDA.F90 xc_GGA.F90
hamiltonian.o: hamiltonian_inc.F90 h_external_pot.F90
eigen.o: eigen_cg1.F90 eigen_cg2.F90
states.o: states_inc.F90
td.o: td_init.F90 td_rti.F90
mix.o: mix_inc.F90
