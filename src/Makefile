###############################################################################
#                                                                             #
#                         Makefile for TDDFT                                  #
#                                                                             #
###############################################################################
PROG_NAME = octopus.x

#############################################################################
# modules                                                                   # 
#############################################################################

F90MOD_SRC = fdf.F90 global.F90 pulpo.F90 units.F90 math.F90 spline.F90 \
        fft.F90 ps_kb.F90 ps.F90 specie.F90 atom.F90 mesh3D.F90 \
	hartree.F90 xc.F90\
        states.F90 systm.F90 hamiltonian.F90 run.F90
F90MOD_OBJ= $(F90MOD_SRC:.F90=.o)

SRCS_90 = ps_xc.F90 main.F90
OBJS_90=$(SRCS_90:.F90=.o)

SPECTRUM_SRC = global.F90 fdf.F90 spectrum.F90 calc_spectrum.F90
SPECTRUM_OBJ = $(SPECTRUM_SRC:.F90=.o)

PLAN_SRC = global.F90 fdf.F90 psm_consts.F90 geom.F90 mesh.F90 \
	math.F90 fermi.F90 \
	charge_plane.F90
PLAN_OBJ = $(PLAN_SRC:.F90=.o)

###############################################################################
# FORTRAN 77/Fortran 90 files                                                 #
###############################################################################

ALL_SRCS=$(F90MOD_SRC) $(SRCS_90)
ALL_OBJS=$(F90MOD_OBJ) $(OBJS_90)

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
all: what tddft

prep:
	for i in $(ALL_SRCS) $(SPECTRUM_SRC); do \
	  $(CPP) $(CPPFLAGS) $$i >tmp; \
	  mv -f tmp ../src.$(DEST_ARCH)/`echo $$i | sed 's/F90/f90/'` ; \
	done

what:
	@echo
	@echo " TDDFT compilation for $(ARCH) architecture......"
	@echo


tddft: $(ALL_OBJS)
	$(FC) $(LDFLAGS) $(ALL_OBJS) \
	  $(TDDFTLIBS) -o $(EXEC_DIR)/$(PROG_NAME) 

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
