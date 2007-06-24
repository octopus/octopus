# Define empty rule to override the native rule of automake (otherwise the native
# rule will have precedence over the pattern based rules below). Could not figure out
# how to remove .F90 from .SUFFIXES, which would be an alternative route.
.F90.o:

# compilation is a two step process: first we preprocess F90 files to generate _oct.f90 files
%_oct.f90: %.F90
	@FCCPP@ @CPPFLAGS@ $(AM_CPPFLAGS) -I. $< > $*_oct.f90
	@if [ "@DEBUG@" = "no" ]; then \
		cat $*_oct.f90 | grep -v pop_sub | grep -v push_sub >$*_oct.f91; \
		mv -f $*_oct.f91 $*_oct.f90; \
	fi
	@perl -pi -e 's/\\newline/\n/g; s/\\cardinal/#/g' $*_oct.f90

# then we actually compile _oct.f90 files and cleanup
%.o: %_oct.f90
	@FC@ @FCFLAGS@ @FCFLAGS_NETCDF@ $(AM_FCFLAGS) -c @FCFLAGS_f90@ -o $@ $*_oct.f90
	@rm -f $*_oct.f90

# ctags
CTAGS=ctags-exuberant -e

# cleaning
CLEANFILES = *~ *.bak *.mod *.il *.d *.pc* ifc* *_oct.f90 config_F90.h

!! Local Variables:
!! mode: Makefile
!! coding: utf-8
!! End:
