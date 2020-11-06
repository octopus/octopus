# Copyright 2015 Lorenz Hüdepohl
#
# This file is part of fdep and licensed under the MIT license
# see the file LICENSE for more information
#

define translate_name
$(subst -,_,$(subst .,_,$1))
endef

_f90_verbose = $(_f90_verbose_$(V))
_f90_verbose_ = $(_f90_verbose_$(AM_DEFAULT_VERBOSITY))
_f90_verbose_0 = @echo "  $1";
_f90_only_verbose = $(_f90_only_verbose_$(V))
_f90_only_verbose_ = @
_f90_only_verbose_0 = @
_f90_targets = $(call translate_name,$(PROGRAMS) $(LTLIBRARIES))

FORTRAN_CPP ?= cpp -P -traditional -Wall -Werror

# $1 source files
#
# returns: file without any .F90 .f90 .F .f extension
define strip_fortran_ext
$(patsubst %.F90,%,$(patsubst %.f90,%,$(patsubst %.F,%,$(patsubst %.f,%,$1))))
endef

# $1 program
#
# returns:
#  '1' if object files for target $1 are prefixed due to 'per-target' flags,
#  '' (the empty string) otherwise. See the automake manual for 'per-target'
#  compilation
#
define is_per_target
$(if $(filter $(call strip_fortran_ext,$(firstword $(call fortran_sources,$1))),$(patsubst %.o,%,$(patsubst %.lo,%,$($1_OBJECTS)))),,1)
endef

# $1 top-level target name (i.e. an entry of _f90_targets)
#
# returns: all target source files matching *.F90 *.f90 *.F *.f
define fortran_sources
$(filter %.F90 %.f90 %.F %.f,$($1_SOURCES))
endef

# $1 top-level target name
#
# returns: the appropriate extension (i.e. 'o' for normal programs, '.lo' for libraries)
define object_extension
$(if $(filter $1,$(subst -,_,$(PROGRAMS))),o,lo)
endef

_f90_depdir=$(builddir)/.fortran_dependencies
_f90_depmodfile = $(_f90_depdir)/dependencies_modules.mk
_f90_depincfile = $(_f90_depdir)/dependencies_includes.mk

# $1 source file
# $2 stem
# $3 program
# $4 kind of file ('use' or 'def')
define modinfo_name
$(_f90_depdir)/$(dir $1)$(2)$(call strip_fortran_ext,$(notdir $1)).$4_mods.$(call object_extension,$3)
endef

# $1 source_file
# $2 stem
# $3 program
define module_use_target
$(eval _use_mods += $(call modinfo_name,$1,$2,$3,use))
$(call modinfo_name,$1,$2,$3,use): $1 $(dir $1)$(am__dirstamp) ${HEADERS_DEP} | $(_f90_depdir)/$(dir $1)
	$(call _f90_verbose,F90 USE  [$3] $$<)${FORTRAN_CPP} $$< | $(top_srcdir)/build/preprocess.pl - ${PERL_ARGS} | \
		grep -i -o '^ *use [^ ,!:]*' | sed 's/^[[:space:]]*//;' | tr '[:upper:]' '[:lower:]' | sort -u > $$@

endef

# $1 source_file
# $2 stem
# $3 program
define module_def_target
$(eval _def_mods += $(call modinfo_name,$1,$2,$3,def))
$(call modinfo_name,$1,$2,$3,def): $1 $(dir $1)$(am__dirstamp) ${HEADERS_DEP} | $(_f90_depdir)/$(dir $1)
	$(call _f90_verbose,F90 MOD  [$3] $$<)${FORTRAN_CPP} $$< | $(top_srcdir)/build/preprocess.pl - ${PERL_ARGS} | \
		grep -i -o '^ *module [^!]*' | sed 's/^[[:space:]]*//;' | tr '[:upper:]' '[:lower:]' | grep -v "\<procedure\>\|\<intrinsic\>" > $$@ || true

endef

# $1 source_file
# $2 stem
# $3 program
define module_inc_target
$(eval _inc_mods += $(call modinfo_name,$1,$2,$3,inc))
$(call modinfo_name,$1,$2,$3,inc): $1 $(dir $1)$(am__dirstamp) ${HEADERS_DEP} | $(_f90_depdir)/$(dir $1)
	$(call _f90_verbose,F90 INC  [$3] $$<) cat $$< | \
	  grep -i -o '^ *# *include [^!]*' | \
	  sed 's-^[[:space:]]*#[[:space:]]*include[[:space:]]*-$(top_srcdir)/src/*/-;' | \
	  tr -d \" | grep -v \< | sort | uniq > $$@ || true

endef

# $1 source_file
# $2 stem
# $3 program
# only call the use and def targets once, not for each program
define module_targets
$(eval _$(3)_use_mods += $(call modinfo_name,$1,$2,$3,use))
$(if $(findstring $(call modinfo_name,$1,$2,$3,use),$(_use_mods)),,$(call module_use_target,$1,$2,$3))

$(eval _$(3)_def_mods += $(call modinfo_name,$1,$2,$3,def))
$(if $(findstring $(call modinfo_name,$1,$2,$3,def),$(_def_mods)),,$(call module_def_target,$1,$2,$3))

$(eval _$(3)_inc_mods += $(call modinfo_name,$1,$2,$3,inc))
$(if $(findstring $(call modinfo_name,$1,$2,$3,inc),$(_inc_mods)),,$(call module_inc_target,$1,$2,$3))
endef
$(foreach p,$(_f90_targets),$(if $(call is_per_target,$p),$(foreach s,$(call fortran_sources,$p),$(eval $(call module_targets,$s,$p-,$p))),$(foreach s,$(call fortran_sources,$p),$(eval $(call module_targets,$s,,$p)))))


# $1 target-name
define recursive_lib_deps
$(foreach l,$(call translate_name,$($1_LDADD) $($1_LIBADD)),$l $(call recursive_lib_deps,$l))
endef

# check also if make goals are empty
define is_clean
$(if $(strip $(MAKECMDGOALS)),$(if $(filter-out mostlyclean clean distclean maintainer-clean am--depfiles,$(MAKECMDGOALS)),0,1),0)
endef

define newline


endef


# ignore the dependencies for cleaning
ifneq ($(call is_clean),1)
# always include the dependencies for include statements
include $(_f90_depincfile)

ifeq ($(NODEP),1)
# If NODEP is set to 1, only include module dependencies if none of the targets
# in PROGRAMS or LTLIBRARIES is available. This ensures that a complete build
# happened before ignoring dependencies due to modules.
ifeq ($(wildcard $(PROGRAMS) $(LTLIBRARIES)),)
include $(_f90_depmodfile)
else
$(warning Ignoring Fortran module dependencies. If the compilation fails, please recompile without NODEP)
endif
else
# if NODEP is not set to 1, always include module dependencies
include $(_f90_depmodfile)
endif
endif

# $1 program
define program_module_dependencies
	$(_f90_only_verbose){ $(foreach argument,$(_$p_use_mods) $(_$p_def_mods) $(foreach l,$(call recursive_lib_deps,$p),$(_$l_use_mods) $(_$l_def_mods)),echo $(argument); ) true; } | \
	$(top_srcdir)/src/fdep/fortran_dependencies.pl mod $p >> $@ || { rm $@; exit 1; }

endef

# $1 program
define program_include_dependencies
	$(_f90_only_verbose){ $(foreach argument,$(_$p_inc_mods) $(foreach l,$(call recursive_lib_deps,$p),$(_$l_inc_mods)),echo $(argument); ) true; } | \
	$(top_srcdir)/src/fdep/fortran_dependencies.pl inc $p >> $@ || { rm $@; exit 1; }

endef

# gather dependencies for module files
$(_f90_depmodfile): $(top_srcdir)/src/fdep/fortran_dependencies.pl $(foreach p,$(_f90_targets),$(_$p_use_mods) $(_$p_def_mods))
	$(call _f90_verbose,F90 DEPS $@)echo > $@;
	$(foreach p,$(_f90_targets),$(call program_module_dependencies,$p))

# gather dependencies for includes
$(_f90_depincfile): $(top_srcdir)/src/fdep/fortran_dependencies.pl $(foreach p,$(_f90_targets),$(_$p_inc_mods))
	$(call _f90_verbose,F90 INCS $@)echo > $@;
	$(foreach p,$(_f90_targets),$(call program_include_dependencies,$p))

$(_f90_depdir):
	@mkdir $@

$(_f90_depdir)/%/:
	@mkdir -p $@


.PHONY: clean-fdep
clean-fdep:
	rm -rf $(_f90_depdir)
