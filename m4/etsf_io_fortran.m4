# -*- Autoconf -*-
#
# Copyright (c) 2005-2008 ABINIT Group (Yann Pouillon)
# All rights reserved.
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Fortran compilers support
#



# _ABI_CHECK_FC_ABSOFT(COMPILER)
# ------------------------------
#
# Checks whether the specified Fortran compiler is the ABSoft Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_ABSOFT],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the ABSoft Fortran compiler])
 fc_info_string=`$1 -V 2> /dev/null`
 abi_result=`echo "${fc_info_string}" | grep '^Pro Fortran'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([FC_ABSOFT],1,[Define to 1 if you are using the ABSOFT Fortran compiler])
  fc_type="absoft"
  fc_version=`echo "${abi_result}" | sed -e 's/Pro Fortran //'`
  if test "${fc_version}" = "${abi_result}"; then
   fc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_ABSOFT



# _ABI_CHECK_FC_COMPAQ(COMPILER)
# ------------------------------
#
# Checks whether the specified Fortran compiler is the COMPAQ Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_COMPAQ],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the Compaq Fortran compiler])
 fc_info_string=`$1 -version 2>&1 | sed -e 's/^	//' | grep '^Compaq Fortran Compiler'`
 abi_result="${fc_info_string}"
 if test "${abi_result}" = ""; then
  fc_info_string=`$1 -version 2>&1 | sed -e 's/^	//' | grep '^HP Fortran Compiler'`
  abi_result="${fc_info_string}"
 fi
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([FC_COMPAQ],1,[Define to 1 if you are using the COMPAQ Fortran compiler])
  fc_type="compaq"
  fc_version=`echo "${abi_result}" | sed -e 's/.* V//;s/-.*//'`
  if test "${fc_version}" = "${abi_result}"; then
   fc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_COMPAQ



# _ABI_CHECK_FC_FUJITSU(COMPILER)
# -------------------------------
#
# Checks whether the specified Fortran compiler is the Fujitsu Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_FUJITSU],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the Fujitsu Fortran compiler])
 fc_info_string=`$1 -V 2> /dev/null`
 abi_result=`echo "${fc_info_string}" | grep '^Fujitsu Fortran'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([FC_FUJITSU],1,[Define to 1 if you are using the Fujitsu Fortran compiler])
  fc_type="fujitsu"
  fc_version=`echo "${abi_result}" | sed -e 's/.*Driver //;s/ .*//'`
  if test "${fc_version}" = "${abi_result}"; then
   fc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_FUJITSU



# _ABI_CHECK_FC_G95(COMPILER)
# ---------------------------
#
# Checks whether the specified Fortran compiler is the G95 Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_G95],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the G95 Fortran compiler])
 fc_info_string=`$1 --version 2>&1`
 abi_result=`echo "${fc_info_string}" | grep '^G95'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([FC_G95],1,[Define to 1 if you are using the G95 Fortran compiler])
  fc_type="g95"
  fc_version=`echo ${abi_result} | sed -e 's/.*GCC //; s/ .*//'`
  if test "${fc_version}" = "${abi_result}"; then
   fc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_G95



# _ABI_CHECK_FC_GCC(COMPILER)
# ---------------------------
#
# Checks whether the specified Fortran compiler is the GCC Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_GCC],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the GCC Fortran compiler])
 fc_info_string=`$1 --version 2>&1`
 abi_result=`echo "${fc_info_string}" | grep '^GNU Fortran'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([FC_GCC],1,[Define to 1 if you are using the GNU Fortran compiler])
  fc_type="gcc"
  fc_version=`echo ${abi_result} | sed -e 's/.*(GCC) //; s/.*GCC //; s/ .*//'`
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_GCC



# _ABI_CHECK_FC_HITACHI(COMPILER)
# -------------------------------
#
# Checks whether the specified Fortran compiler is the Hitachi Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_HITACHI],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the Hitachi Fortran compiler])
 fc_info_string=`$1 -V 2> /dev/null`
 abi_result=`echo "${fc_info_string}" | grep '^Hitachi Fortran'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([FC_HITACHI],1,[Define to 1 if you are using the Hitachi Fortran compiler])
  fc_type="hitachi"
  fc_version=`echo "${abi_result}" | sed -e 's/.*Driver //;s/ .*//'`
  if test "${fc_version}" = "${abi_result}"; then
   fc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_HITACHI



# _ABI_CHECK_FC_IBM(COMPILER)
# ---------------------------
#
# Checks whether the specified Fortran compiler is the IBM XL Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_IBM],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the IBM XL Fortran compiler])
 fc_info_string=`$1 -qversion 2>&1`
 fc_garbage=`$1 -qversion 2>&1 | wc -l | sed -e 's/ //g'`
 abi_result=`echo "${fc_info_string}" | grep 'IBM(R) XL Fortran'`
 if test "${abi_result}" = ""; then
  abi_result=`echo "${fc_info_string}" | grep 'IBM XL Fortran'`
 fi
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
  if test "${fc_garbage}" -gt 50; then
   AC_DEFINE([FC_IBM],1,[Define to 1 if you are using the IBM XL Fortran compiler])
   fc_type="ibm"
   fc_version="UNKNOWN"
   abi_result="yes"
  fi
 else
  AC_DEFINE([FC_IBM],1,[Define to 1 if you are using the IBM XL Fortran compiler])
  fc_type="ibm"
  fc_version=`echo "${abi_result}" | sed -e 's/.* V\([[0-9\.]]*\) .*/\1/'`
  if test "${fc_version}" = "${abi_result}"; then
   fc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_IBM



# _ABI_CHECK_FC_INTEL(COMPILER)
# -----------------------------
#
# Checks whether the specified Fortran compiler is the Intel Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_INTEL],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the Intel Fortran compiler])
 fc_info_string=`$1 -v -V 2>&1 | sed -e '/^ifc: warning/d'`
 abi_result=`echo "${fc_info_string}" | grep '^Intel(R) Fortran'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([FC_INTEL],1,[Define to 1 if you are using the Intel Fortran compiler])
  fc_type="intel"
  fc_version=`echo "${fc_info_string}" | grep '^Version' | sed -e 's/Version //;s/ .*//;s/ //g' | head -n 1`
  if test "${fc_version}" = ""; then
   fc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_INTEL



# _ABI_CHECK_FC_MIPSPRO(COMPILER)
# -------------------------------
#
# Checks whether the specified Fortran compiler is the MIPSpro Fortran
# compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_MIPSPRO],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the MIPSpro Fortran compiler])
 fc_info_string=`$1 -version 2>&1 | sed -e '/^$/d'`
 abi_result=`echo "${fc_info_string}" | grep '^MIPSpro'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([FC_MIPSPRO],1,[Define to 1 if you are using the MIPSpro Fortran compiler])
  fc_type="mipspro"
  fc_version=`echo "${abi_result}" | sed -e 's/.*Version //'`
  if test "${fc_version}" = "${abi_result}"; then
   fc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_MIPSPRO



# _ABI_CHECK_FC_OPEN64(COMPILER)
# ------------------------------
#
# Checks whether the specified Fortran compiler is the Open64
# Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_OPEN64],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the PathScale Fortran compiler])
 fc_info_string=`$1 --version 2>&1`
 abi_result=`echo "${fc_info_string}" | grep '^Open64'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([FC_OPEN64],1,[Define to 1 if you are using the Open64 Fortran compiler])
  fc_type="open64"
  fc_version=`echo "${abi_result}" | sed -e 's/.* Version //; s/ .*//'`
  if test "${fc_version}" = "${abi_result}"; then
   fc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_OPEN64



# _ABI_CHECK_FC_PATHSCALE(COMPILER)
# ---------------------------------
#
# Checks whether the specified Fortran compiler is the PathScale
# Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_PATHSCALE],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the PathScale Fortran compiler])
 fc_info_string=`$1 -version 2>&1`
 abi_result=`echo "${fc_info_string}" | grep '^PathScale'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([FC_PATHSCALE],1,[Define to 1 if you are using the PathScale Fortran compiler])
  fc_type="pathscale"
  fc_version=`echo "${abi_result}" | sed -e 's/.* Version //; s/ .*//'`
  if test "${fc_version}" = "${abi_result}"; then
   fc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_PATHSCALE



# _ABI_CHECK_FC_PGI(COMPILER)
# ---------------------------
#
# Checks whether the specified Fortran compiler is the Portland Group
# Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_PGI],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the Portland Group Fortran compiler])
 fc_info_string=`$1 -V 2>&1 | sed -e '/^$/d'`
 abi_result=`echo "${fc_info_string}" | grep '^pgf9[[05]]' | grep -v 'No files to process'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([FC_PGI],1,[Define to 1 if you are using the Portland Group Fortran compiler])
  fc_type="pgi"
  fc_version=`echo "${abi_result}" | sed -e 's/^pgf9[[05]] //' | sed -e 's/-.*//'`
  if test "${fc_version}" = "${abi_result}"; then
   fc_version="UNKNOWN"
  else
   if test "${fc_version}" = "6.0"; then
        AC_DEFINE([FC_PGI6],1,[Define to 1 if you are using the Portland Group Fortran compiler version 6])
   fi
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_PGI



# _ABI_CHECK_FC_SUN(COMPILER)
# ---------------------------
#
# Checks whether the specified Fortran compiler is the Sun WorkShop Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_SUN],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the Sun WorkShop Fortran compiler])
 fc_info_string=`$1 -V 2>&1 | head -n 1`
 abi_result=`echo "${fc_info_string}" | grep 'Sun' | grep 'Fortran 95'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([FC_SUN],1,[Define to 1 if you are using the Sun WorkShop])
  fc_type="sun"
  fc_version=`echo "${abi_result}" | sed -e 's/.* Fortran 95 //;s/ .*//'`
  if test "${fc_version}" = "${abi_result}" -o "${fc_version}" = ""; then
   fc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_SUN



 ##############################################################################



# _ABI_CHECK_FC_EXIT()
# --------------------
#
# Checks whether the Fortran compiler supports the exit() subroutine.
#
AC_DEFUN([_ABI_CHECK_FC_EXIT],
[dnl Init
 fc_has_exit="no"

 AC_MSG_CHECKING([whether the Fortran compiler accepts exit()])

 dnl Try to compile a program calling exit()
 AC_LANG_PUSH([Fortran])
 AC_LINK_IFELSE([AC_LANG_PROGRAM([],
  [[
      call exit(1)
  ]])], [fc_has_exit="yes"])
 AC_LANG_POP()

 if test "${fc_has_exit}" = "yes"; then
  AC_DEFINE([HAVE_FC_EXIT],1,
   [Define to 1 if your Fortran compiler supports exit()])
 fi

 AC_MSG_RESULT(${fc_has_exit})
]) # _ABI_CHECK_FC_EXIT



# _ABI_CHECK_FC_FLUSH()
# ---------------------
#
# Checks whether the Fortran compiler supports the flush() subroutine.
#
AC_DEFUN([_ABI_CHECK_FC_FLUSH],
[dnl Init
 fc_has_flush="no"

 AC_MSG_CHECKING([whether the Fortran compiler accepts flush()])

 dnl Try to compile a program calling flush()
 AC_LANG_PUSH([Fortran])
 AC_LINK_IFELSE([AC_LANG_PROGRAM([],
  [[
      call flush()
  ]])], [fc_has_flush="yes"])
 AC_LANG_POP()

 if test "${fc_has_flush}" = "yes"; then
  AC_DEFINE([HAVE_FC_FLUSH],1,
   [Define to 1 if your Fortran compiler supports flush()])
 fi

 AC_MSG_RESULT(${fc_has_flush})
]) # _ABI_CHECK_FC_FLUSH



# ABI_PROG_FC()
# -------------
#
# Tries to determine which type of Fortran compiler is installed.
#
AC_DEFUN([ABI_PROG_FC],
[dnl Init
 if test "${fc_type}" = ""; then
  fc_type="UNKNOWN"
 fi
 if test "${fc_version}" = ""; then
  fc_version="UNKNOWN"
 fi
 fc_wrap="no"

 dnl Determine Fortran compiler type (the order is important)
 AC_MSG_CHECKING([which type of Fortran compiler we have])

 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_G95(${FC})
 fi
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_GCC(${FC})
 fi
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_INTEL(${FC})
 fi
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_PATHSCALE(${FC})
 fi
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_PGI(${FC})
 fi
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_COMPAQ(${FC})
 fi
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_ABSOFT(${FC})
 fi
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_MIPSPRO(${FC})
 fi
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_OPEN64(${FC})
 fi
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_FUJITSU(${FC})
 fi
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_SUN(${FC})
 fi
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_HITACHI(${FC})
 fi
 dnl Always keep that one at the end
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_IBM(${FC})
 fi

 dnl Fall back to generic when detection fails
 if test "${fc_type}" = "UNKNOWN"; then
  fc_type="generic"
  fc_version="0.0"
 fi

 dnl Normalise Fortran compiler version
 fc_version=`echo ${fc_version} | cut -d. -f1-2`

 dnl Display final result
 AC_MSG_RESULT([${fc_type} ${fc_version}])

 dnl Schedule compiler info for substitution
 AC_SUBST(fc_type)
 AC_SUBST(fc_version)
 AC_SUBST(fc_wrap)

 dnl Further explore compiler peculiarities
 dnl Not needed here, commented on adding to Octopus -- DAS
 dnl _ABI_CHECK_FC_EXIT
 dnl _ABI_CHECK_FC_FLUSH
]) # ABI_PROG_FC
