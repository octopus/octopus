!! Copyright (C) 2003-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$

#include "config_F90.h"

#define MAX_SPIN 4

#if defined(LONG_LINES)
#  define _newline_
#  define _anl_
#else
#  define _newline_    \newline
#  define _anl_      & \newline
#endif

#if defined(LONG_LINES)
#  define CARDINAL
#else
#  if defined(F90_ACCEPTS_LINE_NUMBERS)
#    define CARDINAL _newline_\cardinal __LINE__ __FILE__
#  else
#    define CARDINAL _newline_
#  endif
#endif

#define __STRING(x)     #x

#if !defined(NDEBUG)
#  define ASSERT(expr)  \
  if(.not.(expr)) _anl_ \
     call assert_die(__STRING(expr), _anl_ __FILE__, _anl_  __LINE__) \
  CARDINAL
#else
#  define ASSERT(expr)
#endif

#ifdef HAVE_FC_SIZEOF
#  define SIZEOF(x) sizeof(x)
#else
#  define SIZEOF(x) 1
#endif

#if defined(NDEBUG)
#  define SAFE_ALLOCATE(x) allocate(x)
#  define SAFE_DEALLOCATE_P(x) if(associated(x)) then; deallocate(x); nullify(x); end if
#  define SAFE_DEALLOCATE_A(x) if(allocated(x)) then; deallocate(x); end if
#else
#  define SAFE_ALLOCATE(x)			\
  allocate(x, stat=global_alloc_err); _newline_ \
  if(iand(prof_vars%mode, PROFILING_MEMORY).ne.0 .or. global_alloc_err.ne.0) _anl_ \
  global_sizeof = SIZEOF(x); _newline_	\
  if(iand(prof_vars%mode, PROFILING_MEMORY).ne.0) _anl_ \
    call profiling_memory_allocate(_anl_ #x, _anl_ __FILE__, _anl_ __LINE__, _anl_ global_sizeof); _newline_ \
  if(global_alloc_err.ne.0) _anl_ \
    call alloc_error(global_sizeof, _anl_ __FILE__, _anl_ __LINE__); \
  CARDINAL

#  define MY_DEALLOCATE(x) \
  global_sizeof = SIZEOF(x); _newline_ \
  deallocate(x, stat=global_alloc_err); _newline_ \
  if(iand(prof_vars%mode, PROFILING_MEMORY).ne.0) _anl_ \
    call profiling_memory_deallocate(#x, _anl_ __FILE__, _anl_ __LINE__, _anl_ global_sizeof); _newline_ \
  if(global_alloc_err.ne.0) _anl_ \
    call dealloc_error(global_sizeof, _anl_ __FILE__, _anl_ __LINE__); \
  CARDINAL

#  define SAFE_DEALLOCATE_P(x) \
  if(associated(x)) then; _newline_ \
    MY_DEALLOCATE(x);     _newline_ \
    nullify(x);           _newline_ \
  end if
#  define SAFE_DEALLOCATE_A(x) \
  if(allocated(x)) then;  _newline_ \
    MY_DEALLOCATE(x);     _newline_ \
  end if

#endif

#define ALLOCATE(a,b) _DEPRECATED_PLEASE_USE_SAFE_ALLOCATE_

#define REAL_DOUBLE real(8)
#define REAL_SINGLE real(4)

#if defined(SINGLE_PRECISION)
#  define REAL_PRECISION 4
#  define FLOAT     	 real(4)
#  define MPI_FLOAT 	 MPI_REAL
#  define CMPLX     	 complex(4)
#  define MPI_CMPLX 	 MPI_COMPLEX
#  define PREC(x)   	 s ## x
#  define ZPREC(x)   	 c ## x
#  define CNST(x)   	 x ## _4
#  define XC_F90(x)      xc_s_f90_ ## x
#else
#  define REAL_PRECISION 8
#  define FLOAT          real(8)
#  define MPI_FLOAT 	 MPI_DOUBLE_PRECISION
#  define CMPLX     	 complex(8)
#  define MPI_CMPLX 	 MPI_DOUBLE_COMPLEX
#  define PREC(x)   	 d ## x
#  define ZPREC(x)   	 z ## x
#  define CNST(x)   	 x ## _8
#  define XC_F90(x)      xc_f90_ ## x
#endif

#define   TOFLOAT(x) real(x, REAL_PRECISION)
#define   TOCMPLX(x, y) cmplx(x, y, REAL_PRECISION)

#define M_ONE CNST(1.0)
#define M_ZERO CNST(0.0)

#define GS_DIR "gs/"
#define STATIC_DIR "static/"
#define EM_RESP_DIR "em_resp/"
#define EM_RESP_FD_DIR "em_resp_fd/"
#define KDOTP_DIR "kdotp/"
#define VIB_MODES_DIR "vib_modes/"
#define VDW_DIR "vdw/"
#define CASIDA_DIR "casida/"
#define OCT_DIR "opt-control/"

#ifdef HAVE_MPI
#ifdef HAVE_MPI2
#define MPI_IN_PLACE_OR(x) MPI_IN_PLACE
#else
#define MPI_IN_PLACE_OR(x) x
#endif
#endif

! the TOSTRING macro converts a macro into a string
! do not use the STRINGIFY macro
#define STRINGIFY(x) #x
#define TOSTRING(x)  STRINGIFY(x)

#define INCR(x, y) x = (x) + (y)

#define PUSH_SUB(routine) \
  if(in_debug_mode) then; _newline_ \
    call push_sub(__FILE__//"."//TOSTRING(routine)); _newline_ \
  endif
#define POP_SUB(routine) \
  if(in_debug_mode) then; _newline_ \
    call pop_sub(__FILE__//"."//TOSTRING(routine)); _newline_ \
  endif

! the leading dimension of the array
#define LD(a) ubound(a,dim=1)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
