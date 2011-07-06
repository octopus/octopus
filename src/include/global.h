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

! If the compiler accepts line number markers, then "CARDINAL" will
! put them.  Otherwise, just a new line. Note that the "cardinal" and
! "newline" words are substituted by the program preprocess.pl by the
! ampersand and by a real new line just before compilation.
#if defined(F90_ACCEPTS_LINE_NUMBERS)
#    define CARDINAL \newline\cardinal __LINE__ __FILE__
#else
#    define CARDINAL \newline
#endif


! The assertions are ignored if the code is compiled in not-debug mode (NDEBUG
! is defined). Otherwise it is merely a logical assertion that, when fails,
! prints out the assertion string, the file, and the line. The subroutine
! aassert_die is in the global_m module.
#define __STRING(x)     #x
#if !defined(NDEBUG)
#  define ASSERT(expr)  \
  if(.not.(expr)) & \newline \
     call assert_die(__STRING(expr), & \newline __FILE__, & \newline  __LINE__) \
  CARDINAL
#else
#  define ASSERT(expr)
#endif


! Some compilers will not have the sizeof intrinsic.
#ifdef HAVE_FC_SIZEOF
#  define SIZEOF(x) sizeof(x)
#else
#  define SIZEOF(x) 1
#endif


! In octopus, one should normally use the SAFE_(DE)ALLOCATE macros below, which emit
! a helpful error if the the allocation or deallocation fails. The "MY_DEALLOCATE" macro
! is only used in this file; in the code, one should use SAFE_DEALLOCATE_P for pointers
! and SAFE_DEALLOCATE_A for arrays.
#if defined(NDEBUG)
#  define SAFE_ALLOCATE(x) allocate(x)
#  define SAFE_DEALLOCATE_P(x) if(associated(x)) then; deallocate(x); nullify(x); end if
#  define SAFE_DEALLOCATE_A(x) if(allocated(x)) then; deallocate(x); end if
#else
#  define SAFE_ALLOCATE(x)			\
  allocate(x, stat=global_alloc_err); \newline \
  if(iand(prof_vars%mode, PROFILING_MEMORY).ne.0 .or. global_alloc_err.ne.0) & \newline \
  global_sizeof = SIZEOF(x); \newline	\
  if(iand(prof_vars%mode, PROFILING_MEMORY).ne.0) & \newline \
    call profiling_memory_allocate(& \newline #x, & \newline __FILE__, & \newline __LINE__, & \newline global_sizeof); \newline \
  if(global_alloc_err.ne.0) & \newline \
    call alloc_error(global_sizeof, & \newline __FILE__, & \newline __LINE__); \
  CARDINAL

#  define MY_DEALLOCATE(x) \
  global_sizeof = SIZEOF(x); \newline \
  deallocate(x, stat=global_alloc_err); \newline \
  if(iand(prof_vars%mode, PROFILING_MEMORY).ne.0) & \newline \
    call profiling_memory_deallocate(#x, & \newline __FILE__, & \newline __LINE__, & \newline global_sizeof); \newline \
  if(global_alloc_err.ne.0) & \newline \
    call dealloc_error(global_sizeof, & \newline __FILE__, & \newline __LINE__); \
  CARDINAL

#  define SAFE_DEALLOCATE_P(x) \
  if(associated(x)) then; \newline \
    MY_DEALLOCATE(x);     \newline \
    nullify(x);           \newline \
  end if
#  define SAFE_DEALLOCATE_A(x) \
  if(allocated(x)) then;  \newline \
    MY_DEALLOCATE(x);     \newline \
  end if

#endif

! This was used in the past and should not be used any more.
#define ALLOCATE(a,b) _DEPRECATED_PLEASE_USE_SAFE_ALLOCATE_


! The following macros facilitate the use of real or complex variables,
! and the possibility of compiling the code in single or double precision.
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


! The code directories should be defined here, and not hard coded in the Fortran files.
#define GS_DIR "gs/"
#define STATIC_DIR "static/"
#define EM_RESP_DIR "em_resp/"
#define EM_RESP_FD_DIR "em_resp_fd/"
#define KDOTP_DIR "kdotp/"
#define VIB_MODES_DIR "vib_modes/"
#define VDW_DIR "vdw/"
#define CASIDA_DIR "casida/"
#define OCT_DIR "opt-control/"


! The MPI1 and MPI2 standards are different regarding the MPI_IN_PLACE constant. In
! the code, just use the MPI_IN_PLACE_OR defined here.
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


! Whenever a procedure is not called too many times, one should start it
! and finish it with the PUSH_SUB and POP_SUB macros, which are these
! pieces of code that call the push_sub and pop_sub routines defined
! in the messages_m module.
#define PUSH_SUB(routine) \
  if(in_debug_mode) then; \newline \
    call push_sub(__FILE__//"."//TOSTRING(routine)); \newline \
  endif
#define POP_SUB(routine) \
  if(in_debug_mode) then; \newline \
    call pop_sub(__FILE__//"."//TOSTRING(routine)); \newline \
  endif

! The leading dimension of the array
#define LD(a) ubound(a,dim=1)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
