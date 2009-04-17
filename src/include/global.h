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

!#ifndef HAVE_FC_SIZEOF
#define sizeof(x) 8.0_8
!#endif

#define MAX_SPIN 4

#if defined(F90_ACCEPTS_LINE_NUMBERS)
#  define CARDINAL \newline\cardinal __LINE__ __FILE__
#else
#  define CARDINAL \newline
#endif

#define __STRING(x)     #x

#if !defined(NDEBUG)
#  if defined(LONG_LINES)
#    define ASSERT(expr) if(.not.(expr)) \
       call assert_die(__STRING(expr), __FILE__, __LINE__)
#  else
#    define ASSERT(expr) \
       if(.not.(expr)) then                    \newline \
         call assert_die (__STRING(expr), &    \newline \
                          __FILE__, __LINE__)  \newline \
       end if                                  \
       CARDINAL
#  endif
#else
#  define ASSERT(expr)
#endif

#if defined(NDEBUG)
#  define ALLOCATE(x, size) allocate(x)
#  define DEALLOCATE(x) deallocate(x)
#  define SAFE_DEALLOCATE_P(x) if(associated(x)) then; DEALLOCATE(x); nullify(x); end if
#  define SAFE_DEALLOCATE_A(x) if(allocated(x)) then; DEALLOCATE(x); end if
#else

#if defined(LONG_LINES)
#  define ALLOCATE(x, size) allocate(x, stat=global_alloc_err); \
     if(profiling_space .and. size > 0) call profiling_memory_allocate(#x, __FILE__, __LINE__, (size)*dble(sizeof(x))); \
     if(global_alloc_err.ne.0) call alloc_error((size), __FILE__, __LINE__)

#  define DEALLOCATE(x) deallocate(x, stat=global_alloc_err); \
     if(profiling_space) call profiling_memory_deallocate(#x, __FILE__, __LINE__); \
     if(global_alloc_err.ne.0) call alloc_error(-1, __FILE__, __LINE__)

#  define SAFE_DEALLOCATE_P(x) if(associated(x)) then; DEALLOCATE(x); nullify(x); end if
#  define SAFE_DEALLOCATE_A(x) if(allocated(x)) then; DEALLOCATE(x); end if

#else
#  define ALLOCATE(x, size) \
     allocate(x, stat=global_alloc_err)                 \newline \
       if(profiling_space .and. &                       \newline \
       size > 0) then                                   \newline \
       call profiling_memory_allocate(&                 \newline \
  #x,       &                                           \newline \
  __FILE__, &                                           \newline \
  __LINE__, &                                           \newline \
     (size)*&                                             \newline \
     dble(sizeof(x)))                                   \newline \
     end if                                             \newline \
     if(global_alloc_err.ne.0) then                     \newline \
       call alloc_error((size),  &                      \newline \
  __FILE__, &                                           \newline \
  __LINE__)                                             \newline \
     end if                                             \
     CARDINAL

#  define DEALLOCATE(x) \
     deallocate(x, stat=global_alloc_err)               \newline \
     if(profiling_space) then                           \newline \
       call profiling_memory_deallocate(&               \newline \
  #x,       &                                           \newline \
  __FILE__, &                                           \newline \
  __LINE__)                                             \newline \
     end if                                             \newline \
     if(global_alloc_err.ne.0) then                     \newline \
       call alloc_error(-1, &                           \newline \
  __FILE__, &                                           \newline \
  __LINE__)                                             \newline \
     end if                                             \
     CARDINAL

#  define SAFE_DEALLOCATE_P(x) if(associated(x)) then \newline DEALLOCATE(x) \newline nullify(x) \newline end if
#  define SAFE_DEALLOCATE_A(x) if(allocated(x)) then \newline DEALLOCATE(x) \newline end if

#endif
#endif


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

#define STATIC_DIR "static"
#define EM_RESP_RESTART_DIR "em_resp/"
#define KDOTP_RESTART_DIR "kdotp/"
#define PHONONS_RESTART_DIR "vib_modes/"
#define VDW_RESTART_DIR "pol_lr/"
