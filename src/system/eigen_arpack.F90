!! Copyright (C) 2008 X. Andrade
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
!! $Id: eigen_rmmdiis.F90 5954 2009-10-17 20:53:52Z xavier $

#include "global.h"


module eigen_arpack_m


use batch_m
  use comm_m
  use datasets_m
  use global_m
  use grid_m
  use hamiltonian_m
  use io_m
  use lalg_basic_m
  use lalg_adv_m
  use loct_m
  use math_m
  use mesh_m
  use mesh_batch_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use mpi_lib_m
  use parser_m
  use profiling_m
  use states_m
  use states_calc_m
  use varinfo_m

  implicit none

  private
  public ::                     &
    eigen_arpack_t,             &
    deigen_solver_arpack,       &
    zeigen_solver_arpack,       &
    arpack_init

    type eigen_arpack_t
      integer          :: arnoldi_vectors !< number of Arnoldi vectors
      character(len=2) :: sort            !< which eigenvalue sorting
      integer          :: init_resid      !< inital residual strategy
      CMPLX            :: rotation        !< rotate spectrum by complex number before determining order
      logical          :: use_parpack
    end type eigen_arpack_t

  contains

  subroutine arpack_init(this, gr, nst)
    type(eigen_arpack_t),  intent(inout) :: this
    type(grid_t),          intent(in)    :: gr
    integer,               intent(in)    :: nst

    FLOAT   :: rotate_spectrum_angle
    logical :: use_parpack

    PUSH_SUB(arpack_init)
#if defined(HAVE_ARPACK)

    ! XXX yuck, this parameter is also given in Hamiltonian.  How should it be transferred from there to here
    ! without parsing the parameter again?  For now we'll just parse it again
    rotate_spectrum_angle = M_ZERO

    call parse_float(datasets_check('ComplexScalingRotateSpectrum'), M_ZERO, rotate_spectrum_angle)
    call messages_print_var_value(stdout, "ComplexScalingRotateSpectrum", rotate_spectrum_angle)

    this%rotation = exp(M_zI * rotate_spectrum_angle)
     
    use_parpack = .false.
#if defined(HAVE_PARPACK)    
    use_parpack = gr%mesh%parallel_in_domains
    this%use_parpack = use_parpack
    
    !%Variable EigensolverParpack 
    !%Type logical 
    !%Section SCF::Eigensolver 
    !%Description 
    !% Use PARPACK.
    !%End 
    call parse_logical(datasets_check('EigensolverParpack'), use_parpack, this%use_parpack) 
    call messages_print_var_value(stdout, "EigensolverParpack", this%use_parpack)
    
#endif


    
    
    !%Variable EigensolverArnoldiVectors 
    !%Type integer 
    !%Section SCF::Eigensolver 
    !%Description 
    !% This indicates how many Arnoldi vectors are generated 
    !% It must satisfy EigenSolverArnoldiVectors - Number Of Eigenvectors >= 2. 
    !% See the ARPACK documentation for more details. It will default to  
    !% twice the number of eigenvectors (which is the number of states) 
    !%End 
    call parse_integer(datasets_check('EigensolverArnoldiVectors'), 2*nst, this%arnoldi_vectors) 
    if(this%arnoldi_vectors - nst < M_TWO) call input_error('EigensolverArnoldiVectors') 
    call messages_print_var_value(stdout, "EigensolverArnoldiVectors", this%arnoldi_vectors)
    
    !%Variable EigensolverArpackSort
    !%Type string 
    !%Default SR 
    !%Section SCF::Eigensolver 
    !%Description 
    !% Eigenvalues sorting strategy (case sensitive).
    !% From ARPACK documentation: 
    !% 'LM' -> want eigenvalues of largest magnitude.
    !% 'SM' -> want eigenvalues of smallest magnitude.
    !% 'LR' -> want eigenvalues of largest real part.
    !% 'SR' -> want eigenvalues of smallest real part.
    !% 'LI' -> want eigenvalues of largest imaginary part.
    !% 'SI' -> want eigenvalues of smallest imaginary part.
    !%End 
    call parse_string(datasets_check('EigensolverArpackSort'), "SR", this%sort)
    if(this%sort /= "LM"  .and. &
       this%sort /= "SM"  .and. &
       this%sort /= "LR"  .and. &
       this%sort /= "SR"  .and. &
       this%sort /= "LI"  .and. &
       this%sort /= "SI") call input_error('EigensolverArpackSort')
    call messages_print_var_value(stdout, "EigensolverArpackSort", this%sort)
    

 
    !%Variable EigensolverArpackIntialResid
    !%Type integer
    !%Default constant 
    !%Section SCF::Eigensolver
    !%Description
    !% Initial residual vector.
    !%Option constant 2
    !% Initial residual vector constant = 1.
    !%Option rand 0
    !% Random residual vector.
    !%Option calc 1
    !% resid = H*psi - epsilon*psi.
    !%End
    call parse_integer(datasets_check('EigensolverArpackIntialResid'), 2, this%init_resid)
    if(.not.varinfo_valid_option('EigensolverArpackIntialResid', this%init_resid))&
       call input_error('EigensolverArpackIntialResid')
    call messages_print_var_option(stdout, "EigensolverArpackIntialResid", this%init_resid)

#else 

    write(message(1), '(a)') 'Eigensolver = arpack requires arpack or parpack libaries.' 
    write(message(2), '(a)') 'Provide a different EigenSolver or recompile with p/arpack support.' 
    call messages_fatal(2)     

#endif
!HAVE_ARPACK
  
    POP_SUB(arpack_init)    
  end subroutine arpack_init
  
#if defined(HAVE_ARPACK)  
  !--------------------------------------------    
  subroutine arpack_debug(debug_level)
    integer, intent(in) :: debug_level

! Modified from ARPACK debug.h 
! I don't think this is going to change too much.. well at least it didn't since 1997 :)
!
!\SCCS Information: @(#) 
! FILE: debug.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2 
!
!     %---------------------------------%
!     | See debug.doc for documentation |
!     %---------------------------------%
    integer  logfil, ndigit, mgetv0, &
             msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd, &
             mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd, &
             mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
    common /debug/ &
             logfil, ndigit, mgetv0, &
             msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd, &
             mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd, &
             mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd

    PUSH_SUB(arpack_debug)
      
    ndigit = -3
    logfil = 6
    mnaitr = 0
    mnapps = 0
    mnaupd = debug_level + 1
    mnaup2 = debug_level + 1
    mneigh = 0
    mneupd = debug_level + 1
    mngets = debug_level + 1
    
    mcaupd = debug_level + 1
    mcaup2 = debug_level + 1
    mceigh = 0
    mceupd = debug_level + 1
    mcgets = debug_level + 1
    
    POP_SUB(arpack_debug)
  end subroutine arpack_debug
  
    
  !----------------------------------------------------
  subroutine arpack_check_error(sub, info)
    integer,           intent(in) :: info
    character(len= *), intent(in) :: sub
    
    integer :: msg_lines
    logical :: OK
    
    PUSH_SUB(arpack_check_error)
    
    msg_lines = 1
    OK = .false.
  
    if (sub == 'neupd') then
      
      select case (info)
        case (0)
          OK = .true.
          
        case (1)
          write(message(2),'(a)'), 'The Schur form computed by LAPACK routine csheqr'
          write(message(3),'(a)'), 'could not be reordered by LAPACK routine ztrsen.'
          write(message(4),'(a)'), 'Re-enter subroutine pzneupd with IPARAM(5)=NCV and'
          write(message(5),'(a)'), 'increase the size of the array D to have'
          write(message(6),'(a)'), 'dimension at least dimension NCV and allocate at least NCV'
          write(message(7),'(a)'), 'columns for Z. NOTE: Not necessary if Z and V share'
          write(message(8),'(a)'), 'the same space. Please notify the authors if this error'
          write(message(9),'(a)'), 'occurs.'
          msg_lines = 9
         
        case(-1) 
          write(message(2),'(a)'), 'N must be positive.'
          msg_lines = 2
                         
        case(-2)
          write(message(2),'(a)'), 'NEV must be positive.'
          msg_lines = 2
          
        case(-3)
          write(message(2),'(a)'), 'NCV-NEV >= 2 and less than or equal to N.'
          msg_lines = 2
          
        case(-5)
          write(message(2),'(a)'), 'WHICH must be one of "LM", "SM", "LR", "SR", "LI", "SI"'
          msg_lines = 2
          
        case(-6)
          write(message(2),'(a)'), 'BMAT must be one of "I" or "G".'
          msg_lines = 2
          
        case(-7)
          write(message(2),'(a)'), 'Length of private work WORKL array is not sufficient.'
          msg_lines = 2
          
        case(-8)
          write(message(2),'(a)'), 'Error return from LAPACK eigenvalue calculation.'
          write(message(3),'(a)'), 'This should never happened.'
          msg_lines = 3
          
        case(-9)
          write(message(2),'(a)'), 'Error return from calculation of eigenvectors.'
          write(message(3),'(a)'), 'Informational error from LAPACK routine ztrevc.'
          msg_lines = 3
                               
        case(-10)
          write(message(2),'(a)'), 'IPARAM(7) must be 1,2,3'
          msg_lines = 2
                        
        case(-11)
          write(message(2),'(a)'), 'PARAM(7) = 1 and BMAT = "G" are incompatible.'
          msg_lines = 2       
                            
        case(-12)
          write(message(2),'(a)'), 'HOWMNY = "S" not yet implemented'
          msg_lines = 2
          
        case(-13)
          write(message(2),'(a)'), 'OWMNY must be one of "A" or "P" if RVEC = .true.'
          msg_lines = 2
          
        case(-14)
          write(message(2),'(a)'), 'PZNAUPD did not find any eigenvalues to sufficient'
          write(message(3),'(a)'), 'accuracy.'
          msg_lines = 3
          
        case(-15)
          write(message(2),'(a)'), 'ZNEUPD got a different count of the number of converged'
          write(message(3),'(a)'), 'Ritz values than ZNAUPD got.  This indicates the user'
          write(message(4),'(a)'), 'probably made an error in passing data from ZNAUPD to'
          write(message(5),'(a)'), 'ZNEUPD or that the data was modified before entering'
          write(message(6),'(a)'), 'ZNEUPD.'
          msg_lines = 6          
          
      end select
      
    else if( sub == 'naupd') then
  
      select case (info)
        case (0)
          OK = .true.
          
        case (1)
          write(message(2),'(a)'), 'Maximum number of iterations taken.'
          write(message(3),'(a)'), 'All possible eigenvalues of OP has been found. IPARAM(5)'
          write(message(4),'(a)'), 'returns the number of wanted converged Ritz values.'
          msg_lines = 4
          OK = .true.
          
        case (2)        
          write(message(2),'(a)'), 'No longer an informational error. Deprecated starting'
          write(message(3),'(a)'), 'with release 2 of ARPACK.'
          msg_lines = 3
          OK = .true.
          
        case (3)
          write(message(2),'(a)'), 'No shifts could be applied during a cycle of the'
          write(message(3),'(a)'), 'Implicitly restarted Arnoldi iteration. One possibility'
          write(message(4),'(a)'), 'is to increase the size of NCV relative to NEV.'
          write(message(5),'(a)'), 'See remark 4 below.'
          msg_lines = 5
          OK = .true.
                
        case (-1)
           write(message(2),'(a)'), 'N must be positive.'
           msg_lines = 2
           
        case (-2)       
           write(message(2),'(a)'), 'NEV must be positive.'
           msg_lines = 2
           
        case (-3)
           write(message(2),'(a)'), 'NCV-NEV >= 2 and less than or equal to N.'
           msg_lines = 2
           
        case (-4)
           write(message(2),'(a)'), 'The maximum number of Arnoldi update iteration'          
           write(message(3),'(a)'), 'must be greater than zero.'
           msg_lines = 3
                
        case (-5)
           write(message(2),'(a)'), 'WHICH must be one of "LM", "SM", "LR", "SR", "LI", "SI"'
           msg_lines = 2
           
        case (-6)
           write(message(2),'(a)'), 'BMAT must be one of "I" or "G".'
           msg_lines = 2
           
        case (-7)
           write(message(2),'(a)'), 'Length of private work array is not sufficient.'
           msg_lines = 2
           
        case (-8)
           write(message(2),'(a)'), 'Error return from LAPACK eigenvalue calculation;'
           msg_lines = 2
           
        case (-9)
           write(message(2),'(a)'), 'Starting vector is zero.'
           msg_lines = 2
           
        case (-10)
           write(message(2),'(a)'), 'IPARAM(7) must be 1,2,3.'
           msg_lines = 2
           
        case (-11)
           write(message(2),'(a)'), 'IPARAM(7) = 1 and BMAT = "G" are incompatable.'
           msg_lines = 2
           
        case (-12)
           write(message(2),'(a)'), 'IPARAM(1) must be equal to 0 or 1.'
           msg_lines = 2
           
        case (-9999)
           write(message(2),'(a)'), 'Could not build an Arnoldi factorization.'
           write(message(3),'(a)'), 'User input error highly likely.  Please'
           write(message(4),'(a)'), 'check actual array dimensions and layout.'
           write(message(5),'(a)'), 'IPARAM(5) returns the size of the current Arnoldi'
           write(message(6),'(a)'), 'factorization.'
           msg_lines = 6
          
      end select
        
        
    else
      write(message(1),'(a)') 'Unrecognized arpack subroutine '
      call messages_fatal(1)  
      
    end if
 
    if(.not. OK) then
      write(message(1),'(a,a,a,i5)') 'Error with P/ARPACK ', sub, ', info = ', info
      call messages_fatal(msg_lines)
    else if(msg_lines >= 2) then      
      write(message(1),'(a,a,a,i5)') 'P/ARPACK ',sub, ', info = ', info
      call messages_warning(msg_lines)
    end if
    
    POP_SUB(arpack_check_error)
  end subroutine arpack_check_error

#endif
!HAVE_ARPACK  

#include "real.F90" 
#include "eigen_arpack_inc.F90" 
#include "undef.F90" 

#include "complex.F90" 
#include "eigen_arpack_inc.F90" 
#include "undef.F90" 


  end module eigen_arpack_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
