!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

#include "global.h"

module ode_solver
  use global
  use messages
  use datasets_mod
  use lib_oct_parser

  implicit none

  private
  public ::                     &
    ode_solver_type,            &
    dode_solver_init,           &
    dode_solver_create,         &
    dode_solver_run,            &
    dode_solver_end,            &
    zode_solver_init,           &
    zode_solver_create,         &
    zode_solver_run,            &
    zode_solver_end

  integer, public, parameter :: &
    ODE_RK4    =  1,            &
    ODE_FB78   =  2,            &
    ODE_VR89   =  3,            &
    ODE_PD89   =  4,            &
    ODE_MINVAL =  ODE_RK4,      &
    ODE_MAXVAL =  ODE_PD89

  type ode_solver_type
    integer :: solver_type     ! what solver to use (see ODE_* variables above)
    integer :: nsteps          ! how many steps to use
    integer :: nsize           ! how many odes to solve simultaneously
    integer :: vsize           ! vector size of ode method (used internally)
    FLOAT   :: tmax, tmin      ! integrate ODE from tmin to tmax
    logical :: adaptive_steps  ! should we use adaptive steps?
    logical :: full_solution   ! if true the solution will be returned for all t, otherwise only at the endpoint t=tmax.
    FLOAT, pointer :: a(:,:), b(:), c(:), e(:) ! coefficients dor the ode solver
  end type ode_solver_type


contains

  ! ---------------------------------------------------------
  ! coefficients for standard Runge-Kutta 4.th order
  subroutine ode_rk4_coeff(os)
    type(ode_solver_type), intent(inout)  :: os

    call push_sub('ode_solver.ode_rk4_coeff')

    os%b = (/ CNST(1.0) / CNST(6.0), &
      CNST(  1.0) / CNST(3.0),       &
      CNST(  1.0) / CNST(3.0),       &
      CNST(  1.0) / CNST(6.0) /)

    os%c = (/ CNST(0.0),             &
      CNST(  1.0) / CNST(2.0),       &
      CNST(  1.0) / CNST(2.0),       &
      CNST(  1.0) /)

    os%a = M_ZERO

    os%a(2,1:1) = (/ CNST(1.0) / CNST(2.0) /)

    os%a(3,1:2) = (/ CNST(0.0),  CNST(1.0) / CNST(2.0) /)

    os%a(4,1:3) = (/ CNST(0.0),  CNST(0.0),  CNST(1.0) /)

    call pop_sub()
  end subroutine ode_rk4_coeff


  ! ---------------------------------------------------------
  ! coefficients for Fehlberg 7/8.th order
  subroutine ode_fb78_coeff(os)
    type(ode_solver_type), intent(inout)  :: os

    call push_sub('ode_solver.ode_fb78_coeff')

    os%b = (/ CNST(41.0) / CNST( 840.0), &
      CNST(   0.0),                      &
      CNST(   0.0),                      &
      CNST(   0.0),                      &
      CNST(   0.0),                      &
      CNST(  34.0) / CNST( 105.0),       &
      CNST(   9.0) / CNST(  35.0),       &
      CNST(   9.0) / CNST(  35.0),       &
      CNST(   9.0) / CNST( 280.0),       &
      CNST(   9.0) / CNST( 280.0),       &
      CNST(  41.0) / CNST( 840.0),       &
      CNST(   0.0),                      &
      CNST(   0.0) /)

    os%e = (/ CNST(0.0),                 &
      CNST(   0.0),                      &
      CNST(   0.0),                      &
      CNST(   0.0),                      &
      CNST(   0.0),                      &
      CNST(  34.0) / CNST( 105.0),       &
      CNST(   9.0) / CNST(  35.0),       &
      CNST(   9.0) / CNST(  35.0),       &
      CNST(   9.0) / CNST( 280.0),       &
      CNST(   9.0) / CNST( 280.0),       &
      CNST(   0.0),                      &
      CNST(  41.0) / CNST( 840.0),       &
      CNST(  41.0) / CNST( 840.0) /)

    os%c = (/ CNST( 0.0),                &
      CNST(   2.0) / CNST(  27.0),       &
      CNST(   1.0) / CNST(   9.0),       &
      CNST(   1.0) / CNST(   6.0),       &
      CNST(   5.0) / CNST(  12.0),       &
      CNST(   1.0) / CNST(   2.0),       &
      CNST(   5.0) / CNST(   6.0),       &
      CNST(   1.0) / CNST(   6.0),       &
      CNST(   2.0) / CNST(   3.0),       &
      CNST(   1.0) / CNST(   3.0),       &
      CNST(   1.0),                      &
      CNST(   0.0),                      &
      CNST(   1.0) /)

    os%a = M_ZERO

    os%a( 2,   1) = CNST(2.0) / CNST(27.0)

    os%a( 3, 1:2) = (/ CNST(1.0) / CNST(36.0),    &
      CNST(   1.0) / CNST(  12.0) /)

    os%a( 4, 1:3) = (/ CNST(1.0) / CNST(24.0),    &
      CNST(   0.0),                               &
      CNST(   1.0) / CNST(   8.0) /)


    os%a( 5, 1:4) = (/ CNST(5.0) / CNST(12.0),    &
      CNST(   0.0),                               &
      CNST( -25.0) / CNST(  16.0),                &
      CNST(  25.0) / CNST(  16.0) /)

    os%a( 6, 1:5) = (/ CNST(1.0) / CNST(20.0),    &
      CNST(   0.0),                               &
      CNST(   0.0),                               &
      CNST(   1.0) / CNST(   4.0),                &
      CNST(   1.0) / CNST(   5.0) /)

    os%a( 7, 1:6) = (/ CNST(-25.0) / CNST(108.0), &
      CNST(   0.0),                               &
      CNST(   0.0),                               &
      CNST( 125.0) / CNST( 108.0),                &
      CNST( -65.0) / CNST(  27.0),                &
      CNST( 125.0) / CNST(  54.0) /)

    os%a( 8, 1:7) = (/ CNST( 31.0) / CNST(300.0), &
      CNST(   0.0),                               &
      CNST(   0.0),                               &
      CNST(   0.0),                               &
      CNST(  61.0) / CNST( 225.0),                &
      CNST(  -2.0) / CNST(   9.0),                &
      CNST(  13.0) / CNST( 900.0) /)

    os%a( 9, 1:8) = (/ CNST(  2.0),               &
      CNST(   0.0),                               &
      CNST(   0.0),                               &
      CNST( -53.0) / CNST(   6.0),                &
      CNST( 704.0) / CNST(  45.0),                &
      CNST(-107.0) / CNST(   9.0),                &
      CNST(  67.0) / CNST(  90.0),                &
      CNST(   3.0) /)

    os%a(10, 1:9) = (/ CNST(-91.0) / CNST(108.0), &
      CNST(   0.0),                               &
      CNST(   0.0),                               &
      CNST(  23.0) / CNST( 108.0),                &
      CNST(-976.0) / CNST( 135.0),                &
      CNST( 311.0) / CNST(  54.0),                &
      CNST( -19.0) / CNST(  60.0),                &
      CNST(  17.0) / CNST(   6.0),                &
      CNST(  -1.0) / CNST(  12.0) /)

    os%a(11,1:10) = (/ CNST(2383.0) / CNST(4100.0),  &
      CNST(   0.0),                                  &
      CNST(   0.0),                                  &
      CNST(-341.0) / CNST( 164.0),                   &
      CNST(4496.0) / CNST(1025.0),                   &
      CNST(-301.0) / CNST(  82.0),                   &
      CNST(2133.0) / CNST(4100.0),                   &
      CNST(  45.0) / CNST(  82.0),                   &
      CNST(  45.0) / CNST( 164.0),                   &
      CNST(  18.0) / CNST(  41.0) /)

    os%a(12,1:11) = (/ CNST(   3.0) / CNST( 205.0),  &
      CNST(   0.0),                                  &
      CNST(   0.0),                                  &
      CNST(   0.0),                                  &
      CNST(   0.0),                                  &
      CNST(  -6.0) / CNST(  41.0),                   &
      CNST(  -3.0) / CNST( 205.0),                   &
      CNST(  -3.0) / CNST(  41.0),                   &
      CNST(   3.0) / CNST(  41.0),                   &
      CNST(   6.0) / CNST(  41.0),                   &
      CNST(   0.0) /)

    os%a(13,1:12) = (/ CNST(-1777.0) / CNST(4100.0), &
      CNST(   0.0),                                  &
      CNST(   0.0),                                  &
      CNST(-341.0) / CNST(     164.0),               &
      CNST(4496.0) / CNST(    1025.0),               &
      CNST(-289.0) / CNST(      82.0),               &
      CNST(2193.0) / CNST(    4100.0),               &
      CNST(  51.0) / CNST(      82.0),               &
      CNST(  33.0) / CNST(     164.0),               &
      CNST(  12.0) / CNST(      41.0),               &
      CNST(   0.0),                                  &
      CNST(   1.0) /)

    call pop_sub()
  end subroutine ode_fb78_coeff


  ! ---------------------------------------------------------
  ! coefficients for Verner 8/9.th order
  !
  subroutine ode_vr89_coeff(os)
    type(ode_solver_type), intent(inout)  :: os

    FLOAT :: SQRT6

    call push_sub('ode_solver.ode_vr89_coeff')


    SQRT6 = sqrt(M_SIX)

    os%b = (/ CNST(   103.0) / CNST( 1680.0), &
      CNST(   0.0),                           &
      CNST(   0.0),                           &
      CNST(   0.0),                           &
      CNST(   0.0),                           &
      CNST(   0.0),                           &
      CNST(   0.0),                           &
      CNST( -27.0) / CNST(  140.0),           &
      CNST(  76.0) / CNST(  105.0),           &
      CNST(-201.0) / CNST(  280.0),           &
      CNST(1024.0) / CNST( 1365.0),           &
      CNST(   3.0) / CNST( 7280.0),           &
      CNST(  12.0) / CNST(   35.0),           &
      CNST(   9.0) / CNST(  280.0),           &
      CNST(   0.0),                           &
      CNST(   0.0)                            &
      /)

    os%e = (/ CNST( -1911.0) / CNST( 109200.0), &
      CNST(   0.0),                             &
      CNST(   0.0),                             &
      CNST(   0.0),                             &
      CNST(   0.0),                             &
      CNST(   0.0),                             &
      CNST(   0.0),                             &
      CNST(  34398.0) / CNST(109200.0),         &
      CNST( -61152.0) / CNST(109200.0),         &
      CNST( 114660.0) / CNST(109200.0),         &
      CNST(-114688.0) / CNST(109200.0),         &
      CNST(    -63.0) / CNST(109200.0),         &
      CNST( -13104.0) / CNST(109200.0),         &
      CNST(  -3510.0) / CNST(109200.0),         &
      CNST(  39312.0) / CNST(109200.0),         &
      CNST(   6058.0) / CNST(109200.0)          &
      /)

    os%c = (/ CNST(0.0),                        &
      ( CNST(1.0) / CNST(12.0) ),               &
      ( CNST(1.0) / CNST( 9.0) ),               &
      ( CNST(1.0) / CNST( 6.0) ),               &
      ( CNST(2.0) * ( CNST( 1.0) + SQRT6 )) / CNST(15.0), &
      ( CNST(6.0) + SQRT6 ) / CNST(15.0),       &
      ( CNST(6.0) - SQRT6 ) / CNST(15.0),       &
      ( CNST(2.0) / CNST( 3.0) ),               &
      ( CNST(1.0) / CNST( 2.0) ),               &
      ( CNST(1.0) / CNST( 3.0) ),               &
      ( CNST(1.0) / CNST( 4.0) ),               &
      ( CNST(4.0) / CNST( 3.0) ),               &
      ( CNST(5.0) / CNST( 6.0) ),               &
      ( CNST(1.0) ),                            &
      ( CNST(1.0) / CNST( 6.0) ),               &
      ( CNST(1.0) )                             &
      /)


    os%a = M_ZERO

    os%a( 2,   1) = CNST(1.0) / CNST(12.0)

    os%a( 3, 1:2) = (/ CNST(1.0) / CNST(27.0),  &
      ( CNST(      2.0) / CNST(27.0) ) /)

    os%a( 4, 1:3) = (/ CNST(1.0) / CNST(24.0),  &
      ( CNST(      0.0) ),                      &
      ( CNST(      3.0) / CNST(  24.0) ) /)

    os%a( 5, 1:4) = (/ (CNST(4.0) + CNST(94.0)*SQRT6) / CNST(375.0),    &
      ( CNST(      0.0) ),                                              &
      (-CNST(    282.0) - CNST( 252.0) * SQRT6 ) / CNST( 375.0),        &
      ( CNST(    328.0) + CNST( 208.0) * SQRT6 ) / CNST( 375.0) /)

    os%a( 6, 1:5) = (/ (CNST(9.0) - SQRT6) / CNST(150.0),               &
      ( CNST(      0.0) ),                                              &
      ( CNST(      0.0) ),                                              &
      ( CNST(    312.0) + CNST(  32.0) * SQRT6 ) / CNST(1425.0),        &
      ( CNST(     69.0) + CNST(  29.0) * SQRT6 ) / CNST( 570.0) /)

    os%a( 7, 1:6) = (/ ( CNST(927.0) - CNST(347.0) * SQRT6 ) / CNST(1250.0), &
      ( CNST(      0.0) ),                                              &
      ( CNST(      0.0) ),                                              &
      ( CNST( -16248.0) + CNST(7328.0) * SQRT6 ) / CNST(9375.0),        &
      ( CNST(   -489.0) + CNST( 179.0) * SQRT6 ) / CNST(3750.0),        &
      ( CNST(  14268.0) - CNST(5798.0) * SQRT6 ) / CNST(9375.0) /)

    os%a( 8, 1:7) = (/ CNST( 4.0) / CNST(   54.0),                      &
      ( CNST(      0.0) ),                                              &
      ( CNST(      0.0) ),                                              &
      ( CNST(      0.0) ),                                              &
      ( CNST(      0.0) ),                                              &
      ( CNST(     16.0) - SQRT6 ) / CNST( 54.0),                        &
      ( CNST(     16.0) + SQRT6 ) / CNST( 54.0) /)

    os%a( 9, 1:8) = (/ CNST(38.0) / CNST(  512.0),                      &
      ( CNST(      0.0) ),                                              &
      ( CNST(      0.0) ),                                              &
      ( CNST(      0.0) ),                                              &
      ( CNST(      0.0) ),                                              &
      ( CNST(    118.0) - CNST(  23.0) * SQRT6 ) / CNST( 512.0),        &
      ( CNST(    118.0) + CNST(  23.0) * SQRT6 ) / CNST( 512.0),        &
      ( CNST(    -18.0) / CNST( 512.0) ) /)

    os%a(10, 1:9) = (/ CNST( 11.0) / CNST( 144.0),                      &
      ( CNST(      0.0) ),                                              &
      ( CNST(      0.0) ),                                              &
      ( CNST(      0.0) ),                                              &
      ( CNST(      0.0) ),                                              &
      ( CNST(    266.0) - SQRT6 ) / CNST(864.0),                        &
      ( CNST(    266.0) + SQRT6 ) / CNST(864.0),                        &
      ( CNST(     -1.0) / CNST(  16.0) ),                               &
      ( CNST(     -8.0) / CNST(  27.0) ) /)

    os%a(11,1:10) = (/ ( CNST(5034.0) - CNST( 271.0) * SQRT6 ) / CNST(61440.0), &
      ( CNST(      0.0) ),                                              &
      ( CNST(      0.0) ),                                              &
      ( CNST(      0.0) ),                                              &
      ( CNST(      0.0) ),                                              &
      ( CNST(      0.0) ),                                              &
      ( CNST(   7859.0) - CNST(1626.0) * SQRT6 ) / CNST(10240.0),       &
      ( CNST(  -2232.0) + CNST( 813.0) * SQRT6 ) / CNST(20480.0),       &
      ( CNST(   -594.0) + CNST( 271.0) * SQRT6 ) / CNST(  960.0),       &
      ( CNST(    657.0) - CNST( 813.0) * SQRT6 ) / CNST( 5120.0) /)

    os%a(12,1:11) = (/ ( CNST(5996.0) - CNST(3794.0) * SQRT6 ) / CNST(405.0), &
      ( CNST(      0.0) ),                                              &
      ( CNST(      0.0) ),                                              &
      ( CNST(      0.0) ),                                              &
      ( CNST(      0.0) ),                                              &
      ( CNST(  -4342.0) - CNST(   338.0) * SQRT6 ) / CNST(  9.0),       &
      ( CNST( 154922.0) - CNST( 40458.0) * SQRT6 ) / CNST(135.0),       &
      ( CNST(  -4176.0) + CNST(  3794.0) * SQRT6 ) / CNST( 45.0),       &
      ( CNST(-340864.0) + CNST(242816.0) * SQRT6 ) / CNST(405.0),       &
      ( CNST(  26304.0) - CNST( 15176.0) * SQRT6 ) / CNST( 45.0),       &
      ( CNST( -26624.0) / CNST(    81.0) ) /)

    os%a(13,1:12) = (/ ( CNST( 3793.0) + CNST(2168.0) * SQRT6 ) / CNST(103680.0), &
      ( CNST(      0.0) ),                                              &
      ( CNST(      0.0) ),                                              &
      ( CNST(      0.0) ),                                              &
      ( CNST(      0.0) ),                                              &
      ( CNST(   4042.0) + CNST(  2263.0) * SQRT6 ) / CNST(13824.0),     &
      ( CNST(-231278.0) + CNST( 40717.0) * SQRT6 ) / CNST(69120.0),     &
      ( CNST(   7947.0) - CNST(  2168.0) * SQRT6 ) / CNST(11520.0),     &
      ( CNST(   1048.0) - CNST(   542.0) * SQRT6 ) / CNST(  405.0),     &
      ( CNST(  -1383.0) + CNST(   542.0) * SQRT6 ) / CNST(  720.0),     &
      ( CNST(   2624.0) / CNST(  1053.0) ) ,                            &
      ( CNST(      3.0) /  CNST( 1664.0) )  /)

    os%a(14,1:13) = (/ CNST( -137.0) / CNST(1296.0),                    &
      ( CNST(      0.0) ),                                              &
      ( CNST(      0.0) ),                                              &
      ( CNST(      0.0) ),                                              &
      ( CNST(      0.0) ),                                              &
      ( CNST(   5642.0) - CNST(   337.0) * SQRT6 ) / CNST(864.0),       &
      ( CNST(   5642.0) + CNST(   337.0) * SQRT6 ) / CNST(864.0),       &
      ( CNST(   -299.0) / CNST(    48.0) ),                             &
      ( CNST(    184.0) / CNST(    81.0) ),                             &
      ( CNST(    -44.0) / CNST(     9.0) ),                             &
      ( CNST(  -5120.0) / CNST(  1053.0) ),                             &
      ( CNST(    -11.0) / CNST(   468.0) ),                             &
      ( CNST(     16.0) / CNST(     9.0) ) /)

    os%a(15,1:14) = (/ ( CNST( 33617.0) - CNST(2168.0) * SQRT6 ) / CNST(518400.0), &
      ( CNST(      0.0) ),                                              &
      ( CNST(      0.0) ),                                              &
      ( CNST(      0.0) ),                                              &
      ( CNST(      0.0) ),                                              &
      ( CNST(  -3846.0) + CNST(    31.0) * SQRT6 ) / CNST( 13824.0),    &
      ( CNST( 155338.0) - CNST( 52807.0) * SQRT6 ) / CNST(345600.0),    &
      ( CNST( -12537.0) + CNST(  2168.0) * SQRT6 ) / CNST( 57600.0),    &
      ( CNST(     92.0) + CNST(   542.0) * SQRT6 ) / CNST(  2025.0),    &
      ( CNST(  -1797.0) - CNST(   542.0) * SQRT6 ) / CNST(  3600.0),    &
      ( CNST(    320.0) / CNST(   567.0) ),                             &
      ( CNST(     -1.0) / CNST(  1920.0) ),                             &
      ( CNST(      4.0) / CNST(   105.0) ),                             &
      ( CNST(      0.0) ) /)

    os%a(16,1:15) = (/ ( CNST( -36487.0) - CNST(30352.0) * SQRT6 ) / CNST(279600.0), &
      ( CNST(      0.0) ),                                              &
      ( CNST(      0.0) ),                                              &
      ( CNST(      0.0) ),                                              &
      ( CNST(      0.0) ),                                              &
      ( CNST( -29666.0) - CNST(  4499.0) * SQRT6 ) / CNST(  7456.0),    &
      ( CNST(2779182.0) - CNST(615973.0) * SQRT6 ) / CNST(186400.0),    &
      ( CNST( -94329.0) + CNST( 91056.0) * SQRT6 ) / CNST( 93200.0),    &
      ( CNST(-232192.0) + CNST(121408.0) * SQRT6 ) / CNST( 17475.0),    &
      ( CNST( 101226.0) - CNST( 22764.0) * SQRT6 ) / CNST(  5825.0),    &
      ( CNST(-169984.0) / CNST(  9087.0) ),                             &
      ( CNST(    -87.0) / CNST( 30290.0) ),                             &
      ( CNST(    492.0) / CNST(  1165.0) ),                             &
      ( CNST(      0.0) ),                                              &
      ( CNST(   1260.0) / CNST(   233.0) ) /)

    call pop_sub()
  end subroutine ode_vr89_coeff


  ! ---------------------------------------------------------
  ! coefficients for Prince-Dormand 8/9.th order
  !
  subroutine ode_pd89_coeff(os)
    type(ode_solver_type), intent(inout)  :: os

    call push_sub('ode_solver.ode_pd89_coeff')


    os%b = (/ CNST( 14005451.0) / CNST( 335480064.0), &
      CNST(          0.0),                            &
      CNST(          0.0),                            &
      CNST(          0.0),                            &
      CNST(          0.0),                            &
      CNST(  -59238493.0) / CNST(1068277825.0),       &
      CNST(  181606767.0) / CNST( 758867731.0),       &
      CNST(  561292985.0) / CNST( 797845732.0),       &
      CNST(-1041891430.0) / CNST(1371343529.0),       &
      CNST(  760417239.0) / CNST(1151165299.0),       &
      CNST(  118820643.0) / CNST( 751138087.0),       &
      CNST( -528747749.0) / CNST(2220607170.0),       &
      CNST(          1.0) / CNST(         4.0)        &
      /)

    os%e = (/ CNST( 13451932.0) / CNST( 455176623.0), &
      CNST(          0.0),                            &
      CNST(          0.0),                            &
      CNST(          0.0),                            &
      CNST(          0.0),                            &
      CNST(  -808719846.0) / CNST(  976000145.0),     &
      CNST(  1757004468.0) / CNST( 5645159321.0),     &
      CNST(   656045339.0) / CNST(  265891186.0),     &
      CNST( -3867574721.0) / CNST( 1518517206.0),     &
      CNST(   465885868.0) / CNST(  322736535.0),     &
      CNST(    53011238.0) / CNST(  667516719.0),     &
      CNST(           2.0) / CNST(         45.0),     &
      CNST(           0.0)                            &
      /)

    os%c = (/ CNST(        0.0),                      &
      CNST(          1.0) / CNST(        18.0),       &
      CNST(          1.0) / CNST(        12.0),       &
      CNST(          1.0) / CNST(         8.0),       &
      CNST(          5.0) / CNST(        16.0),       &
      CNST(          3.0) / CNST(         8.0),       &
      CNST(         59.0) / CNST(       400.0),       &
      CNST(         93.0) / CNST(       200.0),       &
      CNST( 5490023248.0) / CNST(9719169821.0),       &
      CNST(         13.0) / CNST(        20.0),       &
      CNST( 1201146811.0) / CNST(1299019798.0)        &
      /)

    os%a = M_ZERO

    os%a( 2,   1) = CNST(1.0) / CNST(18.0)

    os%a( 3, 1:2) = (/ CNST(1.0) / CNST(48.0),        &
      CNST(  1.0) / CNST(16.0) /)

    os%a( 4, 1:3) = (/ CNST(1.0) / CNST(32.0),        &
      CNST(  0.0),                                    &
      CNST(  3.0) / CNST(32.0) /)


    os%a( 5, 1:4) = (/ CNST(5.0) / CNST(16.0),        &
      CNST(  0.0),                                    &
      CNST(-75.0) / CNST(64.0),                       &
      CNST( 75.0) / CNST(64.0) /)

    os%a( 6, 1:5) = (/ CNST(3.0) / CNST(80.0),        &
      CNST(  0.0),                                    &
      CNST(  0.0),                                    &
      CNST(  3.0) / CNST(16.0),                       &
      CNST(  3.0) / CNST(20.0) /)

    os%a( 7, 1:6) = (/ CNST(29443841.0) / CNST(614563906.0), &
      CNST(         0.0),                                    &
      CNST(         0.0),                                    &
      CNST(  77736538.0) / CNST( 692538347.0),               &
      CNST( -28693883.0) / CNST(1125000000.0),               &
      CNST(  23124283.0) / CNST(1800000000.0) /)

    os%a( 8, 1:7) = (/ CNST(16016141.0) / CNST(946692911.0), &
      CNST(         0.0),                                    &
      CNST(         0.0),                                    &
      CNST(  61564180.0) / CNST( 158732637.0),               &
      CNST(  22789713.0) / CNST( 633445777.0),               &
      CNST( 545815736.0) / CNST(2771057229.0),               &
      CNST(-180193667.0) / CNST(1043307555.0) /)

    os%a( 9, 1:8) = (/ CNST(39632708.0) / CNST(573591083.0), &
      CNST(         0.0),                                    &
      CNST(         0.0),                                    &
      CNST(-433636366.0) / CNST( 683701615.0),               &
      CNST(-421739975.0) / CNST(2616292301.0),               &
      CNST( 100302831.0) / CNST( 723423059.0),               &
      CNST( 790204164.0) / CNST( 839813087.0),               &
      CNST( 800635310.0) / CNST(3783071287.0) /)

    os%a(10, 1:9) = (/ CNST( 246121993.0) / CNST(1340847787.0), &
      CNST(           0.0),                                     &
      CNST(           0.0),                                     &
      CNST(-37695042795.0) / CNST(15268766246.0),               &
      CNST(  -309121744.0) / CNST( 1061227803.0),               &
      CNST(   -12992083.0) / CNST(  490766935.0),               &
      CNST(  6005943493.0) / CNST( 2108947869.0),               &
      CNST(   393006217.0) / CNST( 1396673457.0),               &
      CNST(   123872331.0) / CNST( 1001029789.0) /)

    os%a(11,1:10) = (/ CNST(-1028468189.0) / CNST(846180014.0), &
      CNST(           0.0),                                     &
      CNST(           0.0),                                     &
      CNST(  8478235783.0) / CNST(  508512852.0),               &
      CNST(  1311729495.0) / CNST( 1432422823.0),               &
      CNST(-10304129995.0) / CNST( 1701304382.0),               &
      CNST(-48777925059.0) / CNST( 3047939560.0),               &
      CNST( 15336726248.0) / CNST( 1032824649.0),               &
      CNST(-45442868181.0) / CNST( 3398467696.0),               &
      CNST(  3065993473.0) / CNST(  597172653.0) /)

    os%a(12,1:11) = (/ CNST(  185892177.0) / CNST(718116043.0), &
      CNST(           0.0),                                     &
      CNST(           0.0),                                     &
      CNST( -3185094517.0) / CNST(  667107341.0),               &
      CNST(  -477755414.0) / CNST( 1098053517.0),               &
      CNST(  -703635378.0) / CNST(  230739211.0),               &
      CNST(  5731566787.0) / CNST( 1027545527.0),               &
      CNST(  5232866602.0) / CNST(  850066563.0),               &
      CNST( -4093664535.0) / CNST(  808688257.0),               &
      CNST(  3962137247.0) / CNST( 1805957418.0),               &
      CNST(    65686358.0) / CNST(  487910083.0)  /)

    os%a(13,1:12) = (/ CNST(  403863854.0) / CNST(491063109.0), &
      CNST(           0.0),                                     &
      CNST(           0.0),                                     &
      CNST( -5068492393.0) / CNST(  434740067.0),               &
      CNST(  -411421997.0) / CNST(  543043805.0),               &
      CNST(   652783627.0) / CNST(  914296604.0),               &
      CNST( 11173962825.0) / CNST(  925320556.0),               &
      CNST(-13158990841.0) / CNST( 6184727034.0),               &
      CNST(  3936647629.0) / CNST( 1978049680.0),               &
      CNST(  -160528059.0) / CNST(  685178525.0),               &
      CNST(   248638103.0) / CNST( 1413531060.0),               &
      CNST(           0.0) /)

    call pop_sub()
  end subroutine ode_pd89_coeff


#include "undef.F90"
#include "complex.F90"
#include "ode_solver_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "ode_solver_inc.F90"

end module ode_solver
