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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! /* Calculates exp(factor*in), and places it into out, where in and out
! are order x order matrices. Useful only for small matrices.
!
! It uses two methods; (i) decomposition of A=in into U(+)DU [where U is 
! unitary and D is diagonal] for hermitian matrices, and making use of the
! identity exp(A) = U(+)exp(D)exp(U), and (ii) polynomial expansion after
! proper scaling and squaring. This is simple minded, but probably enough
! for the cases octopus needs. I have learnt from C. Moler and C. Van Loan,
! SIAM Rev. 45, 3 (2003).
!
! This subroutine is just a driver to select the method, by calling
! zmatexp_scaleandsquare or zmatexp_decomposition.
!
! order : the order of the matrix.
! in    : input matrix.
! out   : output matrix, out = exp(factor*in) if subroutine.
! factor: the factor to multiply in before exponentiation.
! norm  : (optional) The euclidean norm of the matrix, or at least a reasonable
!         approximation. It is only needed by the scale and square method. If
!         this argument is not passed, the estimation sum(abs(in(:, :))) will
!         be used.
! method: on input, should be 1 for the scale and square method, and 2 for the
!         decomposition method. If this argument is not passed, it defaults to
!         1. */
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine X(matexp) (order, in, out, factor, norm, method)
  integer, intent(in)            :: order
  R_TYPE, intent(in)        :: in(order, order)
  R_TYPE, intent(out)       :: out(order, order)
  R_TYPE, intent(in)        :: factor
  FLOAT, intent(in), optional :: norm
  integer, intent(in), optional  :: method

  integer, parameter :: SCALEANDSQUARE = 1, &
                        DECOMPOSITION  = 2
  integer :: lmethod
  FLOAT :: lnorm

  lmethod = SCALEANDSQUARE
  if(present(method)) lmethod = method

  select case(lmethod)
  case(SCALEANDSQUARE)
     if(present(norm)) then
        lnorm = abs(factor)*norm
     else
        lnorm = abs(factor)*sum(abs(in(:, :)))
     endif
     if(lnorm <= M_ZERO) then
        write(message(1),'(a)') 'Internal [matexp]: invalid "norm" variable value'
        call write_fatal(1)
     endif
     call X(matexp_scaleandsquare)(order, in, out, factor, lnorm)
  case(DECOMPOSITION)
     call X(matexp_decomposition)(order, in, out, factor)
  case default
    write(message(1),'(a,i5,a)') 'Internal [matexp]: "method" input argument,', lmethod,', not valid.'
    call write_fatal(1)
  end select

end subroutine X(matexp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! /* Calculates exp(factor*in), and places it into out, where in and out
! are order x order complex matrices. Useful only for small matrices.
!
! It just computes the Taylor expansion approximation to the exponential,
! up to order exporder. This is very unreliable and ineffective as soon as
! the norm of the matrix is not very small. So it should be used only after
! proper scaling and squaring.
!
! order    : the order of the matrix.
! in       : input matrix.
! out      : output matrix, out = exp(factor*in) if subroutine.
! factor   : the factor to multiply in before exponentiation.
! exporder : the order of the Taylor expansion to be used. */ 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine X(matexp_polynomial)(order, in, out, factor, exporder)
  integer, intent(in)      :: order, exporder
  R_TYPE, intent(in)  :: in(:, :)
  R_TYPE, intent(out) :: out(:, :)
  R_TYPE, intent(in)  :: factor 

  integer :: n
  R_TYPE :: zfact
  R_TYPE, allocatable :: aux(:, :), dd(:, :)

  allocate(aux(order, order), dd(order, order))
  out = R_TOTYPE(M_ZERO)
  do n = 1, order
     out(n, n) = M_ONE
  enddo
  zfact = M_ONE
  call X(copy) (order**2, out(1, 1), 1, aux(1, 1), 1)
  do n = 1, exporder
     call X(copy) (order**2, aux(1, 1), 1, dd(1, 1), 1)
     call X(gemm) ('n', 'n', order, order, order, R_TOTYPE(M_ONE), dd(1, 1), &
                   order, in(1, 1), order, R_TOTYPE(M_ZERO), aux(1 ,1), order)
     zfact = zfact*factor/n
     call X(axpy)(order**2, zfact, aux(1, 1), 1, out(1, 1), 1)
  enddo
  deallocate(aux, dd)

end subroutine X(matexp_polynomial)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! /* Calculates exp(factor*in), and places it into out, where in and out
! are order x order matrices. Useful only for small matrices.
!
! It uses the scale and square method (see C. Moler and C. Van Loan, SIAM Rev.
! 45, 3 (2003). The exponentiation of the scaled matrix is done through a
! 12th order Taylor expansion, via zmatexp_polynomial.
!
! order : the order of the matrix.
! in    : input matrix.
! out   : output matrix, out = exp(factor*in) if subroutine.
! factor: the factor to multiply in before exponentiation.
! norm  : (optional) The euclidean norm of the matrix, (times factor), 
!         or at least a reasonable approximation */
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine X(matexp_scaleandsquare)(order, in, out, factor, norm)
  integer, intent(in)  :: order
  R_TYPE, intent(in)   :: in(:, :)
  R_TYPE, intent(out)  :: out(:, :)
  R_TYPE, intent(in)   :: factor 
  FLOAT, intent(in) :: norm

  integer :: i, j
  R_TYPE, allocatable :: aux(:, :)

  j = max(int(log(norm)/log(M_TWO)) + 1, 0)

  allocate(aux(order, order))

  call X(matexp_polynomial)(order, in, out, factor/2**j, 12)
  
  do i = 1, j
     call X(copy)(order**2, out(1, 1), 1, aux(1, 1), 1)
     call X(gemm)('n', 'n', order, order, order, R_TOTYPE(M_ONE), aux(1, 1), &
                   order, aux(1, 1), order, R_TOTYPE(M_ZERO), out(1, 1), order)
  enddo

  deallocate(aux)  
end subroutine X(matexp_scaleandsquare)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! /* Calculates exp(factor*in), and places it into out, where in and out
! are order x order complex matrices. Useful only for small *hermitian* matrices.
!
! It makes a decomposition of A=in into A = U(+)DU [where U is 
! unitary and D is diagonal] and makes use of the
! identity exp(A) = U(+)exp(D)exp(U). This is only meaningful if A is hermitian.
!
! order : the order of the matrix.
! in    : input matrix.
! out   : output matrix, out = exp(factor*in) if subroutine.
! factor: the factor to multiply in before exponentiation. */
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine X(matexp_decomposition)(order, in, out, factor)
  integer, intent(in)  :: order
  R_TYPE,  intent(in)  :: in(order, order)
  R_TYPE,  intent(out) :: out(order, order)
  R_TYPE,  intent(in)  :: factor

  integer ::  n, info, lwork
  R_TYPE :: zfact
  R_TYPE, allocatable :: aux(:,:), dd(:,:)
  FLOAT, allocatable :: w(:)

  allocate(aux(order, order), dd(order, order), w(order))
  call X(iagonalise) (order, in, aux, w)

  dd = M_z0
  do n = 1, order
     dd(n, n) = exp(factor*w(n))
  enddo

  call X(gemm)('n', 'c', order, order, order, R_TOTYPE(M_ONE), &
       dd(1, 1),  order, aux(1, 1), order, R_TOTYPE(M_ZERO), out(1, 1), order)
  call X(copy)(order**2, out(1, 1), 1, dd(1, 1), 1)
  call X(gemm)('n', 'n', order, order, order, R_TOTYPE(M_ONE), &
       aux(1, 1), order, dd(1, 1),  order, R_TOTYPE(M_ZERO), out(1, 1), order) 

  deallocate(aux, dd, w)
end subroutine X(matexp_decomposition)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine performs polynomial extrapolation (interpolation) of orders
! (1-4). It assumes:
!                      f(t) = sum_{j=0}^{order} a_j * t^k,
! and the inputs are f(0), f(-dt), f(-2*dt), ..., f(-order*dt).
!
! (I know there is a smarter and more general way to set the coefficients,
! but for the moment this works. If someday higuer orders are needed, or I am
! bored, I will put a more general formula)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine X(extrapolate)(order, n, v, vex, dt, t)
  integer, intent(in)  :: order, n
  R_TYPE, intent(in)   :: v(*)
  R_TYPE, intent(out)  :: vex(*)
  FLOAT, intent(in) :: dt, t

  integer :: j
  FLOAT :: x
  R_TYPE, allocatable :: c(:)

  x = (t/dt)
  allocate(c(0:order))
  ! I got this coefficients from mathematica...
  select case(order)
  case(1)
    c(0) = 1.0 + x
    c(1) =     - x 
  case(2)
    c(0) = 1.0 + (3.0/2.0)*x + (1.0/2.0)*x**2
    c(1) =     -      2.0 *x -           x**2
    c(2) =       (1.0/2.0)*x + (1.0/2.0)*x**2
  case(3)
    c(0) = 1.0 + (11.0/6.0)*x +           x**2 + (1.0/6.0)*x**3
    c(1) =            -3.0 *x - (5.0/2.0)*x**2 - (1.0/2.0)*x**3
    c(2) =       ( 3.0/2.0)*x +      2.0 *x**2 + (1.0/2.0)*x**3
    c(3) =     - ( 1.0/3.0)*x - (1.0/2.0)*x**2 - (1.0/6.0)*x**3
  case(4)
    c(0) = 1.0 + (25.0/12.0)*x + (35.0/24.0)*x**2 + ( 5.0/12.0)*x**3 + ( 1.0/24.0)*x**4
    c(1) =     -        4.0 *x - (13.0/ 3.0)*x**2 - ( 3.0/ 2.0)*x**3 - ( 1.0/ 6.0)*x**4
    c(2) =              3.0 *x + (19.0/ 4.0)*x**2 +        2.0 *x**3 + ( 1.0/ 4.0)*x**4
    c(3) =     - ( 4.0/ 3.0)*x - ( 7.0/ 3.0)*x**2 - ( 7.0/ 6.0)*x**3 - ( 1.0/ 6.0)*x**4
    c(4) =       ( 1.0/ 4.0)*x + (11.0/24.0)*x**2 + ( 1.0/ 4.0)*x**3 + ( 1.0/24.0)*x**4
  case default
    message(1) = 'extrapolate: Unsupported order.'
    call write_fatal(1)
  end select

  vex(1:n) = R_TOTYPE(M_ZERO)
  do j = 0, order
    call X(axpy)(n, c(j), v(n*j+1), 1, vex(1), 1)
  enddo

  deallocate(c)
end subroutine X(extrapolate)
