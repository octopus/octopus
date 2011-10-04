!! Copyright (C) 2011 X. Andrade
!! Copyright (C) 2007 Ewout van den Berg and Michael P. Friedlander,
!! Department of Computer Science, University of British Columbia, Canada.
!! 
!! This program is a free translation to Fortran of part of the spgl1
!! library, version 1.0.
!!
!! This program is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!
!! $Id: bpdn.f90 3550 2007-11-19 14:32:49Z marques $

module bpdn_m
  implicit none
  
  private

  public ::    &
    bpdn

  interface 
    subroutine spgl1_projector(xx, cc, tau, nn)
      real(8), intent(inout) :: xx(1:nn)
      real(8), intent(in)    :: cc(1:nn)
      real(8), intent(in)    :: tau
      integer, intent(in)    :: nn
    end subroutine spgl1_projector
  end interface

  integer, parameter ::     &
    EXIT_CONVERGED  = 0,    &
    EXIT_ITERATIONS = 1,    &
    EXIT_NODESCENT  = 2

  logical, parameter :: debug = .true.

contains
  
  subroutine bpdn(nn, mm, aa, bb, sigma, xx)
    integer, intent(in)    :: nn
    integer, intent(in)    :: mm
    real(8), intent(in)    :: aa(:, :) !(1:nn, 1:mm)
    real(8), intent(in)    :: bb(:)    !(1:nn)
    real(8), intent(in)    :: sigma
    real(8), intent(out)   :: xx(:)    !(1:mm)
    integer, parameter :: nprevvals = 3
    real(8), parameter :: stepmin = 1.0e-16_8, stepmax = 1.0e+5_8
    real(8), parameter :: opttol = 1.0e-6_8
    real(8), parameter :: bptol = 1.0e-6_8
    real(8), parameter :: dectol = 1.0e-4_8

    real(8) :: lastffv(1:nprevvals)
    real(8) :: tau, tauold
    real(8) :: bbnorm, ff, ffbest, ffold, dxxnorm, gstep, gradnorm, gap, rgap, resnorm
    real(8), allocatable :: res(:), grad(:), xxbest(:), dxx(:), tmp(:), xxnew(:), ss(:), gradnew(:), yy(:)
    real(8) :: aerror, rerror, sts, sty
    integer :: nnzxold, nnzxiter, iter, maxits, lserr
    logical :: done, testsubopt, testrelchange1, testrelchange2, testupdatetau

    allocate(res(1:nn))
    allocate(grad(1:mm))
    allocate(xxbest(1:mm))
    allocate(dxx(1:mm))
    allocate(tmp(1:mm))
    allocate(xxnew(1:mm))
    allocate(ss(1:mm))
    allocate(gradnew(1:mm))
    allocate(yy(1:mm))

    bbnorm = norm_2(nn, bb)

    print*, 'Number of columns = ', nn
    print*, 'Number of rows    = ', mm
    print*, 'b 2-norm          = ', bbnorm
    
    tau = 0.0_8

    print*, "MATRIX ", ubound(aa), nn, mm

    ! Project the starting point and evaluate function and gradient.
    tmp = 1.0_8 ! some starting point
    call spgl1_projector(xx, tmp, tau, mm)
    res(1:nn) = bb(1:nn) - matmul(aa(1:nn, 1:mm), xx(1:mm))
    grad(1:mm) = -matmul(transpose(aa(1:nn, 1:mm)), res(1:nn))
    ff = dot_product(res(1:nn), res(1:nn))/2.0_8
    
    ! Required for nonmonotone strategy.
    lastffv(1:nprevvals) = ff
    ffbest               = ff
    xxbest(1:mm)         = xx(1:mm)
    ffold                = ff
    nnzxold              = -1
    nnzxiter             = 0
    
    ! Compute projected gradient direction and initial steplength.
    call spgl1_projector(dxx, xx(1:mm) - grad(1:mm), tau, mm)
    dxx(1:mm) = dxx(1:mm) - xx(1:mm)
    dxxnorm = maxval(abs(dxx(1:mm)))

    if(dxxnorm < (1.0_8 / stepmax)) then
      gstep = stepmax
    else
      gstep = min(stepmax, max(stepmin, 1.0_8/dxxnorm) );
    end if
    
    maxits = 2*nn
    iter = 0
    done = .false.
    do
      if(debug) then
        print*, 'Iter ', iter
      end if

      gradnorm = maxval(abs(grad))
      resnorm = norm_2(nn, res)
      gap = sum(res(1:nn)*(res(1:nn) - bb(1:nn))) + tau*gradnorm
      rgap = abs(gap)/max(1.0_8, ff)
      aerror = resnorm - sigma
      rerror = abs(ff - sigma**2/2.0_8)/max(1.0_8, ff)
      
      testsubopt = rgap <= max(opttol, rerror)

      if(testsubopt) then
        if(resnorm <= bptol*bbnorm) done = .true.
        if(gradnorm <= bptol*resnorm) done =  .true.
        if(abs(aerror) <= opttol*max(1.0_8, resnorm)) done = .true.
        if(resnorm <= sigma) done = .true.
      end if

      testrelchange1 = abs(ff - ffold) <= dectol*ff
      testrelchange2 = abs(ff - ffold) <= 1.0e-1_8*ff*abs(resnorm - sigma)
      testupdatetau = (testrelchange1 .and. resnorm > 2.0_8*sigma) &
        .or. (testrelchange2 .and. resnorm <= 2.0_8*sigma) .and. .not. done

      if(testupdatetau) then
        tauold = tau
        tau = tau + resnorm*aerror/gradnorm
        if(tau < tauold) then
          ! The one-norm ball has decreased.  Need to make sure that the
          ! next iterate if feasible, which we do by projecting it.
          tmp(1:mm) = xx(1:mm)
          call spgl1_projector(xx, tmp, tau, mm)
        end if
      end if

      ! Too many its and not converged.
      if (.not. done .and. iter >= maxits) then
        done = .true. 
      end if
    
      if(done) exit

      ! Iterations begin here.
      iter = iter + 1

      ffold = ff

      call spglinecurvy(mm, nn, aa, bb, xx, grad, maxval(lastffv), gstep, tau, xxnew, res, ff, lserr)

      if(ff > 0.5_8*sigma**2) then
        lastffv(mod(iter, nprevvals) + 1) = ff
        if(ffbest > ff) then
          ffbest = ff
          xxbest(1:mm) = xxnew(1:mm)
        end if
      end if
      
      ! Update gradient and compute new step length.
      gradnew(1:mm) = -matmul(transpose(aa(1:nn, 1:mm)), res(1:nn))
      ss(1:mm) = xxnew(1:mm) - xx(1:mm)
      yy(1:mm) = gradnew(1:mm) - grad(1:mm)
      xx(1:mm) = xxnew(1:mm)
      grad(1:mm) = gradnew(1:mm)
      sts = dot_product(ss(1:mm), ss(1:mm))
      sty = dot_product(ss(1:mm), yy(1:mm))
      if(sty <= 0.0_8) then
        gstep = stepmax
      else
        gstep = min(stepmax, max(stepmin, sts/sty) );
      end if

    end do

    deallocate(res)
    deallocate(grad)
    deallocate(xxbest)
    deallocate(dxx)
    deallocate(tmp)
    deallocate(xxnew)
    deallocate(ss)
    deallocate(gradnew)
    deallocate(yy)

  end subroutine bpdn

  !----------------------------------------------------------

  subroutine spglinecurvy(mm, nn, aa, bb, xx, gg, fmax, stepmax, tau, xxnew, resnew, fnew, ierr)
    integer, intent(in)    :: mm
    integer, intent(in)    :: nn
    real(8), intent(in)    :: aa(:, :) !(1:nn, 1:mm)
    real(8), intent(in)    :: bb(:)    !(1:nn)
    real(8), intent(in)    :: xx(:)    !(1:mm)
    real(8), intent(in)    :: gg(:)    !(1:mm)
    real(8), intent(in)    :: fmax
    real(8), intent(in)    :: stepmax
    real(8), intent(in)    :: tau
    real(8), intent(out)   :: xxnew(:)  !(1:mm)
    real(8), intent(out)   :: resnew(:) !(1:nn)
    real(8), intent(out)   :: fnew
    integer, intent(out)   :: ierr

    real(8), parameter :: gamma  = 1.0e-4_8
    integer, parameter :: maxiter = 10
    integer :: iter, nsafe
    real(8) :: scale, gts, ssnorm, ssnormold, ggnorm, step
    real(8), allocatable :: ss(:)

    allocate(ss(1:mm))

    step = stepmax
    ! safeguard scaling variables
    scale = 1.0_8 
    ssnorm = 0.0_8
    nsafe = 0

    iter = 0
    do

      call spgl1_projector(xxnew, xx(1:mm) - step*scale*gg(1:mm), tau, mm)
      resnew(1:nn) = bb(1:nn) - matmul(aa(1:nn, 1:mm), xxnew(1:mm))
      fnew = dot_product(resnew(1:nn), resnew(1:nn))/2.0_8
      ss(1:mm) = xxnew(1:mm) - xx(1:mm)
      gts = scale*dot_product(gg(1:mm), ss(1:mm))

      if(debug) then
        print*, ' LS ', iter, fNew, step, gts, scale
      end if
             
      if(gts >= 0.0_8) then
        ierr = EXIT_NODESCENT
        exit
      end if
      
      if(fnew < fmax - gamma*step*abs(gts)) then
        ierr = EXIT_CONVERGED
        exit
      else if(iter >= maxiter) then
        ierr = EXIT_ITERATIONS
        exit
      end if

      iter = iter + 1
      step = step/2.0_8

      ! Safeguard: If stepMax is huge, then even damped search
      ! directions can give exactly the same point after projection.  If
      ! we observe this in adjacent iterations, we drastically damp the
      ! next search direction.
      ssnormold = ssnorm
      ssnorm = norm_2(mm, ss)/sqrt(dble(mm))
      if(abs(ssnorm - ssnormold) < 1.0e-6*ssnorm) then
        ggnorm = norm_2(mm, gg)/sqrt(dble(mm))
        scale = ssnorm/ggnorm/(2.0_8**nsafe)
        nsafe = nsafe + 1
      end if
    end do

    deallocate(ss)

  end subroutine spglinecurvy

  ! --------------------------------------------
  
  real(8) pure function norm_2(nn, vv)
    integer, intent(in) :: nn
    real(8), intent(in) :: vv(:)

    norm_2 = sqrt(dot_product(vv(1:nn), vv(1:nn)))
  end function norm_2

end module bpdn_m
