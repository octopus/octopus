!! Copyright (C) 2011 X. Andrade
!! Copyright (C) 2007 Ewout van den Berg and Michael P. Friedlander,
!! Department of Computer Science, University of British Columbia, Canada.
!! 
!! This program is a free translation to Fortran of part of the spgl1
!! library, version 1.0.
!!
!! This program is free software: you can redistribute it and/or modify
!! it under the terms of the GNU Lesser General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU Lesser General Public License for more details.
!!
!! You should have received a copy of the GNU Lesser General Public License
!! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!
!! $Id: bpdn.f90 3550 2007-11-19 14:32:49Z marques $

module bpdn_m
  implicit none
  
  private

  public ::                 &
    bpdn,                   &
    bpdn_matrix,            &
    bpdn_matrix_init,       &
    bpdn_matrix_end,        &
    bpdn_matrix_set_delta

  interface
    subroutine spgl1_projector(xx, cc, tau, nn)
      implicit none

      integer, intent(in)    :: nn
      real(8), intent(inout) :: xx !(1:nn)
      real(8), intent(in)    :: cc !(1:nn)
      real(8), intent(in)    :: tau
    end subroutine spgl1_projector
  end interface

  integer, parameter ::        &
    EXIT_LINE_ERROR    =-2,    &
    EXIT_ITERATIONS    =-1,    &
    EXIT_CONVERGED     = 0,    &
    EXIT_ROOT_FOUND    = 2,    &
    EXIT_BPSOL1_FOUND  = 3,    &
    EXIT_BPSOL2_FOUND  = 4,    &
    EXIT_OPTIMAL       = 5,    &
    EXIT_SUBOPTIMAL_BP = 7,    &
    EXIT_ACTIVE_SET    = 9,    &
    EXIT_NODESCENT     = 10

  logical, parameter :: debug = .false.

  interface
    subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      implicit none

      character, intent(in)    :: trans
      integer,   intent(in)    :: m
      integer,   intent(in)    :: n
      real(8),   intent(in)    :: alpha
      real(8),   intent(in)    :: a !(1:lda, :)
      integer,   intent(in)    :: lda
      real(8),   intent(in)    :: x !(:)
      integer,   intent(in)    :: incx
      real(8),   intent(in)    :: beta
      real(8),   intent(in)    :: y !(:)
      integer,   intent(in)    :: incy
    end SUBROUTINE DGEMV
  end interface

  real(8), parameter :: opttol = 1.0e-6_8
  
  integer, parameter, public ::    &
    EXPLICIT_MATRIX = 1,           &
    COS_MATRIX      = 2,           &
    SIN_MATRIX      = 3

  type bpdn_matrix
    integer          :: type
    integer          :: nn
    integer          :: mm
    real(8), pointer :: matrix(:, :) !(1:nn, 1:mm)
    real(8)          :: dnn
    real(8)          :: dmm
  end type bpdn_matrix

contains

  subroutine bpdn_matrix_init(this, nn, mm, type)
    type(bpdn_matrix), intent(out)   :: this
    integer,           intent(in)    :: nn
    integer,           intent(in)    :: mm
    integer,           intent(in)    :: type

    this%nn = nn
    this%mm = mm
    this%type = type

    nullify(this%matrix)
    
    if(this%type == EXPLICIT_MATRIX) then
      allocate(this%matrix(1:this%nn, 1:this%mm))
    end if

  end subroutine bpdn_matrix_init

  subroutine bpdn_matrix_end(this)
    type(bpdn_matrix), intent(inout)   :: this
    
    if(this%type == EXPLICIT_MATRIX) then
      deallocate(this%matrix)
    end if
    this%type = 0

  end subroutine bpdn_matrix_end
  
  subroutine bpdn_matrix_set_delta(this, dnn, dmm)
    type(bpdn_matrix), intent(out)   :: this
    real(8),           intent(in)    :: dnn
    real(8),           intent(in)    :: dmm

    this%dnn = dnn
    this%dmm = dmm
  end subroutine bpdn_matrix_set_delta

  subroutine bpdn(nn, mm, aa, bb, sigma, xx, ierr, activesetit)
    type(bpdn_matrix), intent(in)    :: aa
    real(8),           intent(in)    :: bb(:)    !(1:nn)
    real(8),           intent(in)    :: sigma
    real(8),           intent(out)   :: xx(:)    !(1:mm)
    integer, optional, intent(out)   :: ierr     ! < 0 : some error occurred. >= 0 : solution found
    integer, optional, intent(in)    :: activesetit

    integer, parameter :: nprevvals = 3
    real(8), parameter :: stepmin = 1.0e-16_8, stepmax = 1.0e+5_8
    real(8), parameter :: bptol = 1.0e-6_8
    real(8), parameter :: dectol = 1.0e-4_8

    real(8) :: lastffv(1:nprevvals)
    real(8) :: tau, tauold
    real(8) :: bbnorm, ff, ffbest, ffold, dxxnorm, gstep, gradnorm, gap, rgap, resnorm
    real(8), allocatable :: res(:), grad(:), xxbest(:), dxx(:), tmp(:), xxnew(:), ss(:), gradnew(:), yy(:)
    real(8) :: aerror, rerror1, rerror2, sts, sty
    integer :: nnziter, iter, maxits, lserr, status, maxlineerrors, nnzdiff, nnzg, nnzx
    logical :: done, testrelchange1, testrelchange2, testupdatetau, singletau
    logical, allocatable :: nnzidx(:)
    integer :: nn, mm

    nn = aa%nn
    mm = aa%mm

    allocate(res(1:nn))
    allocate(grad(1:mm))
    allocate(xxbest(1:mm))
    allocate(dxx(1:mm))
    allocate(tmp(1:mm))
    allocate(xxnew(1:mm))
    allocate(ss(1:mm))
    allocate(gradnew(1:mm))
    allocate(yy(1:mm))
    allocate(nnzidx(1:mm))

    bbnorm = norm_2(nn, bb)

    print*, 'Number of columns = ', nn
    print*, 'Number of rows    = ', mm
    print*, 'b 2-norm          = ', bbnorm

    singletau = .false.

    ! Quick exit sigma is too large: tau = 0 will short-circuit the loop.
    if(bbnorm <= sigma) then
      tau = 0.0_8
      singletau = .true.
    end if
  
    tau = 0.0_8

    print*, 'Sigma             = ', sigma
    print*, 'Initial tau       = ', tau

    ! Project the starting point and evaluate function and gradient.
    tmp(1:mm) = 1.0_8 ! some starting point
    call spgl1_projector(xx(1), tmp(1), tau, mm)
    call residual(aa, bb, xx, res)
    call calc_grad(aa, res, grad)
    ff = dotp(nn, res, res)/2.0_8

    ! Required for nonmonotone strategy.
    lastffv(1:nprevvals) = ff
    ffbest               = ff
    xxbest(1:mm)         = xx(1:mm)
    ffold                = ff
    nnziter              = 0
    nnzidx(1:mm)         = .false.
    testupdatetau        = .false.

    ! Compute projected gradient direction and initial steplength.
    tmp(1:mm) = xx(1:mm) - grad(1:mm)
    call spgl1_projector(dxx(1), tmp(1), tau, mm)
    dxx(1:mm) = dxx(1:mm) - xx(1:mm)
    dxxnorm = maxval(abs(dxx(1:mm)))

    if(dxxnorm < (1.0_8 / stepmax)) then
      gstep = stepmax
    else
      gstep = min(stepmax, max(stepmin, 1.0_8/dxxnorm) );
    end if

    maxlineerrors = 10
    maxits = max(2000, 2*nn)
    iter = 0
    done = .false.
    status = -1
    do
      if(debug) then
        print*, 'Iter ', iter
      end if

      gradnorm = maxval(abs(grad))
      resnorm = norm_2(nn, res)
      gap = sum(res(1:nn)*(res(1:nn) - bb(1:nn))) + tau*gradnorm
      rgap = abs(gap)/max(1.0_8, ff)
      aerror = resnorm - sigma
      rerror1 = abs(aerror)/max(1.0_8, resnorm)
      rerror2 = abs(ff - sigma**2/2.0_8)/max(1.0_8, ff)

      ! count number of consecutive iterations with identical support.
      call active_vars(mm, xx, grad, nnzidx, nnzx, nnzg, nnzdiff)

      if(nnzDiff > 0) then
        nnziter = 0;
      else
        nnziter = nnziter + 1;
        if(present(activesetit)) then
          if(nnziter >= activesetit) then
            done = .true.
            status = EXIT_ACTIVE_SET
          end if
        end if
      end if

      if(singletau) then

        if(rgap <= opttol .or. resnorm < opttol*bbnorm) then
          status = EXIT_OPTIMAL
          done = .true.
        end if

      else

        if(rgap <= max(opttol, rerror2) .or. rerror1 <= opttol) then
          if(resnorm <= sigma) then
            status = EXIT_SUBOPTIMAL_BP
            done = .true.
          end if
          if(abs(aerror) <= opttol*max(1.0_8, resnorm)) then
            status = EXIT_ROOT_FOUND
            done = .true.
          end if
          if(gradnorm <= bptol*resnorm) then
            status = EXIT_BPSOL2_FOUND
            done = .true.
          end if
          if(resnorm <= bptol*bbnorm) then
            status = EXIT_BPSOL1_FOUND
            done = .true.
          end if
        end if

        testrelchange1 = abs(ff - ffold) <= dectol*ff
        testrelchange2 = abs(ff - ffold) <= 1.0e-1_8*ff*abs(resnorm - sigma)
        testupdatetau = ((testrelchange1 .and. resnorm > 2.0_8*sigma) &
          .or. (testrelchange2 .and. resnorm <= 2.0_8*sigma)) &
          .and. .not. done .and. .not. testupdatetau

        if(testupdatetau) then
          tauold = tau
          tau = max(0.0_8, tau + resnorm*aerror/gradnorm)
          if(tau < tauold) then
            ! The one-norm ball has decreased.  Need to make sure that the
            ! next iterate if feasible, which we do by projecting it.
            tmp(1:mm) = xx(1:mm)
            call spgl1_projector(xx(1), tmp(1), tau, mm)
          end if
        end if
        
      end if

      ! Too many its and not converged.
      if (.not. done .and. iter >= maxits) then
        status = EXIT_ITERATIONS
        done = .true.
      end if

      print*, 'iter = ', iter, 'tau =', tau, 'sigma(tau) =', resnorm


      if(done) exit

      ! Iterations begin here.
      iter = iter + 1

      ffold = ff

      call spg_line_curvy(aa, bb, xx, grad, maxval(lastffv), gstep, tau, xxnew, res, ff, lserr)

      if(lserr /= EXIT_CONVERGED) then
        if(maxlineerrors <= 0) then
          status = EXIT_LINE_ERROR
          done = .true.
        else
          maxlineerrors = maxLineErrors - 1
        end if
      end if

      if(singletau .or. ff > 0.5_8*sigma**2) then
        lastffv(mod(iter, nprevvals) + 1) = ff
        if(ffbest > ff) then
          ffbest = ff
          xxbest(1:mm) = xxnew(1:mm)
        end if
      end if

      ! Update gradient and compute new step length.
      call calc_grad(aa, res, gradnew)
      ss(1:mm) = xxnew(1:mm) - xx(1:mm)
      yy(1:mm) = gradnew(1:mm) - grad(1:mm)
      xx(1:mm) = xxnew(1:mm)
      grad(1:mm) = gradnew(1:mm)
      sts = dotp(mm, ss, ss)
      sty = dotp(mm, ss, yy)
      if(sty <= 0.0_8) then
        gstep = stepmax
      else
        gstep = min(stepmax, max(stepmin, sts/sty))
      end if

    end do

    if(singletau .and. ff > ffbest) then
      resnorm = sqrt(2.0_8*ffbest)
      print*, 'Restoring best iterate to objective ', resnorm
      xx(1:mm) = xxbest(1:mm)
      call residual(aa, bb, xx, res)
      call calc_grad(aa, res, grad)
    end if

    print*, 'Iterations        = ', iter
    print*, 'Final tau         = ', tau
    print*, 'Final sigma(tau)  = ', resnorm

    ! Print final output.
    select case (status)
    case(EXIT_OPTIMAL)
      print*, 'EXIT -- Optimal solution found'
    case(EXIT_ITERATIONS)
      print*, 'ERROR EXIT -- Too many iterations'
    case(EXIT_ROOT_FOUND)
      print*, 'EXIT -- Found a root'
    case(EXIT_BPSOL1_FOUND, EXIT_BPSOL2_FOUND)
      print*, 'EXIT -- Found a BP solution'
    case(EXIT_LINE_ERROR)
      print*, 'ERROR EXIT -- Linesearch error'
    case(EXIT_SUBOPTIMAL_BP)
      print*, 'EXIT -- Found a suboptimal BP solution'
    case(EXIT_ACTIVE_SET)
      print*, 'EXIT -- Found a possible active set'
    case default 
      stop 'Unknown termination condition'
    end select

    if(present(ierr)) then
      ierr = status
    end if

    deallocate(res)
    deallocate(grad)
    deallocate(xxbest)
    deallocate(dxx)
    deallocate(tmp)
    deallocate(xxnew)
    deallocate(ss)
    deallocate(gradnew)
    deallocate(yy)
    deallocate(nnzidx)

  end subroutine bpdn

  !----------------------------------------------------------

  subroutine spg_line_curvy(aa, bb, xx, gg, fmax, stepmax, tau, xxnew, resnew, fnew, ierr)
    type(bpdn_matrix), intent(in) :: aa
    real(8),           intent(in)    :: bb(:)    !(1:nn)
    real(8),           intent(in)    :: xx(:)    !(1:mm)
    real(8),           intent(in)    :: gg(:)    !(1:mm)
    real(8),           intent(in)    :: fmax
    real(8),           intent(in)    :: stepmax
    real(8),           intent(in)    :: tau
    real(8),           intent(out)   :: xxnew(:)  !(1:mm)
    real(8),           intent(out)   :: resnew(:) !(1:nn)
    real(8),           intent(out)   :: fnew
    integer,           intent(out)   :: ierr

    real(8), parameter :: gamma  = 1.0e-4_8
    integer, parameter :: maxiter = 10
    integer :: iter, nsafe
    real(8) ::  gts, ssnorm, ssnormold, step, scale, ggnorm
    real(8), allocatable :: ss(:)

    integer :: nn, mm
    
    nn = aa%nn
    mm = aa%mm

    allocate(ss(1:mm))

    step = stepmax
    scale = 1.0_8

    ! safeguard scaling variables
    ssnorm = 0.0_8
    nsafe = 0

    iter = 0
    do
      ss(1:mm) = xx(1:mm) - step*scale*gg(1:mm)
      call spgl1_projector(xxnew(1), ss(1), tau, mm)
      call residual(aa, bb, xxnew, resnew)
      fnew = dotp(nn, resnew, resnew)/2.0_8
      ss(1:mm) = xxnew(1:mm) - xx(1:mm)
      gts = scale*dotp(mm, gg, ss)

      if(gts >= 0) then
        ierr = EXIT_NODESCENT;
        exit
      end if

      if(debug) then
        print*, ' LS ', iter, fNew, step, gts
      end if
             
      if(fnew < fmax + gamma*gts) then
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

  end subroutine spg_line_curvy

  ! -----------------------------------------

  ! Find the current active set.
  subroutine active_vars(mm, xx, gg, nnzidx, nnzx, nnzg, nnzdiff)
    integer, intent(in)    :: mm
    real(8), intent(in)    :: xx(:)
    real(8), intent(in)    :: gg(:)
    logical, intent(inout) :: nnzidx(:) !< vector of primal/dual indicators.
    integer, intent(out)   :: nnzx      !< the number of nonzero x.
    integer, intent(out)   :: nnzg      !< is the number of elements in nnzIdx.
    integer, intent(out)   :: nnzdiff   !< the no. of elements that changed in the support

    real(8), parameter :: xtol  = 10.0_8*opttol
    real(8), parameter :: gtol  = 10.0_8*opttol
    real(8) :: ggnorm, z1, z2
    logical, allocatable :: nnzold(:)
    integer :: im
    logical :: xpos, xneg

    ggnorm = maxval(abs(gg(1:mm)))

    allocate(nnzold(1:mm))
    nnzold(1:mm) = nnzidx(1:mm)

    nnzg = 0
    nnzx = 0
    nnzdiff = 0
    do im = 1, mm
      ! Reduced costs for postive & negative parts of x.
      z1 = ggnorm + gg(im)
      z2 = ggnorm - gg(im)
      
      ! Primal/dual based indicators.
      xpos = xx(im) >  xtol .and. z1 < gtol
      xneg = xx(im) < -xtol .and. z2 < gtol
      nnzidx(im) = xpos .or. xneg
      if(nnzidx(im))  nnzg = nnzg + 1
      if(abs(xx(im)) >= xtol) nnzx = nnzx + 1
      if(.not. nnzidx(im) .eqv. nnzold(im)) nnzdiff = nnzdiff + 1 
    end do

  end subroutine active_vars


  ! --------------------------------------------
  
  real(8) function dotp(nn, vv1, vv2)
    integer, intent(in) :: nn
    real(8), intent(in) :: vv1(:)
    real(8), intent(in) :: vv2(:)

    interface 
      real(8) function ddot(n, x, incx, y, incy)
        implicit none
        
        integer, intent(in) :: n
        real(8), intent(in) :: x
        integer, intent(in) :: incx
        real(8), intent(in) :: y
        integer, intent(in) :: incy
      end function ddot
    end interface

    dotp = ddot(nn, vv1(1), 1, vv2(1), 1)
  end function dotp

  ! --------------------------------------------
  
  real(8) function norm_2(nn, vv)
    integer, intent(in) :: nn
    real(8), intent(in) :: vv(:)

    interface 
      real(8) function dnrm2(n, x, incx)
        implicit none
        
        integer, intent(in) :: n
        real(8), intent(in) :: x
        integer, intent(in) :: incx      
      end function dnrm2
    end interface

    norm_2 = dnrm2(nn, vv(1), 1)
  end function norm_2

  ! -------------------------------------------
  
  subroutine residual(aa, bb, xx, res)
    type(bpdn_matrix), intent(in)    :: aa
    real(8),           intent(in)    :: bb(:)    !(1:aa%nn)
    real(8),           intent(in)    :: xx(:)    !(1:aa%mm)
    real(8),           intent(out)   :: res(:)   !(1:nn)

    integer :: inn, imm, comp

    select case(aa%type)
    case(EXPLICIT_MATRIX)
      res(1:aa%nn) = bb(1:aa%nn)
      call dgemv(trans = 'n', m = aa%nn, n = aa%mm, alpha = -1.0_8, a = aa%matrix(1, 1), lda = ubound(aa%matrix, dim = 1), &
        x = xx(1), incx = 1, beta = 1.0_8, y = res(1), incy = 1)

    case(SIN_MATRIX, COS_MATRIX)

      comp = 0
      if(aa%type == SIN_MATRIX) comp = 1

      call expmm(aa%mm, aa%nn, xx(1), res(1), aa%dmm, aa%dnn,  1)
        
      do inn = 1, aa%nn
        res(inn) = bb(inn) + res(inn)
      end do

    end select

  end subroutine residual

  ! -------------------------------------------

  subroutine calc_grad(aa, res, grad)
    type(bpdn_matrix), intent(in)    :: aa
    real(8),           intent(in)    :: res(:)   !(1:aa%nn)
    real(8),           intent(out)   :: grad(:)  !(1:aa%mm)

    integer :: inn, imm, comp

    select case(aa%type)
    case(EXPLICIT_MATRIX)
      grad = 0.0_8
      call dgemv(trans = 't', m = aa%nn, n = aa%mm, alpha = -1.0_8, a = aa%matrix(1, 1), lda = ubound(aa%matrix, dim = 1), &
        x = res(1), incx = 1, beta = 0.0_8, y = grad(1), incy = 1)

    case(SIN_MATRIX, COS_MATRIX)
      comp = 0
      if(aa%type == SIN_MATRIX) comp = 1

      call expmm(aa%nn, aa%mm, res(1), grad(1), aa%dnn, aa%dmm, comp)
      
    end select

  end subroutine calc_grad

end module bpdn_m
