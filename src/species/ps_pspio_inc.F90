!! Copyright (C) 2012 M. Oliveira
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
!! $Id: ps_upf.F90 8448 2011-11-01 13:35:20Z xavier $

  subroutine ps_pspio_init(ps, z, lmax, lloc, ispin, filename)
    type(ps_t),        intent(out)   :: ps
    integer,           intent(inout) :: lmax
    integer,           intent(in)    :: lloc, ispin
    FLOAT,             intent(in)    :: z
    character(len=*),  intent(in)    :: filename

#if HAVE_PSPIO

    logical :: found, has_nlcc
    integer :: ierr, idir, n_kbproj, wave_eq, ip, ir, nrc, ist, is, ikb, ikbc, l
    FLOAT :: j, x
    character(len=3) :: psp_dir(3) = (/"PSF", "FHI", "UPF"/)
    character(len=256) :: filename2
    FLOAT, parameter :: threshold = CNST(0.5e-7)
    logical, allocatable :: was_init(:,:)
    FLOAT, allocatable :: r_tmp(:), mesh_func(:)
    type(pspio_f90_pspdata_t)   :: pspdata
    type(pspio_f90_mesh_t)      :: mesh
    type(pspio_f90_state_t)     :: state
    type(pspio_f90_potential_t) :: vlocal
    type(pspio_f90_projector_t) :: kb_projector
    type(pspio_f90_xc_t)        :: xc

    PUSH_SUB(ps_pspio_init)

    call messages_experimental("Reading pseudopotential file using PSPIO library")

    ! Find out the file
    filename2 = trim(filename)
    inquire(file=filename2, exist=found)

    if(.not. found) then
      do idir = 1, size(psp_dir)
        filename2 = trim(conf%share) // "/PP/" // trim(psp_dir(idir)) // "/" // trim(filename)
        inquire(file=filename2, exist=found)

        if(found) exit

        if (idir == size(psp_dir)) then
          message(1) = "Pseudopotential file '" // trim(filename) // " not found"
          call messages_fatal(1)
        end if
      end do
    end if

    message(1) = "Reading pseudopotential from file:"
    write(message(2), '(6x,3a)') "'", trim(filename2), "'"
    call messages_info(2)

    ! Init pspio data structure and parse file
    ierr = pspio_f90_pspdata_init(pspdata)
    call check_error()
    ierr = pspio_f90_pspdata_read(pspdata, PSPIO_UNKNOWN, filename2)
    call check_error()

    ! General info
    ps%ispin = ispin
    ps%z = z
    ps%conf%z = nint(z)
    ierr = pspio_f90_pspdata_get_symbol(pspdata, ps%conf%symbol)
    call check_error()
    ierr = pspio_f90_pspdata_get_l_max(pspdata, ps%l_max)
    call check_error()
    lmax = ps%l_max
    ierr = pspio_f90_pspdata_get_wave_eq(pspdata, wave_eq)
    call check_error()
    ierr = pspio_f90_pspdata_get_zvalence(pspdata, ps%z_val)
    call check_error()

    write(message(1), '(a,i2,a)') "Info: l = ", ps%l_max, " is maximum angular momentum considered."
    call messages_info(1)

    ! Mesh
    ierr = pspio_f90_pspdata_get_mesh(pspdata, mesh)
    call check_error()
    ierr = pspio_f90_mesh_get_np(mesh, ps%g%nrval)
    call check_error()

    SAFE_ALLOCATE(r_tmp(1:ps%g%nrval))
    ierr = pspio_f90_mesh_get_r(mesh, r_tmp(1))
    call check_error()
    if(any(abs(r_tmp(2:ps%g%nrval)) < M_EPSILON)) then
      ! only the first point is allowed to be zero
      message(1) = "Illegal zero values in PSPIO radial grid"
      call messages_fatal(1)
    endif
    if (r_tmp(1) == M_ZERO) then
      ip = 1
    else
      ip = 2
      ps%g%nrval = ps%g%nrval + 1
    end if
    nullify(ps%g%drdi, ps%g%s)
    SAFE_ALLOCATE(ps%g%rofi(1:ps%g%nrval))
    SAFE_ALLOCATE(ps%g%r2ofi(1:ps%g%nrval))
    ps%g%rofi(1) = M_ZERO
    ps%g%rofi(ip:ps%g%nrval) = r_tmp(1:ps%g%nrval-ip+1)
    ps%g%r2ofi = ps%g%rofi**2
    SAFE_DEALLOCATE_A(r_tmp)
    SAFE_ALLOCATE(mesh_func(ps%g%nrval))

    ! Check if we have KB projectors
    ierr = pspio_f90_pspdata_get_n_kbproj(pspdata, n_kbproj)
    call check_error()
    if (n_kbproj == 0) then
      message(1) = "Not yet implemented"
      call messages_fatal(1)
    end if

    ! Local potential
    ierr = pspio_f90_pspdata_get_l_local(pspdata, ps%l_loc)
    call check_error()
    ierr = pspio_f90_pspdata_get_vlocal(pspdata, vlocal)
    call check_error()
    do ir = 0, ps%g%nrval - ip
      ierr = pspio_f90_potential_eval(vlocal, ps%g%rofi(ip+ir), mesh_func(ip+ir))
      call check_error()
    end do
    if (ip == 2) mesh_func(1) = linear_extrapolate(ps%g%rofi(1), ps%g%rofi(2), &
         ps%g%rofi(3), mesh_func(2), mesh_func(3))
    do ir = ps%g%nrval-1, 2, -1
      if(abs(mesh_func(ir)*ps%g%rofi(ir) + ps%z_val) > threshold) exit
    end do
    ps%rc_max = ps%g%rofi(ir + 1)
    call spline_init(ps%vl)
    call spline_fit(ps%g%nrval, ps%g%rofi, mesh_func, ps%vl)

    ! KB projectors    
    if (wave_eq == PSPIO_DIRAC) then
      ps%kbc = 2
    else
      ps%kbc = 1
    end if

    SAFE_ALLOCATE(ps%kb (0:ps%l_max, 1:ps%kbc))
    SAFE_ALLOCATE(ps%dkb(0:ps%l_max, 1:ps%kbc))
    SAFE_ALLOCATE(ps%h  (0:ps%l_max, 1:ps%kbc, 1:ps%kbc))
    SAFE_ALLOCATE(ps%k  (0:ps%l_max, 1:ps%kbc, 1:ps%kbc))
    SAFE_ALLOCATE(was_init(0:ps%l_max, 1:ps%kbc))
    call spline_init(ps%kb)
    call spline_init(ps%dkb)

    was_init = .false.
    ps%h = M_ZERO
    ps%k = M_ZERO
    do ikb = 1, n_kbproj
      ierr = pspio_f90_pspdata_get_kb_projector(pspdata, ikb, kb_projector)
      call check_error()

      ierr = pspio_f90_projector_get_l(kb_projector, l)
      call check_error()
      ierr = pspio_f90_projector_get_j(kb_projector, j)
      call check_error()
      ikbc = 1
      if (j == real(l, REAL_PRECISION) - M_HALF) ikbc = 2

      ierr = pspio_f90_projector_get_energy(kb_projector, ps%h(l, ikbc, ikbc))
      call check_error()

      do ir = 0, ps%g%nrval - ip
        ierr = pspio_f90_projector_eval(kb_projector, ps%g%rofi(ip+ir), mesh_func(ip+ir))
        call check_error()
      end do
      if (ip == 2) mesh_func(1) = linear_extrapolate(ps%g%rofi(1), ps%g%rofi(2), &
           ps%g%rofi(3), mesh_func(2), mesh_func(3))
      do ir = ps%g%nrval-1, 2, -1
        if(abs(mesh_func(ir)) > threshold) exit
      end do
      ps%rc_max = max(ps%g%rofi(ir + 1), ps%rc_max)
      mesh_func(ir+1:ps%g%nrval) = M_ZERO
      call spline_fit(ps%g%nrval, ps%g%rofi, mesh_func, ps%kb(l, ikbc))
      was_init(l, ikbc) = .true.
    end do

    !Make sure all projector splines were initialized
    mesh_func = M_ZERO
    do l = 0, ps%l_max
      do ikbc = 1, ps%kbc
        if (.not. was_init(l, ikbc)) then
          call spline_fit(ps%g%nrval, ps%g%rofi, mesh_func, ps%kb(l, ikbc))
        end if
      end do
    end do
    SAFE_DEALLOCATE_A(was_init)

    ! Increasing radius a little, just in case.
    ! I have hard-coded a larger increase of the cutoff for the filtering.
    ps%rc_max = ps%rc_max*CNST(1.5)

    ! States
    ierr = pspio_f90_pspdata_get_n_states(pspdata, ps%conf%p)
    call check_error()
    SAFE_ALLOCATE(ps%ur   (1:ps%conf%p, 1:ps%ispin))
    SAFE_ALLOCATE(ps%ur_sq(1:ps%conf%p, 1:ps%ispin))
    do ist = 1, ps%conf%p
      ierr = pspio_f90_pspdata_get_state(pspdata, ist, state)
      call check_error()

      !Quantum numbers
      ierr = pspio_f90_state_get_qn(state, ps%conf%n(ist), ps%conf%l(ist), j)
      call check_error()

      !Occupations
      ierr = pspio_f90_state_get_occ(state, ps%conf%occ(ist, 1))
      call check_error()
      if(ps%ispin == 2) then
        ! Spin-dependent pseudopotentials are not supported, so we need to fix the occupations
        ! if we want to have a spin-dependent atomic density.      
        x = ps%conf%occ(l, 1)
        ps%conf%occ(ist, 1) = min(x, real(2*ps%conf%l(ist)+1, REAL_PRECISION))
        ps%conf%occ(ist, 2) = x - ps%conf%occ(l, 1)
      end if

      !Wavefunctions
      do ir = 0, ps%g%nrval - ip
        ierr = pspio_f90_state_wf_eval(state, ps%g%rofi(ip+ir), mesh_func(ip+ir))
        call check_error()
      end do
      if (ip == 2) mesh_func(1) = linear_extrapolate(ps%g%rofi(1), ps%g%rofi(2), &
           ps%g%rofi(3), mesh_func(2), mesh_func(3))
      do is = 1, ps%ispin
        call spline_fit(ps%g%nrval, ps%g%rofi, mesh_func, ps%ur(ist, is))
      end do
      do is = 1, ps%ispin
        call spline_fit(ps%g%nrval, ps%g%r2ofi, mesh_func, ps%ur_sq(ist, is))
      end do

    end do

    !Non-linear core-corrections
    ierr = pspio_f90_pspdata_get_xc(pspdata, xc)
    call check_error()

    call spline_init(ps%core)
    ierr = pspio_f90_xc_has_nlcc(xc, has_nlcc)
    call check_error()
    if (has_nlcc) then
      ps%icore=''

      ! get core density
      do ir = 0, ps%g%nrval - ip
        ierr = pspio_f90_xc_nlcc_eval(xc, ps%g%rofi(ip+ir), mesh_func(ip+ir))
        call check_error()
      end do
      if (ip == 2) mesh_func(1) = linear_extrapolate(ps%g%rofi(1), ps%g%rofi(2), &
           ps%g%rofi(3), mesh_func(2), mesh_func(3))

      ! find cutoff radius
      do ir = ps%g%nrval-1, 1, -1
        if(mesh_func(ir) > eps) then
          nrc = ir + 1
          exit
        end if
      end do
      mesh_func(nrc:ps%g%nrval) = M_ZERO

      call spline_fit(ps%g%nrval, ps%g%rofi, mesh_func, ps%core)
    else
      ps%icore = 'nc'
    end if

    ! Fix the threshold to calculate the radius of the projector-function localization spheres:
    call messages_obsolete_variable('SpecieProjectorSphereThreshold', 'SpeciesProjectorSphereThreshold')

    !No variable description, as it is already in ps.F90
    call parse_float(datasets_check('SpeciesProjectorSphereThreshold'), &
      CNST(0.001), ps%projectors_sphere_threshold)
    if(ps%projectors_sphere_threshold <= M_ZERO) call input_error('SpeciesProjectorSphereThreshold')
    ps%has_long_range = .true.
    ps%is_separated = .false.

    !Free memory
    ierr = pspio_f90_pspdata_free(pspdata)
    call check_error()
    SAFE_DEALLOCATE_A(mesh_func)

#else
    message(1) = 'PSPIO selected for pseudopotential parsing, but the code was compiled witout PSPIO support.'
    call messages_fatal(1)
#endif

    POP_SUB(ps_pspio_init)

#if HAVE_PSPIO
  contains

    subroutine check_error()
      if (ierr /= PSPIO_SUCCESS) then
        ierr = pspio_f90_error_flush()
        message(1) = "PSPIO error"
        call messages_fatal(1)
      end if
    end subroutine check_error

#endif
  end subroutine ps_pspio_init

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
