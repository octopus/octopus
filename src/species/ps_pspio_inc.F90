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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

  ! ---------------------------------------------------------
  subroutine ps_pspio_init(ps, namespace, label, z, lmax, lloc, ispin, filename)
    type(ps_t),        intent(out)   :: ps
    type(namespace_t), intent(in)    :: namespace
    character(len=10), intent(in)    :: label
    integer,           intent(inout) :: lmax
    integer,           intent(in)    :: lloc, ispin
    FLOAT,             intent(in)    :: z
    character(len=*),  intent(in)    :: filename

#ifdef HAVE_PSPIO

    logical :: found, has_kb
    integer :: idir
    character(len=3) :: psp_dir(3) = (/"PSF", "FHI", "UPF"/)
    character(len=MAX_PATH_LEN) :: filename2
    type(fpspio_pspdata_t)   :: pspdata

    PUSH_SUB(ps_pspio_init)

    call messages_experimental("Reading pseudopotential file using PSPIO library")

    ! Find out the file
    filename2 = trim(filename)
    inquire(file=filename2, exist=found)

    if(.not. found) then
      do idir = 1, size(psp_dir)
        filename2 = trim(conf%share) // "/pseudopotentials/" // trim(psp_dir(idir)) // "/" // trim(filename)
        inquire(file=filename2, exist=found)

        if(found) exit

        if (idir == size(psp_dir)) then
          message(1) = "Pseudopotential file '" // trim(filename) // " not found"
          call messages_fatal(1, namespace=namespace)
        end if
      end do
    end if

    message(1) = "Reading pseudopotential from file:"
    write(message(2), '(6x,3a)') "'", trim(filename2), "'"
    call messages_info(2)

    ! Init pspio data structure and parse file
    call check_error(fpspio_pspdata_alloc(pspdata), namespace)
    call check_error(fpspio_pspdata_read(pspdata, PSPIO_FMT_UNKNOWN, filename2), namespace)

    ! General info
    ps%label = label
    ps%ispin = ispin
    ps%hamann = .false.
    ps%z = z
    ps%conf%z = nint(z)
    call ps_pspio_read_info(ps, pspdata)
    lmax = ps%lmax
    write(message(1), '(a,i2,a)') "Info: l = ", ps%lmax, " is maximum angular momentum considered."
    call messages_info(1)

    ! Mesh
    call ps_pspio_read_mesh(ps, pspdata, namespace)

    ! XC
    call ps_pspio_read_xc(ps, pspdata)

    ! We will first try to read the KB projectors
    call ps_pspio_read_kb_projectors(ps, pspdata, has_kb)

    ! States
    call ps_pspio_read_states(ps, pspdata)

    ! Density
    call ps_pspio_read_density(ps, pspdata)

    ! If we do not have KB projectors, then we read the pseudopotentials
    if (.not. has_kb) then
      call ps_pspio_read_potentials(ps, pspdata, namespace)
    end if

    !No variable description, as it is already in ps.F90
    call parse_variable(namespace, 'SpeciesProjectorSphereThreshold', CNST(0.001), ps%projectors_sphere_threshold)
    if(ps%projectors_sphere_threshold <= M_ZERO) call messages_input_error('SpeciesProjectorSphereThreshold')
    ps%has_long_range = .true.
    ps%is_separated = .false.

    ps%local = ps%lmax == 0 .and. ps%llocal == 0
    
    !Free memory
    call fpspio_pspdata_free(pspdata)

#else
    message(1) = 'PSPIO selected for pseudopotential parsing, but the code was compiled witout PSPIO support.'
    call messages_fatal(1, namespace=namespace)
#endif

    POP_SUB(ps_pspio_init)
  end subroutine ps_pspio_init


#ifdef HAVE_PSPIO

  ! ---------------------------------------------------------
  subroutine ps_pspio_read_info(ps, pspdata)
    type(ps_t),             intent(inout) :: ps
    type(fpspio_pspdata_t), intent(in)    :: pspdata

    PUSH_SUB(ps_pspio_read_info)

    ps%conf%symbol = fpspio_pspdata_get_symbol(pspdata)
    ps%lmax = fpspio_pspdata_get_l_max(pspdata)
    ps%z_val = fpspio_pspdata_get_zvalence(pspdata)

    POP_SUB(ps_pspio_read_info)
  end subroutine ps_pspio_read_info

  ! ---------------------------------------------------------
  subroutine ps_pspio_read_mesh(ps, pspdata, namespace)
    type(ps_t),             intent(inout) :: ps
    type(fpspio_pspdata_t), intent(in)    :: pspdata
    type(namespace_t),      intent(in)    :: namespace

    integer :: ip
    type(fpspio_mesh_t) :: mesh
    FLOAT, pointer :: r_tmp(:)

    PUSH_SUB(ps_pspio_read_mesh)

    mesh = fpspio_pspdata_get_mesh(pspdata)
    ps%g%nrval = fpspio_mesh_get_np(mesh)

    r_tmp => fpspio_mesh_get_r(mesh)
    if(any(abs(r_tmp(2:ps%g%nrval)) < M_EPSILON)) then
      ! only the first point is allowed to be zero
      message(1) = "Illegal zero values in PSPIO radial grid"
      call messages_fatal(1, namespace=namespace)
    end if
    if (abs(r_tmp(1)) <= M_EPSILON) then
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

    POP_SUB(ps_pspio_read_mesh)
  end subroutine ps_pspio_read_mesh

  ! ---------------------------------------------------------
  subroutine ps_pspio_read_states(ps, pspdata)
    type(ps_t),             intent(inout) :: ps
    type(fpspio_pspdata_t), intent(in)    :: pspdata

    integer :: is, ist, ir
    FLOAT :: x
    type(fpspio_state_t) :: state
    FLOAT, allocatable :: wfs(:)

    PUSH_SUB(ps_pspio_read_states)

    ps%conf%p = fpspio_pspdata_get_n_states(pspdata)
    SAFE_ALLOCATE(ps%ur   (1:ps%conf%p, 1:ps%ispin))
    SAFE_ALLOCATE(ps%ur_sq(1:ps%conf%p, 1:ps%ispin))
    do ist = 1, ps%conf%p
      state = fpspio_pspdata_get_state(pspdata, ist)

      !Quantum numbers
      ps%conf%n(ist) = fpspio_qn_get_n(fpspio_state_get_qn(state))
      ps%conf%l(ist) = fpspio_qn_get_l(fpspio_state_get_qn(state))

      !Occupations
      ps%conf%occ(ist, 1) = fpspio_state_get_occ(state)
      if(ps%ispin == 2) then
        ! Spin-dependent pseudopotentials are not supported, so we need to fix the occupations
        ! if we want to have a spin-dependent atomic density.      
        x = ps%conf%occ(ps%conf%l(ist), 1)
        ps%conf%occ(ist, 1) = min(x, real(2*ps%conf%l(ist)+1, REAL_PRECISION))
        ps%conf%occ(ist, 2) = x - ps%conf%occ(ps%conf%l(ist), 1)
      end if

      !Wavefunctions
      SAFE_ALLOCATE(wfs(1:ps%g%nrval))
      do ir = 1, ps%g%nrval
        wfs(ir) = fpspio_state_wf_eval(state, ps%g%rofi(ir))
      end do
      do is = 1, ps%ispin
        call spline_fit(ps%g%nrval, ps%g%rofi, wfs, ps%ur(ist, is))
      end do
      do is = 1, ps%ispin
        call spline_fit(ps%g%nrval, ps%g%r2ofi, wfs, ps%ur_sq(ist, is))
      end do
      SAFE_DEALLOCATE_A(wfs)
    end do

    POP_SUB(ps_pspio_read_states)
  end subroutine ps_pspio_read_states

  ! ---------------------------------------------------------
  subroutine ps_pspio_read_density(ps, pspdata)
    type(ps_t),             intent(inout) :: ps
    type(fpspio_pspdata_t), intent(in)    :: pspdata

    integer :: ir, is
    type(fpspio_meshfunc_t) :: density
    FLOAT, allocatable :: dens(:)
    
    PUSH_SUB(ps_pspio_read_density)

    density = fpspio_pspdata_get_rho_valence(pspdata)
    ps%has_density = fpspio_associated(density)

    SAFE_ALLOCATE(ps%density(1:ps%ispin))
    SAFE_ALLOCATE(ps%density_der(1:ps%ispin))

    call spline_init(ps%density)
    call spline_init(ps%density_der)
    
    if (ps%has_density) then
      SAFE_ALLOCATE(dens(1:ps%g%nrval))

      do ir = 1, ps%g%nrval
        dens(ir) = fpspio_meshfunc_eval(density, ps%g%rofi(ir))
      end do
      
      do is = 1, ps%ispin
        call spline_fit(ps%g%nrval, ps%g%rofi, dens, ps%density(is))
        call spline_der(ps%density(is), ps%density_der(is))
      end do

      SAFE_DEALLOCATE_A(dens)
    end if
    
    POP_SUB(ps_pspio_read_density)
  end subroutine ps_pspio_read_density

  ! ---------------------------------------------------------
  subroutine ps_pspio_read_kb_projectors(ps, pspdata, has_kb)
    type(ps_t),             intent(inout) :: ps
    type(fpspio_pspdata_t), intent(in)    :: pspdata
    logical,                intent(out)   :: has_kb

    integer :: ir, n_kbproj, wave_eq, ikb, ikbc, l
    FLOAT :: j
    FLOAT, parameter :: threshold = CNST(0.5e-7)
    type(fpspio_potential_t) :: vlocal
    type(fpspio_projector_t) :: kb_projector
    logical, allocatable :: was_init(:,:)
    FLOAT, allocatable :: v_local(:), proj(:)

    PUSH_SUB(ps_pspio_read_kb_projectors)

    n_kbproj = fpspio_pspdata_get_n_projectors(pspdata)

    has_kb = .true.
    if (n_kbproj == 0) then
      has_kb = .false.
      POP_SUB(ps_pspio_read_kb_projectors)
      return
    end if

    ! Local potential
    ps%llocal = fpspio_pspdata_get_l_local(pspdata)
    vlocal = fpspio_pspdata_get_vlocal(pspdata)
    SAFE_ALLOCATE(v_local(1:ps%g%nrval))
    do ir = 1, ps%g%nrval
      v_local(ir) = fpspio_potential_eval(vlocal, ps%g%rofi(ir))
    end do
    do ir = ps%g%nrval-1, 2, -1
      if(abs(v_local(ir)*ps%g%rofi(ir) + ps%z_val) > threshold) exit
    end do
    ps%rc_max = ps%g%rofi(ir + 1)
    call spline_init(ps%vl)
    call spline_fit(ps%g%nrval, ps%g%rofi, v_local, ps%vl)
    SAFE_DEALLOCATE_A(v_local)

    ! KB projectors 
    wave_eq = fpspio_pspdata_get_wave_eq(pspdata)
    if (wave_eq == PSPIO_EQN_DIRAC) then
      ps%kbc = 2
    else
      ps%kbc = 1
    end if

    SAFE_ALLOCATE(ps%kb (0:ps%lmax, 1:ps%kbc))
    SAFE_ALLOCATE(ps%dkb(0:ps%lmax, 1:ps%kbc))
    SAFE_ALLOCATE(ps%h  (0:ps%lmax, 1:ps%kbc, 1:ps%kbc))
    SAFE_ALLOCATE(was_init(0:ps%lmax, 1:ps%kbc))
    SAFE_ALLOCATE(proj(1:ps%g%nrval))
    call spline_init(ps%kb)
    call spline_init(ps%dkb)
    nullify(ps%k)

    was_init = .false.
    ps%h = M_ZERO
    do ikb = 1, n_kbproj
      kb_projector = fpspio_pspdata_get_projector(pspdata, ikb)

      l = fpspio_qn_get_l(fpspio_projector_get_qn(kb_projector))
      j = fpspio_qn_get_j(fpspio_projector_get_qn(kb_projector))
      ikbc = 1
      if (j == real(l, REAL_PRECISION) - M_HALF) ikbc = 2

      ps%h(l, ikbc, ikbc) = fpspio_projector_get_energy(kb_projector)

      do ir = 1, ps%g%nrval
        proj(ir) = fpspio_projector_eval(kb_projector, ps%g%rofi(ir))
      end do
      do ir = ps%g%nrval-1, 2, -1
        if(abs(proj(ir)) > threshold) exit
      end do
      ps%rc_max = max(ps%g%rofi(ir + 1), ps%rc_max)
      proj(ir+1:ps%g%nrval) = M_ZERO
      call spline_fit(ps%g%nrval, ps%g%rofi, proj, ps%kb(l, ikbc))
      was_init(l, ikbc) = .true.
    end do

    !Make sure all projector splines were initialized
    proj = M_ZERO
    do l = 0, ps%lmax
      do ikbc = 1, ps%kbc
        if (.not. was_init(l, ikbc)) then
          call spline_fit(ps%g%nrval, ps%g%rofi, proj, ps%kb(l, ikbc))
        end if
      end do
    end do
    SAFE_DEALLOCATE_A(proj)
    SAFE_DEALLOCATE_A(was_init)

    ! Increasing radius a little, just in case.
    ! I have hard-coded a larger increase of the cutoff for the filtering.
    ps%rc_max = ps%rc_max*CNST(1.5)

    POP_SUB(ps_pspio_read_kb_projectors)
  end subroutine ps_pspio_read_kb_projectors

  ! ---------------------------------------------------------
  subroutine ps_pspio_read_potentials(ps, pspdata, namespace)
    type(ps_t),             intent(inout) :: ps
    type(fpspio_pspdata_t), intent(in)    :: pspdata
    type(namespace_t),      intent(in)    :: namespace

    PUSH_SUB(ps_pspio_read_potentials)

    message(1) = "Not yet implemented"
    call messages_fatal(1, namespace=namespace)

    POP_SUB(ps_pspio_read_potentials)
  end subroutine ps_pspio_read_potentials

  ! ---------------------------------------------------------
  subroutine ps_pspio_read_xc(ps, pspdata)
    type(ps_t),             intent(inout) :: ps
    type(fpspio_pspdata_t), intent(in)    :: pspdata

    integer :: ir, nrc
    type(fpspio_xc_t) :: xc
    FLOAT, allocatable :: rho(:)

    PUSH_SUB(ps_pspio_read_xc)

    !Non-linear core-corrections
    xc = fpspio_pspdata_get_xc(pspdata)

    call spline_init(ps%core)
    ps%nlcc = fpspio_xc_get_nlcc_scheme(xc) /= PSPIO_NLCC_NONE
    if (ps%nlcc) then
      ! get core density
      SAFE_ALLOCATE(rho(1:ps%g%nrval))
      do ir = 1, ps%g%nrval
        rho(ir) = fpspio_xc_nlcc_density_eval(xc, ps%g%rofi(ir))
      end do
      
      ! find cutoff radius
      nrc = ps%g%nrval
      do ir = ps%g%nrval-1, 1, -1
        if(rho(ir) > eps) then
          nrc = ir + 1
          exit
        end if
      end do
      rho(nrc:ps%g%nrval) = M_ZERO

      call spline_fit(ps%g%nrval, ps%g%rofi, rho, ps%core)

      SAFE_DEALLOCATE_A(rho)
    end if

    POP_SUB(ps_pspio_read_xc)
  end subroutine ps_pspio_read_xc

  ! ---------------------------------------------------------
  subroutine check_error(ierr, namespace)
    integer,           intent(in) :: ierr
    type(namespace_t), intent(in) :: namespace

    if (ierr /= PSPIO_SUCCESS) then
      call fpspio_error_flush()
      message(1) = "PSPIO error"
      call messages_fatal(1, namespace=namespace)
    end if

  end subroutine check_error

#endif

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
