!! Copyright (C) 2014 Alain Delgado Gran, Carlo Andrea Rozzi, Stefano Corni
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
!! $Id$

#include "global.h"

module pcm_m
  use global_m
  use geometry_m
  use grid_m
  use io_m
  use index_m
  use messages_m
  use mesh_m 
  use parser_m
  use profiling_m
  use simul_box_m
  use species_m

  implicit none

  private

  public :: &
       pcm_t,                &
       pcm_init,             &
       pcm_end,              &
       pcm_charges,          &
       pcm_pot_rs,           &
       pcm_elect_energy,     &
       v_nuclei_cav,         &
       v_electrons_cav_li
  

 !> The cavity hosting the solute molecule is built from a set of 
 !! interlocking spheres with optimized radii centered at the nuclear positions.  
  type, public :: sphere_t  
    FLOAT :: x !
    FLOAT :: y !< center of the sphere
    FLOAT :: z !
    FLOAT :: r !< radius of the sphere (different for each species)
  end type

 !> The resulting cavity is discretized by a set of tesserae.  
  type, public :: tessera_t
    FLOAT :: point(MAX_DIM)  !< representative point of the tessera 
    FLOAT :: area            !< area of the tessera
    FLOAT :: normal(MAX_DIM) !< unitary outgoing vector normal to the tessera surface 
    FLOAT :: r_sphere        !< radius of the sphere to which the tessera belongs
  end type tessera_t

  type pcm_t
    logical                      :: run_pcm       !< If True, PCM calculation is enabled
    integer                      :: n_spheres     !< Number of spheres used to build the VdW cavity
    integer                      :: n_tesserae    !< Total number of tesserae
    type(sphere_t), allocatable  :: spheres(:)    !< See definition for type sphere_t
    type(tessera_t), allocatable :: tess(:)       !< See definition for type tessera_t
    FLOAT                        :: scale_r       !< factor to scale the radii of the spheres used in PCM
    FLOAT, allocatable           :: matrix(:,:)   !< PCM response matrix
    FLOAT, allocatable           :: q_e(:)        !< set of polarization charges due to the solute electrons        
    FLOAT, allocatable           :: q_n(:)        !< set of polarization charges due to the solute nuclei
    FLOAT                        :: qtot_e        !< total polarization charge due to electrons
    FLOAT                        :: qtot_n        !< total polarization charge due to nuclei
    FLOAT, allocatable           :: v_e(:)        !< Electrostatic potential produced by the electrons at each tessera
    FLOAT, allocatable           :: v_n(:)        !< Electrostatic potential produced by nuclei at each tessera
    FLOAT, allocatable           :: v_e_rs(:)     !< PCM real-space potential produced by q_e(:) 
    FLOAT, allocatable           :: v_n_rs(:)     !< PCM real-space potential produced by q_n(:)
    FLOAT, allocatable           :: arg_li(:,:)   !< 
    FLOAT                        :: epsilon_0     !< Static dielectric constant of the solvent 
    FLOAT                        :: epsilon_infty !< Infinite-frequency dielectric constant of the solvent
    FLOAT                        :: gaussian_width!< Parameter to change the width of density of polarization charges  
    integer                      :: n_vertices    !< Number of grid points used to interpolate the Hartree potential
                                                  !! at the tesserae representative points 
    integer, allocatable         :: ind_vh(:,:)   !< Grid points used during interpolation 
    integer                      :: info_unit     !< unit for pcm info file
    integer                      :: counter       !< used to print the number of SCF or TD iterations in energy_calc  
    character(len=80)            :: input_cavity  !< file name containing the geometry of the VdW cavity
  end type pcm_t

  FLOAT, allocatable :: s_mat_act(:,:) !< S_I matrix 
  FLOAT, allocatable :: d_mat_act(:,:) !< D_I matrix
  FLOAT, allocatable :: Sigma(:,:)     !< S_E matrix
  FLOAT, allocatable :: Delta(:,:)     !< D_E matrix in JCP 139, 024105 (2013).

  integer :: nearest_idx_unit   
  integer :: idx_from_coord_unit

  FLOAT, allocatable :: mat_gamess(:,:) !< PCM matrix formatted to be inputed to GAMESS
  FLOAT, allocatable :: sr_dist(:,:)    !< Table storing the distances between tesserae and grid points.

  !> End of variable declaration.
  !! ----------------------------

contains

  !-------------------------------------------------------------------------------------------------------
  !> Initializes the PCM calculation: reads the VdW molecular cavity and generates the PCM response matrix.
  subroutine pcm_init(pcm, geo, grid)
    type(geometry_t), intent(in) :: geo
    type(grid_t), intent(in)     :: grid
    type(pcm_t), intent(out)     :: pcm

    integer :: ia
    integer :: itess
    integer :: jtess
    integer :: cav_unit_test
    integer :: pcmmat_unit
    integer :: pcmmat_gamess_unit
    integer :: iunit
    integer :: vdw_unit
    integer :: grid_unit
    integer :: ip

    integer, parameter :: mxts = 10000

    FLOAT :: rcav_C
    FLOAT :: rcav_O
    FLOAT :: rcav_N
    FLOAT :: rcav_S
    FLOAT :: rcav_F
    FLOAT :: rcav_Na
    FLOAT :: rcav_Cl

    type(tessera_t) :: dum2(1)
    logical :: band

    PUSH_SUB(pcm_init)

    !%Variable Solvation
    !%Type logical
    !%Default no
    !%Section Hamiltonian::PCM
    !%Description
    !% (Experimental) If true, the calculation is performed accounting for solvation effects
    !% in the framework of Integral Equation Formalism Polarizable Continuum Model IEF-PCM
    !% (<i>Chem. Rev.</i> <b>105</b>, 2999 (2005), <i>J. Chem. Phys.</i> <b>107</b>, 3032 (1997),
    !% <i>J. Chem. Phys.</i> <b>139</b>, 024105 (2013)). At the moment, this option is available 
    !% only for ground-state calculations, and only for <tt>TheoryLevel = DFT</tt>.
    !%End

    call parse_variable('Solvation', .false., pcm%run_pcm)
    if (pcm%run_pcm) then
      if (grid%sb%box_shape /= MINIMUM) then
        message(1) = "PCM is only available for BoxShape = minimum"
        call messages_fatal(1)
      else 
        call messages_experimental("polarizable continuum model")
      end if
    else
      POP_SUB(pcm_init)
      return
    end if

    !%Variable RadiusScalingFactor
    !%Type float
    !%Default 1.0
    !%Section Hamiltonian::PCM
    !%Description
    !% Scales the radii of the spheres used to build the solute cavity surface.
    !%End
    call parse_variable('RadiusScalingFactor', M_ONE, pcm%scale_r)

    rcav_C  = CNST(2.4)*P_Ang*pcm%scale_r    ! 
    rcav_O  = CNST(1.8)*P_Ang*pcm%scale_r    !    
    rcav_N  = CNST(1.9)*P_Ang*pcm%scale_r    !
    rcav_S  = CNST(2.0175)*P_Ang*pcm%scale_r ! Angstrom -> Bohr 
    rcav_F  = CNST(1.682)*P_Ang*pcm%scale_r  !
    rcav_Na = CNST(2.772)*P_Ang*pcm%scale_r  !  
    rcav_Cl = CNST(2.172)*P_Ang*pcm%scale_r  !

    !%Variable SolventDielectricConstant
    !%Type float
    !%Default 1.0
    !%Section Hamiltonian::PCM
    !%Description
    !% Static dielectric constant of the solvent (<math>\varepsilon_0</math>). 1.0 indicates gas phase.
    !%End
    call parse_variable('SolventDielectricConstant', M_ONE, pcm%epsilon_0)

    !%Variable PcmSmearingFactor
    !%Type float
    !%Default 1.0
    !%Section Hamiltonian::PCM
    !%Description
    !% Parameter used to control the width (area of each tessera) of the Gaussians used to distribute
    !% the polarization charges on each tessera of the cavity surface. If set to zero, the solvent 
    !% reaction potential in real-space is defined by using point charges.
    !%End
    call parse_variable('PcmSmearingFactor', M_ONE, pcm%gaussian_width)

    if (pcm%gaussian_width == M_ZERO) then
      message(1) = "Info: PCM potential will be defined in terms of polarization point charges"
      call messages_info(1)
    else
      message(1) = "Info: PCM potential is regularized to avoid Coulomb singularity"
      call messages_info(1)        
    end if

    call io_mkdir('pcm')

    !%Variable CavityGeometry
    !%Type string
    !%Section Hamiltonian::PCM
    !%Description
    !% Name of the file containing the geometry of the Van der Waals surface that defines the cavity hosting
    !% the solute molecule in PCM calculations. Tesserae representative points must be in atomic units.
    !%End
    call parse_variable('CavityGeometry', '', pcm%input_cavity)

    if (pcm%input_cavity == '') then
    
     pcm%n_spheres = 0
     band = .false.
     do ia = 1, geo%natoms
       if (geo%atom(ia)%label == 'H') cycle
       pcm%n_spheres = pcm%n_spheres + 1 !counting the number of species different from Hydrogen
     end do
    
     SAFE_ALLOCATE( pcm%spheres(1:pcm%n_spheres) )
    
     pcm%n_spheres = 0
     do ia = 1, geo%natoms
      
      if (geo%atom(ia)%label == 'H') cycle
      pcm%n_spheres = pcm%n_spheres + 1
      
      !> These coordinates are already in atomic units (Bohr)
      pcm%spheres(pcm%n_spheres)%x = geo%atom(ia)%x(1)
      pcm%spheres(pcm%n_spheres)%y = geo%atom(ia)%x(2)
      pcm%spheres(pcm%n_spheres)%z = geo%atom(ia)%x(3)
      
      if (geo%atom(ia)%label == 'C') then
        pcm%spheres(pcm%n_spheres)%r = rcav_C
        band = .true.
      elseif (geo%atom(ia)%label == 'O') then
        pcm%spheres(pcm%n_spheres)%r = rcav_O
        band = .true.
      elseif (geo%atom(ia)%label == 'N') then
        pcm%spheres(pcm%n_spheres)%r = rcav_N
        band = .true.
      elseif (geo%atom(ia)%label == 'S') then
        pcm%spheres(pcm%n_spheres)%r = rcav_S
        band = .true.
      elseif (geo%atom(ia)%label == 'F') then
        pcm%spheres(pcm%n_spheres)%r = rcav_F
        band = .true.
      elseif (geo%atom(ia)%label == 'Na') then
        pcm%spheres(pcm%n_spheres)%r = rcav_Na
        band = .true.
      elseif (geo%atom(ia)%label == 'Cl') then 
        pcm%spheres(pcm%n_spheres)%r = rcav_Cl
        band = .true.
      endif

      if (.not.(band)) then
        write(message(1),'(a,a)') "Missing radius parameter for species", geo%atom(ia)%label
        call messages_fatal(1)
      endif

     end do
    
     pcm%info_unit = io_open(PCM_DIR//'pcm_info.out', action='write')

     write(pcm%info_unit,'(A35)') '# Configuration: Molecule + Solvent'
     write(pcm%info_unit,'(A35)') '# ---------------------------------'
     write(pcm%info_unit,'(A21,F12.3)') '# Epsilon(Solvent) = ', pcm%epsilon_0
     write(pcm%info_unit,'(A1)')'#' 
     write(pcm%info_unit,'(A35,I4)') '# Number of interlocking spheres = ', pcm%n_spheres
     write(pcm%info_unit,'(A1)')'#'  

     write(pcm%info_unit,'(A8,3X,A7,8X,A26,20X,A10)') '# SPHERE', 'ELEMENT', 'CENTER  (X,Y,Z) (A)', 'RADIUS (A)'
     write(pcm%info_unit,'(A8,3X,A7,4X,A43,7X,A10)') '# ------', '-------', &
                               '-------------------------------------------', '----------'  

     pcm%n_spheres = 0
     do ia = 1, geo%natoms
       if (geo%atom(ia)%label == 'H') cycle
       pcm%n_spheres = pcm%n_spheres + 1       
      
       write(pcm%info_unit,'(A1,2X,I3,9X,A2,3X,F14.8,2X,F14.8,2X,F14.8,3X,F14.8)')'#', pcm%n_spheres, &
            geo%atom(ia)%label,       &
            geo%atom(ia)%x*P_a_B,     &
            pcm%spheres(pcm%n_spheres)%r*P_a_B
     end do

     !> Counting the number of tesserae
     call cav_gen(0, 1, pcm%n_spheres, pcm%spheres, pcm%n_tesserae, dum2, pcm%info_unit)

     SAFE_ALLOCATE(pcm%tess(1:pcm%n_tesserae))

     !> Generating the Van der Waals discretized surface of the solute system
     call cav_gen(1, 1, pcm%n_spheres, pcm%spheres, pcm%n_tesserae, pcm%tess, pcm%info_unit)

     message(1) = "Info: van der Waals surface has been calculated"
     call messages_info(1)

    else

     !> The cavity surface will be read from a external file
     iunit = io_open(trim(pcm%input_cavity), status='old', action='read')
      read(iunit,*) pcm%n_tesserae
    
     if (pcm%n_tesserae.gt.mxts) then
       write(message(1),'(a,I5,a,I5)') "total number of tesserae", pcm%n_tesserae, ">",mxts
       call messages_warning(1)     
     endif
    
     SAFE_ALLOCATE(pcm%tess(1:pcm%n_tesserae))
    
     do ia=1, pcm%n_tesserae
       read(iunit,*) pcm%tess(ia)%point(1)
     end do
    
     do ia=1, pcm%n_tesserae
       read(iunit,*) pcm%tess(ia)%point(2)
     end do
    
     do ia=1, pcm%n_tesserae
       read(iunit,*) pcm%tess(ia)%point(3)
     end do
    
     do ia=1, pcm%n_tesserae
       read(iunit,*) pcm%tess(ia)%area
     end do
    
     do ia=1, pcm%n_tesserae
       read(iunit,*) pcm%tess(ia)%r_sphere
     end do
    
     do ia=1, pcm%n_tesserae
       read(iunit,*) pcm%tess(ia)%normal
     end do
    
     call io_close(iunit)
     message(1) = "Info: van der Waals surface has been read from " // trim(pcm%input_cavity)
     call messages_info(1)
    endif

    cav_unit_test = io_open(PCM_DIR//'cavity_mol.xyz', action='write')
    
    write (cav_unit_test,'(2X,I4)') pcm%n_tesserae + geo%natoms
    write (cav_unit_test,'(2X)')
    
    do ia = 1, pcm%n_tesserae
      write(cav_unit_test,'(2X,A2,3X,4f15.8,3X,4f15.8,3X,4f15.8)') 'H', pcm%tess(ia)%point*P_a_B
    end do
    
    do ia = 1, geo%natoms
      write(cav_unit_test,'(2X,A2,3X,4f15.8,3X,4f15.8,3X,4f15.8)') geo%atom(ia)%label,      &
           geo%atom(ia)%x*P_a_B
    end do

    call io_close(cav_unit_test)

    write(pcm%info_unit,'(A1)')'#'  
    write(pcm%info_unit,'(A1,4X,A4,9X,A4,15X,A4,15X,A4,15X,A4,15X,A8,12X,A5,15X,A5)') &
         '#','iter', 'E_ee', 'E_en', 'E_nn', 'E_ne', 'E_M-solv', 'Q_M^e','Q_M^n'
    pcm%counter = 0
    
    pcm%n_vertices = 8
    SAFE_ALLOCATE(pcm%ind_vh(1:pcm%n_tesserae, 1:pcm%n_vertices))
    pcm%ind_vh = INT(M_ZERO)
    
    SAFE_ALLOCATE(pcm%arg_li(1:pcm%n_tesserae, 1:grid%mesh%sb%dim))
    pcm%arg_li = M_ZERO

    !> Creating the list of the nearest grid points to each tessera
    !! to be used to interpolate the Hartree potential at the representative points
    SAFE_ALLOCATE(sr_dist(1:pcm%n_tesserae, 1:grid%mesh%np)) 
    sr_dist = M_ZERO
    do ia = 1, pcm%n_tesserae
      call nearest_cube_vertices(pcm%tess(ia)%point, grid%mesh, pcm%ind_vh(ia,:), pcm%arg_li(ia,:), ia, pcm%n_vertices)
      do ip = 1, grid%mesh%np !running serially np=np_global
        call mesh_r(grid%mesh, ip, sr_dist(ia,ip), origin=pcm%tess(ia)%point)
      end do
    end do
    
    !>Generating the dynamical PCM matrix to be inputed to GAMESS
    pcm%epsilon_infty = CNST(1.7760) 
    SAFE_ALLOCATE( mat_gamess(1:pcm%n_tesserae, 1:pcm%n_tesserae) )
    mat_gamess = M_ZERO
    
    SAFE_ALLOCATE( pcm%matrix(1:pcm%n_tesserae, 1:pcm%n_tesserae) )
    pcm%matrix = M_ZERO
    
    call pcm_matrix(pcm%epsilon_infty, pcm%tess, pcm%n_tesserae, pcm%matrix) 
    
    pcmmat_gamess_unit = io_open(PCM_DIR//'pcm_matrix_gamess_dyn.out', action='write')
    
    do jtess = 1, pcm%n_tesserae
      do itess = 1, pcm%n_tesserae
        write(pcmmat_gamess_unit,*) mat_gamess(itess,jtess)
      end do
    end do

    call io_close(pcmmat_gamess_unit)
    
    pcm%matrix = M_ZERO
    mat_gamess = M_ZERO
    
    call pcm_matrix(pcm%epsilon_0, pcm%tess, pcm%n_tesserae, pcm%matrix) 
    message(1) = "Info: PCM response matrix has been evaluated"
    call messages_info(1)
    
    pcmmat_unit = io_open(PCM_DIR//'pcm_matrix.out', action='write')
    pcmmat_gamess_unit = io_open(PCM_DIR//'pcm_matrix_gamess.out', action='write')
    
    do jtess = 1, pcm%n_tesserae
      do itess = 1, pcm%n_tesserae
        write(pcmmat_unit,*) pcm%matrix(itess,jtess)
        write(pcmmat_gamess_unit,*) mat_gamess(itess,jtess)
      end do
    end do

    call io_close(pcmmat_unit)
    call io_close(pcmmat_gamess_unit)
    
    SAFE_ALLOCATE( pcm%v_n(1:pcm%n_tesserae) )
    SAFE_ALLOCATE( pcm%q_n(1:pcm%n_tesserae) )
    SAFE_ALLOCATE( pcm%v_n_rs(1:grid%mesh%np) )
    pcm%v_n    = M_ZERO
    pcm%q_n    = M_ZERO
    pcm%v_n_rs = M_ZERO

    SAFE_ALLOCATE( pcm%v_e(1:pcm%n_tesserae) )
    SAFE_ALLOCATE( pcm%q_e(1:pcm%n_tesserae) )
    SAFE_ALLOCATE( pcm%v_e_rs(1:grid%mesh%np) )
    pcm%v_e    = M_ZERO
    pcm%q_e    = M_ZERO
    pcm%v_e_rs = M_ZERO
    
    POP_SUB(pcm_init)
  end subroutine pcm_init
      
  !> Calculates the Hartree potential at the tessera representative points by doing 
  !! a 3D linear interpolation. 
  subroutine v_electrons_cav_li(v_e_cav, v_hartree, pcm)
    type(pcm_t), intent(in)  :: pcm
    FLOAT, intent(in)        :: v_hartree(:) !< (1:mesh%np)
    FLOAT, intent(out)       :: v_e_cav(:)   !< (1:n_tess)

    integer :: ia
    
    FLOAT :: C_00
    FLOAT :: C_10
    FLOAT :: C_01
    FLOAT :: C_11
    FLOAT :: C_0
    FLOAT :: C_1
    
    PUSH_SUB(v_electrons_cav_li)    
    
    v_e_cav = M_ZERO
    
    do ia = 1, pcm%n_tesserae

      C_00 = v_hartree( pcm%ind_vh(ia,1) )*( M_ONE - pcm%arg_li(ia,1) ) + &
           v_hartree( pcm%ind_vh(ia,5) )*( pcm%arg_li(ia,1) )
      
      C_10 = v_hartree( pcm%ind_vh(ia,2) )*( M_ONE - pcm%arg_li(ia,1) ) + &
           v_hartree( pcm%ind_vh(ia,6) )*( pcm%arg_li(ia,1) )
      
      C_01 = v_hartree( pcm%ind_vh(ia,4) )*( M_ONE - pcm%arg_li(ia,1) ) + &
           v_hartree( pcm%ind_vh(ia,8) )*( pcm%arg_li(ia,1) )
      
      C_11 = v_hartree( pcm%ind_vh(ia,3) )*( M_ONE - pcm%arg_li(ia,1) ) + &
           v_hartree( pcm%ind_vh(ia,7) )*( pcm%arg_li(ia,1) )
      
      C_0 = C_00*( M_ONE - pcm%arg_li(ia,2) ) + C_10*pcm%arg_li(ia,2)
      
      C_1 = C_01*( M_ONE - pcm%arg_li(ia,2) ) + C_11*pcm%arg_li(ia,2)
      
      v_e_cav(ia) = C_0*( M_ONE - pcm%arg_li(ia,3) ) + C_1*pcm%arg_li(ia,3)
      
    end do
    
    POP_SUB(v_electrons_cav_li)
  end subroutine v_electrons_cav_li

  !> Calculates the classical electrostatic potential geneated by the nuclei at the tesserae.
  !! v_n_cav(ik) = \sum_{I=1}^{natoms} Z_val / |s_{ik} - R_I|
  subroutine v_nuclei_cav(v_n_cav, geo, tess, n_tess)
    FLOAT, intent(out)           :: v_n_cav(:) !< (1:n_tess)
    type(geometry_t), intent(in) :: geo
    type(tessera_t), intent(in)  :: tess(:)    !< (1:n_tess)
    integer, intent(in)          :: n_tess

    FLOAT   :: diff(1:MAX_DIM)
    FLOAT   :: dist
    FLOAT   :: z_ia
    integer :: ik
    integer :: ia

    type(species_t), pointer :: spci 

    PUSH_SUB(v_nuclei_cav)

    v_n_cav = M_ZERO

    do ik = 1, n_tess
      do ia = 1, geo%natoms
        diff = geo%atom(ia)%x - tess(ik)%point

        dist = dot_product( diff, diff )
        dist = sqrt(dist)

        spci => geo%atom(ia)%species
        z_ia = species_zval(spci)

        v_n_cav(ik) = v_n_cav(ik) + z_ia/dist       
      end do
    end do

    v_n_cav = -v_n_cav

    POP_SUB(v_nuclei_cav)
  end subroutine v_nuclei_cav

  !> Calculates the solute-solvent electrostatic interaction energy
  !! E_M-solv = \sum{ik=1}^n_tess { [VHartree(ik) + Vnuclei(ik)]*[q_e(ik) + q_n(ik)] }   
  subroutine pcm_elect_energy(geo, pcm, E_int_ee, E_int_en, E_int_ne, E_int_nn)
    type(geometry_t), intent(in) :: geo
    type(pcm_t), intent(in)      :: pcm
    FLOAT, intent(out)           :: E_int_ee 
    FLOAT, intent(out)           :: E_int_en 
    FLOAT, intent(out)           :: E_int_ne 
    FLOAT, intent(out)           :: E_int_nn 

    FLOAT   :: diff(1:MAX_DIM)
    FLOAT   :: dist
    FLOAT   :: z_ia
    integer :: ik
    integer :: ia
    
    type(species_t), pointer :: spci 
    
    PUSH_SUB(pcm_elect_energy)
     
    E_int_ee = M_ZERO
    E_int_en = M_ZERO
    E_int_ne = M_ZERO
    E_int_nn = M_ZERO
    
    do ik = 1, pcm%n_tesserae
      
      E_int_ee = E_int_ee + pcm%v_e(ik)*pcm%q_e(ik)
      E_int_en = E_int_en + pcm%v_e(ik)*pcm%q_n(ik)
      
      do ia = 1, geo%natoms
        diff = geo%atom(ia)%x - pcm%tess(ik)%point 
        
        dist = dot_product( diff, diff )
        dist = sqrt(dist)

        spci => geo%atom(ia)%species
        z_ia = -species_zval(spci)
        
        E_int_ne = E_int_ne + z_ia*pcm%q_e(ik) / dist
        E_int_nn = E_int_nn + z_ia*pcm%q_n(ik) / dist
        
      end do
    end do
    
    E_int_ee = M_HALF*E_int_ee
    E_int_en = M_HALF*E_int_en
    E_int_ne = M_HALF*E_int_ne
    E_int_nn = M_HALF*E_int_nn
    
    POP_SUB(pcm_elect_energy)
  end subroutine pcm_elect_energy

  !> Creating the list of the nearest 8 cube vertices in real-space 
  !! to calculate the Hartree potential at 'point'
  subroutine nearest_cube_vertices(point, mesh, vert_idx, weight_li, ia, n_vertices)
    FLOAT, intent(in)        :: point(1:MAX_DIM)
    type(mesh_t), intent(in) :: mesh
    integer, intent(out)     :: vert_idx(:)
    FLOAT, intent(out)       :: weight_li(:)
    integer, intent(in)      :: ia
    integer, intent(in)      :: n_vertices

    FLOAT   :: dmin
    integer :: rankmin

    FLOAT   :: coord_0(1:MAX_DIM)
    integer :: sign_x
    integer :: sign_y
    integer :: sign_z
    integer :: point_0(1:MAX_DIM)
    integer :: point_f(1:MAX_DIM)
    integer :: ii

    PUSH_SUB(nearest_cube_vertices)    

    vert_idx(1) = mesh_nearest_point(mesh, point, dmin, rankmin)

    coord_0 = mesh%x(vert_idx(1), :)
    point_0 = NINT(coord_0/mesh%spacing)

    sign_x = INT( sign( CNST(1.0), point(1) - coord_0(1) ) )
    sign_y = INT( sign( CNST(1.0), point(2) - coord_0(2) ) )
    sign_z = INT( sign( CNST(1.0), point(3) - coord_0(3) ) )

    weight_li = ABS(point - coord_0)/mesh%spacing 

    !FRONT CUBE PLANE
    point_f = point_0

    !POINT P2
    point_f(2) = point_f(2) + sign_y
    vert_idx(2) = index_from_coords(mesh%idx, point_f)

    !POINT P3
    point_f(3) = point_f(3) + sign_z
    vert_idx(3) = index_from_coords(mesh%idx, point_f)

    !POINT P4
    point_f(2) = point_f(2) - sign_y
    vert_idx(4) = index_from_coords(mesh%idx, point_f)

    !REAR CUBE PLANE
    point_f = point_0 

    !POINT P5
    point_f(1) = point_f(1) + sign_x
    vert_idx(5) = index_from_coords(mesh%idx, point_f)

    !POINT P6
    point_f(2) = point_f(2) + sign_y
    vert_idx(6) = index_from_coords(mesh%idx, point_f)

    !POINT P7
    point_f(3) = point_f(3) + sign_z
    vert_idx(7) = index_from_coords(mesh%idx, point_f)

    !POINT P8
    point_f(2) = point_f(2) - sign_y
    vert_idx(8) = index_from_coords(mesh%idx, point_f)

    POP_SUB(nearest_cube_vertices)    

  end subroutine nearest_cube_vertices

  !> Calculates the polarization charges at each tessera by using the response matrix 'pcm_mat',
  !! provided the value of the molecular electrostatic potential at 
  !! the tesserae: q_pcm(ia) = \sum_{ib}^{n_tess} pcm_mat(ia,ib)*v_cav(ib).
  subroutine pcm_charges(q_pcm, q_pcm_tot, v_cav, pcm_mat, n_tess)
    FLOAT, intent(out)   :: q_pcm(:)     !< (1:n_tess)
    FLOAT, intent(out)   :: q_pcm_tot
    FLOAT, intent(in)    :: v_cav(:)     !< (1:n_tess)
    FLOAT, intent(in)    :: pcm_mat(:,:) !< (1:n_tess, 1:n_tess)
    integer, intent(in)  :: n_tess

    integer :: ia
    integer :: ib

    PUSH_SUB(pcm_charges)

    q_pcm     = M_ZERO
    q_pcm_tot = M_ZERO

    do ia = 1, n_tess
      do ib = 1, n_tess
        q_pcm(ia) = q_pcm(ia) + pcm_mat(ia,ib)*v_cav(ib) !< transpose matrix might speed up
      end do
      q_pcm_tot = q_pcm_tot + q_pcm(ia)
    end do

    POP_SUB(pcm_charges)
  end subroutine pcm_charges
  
  !> Generates the potential 'v_pcm' in real-space.
  subroutine pcm_pot_rs(v_pcm, q_pcm, tess, n_tess, mesh, width_factor)
    FLOAT,           intent(out) :: v_pcm(:)!< (1:mesh%np) running serially np=np_global
    FLOAT,           intent(in)  :: q_pcm(:)!< (1:n_tess)
    FLOAT,           intent(in)  :: width_factor
    integer,         intent(in)  :: n_tess  
    type(mesh_t),    intent(in)  :: mesh
    type(tessera_t), intent(in)  :: tess(:) !< (1:n_tess)

    FLOAT, parameter :: p_1 = CNST(0.119763)
    FLOAT, parameter :: p_2 = CNST(0.205117)
    FLOAT, parameter :: q_1 = CNST(0.137546)
    FLOAT, parameter :: q_2 = CNST(0.434344)
    FLOAT            :: arg
    FLOAT            :: term
    integer 	     :: ip
    integer          :: ia

    PUSH_SUB(pcm_pot_rs)

    v_pcm = M_ZERO

    if (width_factor /= M_ZERO) then

      do ia = 1, n_tess
        do ip = 1, mesh%np !running serially np=np_global
          arg = sr_dist(ia,ip)/sqrt( tess(ia)%area*width_factor )        
          term = ( 1 + p_1*arg + p_2*arg**2 )/( 1 + q_1*arg + q_2*arg**2 + p_2*arg**3 )
          v_pcm(ip) = v_pcm(ip) + q_pcm(ia)*term/sqrt( tess(ia)%area*width_factor ) !< regularized PCM field
        end do
      end do
      v_pcm = M_TWO*v_pcm/sqrt(M_Pi)

    else

      do ia = 1, n_tess
        do ip = 1, mesh%np !running serially np=np_global
          v_pcm(ip) = v_pcm(ip) + q_pcm(ia)/sr_dist(ia,ip) !< standard PCM field
        end do
      end do

    end if
    POP_SUB(pcm_pot_rs)
  end subroutine pcm_pot_rs

  !> Generates the PCM response matrix. J. Tomassi et al. Chem. Rev. 105, 2999 (2005). 
  subroutine pcm_matrix(eps, tess, n_tess, pcm_mat )
    FLOAT, intent(in)           :: eps
    type(tessera_t), intent(in) :: tess(:)      !< (1:n_tess)
    integer, intent(in)         :: n_tess
    FLOAT, intent(out)          :: pcm_mat(:,:) !< (1:n_tess, 1:n_tess)

    integer :: i
    integer :: info
    integer, allocatable :: iwork(:)

    FLOAT, allocatable :: mat_tmp(:,:)

    PUSH_SUB(pcm_matrix)

    !> Conforming the S_I matrix
    SAFE_ALLOCATE( s_mat_act(1:n_tess, 1:n_tess) )
    call s_i_matrix(n_tess, tess)

    !> Defining the matrix S_E=S_I/eps
    SAFE_ALLOCATE( Sigma(1:n_tess, 1:n_tess) )
    Sigma = s_mat_act/eps

    !> Conforming the D_I matrix
    SAFE_ALLOCATE( d_mat_act(1:n_tess, 1:n_tess) )
    call d_i_matrix(n_tess, tess)

    !> Defining the matrix D_E=D_I 
    SAFE_ALLOCATE( Delta(1:n_tess, 1:n_tess) )
    Delta = d_mat_act

    !> Start conforming the PCM matrix
    pcm_mat = -d_mat_act

    do i=1, n_tess
       pcm_mat(i,i) = pcm_mat(i,i) + M_TWO*M_Pi  
    end do

    SAFE_DEALLOCATE_A(d_mat_act) 
     
    SAFE_ALLOCATE( iwork(1:n_tess) )

    !> Solving for X = S_I^-1*(2*Pi - D_I) 
    ! FIXME: use interface, or routine in lalg_adv_lapack_inc.F90
    call dgesv(n_tess, n_tess, s_mat_act, n_tess, iwork, pcm_mat, n_tess, info)        

    SAFE_DEALLOCATE_A(iwork)

    SAFE_DEALLOCATE_A(s_mat_act)

    !> Computing -S_E*S_I^-1*(2*Pi - D_I)
    pcm_mat = -matmul( Sigma, pcm_mat ) 

    do i=1, n_tess
      pcm_mat(i,i) = pcm_mat(i,i) + M_TWO*M_Pi
    end do

    pcm_mat = pcm_mat - Delta

    SAFE_ALLOCATE( mat_tmp(1:n_tess,1:n_tess) )
    mat_tmp = M_ZERO

    SAFE_ALLOCATE( d_mat_act(1:n_tess,1:n_tess) )  
    call d_i_matrix(n_tess, tess)

    mat_tmp = transpose(d_mat_act)        

    mat_tmp = matmul(Sigma, mat_tmp)

    mat_tmp = mat_tmp + M_TWO*M_Pi*Sigma

    SAFE_DEALLOCATE_A(d_mat_act)

    SAFE_ALLOCATE( s_mat_act(1:n_tess, 1:n_tess) )
    call s_i_matrix(n_tess, tess)

    mat_tmp = mat_tmp + M_TWO*M_Pi*s_mat_act - matmul(Delta, s_mat_act)

    SAFE_DEALLOCATE_A(s_mat_act)
    SAFE_DEALLOCATE_A(Sigma)
    SAFE_DEALLOCATE_A(Delta)

    SAFE_ALLOCATE( iwork(1:n_tess) )

    !> Solving for [(2*pi - D_E)*S_I + S_E*(2*Pi + D_I*)]*X = [(2*Pi - D_E) - S_E*S_I^-1*(2*Pi - D_I)]    
    call dgesv(n_tess, n_tess, mat_tmp, n_tess, iwork, pcm_mat, n_tess, info)		  		  

    SAFE_DEALLOCATE_A(iwork)
    SAFE_DEALLOCATE_A(mat_tmp)

    pcm_mat = -pcm_mat
    
    !   Testing
    do i=1, n_tess
      mat_gamess(i,:) = pcm_mat(i,:)/tess(i)%area
    end do
    
    POP_SUB(pcm_matrix)
    
  end subroutine pcm_matrix

  subroutine s_i_matrix(n_tess, tess)
    integer, intent(in)         :: n_tess
    type(tessera_t), intent(in) :: tess(:)
    
    integer :: ii, jj

    s_mat_act = M_ZERO 

    do ii = 1, n_tess
     do jj = ii, n_tess

        s_mat_act(ii,jj) = s_mat_elem_I(tess(ii), tess(jj))
        if (ii /= jj) s_mat_act(jj,ii) = s_mat_act(ii,jj) !symmetric matrix

     end do
    end do

  end subroutine s_i_matrix

  subroutine d_i_matrix(n_tess, tess)
    integer, intent(in)         :: n_tess
    type(tessera_t), intent(in) :: tess(:)

    integer :: ii, jj

    d_mat_act = M_ZERO 

    do ii = 1, n_tess
      do jj = 1, n_tess !< non-symmetric matrix
        
        d_mat_act(ii,jj) = d_mat_elem_I(tess(ii), tess(jj))
        
      end do
    end do
    
  end subroutine d_i_matrix

  !> electrostatic Green function in vacuo:
  !! G_I(r,r^\prime) = 1 / | r - r^\prime |
  FLOAT function s_mat_elem_I(tessi, tessj)
    type(tessera_t), intent(in) :: tessi
    type(tessera_t), intent(in) :: tessj
    
    FLOAT, parameter :: M_SD_DIAG    = CNST(1.0694)
    FLOAT, parameter :: M_DIST_MIN   = CNST(0.1)
    
    FLOAT :: diff(1:MAX_DIM)
    FLOAT :: dist
    FLOAT :: s_diag
    FLOAT :: s_off_diag

    s_diag = M_ZERO
    s_off_diag = M_ZERO

    diff = tessi%point - tessj%point

    dist = dot_product( diff, diff )
    dist = sqrt(dist)

    if  (dist == M_ZERO) then

      s_diag = M_SD_DIAG*sqrt( M_FOUR*M_Pi/tessi%area )
      s_mat_elem_I = s_diag
      
    else
      
      if ( dist > M_DIST_MIN ) s_off_diag = M_ONE/dist 
      s_mat_elem_I = s_off_diag
      
    end if

  end function s_mat_elem_I

  !> Gradient of the Green function in vacuo GRAD[G_I(r,r^\prime)]
  FLOAT function d_mat_elem_I(tessi, tessj)
    type(tessera_t), intent(in) :: tessi
    type(tessera_t), intent(in) :: tessj
    
    FLOAT, parameter :: M_SD_DIAG    = CNST(1.0694)
    FLOAT, parameter :: M_DIST_MIN   = CNST(0.04)

    FLOAT :: diff(1:MAX_DIM)
    FLOAT :: dist
    FLOAT :: d_diag
    FLOAT :: d_off_diag

    d_diag = M_ZERO
    d_off_diag = M_ZERO

    diff = tessi%point - tessj%point

    dist = dot_product( diff, diff )
    dist = sqrt(dist)

    if (dist == M_ZERO) then
      !> Diagonal matrix elements  
      d_diag = M_SD_DIAG*sqrt(M_FOUR*M_Pi*tessi%area)
      d_diag = -d_diag/(M_TWO*tessi%r_sphere)
      d_mat_elem_I = d_diag
      
    else
      !> off-diagonal matrix elements  
      if (dist > M_DIST_MIN) then
        d_off_diag = dot_product( diff, tessj%normal(:) )	
        d_off_diag = d_off_diag*tessj%area/dist**3
      end if

      d_mat_elem_I = d_off_diag
      
    end if
    
  end function d_mat_elem_I

  !> It builds the solute cavity surface and calculates the vertices,
  !! representative points and areas of the tesserae by using the 
  !! Gauss-Bonnet theorem.
  subroutine cav_gen(i_count, tess_sphere, nesf, sfe, nts, cts, unit_pcminfo)
    integer, intent(in)  :: i_count
    integer, intent(in)  :: tess_sphere
    integer, intent(in)  :: nesf
    integer, intent(out) :: nts
    integer, intent(in)  :: unit_pcminfo
 
    type(sphere_t),   intent(inout) :: sfe(:)
    type(tessera_t),  intent(out)   :: cts(:)

    integer, parameter :: dim_angles = 24
    integer, parameter :: dim_ten = 10
    integer, parameter :: dim_vertices = 122
    integer, parameter :: n_tess_sphere = 60
    integer, parameter :: max_vertices = 6
    integer, parameter :: mxts = 10000

    FLOAT :: thev(1:dim_angles)
    FLOAT :: fiv(1:dim_angles)
    FLOAT :: fir
    FLOAT :: cv(1:dim_vertices, 1:MAX_DIM)
    FLOAT :: th
    FLOAT :: fi
    FLOAT :: cth
    FLOAT :: sth

    FLOAT :: xctst(tess_sphere*n_tess_sphere)
    FLOAT :: yctst(tess_sphere*n_tess_sphere)
    FLOAT :: zctst(tess_sphere*n_tess_sphere)
    FLOAT :: ast(tess_sphere*n_tess_sphere)
    FLOAT :: nctst(MAX_DIM, tess_sphere*n_tess_sphere)
    FLOAT :: pts(1:MAX_DIM, 1:dim_ten)
    FLOAT :: pp(1:MAX_DIM)
    FLOAT :: pp1(1:MAX_DIM)
    FLOAT :: ccc(1:MAX_DIM, 1:dim_ten)
      
    integer :: idum(1:n_tess_sphere*max_vertices)
    integer :: jvt1(1:max_vertices,1:n_tess_sphere)
    integer :: isfet(1:dim_ten*dim_angles)
      
    integer :: ii
    integer :: ia
    integer :: ja
    integer :: nn
    integer :: nsfe
    integer :: its
    integer :: n1
    integer :: n2
    integer :: n3
    integer :: nv
    integer :: i_tes

    FLOAT :: xen
    FLOAT :: yen
    FLOAT :: zen
    FLOAT :: ren
    FLOAT :: area
    FLOAT :: test
    FLOAT :: test2
    FLOAT :: rij
    FLOAT :: dnorm

    FLOAT :: xi
    FLOAT :: yi
    FLOAT :: zi
    FLOAT :: xj
    FLOAT :: yj
    FLOAT :: zj

    FLOAT :: vol
    FLOAT :: stot
    FLOAT :: prod
    FLOAT :: dr  

    logical :: band_iter

    PUSH_SUB(gen_cav)
      
    !> Angles corresponding to the vertices and centres of a polyhedron
    !! within a sphere of unitary radius and centered at the origin
    data thev/ CNST(0.6523581398) , CNST(1.107148718)  , CNST(1.382085796) , &
               CNST(1.759506858)  , CNST(2.034443936)  , CNST(2.489234514) , &
                                    CNST(0.3261790699) , CNST(0.5535743589), &
               CNST(0.8559571251) , CNST(0.8559571251) , CNST(1.017221968) , &
               CNST(1.229116717)  , CNST(1.229116717)  , CNST(1.433327788) , &
               CNST(1.570796327)  , CNST(1.570796327)  , CNST(1.708264866) , &
               CNST(1.912475937)  , CNST(1.912475937)  , CNST(2.124370686) , &
               CNST(2.285635528)  , CNST(2.285635528)  , CNST(2.588018295) , &
               CNST(2.815413584) /
    data fiv/                       CNST(0.6283185307) , M_ZERO            , &
               CNST(0.6283185307) , M_ZERO             , CNST(0.6283185307), &
               M_ZERO             , CNST(0.6283185307) , M_ZERO, 	     &
               CNST(0.2520539002) , CNST(1.004583161)  , CNST(0.6283185307), &
               CNST(0.3293628477) , CNST(0.9272742138) , M_ZERO, 	     &
               CNST(0.3141592654) , CNST(0.9424777961) , CNST(0.6283185307), &
               CNST(0.2989556830) , CNST(0.9576813784) , M_ZERO, 	     &
               CNST(0.3762646305) , CNST(0.8803724309) , CNST(0.6283188307), &
               M_ZERO /
    data fir / CNST(1.256637061)  /

    !> the vector idum, contained in the matrix jvt1, indicates the vertices 
    !! of the tesserae (using less than 19 continuations)
    data (idum(ii),ii = 1, 280) /                                   &
       1, 6, 2, 32, 36, 37, 1, 2, 3, 33, 32, 38, 1, 3, 4, 34,         &
       33, 39, 1, 4, 5, 35, 34, 40, 1, 5, 6, 36, 35, 41, 7, 2, 6, 51, &
       42, 37, 8, 3, 2, 47, 43, 38, 9, 4, 3, 48, 44, 39, 10, 5, 4,    &
       49, 45, 40, 11, 6, 5, 50, 46, 41, 8, 2, 12, 62, 47, 52, 9,     &
       3, 13, 63, 48, 53, 10, 4, 14, 64, 49, 54, 11, 5, 15, 65, 50,   &
       55, 7, 6, 16, 66, 51, 56, 7, 12, 2, 42, 57, 52, 8, 13, 3,      &
       43, 58, 53, 9, 14, 4, 44, 59, 54, 10, 15, 5, 45, 60, 55, 11,   &
       16, 6, 46, 61, 56, 8, 12, 18, 68, 62, 77, 9, 13, 19, 69, 63,   &
       78, 10, 14, 20, 70, 64, 79, 11, 15, 21, 71, 65, 80, 7, 16,     &
       17, 67, 66, 81, 7, 17, 12, 57, 67, 72, 8, 18, 13, 58, 68, 73,  &
       9, 19, 14, 59, 69, 74, 10, 20, 15, 60, 70, 75, 11, 21, 16,     &
       61, 71, 76, 22, 12, 17, 87, 82, 72, 23, 13, 18, 88, 83, 73,    &
       24, 14, 19, 89, 84, 74, 25, 15, 20, 90, 85, 75, 26, 16, 21,    &
       91, 86, 76, 22, 18, 12, 82, 92, 77, 23, 19, 13, 83, 93, 78,    &
       24, 20, 14, 84, 94, 79, 25, 21, 15, 85, 95, 80, 26, 17, 16,    &
       86, 96, 81, 22, 17, 27, 102, 87, 97, 23, 18, 28, 103, 88, 98,  &
       24, 19, 29, 104, 89, 99, 25, 20, 30, 105, 90, 100, 26, 21,     &
       31, 106, 91, 101, 22, 28, 18, 92, 107, 98, 23, 29, 19, 93 /
    data (idum(ii),ii = 281,360) / 				      &
       108, 99, 24, 30, 20, 94, 109, 100, 25, 31, 21, 95, 110, 101,   &
       26, 27, 17, 96, 111, 97, 22, 27, 28, 107, 102, 112, 23, 28,    &
       29, 108, 103, 113, 24, 29, 30, 109, 104, 114, 25, 30, 31,      &
       110, 105, 115, 26, 31, 27, 111, 106, 116, 122, 28, 27, 117,    &
       118, 112, 122, 29, 28, 118, 119, 113, 122, 30, 29, 119, 120,   &
       114, 122, 31, 30, 120, 121, 115, 122, 27, 31, 121, 117, 116 /

    if (i_count == 0) then
     if (tess_sphere == 1) then
         write(unit_pcminfo,'(A1)')  '#' 
         write(unit_pcminfo,'(A34)') '# Number of tesserae / sphere = 60'
         write(unit_pcminfo,'(A1)')  '#' 
     else
         write(unit_pcminfo,'(A1)')  '#' 
         write(unit_pcminfo,'(A35)') '# Number of tesserae / sphere = 240' 
         write(unit_pcminfo,'(A1)')  '#' 
     endif 
    endif

   !> geometrical data are converted to Angstrom and back transformed
   !! to Bohr at the end of the subroutine.
    dr = CNST(0.01)
    dr = dr*P_a_B

    sfe(:)%x = sfe(:)%x*P_a_B
    sfe(:)%y = sfe(:)%y*P_a_B
    sfe(:)%z = sfe(:)%z*P_a_B
    sfe(:)%r = sfe(:)%r*P_a_B

    vol  = M_ZERO
    stot = M_ZERO
    jvt1 = reshape(idum,(/6,60/))

    !> Coordinates of vertices of tesserae in a sphere with unit radius.
    !! the matrix 'cv' and 'xc', 'yc', 'zc' conatin the vertices and 
    !! the centers of 240 tesserae. The matrix 'jvt1(i,j)' denotes the index
    !! of the i-th vertex of the j-th big tessera. On each big tessera
    !! the 6 vertices are ordered as follows:
    !!  
    !!                      1
    !!
    !!                   4     5
    !!
    !!                3     6     2

    cv(1,1)   =  M_ZERO
    cv(1,2)   =  M_ZERO
    cv(1,3)   =  M_ONE

    cv(122,1) =  M_ZERO
    cv(122,2) =  M_ZERO
    cv(122,3) = -M_ONE

    ii = 1
    do ia = 1, dim_angles
       th = thev(ia)
       fi = fiv(ia)
       cth = cos(th)
       sth = sin(th)
     do ja = 1, 5
        fi = fi + fir
        if (ja == 1) fi = fiv(ia)
        ii = ii + 1
        cv(ii,1) = sth*cos(fi)
        cv(ii,2) = sth*sin(fi)
        cv(ii,3) = cth
     enddo
    enddo

    !> Controls whether the tessera is covered or need to be reshaped it
    nn = 0
    do nsfe = 1, nesf
     xen = sfe(nsfe)%x
     yen = sfe(nsfe)%y
     zen = sfe(nsfe)%z
     ren = sfe(nsfe)%r

     xctst(:) = M_ZERO
     yctst(:) = M_ZERO
     zctst(:) = M_ZERO
     ast(:)   = M_ZERO

     do its = 1, n_tess_sphere 
       do i_tes = 1, tess_sphere
        if (tess_sphere == 1) then
         n1 = jvt1(1,its)
         n2 = jvt1(2,its)
         n3 = jvt1(3,its)
        else
         if (i_tes == 1)      then
          n1 = jvt1(1,its)
          n2 = jvt1(5,its)
          n3 = jvt1(4,its)
         elseif (i_tes == 2)  then 
          n1 = jvt1(4,its)
          n2 = jvt1(6,its)
          n3 = jvt1(3,its)
         elseif (i_tes == 3)  then
          n1 = jvt1(4,its)
          n2 = jvt1(5,its)
          n3 = jvt1(6,its)
         elseif (i_tes == 4)  then
          n1 = jvt1(2,its)
          n2 = jvt1(6,its)
          n3 = jvt1(5,its)
         endif
        endif

        pts(1,1) = cv(n1,1)*ren + xen
        pts(2,1) = cv(n1,3)*ren + yen
        pts(3,1) = cv(n1,2)*ren + zen

        pts(1,2) = cv(n2,1)*ren + xen
        pts(2,2) = cv(n2,3)*ren + yen
        pts(3,2) = cv(n2,2)*ren + zen

        pts(1,3) = cv(n3,1)*ren + xen
        pts(2,3) = cv(n3,3)*ren + yen
        pts(3,3) = cv(n3,2)*ren + zen

        pp(:)  = M_ZERO
        pp1(:) = M_ZERO
        nv = 3

        call subtessera(sfe, nsfe, nesf, nv, pts ,ccc, pp, pp1, area)

        if (area == M_ZERO) cycle

        xctst(tess_sphere*(its-1) + i_tes)   = pp(1)
        yctst(tess_sphere*(its-1) + i_tes)   = pp(2)
        zctst(tess_sphere*(its-1) + i_tes)   = pp(3)
        nctst(:,tess_sphere*(its-1) + i_tes) = pp1(:)
        ast(tess_sphere*(its-1) + i_tes)     = area
        isfet(tess_sphere*(its-1) + i_tes)   = nsfe

        enddo
       enddo !> loop through the tesseare on the sphere 'nsfe'

       do its = 1, n_tess_sphere*tess_sphere

        if (ast(its) == M_ZERO) cycle
        nn = nn + 1

        if (nn > mxts) then !> check the total number of tessera
         write(message(1),'(a,I5,a,I5)') "total number of tesserae", nn, ">",mxts
         call messages_warning(1)     
        endif

        if (i_count ==  1) then
          cts(nn)%point(1)  = xctst(its)
          cts(nn)%point(2)  = yctst(its)
          cts(nn)%point(3)  = zctst(its)
          cts(nn)%normal(:) = nctst(:,its)
          cts(nn)%area      = ast(its)
          cts(nn)%r_sphere  = sfe(isfet(its))%r
        endif

       enddo
      enddo !> loop through the spheres

     nts = nn

     if (i_count == 1) then

     !> checks if two tesseare are too close
     test = CNST(0.1)
     test2 = test*test

     band_iter = .false.
     do while (.not.(band_iter))
       band_iter = .true.

       loop_ia: do ia = 1, nts-1
        if (cts(ia)%area == M_ZERO) cycle
         xi = cts(ia)%point(1)
         yi = cts(ia)%point(2)
         zi = cts(ia)%point(3)

        loop_ja: do ja = ia+1, nts
         if (cts(ja)%area == M_ZERO) cycle
          xj = cts(ja)%point(1)
          yj = cts(ja)%point(2)
          zj = cts(ja)%point(3)

          rij = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2

          if (rij > test2) cycle

          write(unit_pcminfo,'(A40,I4,A5,I4,A4,F8.4,A13,F8.4,A3)' ) &
                              '# Warning: The distance between tesserae', &
                               ia,' and ', ja,' is ',sqrt(rij),' A, less than', test,' A.'

         !> calculating the coordinates of the new tessera weighted by the areas
          xi = (xi*cts(ia)%area + xj*cts(ja)%area) / (cts(ia)%area + cts(ja)%area)
          yi = (yi*cts(ia)%area + yj*cts(ja)%area) / (cts(ia)%area + cts(ja)%area)
          zi = (zi*cts(ia)%area + zj*cts(ja)%area) / (cts(ia)%area + cts(ja)%area)

          cts(ia)%point(1) = xi
          cts(ia)%point(2) = yi
          cts(ia)%point(3) = zi

          !> calculating the normal vector of the new tessera weighted by the areas
          cts(ia)%normal = (cts(ia)%normal*cts(ia)%area + cts(ja)%normal*cts(ja)%area)
          dnorm = sqrt( dot_product(cts(ia)%normal, cts(ia)%normal) )
          cts(ia)%normal = cts(ia)%normal/dnorm

          !> calculating the sphere radius of the new tessera weighted by the areas
          cts(ia)%r_sphere = ( cts(ia)%r_sphere*cts(ia)%area + cts(ja)%r_sphere*cts(ja)%area ) / &
                             ( cts(ia)%area + cts(ja)%area )

          !> calculating the area of the new tessera
          cts(ia)%area = cts(ia)%area + cts(ja)%area

          !> deleting tessera ja
          do ii = ja+1, nts
           cts(ii-1) = cts(ii)
          enddo
          nts = nts -1 
          band_iter = .false.
          exit loop_ia

         enddo loop_ja
        enddo loop_ia
       enddo !> while loop

      !> Calculates the cavity volume: vol = \sum_{its=1}^nts A_{its} s*n/3.
      vol = M_ZERO
      do its = 1, nts
       prod = dot_product( cts(its)%point, cts(its)%normal ) 
       vol  = vol + cts(its)%area * prod / M_THREE
       stot = stot + cts(its)%area
      enddo

      write(unit_pcminfo, '(A2)')  '# '
      write(unit_pcminfo, '(A29,I4)')    '# Total number of tesserae = ', nts
      write(unit_pcminfo, '(A30,F12.6)') '# Cavity surface area (A^2) = ' , stot
      write(unit_pcminfo, '(A24,F12.6)') '# Cavity volume (A^3) = '       , vol
      write(unit_pcminfo, '(A2)')  '# '

      !> transforms results into Bohr.
      cts(:)%area     = cts(:)%area*P_Ang*P_Ang
      cts(:)%point(1) = cts(:)%point(1)*P_Ang
      cts(:)%point(2) = cts(:)%point(2)*P_Ang
      cts(:)%point(3) = cts(:)%point(3)*P_Ang
      cts(:)%r_sphere = cts(:)%r_sphere*P_Ang
     endif

     sfe(:)%x=sfe(:)%x*P_Ang
     sfe(:)%y=sfe(:)%y*P_Ang
     sfe(:)%z=sfe(:)%z*P_Ang
     sfe(:)%r=sfe(:)%r*P_Ang

  return
  POP_SUB(cav_gen)
  end subroutine cav_gen

  !> find the uncovered region for each tessera and computes the area,
  !! the representative point (pp) and the unitary normal vector (pp1)
  subroutine subtessera(sfe, ns, nesf, nv, pts, ccc, pp, pp1, area)
    type(sphere_t), intent(in) :: sfe(:) !< (1:nesf)
    integer, intent(in)        :: ns 
    integer, intent(in)        :: nesf
    integer, intent(inout)     :: nv
    FLOAT, intent(inout)       :: pts(:,:) !< (1:MAX_DIM, 1:dim_ten)
    FLOAT, intent(out)         :: ccc(:,:) !< (1:MAX_DIM, 1:dim_ten)
    FLOAT, intent(out)         :: pp(:)    !< (1:MAX_DIM)
    FLOAT, intent(out)         :: pp1(:)   !< (1:MAX_DIM)
    FLOAT, intent(out)         :: area

    FLOAT, parameter   :: tol = -CNST(1.0e-10)
    integer, parameter :: dim_ten = 10

    integer :: intsph(1:dim_ten)
    integer :: nsfe1
    integer :: na
    integer :: icop
    integer :: ll
    integer :: iv1
    integer :: iv2
    integer :: ii
    integer :: icut
    integer :: jj

    FLOAT  :: p1(1:MAX_DIM)
    FLOAT  :: p2(1:MAX_DIM)
    FLOAT  :: p3(1:MAX_DIM)
    FLOAT  :: p4(1:MAX_DIM)
    FLOAT  :: point(1:MAX_DIM)
    FLOAT  :: pscr(1:MAX_DIM,1:dim_ten)
    FLOAT  :: cccp(1:MAX_DIM,1:dim_ten)
    FLOAT  :: pointl(1:MAX_DIM,1:dim_ten)
    FLOAT  :: diff(1:MAX_DIM)

    integer :: ind(1:dim_ten)
    integer :: ltyp(1:dim_ten)
    integer :: intscr(1:dim_ten)

    FLOAT  :: delr
    FLOAT  :: delr2
    FLOAT  :: rc
    FLOAT  :: rc2
    FLOAT  :: dnorm
    FLOAT  :: dist
    FLOAT  :: de2

    PUSH_SUB(subtessera)

    area = M_ZERO
    do jj=1, 3
      ccc(1,jj) = sfe(ns)%x
      ccc(2,jj) = sfe(ns)%y
      ccc(3,jj) = sfe(ns)%z
    enddo

    intsph = ns
    do nsfe1 = 1, nesf 
     if (nsfe1 == ns) cycle
     do jj =1, nv
      intscr(jj) = intsph(jj)
      pscr(:,jj) = pts(:,jj)
      cccp(:,jj) = ccc(:,jj)
     enddo

     icop = 0
     ind = 0
     ltyp = 0

     do ii = 1, nv
      delr2 = ( pts(1,ii) - sfe(nsfe1)%x )**2 + ( pts(2,ii) - sfe(nsfe1)%y )**2 + &
	      ( pts(3,ii) - sfe(nsfe1)%z )**2
      delr = sqrt(delr2)
      if (delr < sfe(nsfe1)%r) then
	ind(ii) = 1
	icop = icop + 1
      endif
     enddo

     if (icop == nv) return

     do ll = 1, nv
      iv1 = ll
      iv2 = ll+1
      if (ll == nv) iv2 = 1
      IF ( (ind(iv1) == 1) .and. (ind(iv2) == 1) ) then
	ltyp(ll) = 0
      else if ( (ind(iv1) == 0) .and. (ind(iv2) == 1) ) then
	ltyp(ll) = 1
      else if ( (ind(iv1) == 1) .and. (ind(iv2) == 0) ) then
	ltyp(ll) = 2
      else if ( (ind(iv1) == 0) .and. (ind(iv2) == 0) ) then
	ltyp(ll) = 4
	diff = ccc(:,ll) - pts(:,ll)
	rc2 = dot_product(diff,diff)
	rc = sqrt(rc2)

	do ii = 1, 11
	 point = pts(:,iv1) + ii * (pts(:,iv2) - pts(:,iv1)) / 11
	 point = point - CCC(:,ll)
	 dnorm = sqrt( dot_product(point, point) )
	 point = point * rc / dnorm + CCC(:,ll)

	 dist = sqrt(  (point(1) - sfe(nsfe1)%x)**2 + ( point(2) - sfe(nsfe1)%y)**2 &
						    + ( point(3) - sfe(nsfe1)%z)**2  )
       
	if ( (dist - sfe(nsfe1)%r) < tol) then
	  ltyp(ll) = 3
	  pointl(:,ll) = point
	  exit
	endif

	enddo
      endif
    enddo

    icut = 0
    do ll = 1, nv
      if ( (ltyp(ll) == 1) .or. (ltyp(ll) == 2) ) icut = icut + 1
      if (ltyp(ll) == 3) icut = icut + 2
    enddo
    icut = icut / 2
    if (icut > 1) return

    na = 1
    do ll = 1, nv

      if (ltyp(ll) == 0) cycle
      iv1 = ll
      iv2 = ll + 1
      if (ll == nv) iv2 = 1

      if (ltyp(ll) == 1) then
       pts(:,na) = pscr(:,iv1)
       ccc(:,na) = cccp(:,iv1)
       intsph(na) = intscr(iv1)
       na = na + 1
       p1 = pscr(:,IV1)
       p2 = pscr(:,IV2)
       p3 = cccp(:,IV1)

       call inter(sfe, p1, p2, p3, p4, nsfe1, 0)
       pts(:,na) = p4

       de2 = ( sfe(nsfe1)%x - sfe(ns)%x )**2 + ( sfe(nsfe1)%y - sfe(ns)%y )**2 + &
	     ( sfe(nsfe1)%z - sfe(ns)%z )**2

       ccc(1,na) = sfe(ns)%x + ( sfe(nsfe1)%x - sfe(ns)%x)* &
			       ( sfe(ns)%r**2 - sfe(nsfe1)%r**2 + de2 ) / (M_TWO*de2)

       ccc(2,na) = sfe(ns)%y + ( sfe(nsfe1)%y - sfe(ns)%y)* &
			       ( sfe(ns)%r**2 - sfe(nsfe1)%r**2 + de2 ) / (M_TWO*de2)

       ccc(3,na) = sfe(ns)%z + ( sfe(nsfe1)%z - sfe(ns)%z)* &
			       ( sfe(ns)%r**2 - sfe(nsfe1)%r**2 + de2 ) / (M_TWO*de2)

       intsph(na) = nsfe1
       na = na + 1
      endif

      if (ltyp(ll) == 2) then
       p1 = pscr(:,iv1)
       p2 = pscr(:,iv2)
       p3 = cccp(:,iv1)

       call inter( sfe, p1, p2, p3, p4, nsfe1, 1 )
       pts(:,na) = p4
       ccc(:,na) = cccp(:,iv1)
       intsph(na) = intscr(iv1)
       na = na + 1
      endif

      if (ltyp(ll) == 3) then
       pts(:,na) = pscr(:,iv1)
       ccc(:,na) = cccp(:,iv1)
       intsph(na) = intscr(iv1)
       na = na + 1
       p1 = pscr(:,iv1)
       p2 = pointl(:,ll)
       p3 = cccp(:,iv1)

       call inter( sfe, p1, p2, p3, p4, nsfe1, 0 )
       pts(:,na) = p4

       de2 = ( sfe(nsfe1)%x - sfe(ns)%x )**2 + ( sfe(nsfe1)%y - sfe(ns)%y )**2 + &
					       ( sfe(nsfe1)%z - sfe(ns)%z )**2

       ccc(1,na) = sfe(ns)%x + ( sfe(nsfe1)%x - sfe(ns)%x )* &
			      ( sfe(ns)%r**2 - sfe(nsfe1)%r**2 + de2 ) / (M_TWO*de2)

       ccc(2,na) = sfe(ns)%y + ( sfe(nsfe1)%y - sfe(ns)%y )* &
			      ( sfe(ns)%r**2 - sfe(nsfe1)%r**2 + de2 ) / (M_TWO*de2)

       ccc(3,na) = sfe(ns)%z + ( sfe(nsfe1)%z - sfe(ns)%z )* &
			      ( sfe(ns)%r**2 - sfe(nsfe1)%r**2 + de2 ) / (M_TWO*de2)

       intsph(na) = nsfe1
       na = na + 1
       p1 = pointl(:,ll)
       p2 = pscr(:,iv2)
       p3 = cccp(:,iv1)

       call inter( sfe, p1, p2, p3, p4, nsfe1, 1 )
       pts(:,na) = p4
       ccc(:,na) = cccp(:,iv1)
       intsph(na) = intscr(iv1)
       na = na + 1
      endif

      if (ltyp(ll) == 4) then
       pts(:,na) = pscr(:,iv1)
       ccc(:,na) = cccp(:,iv1)
       intsph(na) = intscr(iv1)
       na = na + 1
      endif
    enddo

    nv = na - 1
    if (nv > 10) then
     message(1) = "Too many vertices on the tessera"
     call messages_fatal(1)     
    endif
   enddo

   call gaubon( sfe, nv, ns, pts, ccc, pp, pp1, area, intsph)
  return
  POP_SUB(subtessera)
  end subroutine subtessera

!    !> Finds the point 'p4', on the arc 'p1'-'p2' developed from 'p3',
!    !! which is on the surface of sphere 'ns'. p4 is a linear combination 
     !! of p1 and p2 with the 'alpha' parameter optimized iteratively.
  subroutine inter( sfe, p1, p2, p3, p4, ns, ia)
    type(sphere_t), intent(in) :: sfe(:) !< (1:nesf)
    FLOAT, intent(in)          :: p1(:)  !< (1:MAX_DIM)
    FLOAT, intent(in)          :: p2(:)  !< (1:MAX_DIM)
    FLOAT, intent(in)          :: p3(:)  !< (1:MAX_DIM)
    FLOAT, intent(out)         :: p4(:)  !< (1:MAX_DIM)
      
    integer, intent(in) :: ns
    integer, intent(in) :: ia

    FLOAT, parameter :: tol = CNST(1.0e-08)

    integer :: m_iter
    FLOAT  :: r
    FLOAT  :: alpha
    FLOAT  :: delta
    FLOAT  :: dnorm
    FLOAT  :: diff
    FLOAT  :: diff_vec(1:MAX_DIM)

    logical :: band_iter

    diff_vec = p1 - p3
    r = sqrt( dot_product(diff_vec, diff_vec) )

    alpha = M_HALF
    delta = M_ZERO
    m_iter = 1

    band_iter = .false.
    do while(.not.(band_iter))
     if (m_iter > 1000) then
      message(1) = "Too many iterations inside subrotuine inter"
      call messages_fatal(1)     
     endif

     band_iter = .true.

     alpha = alpha + delta
     dnorm = M_ZERO

     p4 = p1 + alpha*(p2-p1)-p3
     dnorm = sqrt( dot_product(p4,p4) )
     p4 = p4*r/dnorm + p3
     diff =( p4(1) - sfe(ns)%x )**2 + ( p4(2) - sfe(ns)%y )**2 + ( p4(3) - sfe(ns)%z )**2
     diff = sqrt(diff) - sfe(ns)%r

     if ( abs(diff) < tol ) return

     if (ia == 0) then
      if (diff > M_ZERO) delta =  M_ONE/(M_TWO**(m_iter+1))
      if (diff < M_ZERO) delta = -M_ONE/(M_TWO**(m_iter+1))
      m_iter = m_iter + 1
      band_iter = .false.
     endif

     if (ia == 1) then
      if (diff > M_ZERO) delta = -M_ONE/(M_TWO**(m_iter+1))
      if (diff < M_ZERO) delta =  M_ONE/(M_TWO**(m_iter+1))
      m_iter = m_iter + 1
      band_iter = .false.
     endif
    enddo

  return
  end subroutine inter

    !> Use the Gauss-Bonnet theorem to calculate the area of the 
    !! tessera with vertices 'pts(3,nv)'. 
    !! Area = R^2 [ 2pi + S(Phi(N)cosT(N)) - S(Beta(N)) ]
    !! Phi(n): length of the arc in radians of the side 'n'. 
    !! T(n): azimuthal angle for the side 'n'
    !! Beta(n): external angle respect to vertex 'n'.
  subroutine gaubon( sfe, nv, ns, pts, ccc, pp, pp1, area, intsph )
    type(sphere_t), intent(in) :: sfe(:)    !< (1:nesf)
    FLOAT, intent(in)          :: pts(:,:)  !< (1:MAX_DIM,1:dim_ten) 
    FLOAT, intent(in)          :: ccc(:,:)  !< (1:MAX_DIM,1:dim_ten)
    FLOAT, intent(inout)       :: pp(:)     !< (1:MAX_DIM)
    FLOAT, intent(inout)       :: pp1(:)    !< (1:MAX_DIM)
    integer, intent(in)        :: intsph(:) !< (1:dim_ten)
    FLOAT, intent(out)         :: area
    integer, intent(in)        :: nv
    integer, intent(in)        :: ns

    FLOAT ::  p1(1:MAX_DIM)
    FLOAT ::  p2(1:MAX_DIM)
    FLOAT ::  p3(1:MAX_DIM)
    FLOAT ::  u1(1:MAX_DIM)
    FLOAT ::  u2(1:MAX_DIM)
    FLOAT ::  point_1(1:MAX_DIM)
    FLOAT ::  point_2(1:MAX_DIM)
    FLOAT :: tpi
    FLOAT :: sum1
    FLOAT :: dnorm
    FLOAT :: dnorm1
    FLOAT :: dnorm2
    FLOAT :: scal
    FLOAT :: cosphin
    FLOAT :: phin
    FLOAT :: costn
    FLOAT :: sum2
    FLOAT :: betan

    integer :: nsfe1
    integer :: ia
    integer :: nn
    integer :: n0
    integer :: n1
    integer :: n2

    PUSH_SUB(gaubon)

    tpi=2*M_Pi
    sum1 = M_ZERO
    point_1 = M_ZERO
    point_2 = M_ZERO
    do nn = 1, nv
     point_1 = pts(:,nn) - ccc(:,nn)
     if (nn < nv) then
      point_2 = pts(:,nn+1) - ccc(:,nn)
     else
      point_2 = pts(:,1) - ccc(:,nn)
     endif

     dnorm1 = sqrt( dot_product(point_1, point_1) )
     dnorm2 = sqrt( dot_product(point_2, point_2) )
     cosphin = dot_product(point_1, point_2) / (dnorm1*dnorm2)

     if (cosphin >  M_ONE) cosphin =  M_ONE
     if (cosphin < -M_ONE) cosphin = -M_ONE

     phin = acos(cosphin)
     nsfe1 = intsph(nn)

     point_1(1) = sfe(nsfe1)%x - sfe(ns)%x
     point_1(2) = sfe(nsfe1)%y - sfe(ns)%y
     point_1(3) = sfe(nsfe1)%z - sfe(ns)%z

     dnorm1 = sqrt( dot_product(point_1, point_1) )

     if (dnorm1 == M_ZERO) dnorm1 = M_ONE
      point_2(1) = pts(1,nn) - sfe(ns)%x
      point_2(2) = pts(2,nn) - sfe(ns)%y
      point_2(3) = pts(3,nn) - sfe(ns)%z

      dnorm2 = sqrt( dot_product(point_2, point_2) )
      costn  = dot_product(point_1, point_2)/(dnorm1*dnorm2)
      sum1 = sum1 + phin * costn
    enddo

    sum2 = M_ZERO
    !> Loop over the vertices
    do nn = 1, nv
      p1 = M_ZERO
      p2 = M_ZERO    
      p3 = M_ZERO  

      n1 = nn
      if (nn > 1)   n0 = nn - 1
      if (nn == 1)  n0 = nv
      if (nn < nv)  n2 = nn + 1
      if (nn == nv) n2 = 1
      
      p1 = pts(:,n1) - ccc(:,n0)
      p2 = pts(:,n0) - ccc(:,n0)
      call vecp(p1, p2, p3, dnorm)
      p2 = p3    

      call vecp(p1, p2, p3, dnorm)
      u1 = p3/dnorm

      p1 = pts(:,n1) - ccc(:,n1)
      p2 = pts(:,n2) - ccc(:,n1)
      call vecp(p1, p2, p3, dnorm)
      p2 = p3

      call vecp(p1, p2, p3, dnorm)
      u2 = p3/dnorm

      betan = acos( dot_product(u1, u2) )
      sum2 = sum2 + (M_Pi - betan)
    enddo

    !> computes the area of the tessera
    area = sfe(ns)%r*sfe(ns)%r*(tpi + sum1 - sum2)

    !> computes the representative point
    pp = M_ZERO

     do ia = 1, nv
      pp(1) = pp(1) + ( pts(1,ia) - sfe(ns)%x )
      pp(2) = pp(2) + ( pts(2,ia) - sfe(ns)%y )
      pp(3) = pp(3) + ( pts(3,ia) - sfe(ns)%z )
     enddo

     dnorm = M_ZERO
     dnorm = sqrt( dot_product(pp,pp) )

     pp(1) = sfe(ns)%x + pp(1) * sfe(ns)%r / dnorm
     pp(2) = sfe(ns)%y + pp(2) * sfe(ns)%r / dnorm
     pp(3) = sfe(ns)%z + pp(3) * sfe(ns)%r / dnorm

     !> finds the internal normal at the representative point
     pp1(1) = (pp(1) - sfe(ns)%x) / sfe(ns)%r
     pp1(2) = (pp(2) - sfe(ns)%y) / sfe(ns)%r
     pp1(3) = (pp(3) - sfe(ns)%z) / sfe(ns)%r

     !> If the area of the tessera is negative (0^-), due to numerical errors, is discarded
     if (area < M_ZERO) area = M_ZERO

  return
  end subroutine gaubon

  !> calculates the vectorial product p3 = p1 x p2
  subroutine vecp(p1, p2, p3, dnorm)
    FLOAT, intent(in)  :: P1(:) !< (1:MAX_DIM)
    FLOAT, intent(in)  :: P2(:) !< (1:MAX_DIM)
    FLOAT, intent(out) :: P3(:) !< (1:MAX_DIM)
    FLOAT, intent(out) :: dnorm

    p3(1) = p1(2)*p2(3) - p1(3)*p2(2)
    p3(2) = p1(3)*p2(1) - p1(1)*p2(3)
    P3(3) = p1(1)*p2(2) - p1(2)*p2(1)

    dnorm = M_ZERO
    dnorm = sqrt( dot_product(p3, p3) )

  return
  POP_SUB(gaubon)
  end subroutine vecp

  subroutine pcm_end(pcm)
    type(pcm_t), intent(inout) :: pcm
    
    PUSH_SUB(pcm_end)
    
    SAFE_DEALLOCATE_A(pcm%spheres)
    SAFE_DEALLOCATE_A(pcm%tess)
    SAFE_DEALLOCATE_A(pcm%matrix)
    SAFE_DEALLOCATE_A(pcm%q_e)
    SAFE_DEALLOCATE_A(pcm%q_n) 
    SAFE_DEALLOCATE_A(pcm%v_e)
    SAFE_DEALLOCATE_A(pcm%v_n)
    SAFE_DEALLOCATE_A(pcm%v_e_rs)
    SAFE_DEALLOCATE_A(pcm%v_n_rs)
    SAFE_DEALLOCATE_A(pcm%ind_vh)
    SAFE_DEALLOCATE_A(pcm%arg_li)
    SAFE_DEALLOCATE_A(sr_dist) 
    
    call io_close(pcm%info_unit)
    
    POP_SUB(pcm_end)
  end subroutine pcm_end
  
end module pcm_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

