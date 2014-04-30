!! Copyright (C) 2014 Alain Delgado
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
  use datasets_m
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

  public :: pcm_t,                &
            pcm_init,             &
            pcm_end,              &
            pcm_charges,          &
            pcm_pot_rs,           &
            pcm_classical_energy, &
            v_nuclei_cav,         &
            v_electrons_cav_li,   &
            v_electrons_cav


 !> The cavity hosting the solute molecule is built from a set of 
 !! interlocking spheres with optimized radii centered at the nuclear positions.  
  type, public :: sphere_t  
    FLOAT :: x(MAX_DIM) !< center of the sphere
    FLOAT :: r          !< radius of the sphere (different for each species)
  end type

 !> The resulting cavity is discretized by a set of tesserae.  
  type, public :: tessera_t
    FLOAT :: x(MAX_DIM) !< representative point of the tessera 
    FLOAT :: area       !< area of the tessera
    FLOAT :: n(MAX_DIM) !< unitary outgoing vector normal to the tessera surface 
    FLOAT :: r_sphere   !< radius of the sphere to which the tessera belongs
  end type tessera_t

  type pcm_t
    logical                      :: run_pcm       !< If True, PCM calculation is enabled
    integer                      :: n_spheres     !< Number of spheres used to build the VdW cavity
    integer                      :: n_tesserae    !< Total number of tesserae
    type(sphere_t), allocatable  :: spheres(:)    !< See definition for type sphere_t
    type(tessera_t), allocatable :: tess(:)       !< See definition for type tessera_t
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
    FLOAT                        :: epsilon_infty !< Infnite-frequency dielectric constant of the solvent
    FLOAT                        :: gaussian_width!< Parameter to change the width of density of polarization charges  
    integer                      :: n_vertices    !< Number of grid points used to interpolate the Hartree potential
                                                  !! at the tesserae representative points 
    integer, allocatable         :: ind_vh(:,:)   !< Grid points used during interpolation 
    integer                      :: info_unit     !< unit for pcm info file 
    character(len=80)            :: input_cavity  !< file name containing the geometry of the VdW cavity
  end type pcm_t

  FLOAT, allocatable :: s_mat_act(:,:) !< S_I matrix 
  FLOAT, allocatable :: d_mat_act(:,:) !< D_I matrix
  FLOAT, allocatable :: Sigma(:,:)     !< S_E matrix
  FLOAT, allocatable :: Delta(:,:)     !< D_E matrix in JCP 139, 024105 (2013).

  integer :: nearest_idx_unit   
  integer :: idx_from_coord_unit

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
    integer :: cav_unit
    integer :: cav_unit_test
    integer :: pcmmat_unit
    integer :: iunit
    integer :: vdw_unit
    integer :: grid_unit
    integer :: pcm_calc_mode

    integer, parameter :: mxts = 10000

    FLOAT :: rcav_C
    FLOAT :: rcav_O
    FLOAT :: rcav_N
    FLOAT :: rcav_S
    FLOAT :: rcav_F
    FLOAT :: rcav_Na
    FLOAT :: rcav_Cl

    PUSH_SUB(pcm_init)

    !%Variable Solvation
    !%Type logical
    !%Default no
    !%Section Hamiltonian::PCM
    !%Description
    !% If true, the calculation is performed accounting for solvation effects
    !% in the framework of Integral Equation Formalism Polarizable Continuum Model IEF-PCM
    !% (Chem. Rev. 105, 2999 (2005), J. Chem. Phys. 107, 3032 (1997),
    !% J. Chem. Phys. 139, 024105 (2013)). At the moment, this option is available 
    !% only for ground state calculations. Experimental.
    !%End

    call parse_logical(datasets_check('Solvation'), .false., pcm%run_pcm)
    if (pcm%run_pcm) then
      if (grid%sb%box_shape /= MINIMUM) then
          message(1) = "PCM is only available for BoxShape = minimum"
          call messages_fatal(1)
      else 
          call messages_experimental("polarizable continuum model")
      endif
    else
      POP_SUB(pcm_init)
      return
    endif

    rcav_C  = CNST(2.4)*P_Ang    ! 
    rcav_O  = CNST(1.8)*P_Ang    !    
    rcav_N  = CNST(1.9)*P_Ang    !
    rcav_S  = CNST(2.0175)*P_Ang ! Angstrom -> Bohr 
    rcav_F  = CNST(1.682)*P_Ang  !
    rcav_Na = CNST(2.772)*P_Ang  !  
    rcav_Cl = CNST(2.172)*P_Ang  !

    !%Variable SolventDielectricConstant
    !%Type float
    !%Default 1.0 (gas phase)
    !%Section Hamiltonian::PCM
    !%Description
    !% Static dielectric constant of the solvent (\epsilon_0).
    !%End
    call parse_float(datasets_check('SolventDielectricConstant'), CNST(1.0), pcm%epsilon_0)

    !%Variable SmearingFactor
    !%Type float
    !%Default 1.0 (The width of the Gaussian density of polarization charges is the area of each tessera)
    !%Section Hamiltonian::PCM
    !%Description
    !% Parameter used to control the width of the Gaussian density of polarization charges.
    !%End
    call parse_float(datasets_check('SmearingFactor'), CNST(1.0), pcm%gaussian_width)

    pcm%n_spheres = 0
    do ia = 1, geo%natoms
       if (geo%atom(ia)%label == 'H') cycle
       pcm%n_spheres = pcm%n_spheres + 1 !counting the number of species different from Hydrogen
    enddo

    SAFE_ALLOCATE( pcm%spheres(1:pcm%n_spheres) )
    
    pcm%n_spheres = 0
    do ia = 1, geo%natoms
       
       if (geo%atom(ia)%label == 'H') cycle
       pcm%n_spheres = pcm%n_spheres + 1
      
       ! These coordinates are already in atomic units (Bohr)
       pcm%spheres(pcm%n_spheres)%x = geo%atom(ia)%x

       if (geo%atom(ia)%label == 'C')  pcm%spheres(pcm%n_spheres)%r = rcav_C
       if (geo%atom(ia)%label == 'O')  pcm%spheres(pcm%n_spheres)%r = rcav_O
       if (geo%atom(ia)%label == 'N')  pcm%spheres(pcm%n_spheres)%r = rcav_N
       if (geo%atom(ia)%label == 'S')  pcm%spheres(pcm%n_spheres)%r = rcav_S
       if (geo%atom(ia)%label == 'F')  pcm%spheres(pcm%n_spheres)%r = rcav_F
       if (geo%atom(ia)%label == 'Na') pcm%spheres(pcm%n_spheres)%r = rcav_Na
       if (geo%atom(ia)%label == 'Cl') pcm%spheres(pcm%n_spheres)%r = rcav_Cl
    enddo

    call io_mkdir('pcm')

    pcm%info_unit = io_open('pcm/pcm_info.out', action='write')

    write(pcm%info_unit,'(2X,A33)') 'Configuration: Molecule + Solvent'
    write(pcm%info_unit,'(2X,A33)') '---------------------------------'
    write(pcm%info_unit,'(2X,A19,F12.3)') 'Epsilon(Solvent) = ', pcm%epsilon_0
    write(pcm%info_unit,'(2X)') 
    write(pcm%info_unit,'(2X,A33,I4)') 'Number of interlocking spheres = ', pcm%n_spheres
    write(pcm%info_unit,'(2X)')  

    write(pcm%info_unit,'(2X,A6,3X,A7,8X,A26,20X,A10)') 'SPHERE', 'ELEMENT', 'CENTER  (X,Y,Z) (A)', 'RADIUS (A)'
    write(pcm%info_unit,'(2X,A6,3X,A7,4X,A43,7X,A10)') '------', '-------', &
                              '-------------------------------------------', '----------'  

    pcm%n_spheres = 0
    do ia = 1, geo%natoms
       if (geo%atom(ia)%label == 'H') cycle
       pcm%n_spheres = pcm%n_spheres + 1       

       write(pcm%info_unit,'(2X,I3,9X,A2,3X,F14.8,2X,F14.8,2X,F14.8,3X,F14.8)') pcm%n_spheres,            &
                              						        geo%atom(ia)%label,       &
                                                                                geo%atom(ia)%x*P_a_B,     &
                                                                                pcm%spheres(pcm%n_spheres)%r*P_a_B
    enddo
    write(pcm%info_unit,'(2X)')  

    !%Variable CavityGeometry
    !%Type string
    !%Section Hamiltonian::PCM
    !%Description
    !% Name of the file containing the geometry of the Van der Waals surface that defines the cavity hosting
    !% the solute molecule in PCM calculations. Tesserae representative points must be in atomic units!.
    !%End
    call parse_string(datasets_check('CavityGeometry'), '', pcm%input_cavity)

    iunit = io_open(trim(pcm%input_cavity), status='old', action='read')
    
    read(iunit,*) pcm%n_tesserae

    if (pcm%n_tesserae.gt.mxts) then
        write(message(1),'(a,I5,a,I5)') "total number of tesserae", pcm%n_tesserae, ">",mxts
        call messages_warning(1)     
    endif

    SAFE_ALLOCATE( pcm%tess(1:pcm%n_tesserae)  )

    do ia=1, pcm%n_tesserae
       read(iunit,*) pcm%tess(ia)%x(1)
    enddo
	
    do ia=1, pcm%n_tesserae
       read(iunit,*) pcm%tess(ia)%x(2)
    enddo
	
    do ia=1, pcm%n_tesserae
       read(iunit,*) pcm%tess(ia)%x(3)
    enddo

    do ia=1, pcm%n_tesserae
       read(iunit,*) pcm%tess(ia)%area
    enddo

    do ia=1, pcm%n_tesserae
       read(iunit,*) pcm%tess(ia)%r_sphere
    enddo

    do ia=1, pcm%n_tesserae
       read(iunit,*) pcm%tess(ia)%n
    enddo
    
    call io_close(iunit)

    message(1) = "Info: van der Waals surface has been read from " // trim(pcm%input_cavity)
    call messages_info(1)

    cav_unit_test = io_open('pcm/cavity_mol.xyz', action='write')

    write (cav_unit_test,'(2X,I4)') pcm%n_tesserae + geo%natoms
    write (cav_unit_test,'(2X)')

    do ia=1, pcm%n_tesserae
       write(cav_unit_test,'(2X,A2,3X,4f15.8,3X,4f15.8,3X,4f15.8)') 'H', pcm%tess(ia)%x*P_a_B
    enddo

    do ia=1, geo%natoms
       write(cav_unit_test,'(2X,A2,3X,4f15.8,3X,4f15.8,3X,4f15.8)') geo%atom(ia)%label,      &
                                                                    geo%atom(ia)%x*P_a_B
    enddo

    call io_close(cav_unit_test)

    pcm%n_vertices = 8
    SAFE_ALLOCATE( pcm%ind_vh(1:pcm%n_tesserae, 1:pcm%n_vertices) )
    pcm%ind_vh = INT(M_ZERO)

    SAFE_ALLOCATE( pcm%arg_li(1:pcm%n_tesserae, 1:MAX_DIM) )
    pcm%arg_li = M_ZERO

    vdw_unit  = io_open('pcm/vdw_cavity.dat', action='write')
    grid_unit = io_open('pcm/grid.dat', action='write')

     do ia=1,pcm%n_tesserae
        write(vdw_unit,*) pcm%tess(ia)%x/P_Ang
     enddo

     do ia=1, grid%mesh%np
        write(grid_unit,*) grid%mesh%x(ia,:)/P_Ang
     enddo

    call io_close(vdw_unit)
    call io_close(grid_unit)
    
    nearest_idx_unit    = io_open('pcm/nearest_index.dat', action='write')
    idx_from_coord_unit = io_open('pcm/index_from_coords.dat', action='write')

    !> Creating the list of the nearest grid points to each tessera
    !! to be used to interpolate the Hartree potential at the representative points
     do ia = 1, pcm%n_tesserae
        call nearest_cube_vertices( pcm%tess(ia)%x, grid%mesh, pcm%ind_vh(ia,:), pcm%arg_li(ia,:), ia, pcm%n_vertices )
     enddo

    call io_close(nearest_idx_unit)
    call io_close(idx_from_coord_unit)

    SAFE_ALLOCATE( pcm%matrix(1:pcm%n_tesserae, 1:pcm%n_tesserae) )
    pcm%matrix = M_ZERO
    call pcm_matrix(pcm%epsilon_0, pcm%tess, pcm%n_tesserae, pcm%matrix) 
    message(1) = "Info: PCM response matrix has been evaluated"
    call messages_info(1)

    pcmmat_unit = io_open('pcm/pcm_matrix.out', action='write')

     do jtess=1, pcm%n_tesserae
      do itess=1, pcm%n_tesserae
         write(pcmmat_unit,*) pcm%matrix(itess,jtess)
      enddo
     enddo
	 
    call io_close(pcmmat_unit)

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
    return
  end subroutine pcm_init
!=======================================================================================
  !> Calculates the Hartree potential at the tessera representative points by taking the 
  !! average of 'v_hartree' over the closest 8 (cube vertices) grid points. 
  subroutine v_electrons_cav(v_e_cav, v_hartree, pcm)
    type(pcm_t), intent(in)  :: pcm
    FLOAT, intent(in)        :: v_hartree(:) !< (1:mesh%np)
    FLOAT, intent(out)       :: v_e_cav(:)   !< (1:n_tess)

    integer :: ia
    integer :: ib

    PUSH_SUB(v_electrons_cav)    

     v_e_cav = M_ZERO

      do ia=1, pcm%n_tesserae
        do ib=1, pcm%n_vertices
           v_e_cav(ia) = v_e_cav(ia) + v_hartree( pcm%ind_vh(ia,ib) )
        enddo

        v_e_cav(ia) = -v_e_cav(ia)/pcm%n_vertices !Notice the explicit minus sign. 

      enddo

    POP_SUB(v_electrons_cav)
  end subroutine v_electrons_cav      
!==================================================================
  !> Calculates the Hartree potential at the tessera representative points by doing 
  !! a 3D linear interpolation. 
  subroutine v_electrons_cav_li(v_e_cav, v_hartree, pcm)
    type(pcm_t), intent(in)  :: pcm
    FLOAT, intent(in)        :: v_hartree(:) !< (1:mesh%np)
    FLOAT, intent(out)       :: v_e_cav(:)   !< (1:n_tess)

    integer :: ia
    integer :: ib

    FLOAT :: C_00
    FLOAT :: C_10
    FLOAT :: C_01
    FLOAT :: C_11
    FLOAT :: C_0
    FLOAT :: C_1

    PUSH_SUB(v_electrons_cav_li)    

     v_e_cav = M_ZERO

     do ia=1, pcm%n_tesserae

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

      enddo
   
    POP_SUB(v_electrons_cav_li)
  end subroutine v_electrons_cav_li      
!==================================================================
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
        diff = geo%atom(ia)%x - tess(ik)%x 
        
        dist = dot_product( diff, diff )
        dist = sqrt(dist)

        spci => geo%atom(ia)%spec
        z_ia = species_zval(spci)

        v_n_cav(ik) = v_n_cav(ik) + z_ia/dist       
     enddo
    enddo

    v_n_cav = -v_n_cav

    POP_SUB(v_nuclei_cav)
  end subroutine v_nuclei_cav
!==================================================================
  !> Calculates the classical electrostatic interaction energy between the polarization
  !! charges and the ions: int_n_pcm = 1/2 \sum_{I=1}^{natoms} \sum_{j=1}^{n_tess} &
  !! Z_val*(q_j^e+q_j^n) / |s_{j} - R_I|
  FLOAT function pcm_classical_energy(geo, q_e, q_n, tess, n_tess)
    type(geometry_t), intent(in) :: geo
    FLOAT                        :: q_e(:)  !< (1:n_tess)
    FLOAT                        :: q_n(:)  !< (1:n_tess)
    type(tessera_t), intent(in)  :: tess(:) !< (1:n_tess)
    integer, intent(in)          :: n_tess

    FLOAT   :: diff(1:MAX_DIM)
    FLOAT   :: dist
    FLOAT   :: z_ia
    integer :: ik
    integer :: ia

    type(species_t), pointer :: spci 

    PUSH_SUB(pcm_classical_energy)
     
    pcm_classical_energy = M_ZERO

    do ik = 1, n_tess
     do ia = 1, geo%natoms
        diff = geo%atom(ia)%x - tess(ik)%x 
        
        dist = dot_product( diff, diff )
        dist = sqrt(dist)

        spci => geo%atom(ia)%spec
        z_ia = species_zval(spci)

        pcm_classical_energy = pcm_classical_energy + z_ia*( q_e(ik) + q_n(ik) ) / dist       
     enddo
    enddo

    pcm_classical_energy = M_HALF*pcm_classical_energy

    POP_SUB(pcm_classical_energy)
  end function pcm_classical_energy
!==================================================================
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
   integer :: tess_unit
   integer :: vertices_unit

   character(LEN=20) :: cc

   PUSH_SUB(nearest_cube_vertices)    

   vert_idx(1) = mesh_nearest_point(mesh, point, dmin, rankmin)
       
   coord_0 = mesh%x(vert_idx(1), :)
   point_0 = NINT(coord_0/mesh%spacing)
   
   write(nearest_idx_unit,*) vert_idx(1)
   write(idx_from_coord_unit,*) index_from_coords(mesh%idx, mesh%sb%dim, point_0)

   sign_x = INT( sign( CNST(1.0), point(1) - coord_0(1) ) )
   sign_y = INT( sign( CNST(1.0), point(2) - coord_0(2) ) )
   sign_z = INT( sign( CNST(1.0), point(3) - coord_0(3) ) )

   weight_li = ABS(point - coord_0)/mesh%spacing 

  !FRONT CUBE PLANE
   point_f = point_0

  !POINT P2
   point_f(2) = point_f(2) + sign_y
   vert_idx(2) = index_from_coords(mesh%idx, mesh%sb%dim, point_f)

  !POINT P3
   point_f(3) = point_f(3) + sign_z
   vert_idx(3) = index_from_coords(mesh%idx, mesh%sb%dim, point_f)

  !POINT P4
   point_f(2) = point_f(2) - sign_y
   vert_idx(4) = index_from_coords(mesh%idx, mesh%sb%dim, point_f)

  !REAR CUBE PLANE
   point_f = point_0 

  !POINT P5
   point_f(1) = point_f(1) + sign_x
   vert_idx(5) = index_from_coords(mesh%idx, mesh%sb%dim, point_f)

  !POINT P6
   point_f(2) = point_f(2) + sign_y
   vert_idx(6) = index_from_coords(mesh%idx, mesh%sb%dim, point_f)

  !POINT P7
   point_f(3) = point_f(3) + sign_z
   vert_idx(7) = index_from_coords(mesh%idx, mesh%sb%dim, point_f)

  !POINT P8
   point_f(2) = point_f(2) - sign_y
   vert_idx(8) = index_from_coords(mesh%idx, mesh%sb%dim, point_f)

   write(cc,'(I3)') ia

   tess_unit     = io_open('pcm/'//trim(cc)//'_tess.dat', action='write')
   vertices_unit = io_open('pcm/'//trim(cc)//'_vertices.dat', action='write')

    write(tess_unit,*) point/P_Ang
   
    do ii=1, n_vertices
       write(vertices_unit,*) mesh%x(vert_idx(ii), :)/P_Ang 
    enddo

   call io_close(tess_unit)
   call io_close(vertices_unit)

   POP_SUB(nearest_cube_vertices)    
  
  end subroutine nearest_cube_vertices
!==============================================================================================
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
       enddo 
       q_pcm_tot = q_pcm_tot + q_pcm(ia)
    enddo
    
    POP_SUB(pcm_charges)
  end subroutine pcm_charges
!==================================================
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
    FLOAT            :: coord_tess(1:MAX_DIM) 
    FLOAT            :: rr
    FLOAT            :: arg
    FLOAT            :: term
    integer 	     :: ip
    integer          :: ia

    PUSH_SUB(pcm_pot_rs)

    v_pcm = M_ZERO

    do ia = 1, n_tess
        coord_tess = tess(ia)%x
       do ip = 1, mesh%np !running serially np=np_global

          call mesh_r(mesh, ip, rr, origin=coord_tess)

          arg = rr/sqrt( tess(ia)%area*width_factor )    

          term = ( 1 + p_1*arg + p_2*arg**2 )/( 1 + q_1*arg + q_2*arg**2 + p_2*arg**3 )
          v_pcm(ip) = v_pcm(ip) + q_pcm(ia)*term/sqrt( tess(ia)%area*width_factor )

       enddo 
    enddo

    v_pcm = M_TWO*v_pcm/sqrt(M_Pi) 

    POP_SUB(pcm_pot_rs)
  end subroutine pcm_pot_rs
!=====================================================================================
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
    enddo

    SAFE_DEALLOCATE_A(d_mat_act) 
     
    SAFE_ALLOCATE( iwork(1:n_tess) )

    !> Solving for X = S_I^-1*(2*Pi - D_I) 
    call dgesv(n_tess, n_tess, s_mat_act, n_tess, iwork, pcm_mat, n_tess, info)        

    SAFE_DEALLOCATE_A(iwork)

    SAFE_DEALLOCATE_A(s_mat_act)

    !> Computing -S_E*S_I^-1*(2*Pi - D_I)
    pcm_mat = -matmul( Sigma, pcm_mat ) 
     
    do i=1, n_tess
       pcm_mat(i,i) = pcm_mat(i,i) + M_TWO*M_Pi
    enddo

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
       
    POP_SUB(pcm_matrix)

  end subroutine pcm_matrix
!==============================================================
  subroutine s_i_matrix(n_tess, tess)

    integer, intent(in)         :: n_tess
    type(tessera_t), intent(in) :: tess(:)
    
    integer :: i
    integer :: j

    s_mat_act = M_ZERO 

    do i = 1, n_tess
     do j = i, n_tess

        s_mat_act(i,j) = s_mat_elem_I( tess(i), tess(j) )
        if (i /= j) s_mat_act(j,i) = s_mat_act(i,j) !symmetric matrix

     enddo
    enddo

  end subroutine s_i_matrix
!==============================================================
  subroutine d_i_matrix(n_tess, tess)

    integer, intent(in)         :: n_tess
    type(tessera_t), intent(in) :: tess(:)

    integer :: i
    integer :: j

    d_mat_act = M_ZERO 

    do i = 1, n_tess
     do j = 1, n_tess !< non-symmetric matrix

        d_mat_act(i,j) = d_mat_elem_I( tess(i), tess(j) )

     enddo
    enddo

  end subroutine d_i_matrix
!===========================================
  !> electrostatic Green function in vacuo:
  !! G_I(r,r^\prime) = 1 / | r - r^\prime |
  FLOAT function s_mat_elem_I( tessi, tessj )
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

    diff = tessi%x - tessj%x

    dist = dot_product( diff, diff )
    dist = sqrt(dist)

    if  ( dist == M_ZERO ) then

        s_diag = M_SD_DIAG*sqrt( M_FOUR*M_Pi/tessi%area )
        s_mat_elem_I = s_diag

    else

        if ( dist > M_DIST_MIN ) s_off_diag = M_ONE/dist 
        s_mat_elem_I = s_off_diag

    endif

    return
    end function s_mat_elem_I
!=================================================================
  !> Gradient of the Green function in vacuo GRAD[G_I(r,r^\prime)]
  FLOAT function d_mat_elem_I( tessi, tessj )
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

    diff = tessi%x - tessj%x

    dist = dot_product( diff, diff )
    dist = sqrt(dist)

    IF (dist == M_ZERO) THEN
        !> Diagonal matrix elements  
        d_diag = M_SD_DIAG*sqrt(M_FOUR*M_Pi*tessi%area)
        d_diag = -d_diag/(M_TWO*tessi%r_sphere)
        d_mat_elem_I = d_diag
   
    else
        !> off-diagonal matrix elements  
        if (dist > M_DIST_MIN) then
            d_off_diag = dot_product( diff, tessj%n(:) )	
            d_off_diag = d_off_diag*tessj%area/dist**3
        endif

        d_mat_elem_I = d_off_diag

    endif
   
    return
    end function d_mat_elem_I
!==============================================================
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

     call io_close(pcm%info_unit)

    POP_SUB(pcm_end)
  end subroutine pcm_end

end module pcm_m
