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
  use io_m
  use index_m
  use messages_m
  use mesh_m 
  use parser_m
  use profiling_m

  implicit none

  private

  public :: cts_act,    &
            ind_vh,     &
            nts_act,    &
            n_vertices, &
            pcm_init,   &
            pcm_mat,    &
            pcm_charges,&
            run_pcm

 !> The cavity hosting the solute molecule is built from a set of 
 !! interlocking spheres with optimized radii centered at the nuclear positions.  
  type sphere_t  
    FLOAT :: x ! 
    FLOAT :: y !< center of the sphere  
    FLOAT :: z !
    FLOAT :: r !< radius of the sphere (different for each species)
  end type

 !> The resulting cavity is discretized by a set of tesserae.  
  type, public :: tess_pcm_t
    FLOAT :: x    !
    FLOAT :: y    !< representative point of the tessera 
    FLOAT :: z    !
    FLOAT :: area !< area of the tessera
    FLOAT :: n(3) !< unitary outgoing vector normal to the tessera surface 
    FLOAT :: rsfe !< radius of the sphere to which the tessera belongs
  end type tess_pcm_t

  integer :: nts_act                 !< total number of tesserae 
  FLOAT, allocatable :: pcm_mat(:,:) !< PCM response matrix 

  type(sphere_t), allocatable   :: sfe_act(:) !< set of spheres used to build the molecular cavity    
  type(tess_pcm_t), allocatable :: cts_act(:) !< tesselation of the Van der Waals surface

  FLOAT, allocatable :: s_mat_act(:,:) !< S_I matrix 
  FLOAT, allocatable :: d_mat_act(:,:) !< D_I matrix
  FLOAT, allocatable :: Sigma(:,:)     !< S_E matrix
  FLOAT, allocatable :: Delta(:,:)     !< D_E matrix in JCP 139, 024105 (2013).

  integer, allocatable :: ind_vh(:,:)
  integer, parameter   :: n_vertices = 8
  integer :: pcminfo_unit
  logical :: run_pcm 

  ! End of variable declaration.
  ! ----------------------------

contains

  !--------------------------------------------------------------------------------------------
  !> Initializes the PCM calculation: generate the molecular cavity and the PCM response matrix
  subroutine pcm_init(geo, mesh)
    type(geometry_t), intent(in) :: geo
    type(mesh_t), intent(in)     :: mesh

!    type(tess_pcm_t) :: dum2(1)
    integer :: ia
    integer :: itess
    integer :: jtess
    integer :: nesf_act
    integer :: cav_unit
    integer :: cav_unit_test
    integer :: pcmmat_unit
    integer :: iunit
    integer, parameter   :: mxts = 10000

    FLOAT   :: epsilon_static
    FLOAT   :: rcav_C
    FLOAT   :: rcav_O
    FLOAT   :: rcav_N
    FLOAT   :: rcav_S
    FLOAT   :: rcav_F
    FLOAT   :: rep_point(1:3)

!    logical :: run_pcm 

    character(len=80) :: str
    integer :: tess_counter

    PUSH_SUB(pcm_init)

    rcav_C = CNST(2.4)*P_Ang    ! 
    rcav_O = CNST(1.8)*P_Ang    !    
    rcav_N = CNST(1.9)*P_Ang    ! Agnstrom -> Bohr 
    rcav_S = CNST(2.0175)*P_Ang !
    rcav_F = CNST(1.682)*P_Ang  ! 

    !%Variable SolventDielectricConstant
    !%Type float
    !%Default 1.0 (gas phase)
    !%Section Hamiltonian::PCM
    !%Description
    !% Static dielectric constant of the solvent (\epsilon_0).
    !%End
    call parse_float(datasets_check('SolventDielectricConstant'), CNST(1.0), epsilon_static)

    !%Variable CavityGeometry
    !%Type string
    !%Section Hamiltonian::PCM
    !%Description
    !% Name of the file containing the geometry of the Van der Waals surface that defines the cavity hosting
    !% the solute molecule in PCM calculations. Tesserae representative points must be in atomic units!.
    !%End

    call parse_string(datasets_check('CavityGeometry'), '', str)

    iunit = io_open(trim(str), status='old', action='read')
    
    read(iunit,*) nts_act

    if (nts_act.gt.mxts) then
        write(message(1),'(a,I5)') "Info: WARNING: total number of tesserae > 10 000 ", nts_act
        call messages_info(1)     
    endif

    SAFE_ALLOCATE(cts_act(1:nts_act))

    do ia=1, nts_act
       read(iunit,*) cts_act(ia)%x
    enddo
	
    do ia=1, nts_act
       read(iunit,*) cts_act(ia)%y
    enddo
	
    do ia=1, nts_act
       read(iunit,*) cts_act(ia)%z
    enddo

    do ia=1, nts_act
       read(iunit,*) cts_act(ia)%area
    enddo

    do ia=1, nts_act
       read(iunit,*) cts_act(ia)%rsfe
    enddo

    do ia=1, nts_act
       read(iunit,*) cts_act(ia)%n
    enddo
    
    call io_close(iunit)

    message(1) = "Info: van der Waals surface has been read from " // trim(str)
    call messages_info(1)

    ! Creating the list of the nearest grid points to each tessera
    ! to be used to "interpolate" the Hartree potential at the tesserae

    SAFE_ALLOCATE( ind_vh(1:nts_act, 1:n_vertices) )
    ind_vh = INT(M_ZERO)

    ! only for testing 
!    OPEN(500, FILE='pcm_cavity.xyz')
!    OPEN(501, FILE='grid.xyz')

!    OPEN(502, FILE='micael_test_1.xyz')
!    OPEN(503, FILE='micael_test_2.xyz')

!    DO ia=1,nts_act
!       WRITE(500,*) cts_act(ia)%x, cts_act(ia)%y, cts_act(ia)%z
!    ENDDO

!    DO ia=1, mesh%np
!   
!       WRITE(501,*) mesh%x(ia,:)
!       WRITE(502,*) mesh%x(ia,:)/P_Ang, NINT( mesh%x(ia,:)/mesh%spacing )
!       WRITE(503,*) mesh%x(ia,:)/P_Ang, mesh%idx%lxyz(ia,:)
    
!    ENDDO
    
!    CLOSE(500)
!    CLOSE(501)

!    CLOSE(502)
!    CLOSE(503)
!    only for testing, will be cleaned soon 
     
!   OPEN(600, FILE='nearest_index.dat')
!   OPEN(601, FILE='index_from_coords.dat')

    do ia = 1, nts_act

       rep_point(1) = cts_act(ia)%x
       rep_point(2) = cts_act(ia)%y
       rep_point(3) = cts_act(ia)%z

      !Creating the list of the nearest grid points to each tessera
      !to "interpolate" the Hartree potential at the tesserae.     
       call nearest_cube_vertices( rep_point, mesh, ind_vh(ia,:), ia )

    enddo

!   CLOSE(600)
!   CLOSE(601)

    nesf_act = 0
    do ia = 1, geo%natoms
       if (geo%atom(ia)%label == 'H') cycle
       nesf_act = nesf_act + 1 !counting the number of species different from Hydrogen
    enddo

    SAFE_ALLOCATE(sfe_act(1:nesf_act))
    
    nesf_act = 0     
    do ia = 1, geo%natoms
       
       if (geo%atom(ia)%label == 'H') cycle
       nesf_act = nesf_act + 1
      
       ! These coordinates are already in atomic units (Bohr)
       sfe_act(nesf_act)%x = geo%atom(ia)%x(1)
       sfe_act(nesf_act)%y = geo%atom(ia)%x(2)
       sfe_act(nesf_act)%z = geo%atom(ia)%x(3)

       if (geo%atom(ia)%label == 'C') sfe_act(nesf_act)%r = rcav_C
       if (geo%atom(ia)%label == 'O') sfe_act(nesf_act)%r = rcav_O
       if (geo%atom(ia)%label == 'N') sfe_act(nesf_act)%r = rcav_N
       if (geo%atom(ia)%label == 'S') sfe_act(nesf_act)%r = rcav_S
       if (geo%atom(ia)%label == 'F') sfe_act(nesf_act)%r = rcav_F                

    enddo

    call io_mkdir('pcm')

    pcminfo_unit = io_open('pcm/pcm_info.out', action='write')

    write(pcminfo_unit,'(2X,A33)') 'Configuration: Molecule + Solvent'
    write(pcminfo_unit,'(2X,A33)') '---------------------------------'
    write(pcminfo_unit,'(2X,A19,F6.3)') 'Epsilon(Solvent) = ', epsilon_static
    write(pcminfo_unit,'(2X)') 
    write(pcminfo_unit,'(2X,A33,I4)') 'Number of interlocking spheres = ', nesf_act
    write(pcminfo_unit,'(2X)')  

    write(pcminfo_unit,'(2X,A6,3X,A7,8X,A26,20X,A10)') 'SPHERE', 'ELEMENT', 'CENTER  (X,Y,Z) (A)', 'RADIUS (A)'
    write(pcminfo_unit,'(2X,A6,3X,A7,4X,A43,7X,A10)') '------', '-------', &
                              '-------------------------------------------', '----------'  

    nesf_act = 0
    do ia = 1, geo%natoms
       if (geo%atom(ia)%label == 'H') cycle
       nesf_act = nesf_act + 1       

       write(pcminfo_unit,'(2X,I3,9X,A2,3X,F14.8,2X,F14.8,2X,F14.8,3X,F14.8)') nesf_act,                 &
                         						       geo%atom(ia)%label,       &
                                                                               geo%atom(ia)%x*P_a_B,     &
                                                                               sfe_act(nesf_act)%r*P_a_B
    enddo
    write(pcminfo_unit,'(2X)')  

    cav_unit_test = io_open('pcm/cavity_mol.xyz', action='write')

    write (cav_unit_test,'(2X,I4)') nts_act+geo%natoms
    write (cav_unit_test,'(2X)')

    do itess=1, nts_act
       write(cav_unit_test,'(2X,A2,3X,4f15.8,3X,4f15.8,3X,4f15.8)') 'H', cts_act(itess)%x*P_a_B, &
                                                                         cts_act(itess)%y*P_a_B, &
                                                                         cts_act(itess)%z*P_a_B
    enddo

    do ia=1, geo%natoms
       write(cav_unit_test,'(2X,A2,3X,4f15.8,3X,4f15.8,3X,4f15.8)') geo%atom(ia)%label,      &
                                                                    geo%atom(ia)%x(1)*P_a_B, &
                                                                    geo%atom(ia)%x(2)*P_a_B, &
                                                                    geo%atom(ia)%x(3)*P_a_B
    enddo

    call io_close(cav_unit_test)

    call pcm_matrix(epsilon_static) ! Calculates the PCM response matrix
    message(1) = "Info: PCM response matrix has been evaluated"
    call messages_info(1)

    pcmmat_unit = io_open('pcm/pcm_matrix.out', action='write')

     do jtess=1, nts_act
      do itess=1, nts_act
         write(pcmmat_unit,*) pcm_mat(itess,jtess)
      enddo
     enddo
	 
    call io_close(pcmmat_unit)
    call io_close(pcminfo_unit)

    POP_SUB(pcm_init)
    return
  end subroutine pcm_init
!==============================================================
  subroutine nearest_cube_vertices(point, mesh, vert_idx, ia)
   FLOAT, intent(in)        :: point(1:3)
   type(mesh_t), intent(in) :: mesh
   integer, intent(out)     :: vert_idx(:)
   integer, intent(in)      :: ia

   FLOAT   :: dmin
   integer :: rankmin

   FLOAT   :: coord_0(1:3)
   integer :: sign_x
   integer :: sign_y
   integer :: sign_z
   integer :: point_0(1:3)
   integer :: point_f(1:3)

!   INTEGER :: ii
!   CHARACTER(LEN=20) :: cc

   PUSH_SUB(nearest_cube_vertices)    

   vert_idx(1) = mesh_nearest_point(mesh, point, dmin, rankmin)
       
   coord_0 = mesh%x(vert_idx(1), :)
   point_0 = NINT(coord_0/mesh%spacing)
   
!   WRITE(600,*) vert_idx(1)
!   WRITE(601,*) index_from_coords(mesh%idx, mesh%sb%dim, point_0)

   sign_x = INT( sign( CNST(1.0), point(1) - coord_0(1) ) )
   sign_y = INT( sign( CNST(1.0), point(2) - coord_0(2) ) )
   sign_z = INT( sign( CNST(1.0), point(3) - coord_0(3) ) ) 

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

!   WRITE(cc,'(I3)') ia

!   OPEN(400, file=trim(cc)//'_tess.dat')
!   OPEN(401, file=trim(cc)//'_points.dat')

!      WRITE(400,*) point/P_Ang
   
!   DO ii=1, n_vertices
!      WRITE(401,*) mesh%x(vert_idx(ii), :)/P_Ang 
!   ENDDO

!   CLOSE(400)
!   CLOSE(401)

   POP_SUB(nearest_cube_vertices)    
  
  end subroutine nearest_cube_vertices
!==============================================================

!==============================================================
  subroutine pcm_charges(q_pcm, q_pcm_tot, v_cav)
!   Calculates the polarization charges at each tessera by using the response matrix 'pcm_mat',
!   defined in 'pcm/pcm.F90', provided the value of the molecule's electrostatic potential at 
!   the tessera: q_pcm(ia) = \sum_{ib}^{nts_act} pcm_mat(ia,ib)*v_cav(ib).

    FLOAT, intent(out) :: q_pcm(1:nts_act)
    FLOAT, intent(in)  :: v_cav(1:nts_act)
    FLOAT, intent(out) :: q_pcm_tot

    integer :: ia
    integer :: ib

    PUSH_SUB(pcm_charges)

    q_pcm     = M_ZERO
    q_pcm_tot = M_ZERO

    do ia = 1, nts_act
       do ib = 1, nts_act
          q_pcm(ia) = q_pcm(ia) + pcm_mat(ia,ib)*v_cav(ib) !< transpose matrix might speed up
       enddo 
       q_pcm_tot = q_pcm_tot + q_pcm(ia)
    enddo
    
    POP_SUB(pcm_charges)
  end subroutine pcm_charges
!==============================================================

  subroutine pcm_matrix(eps)
    FLOAT, intent(in) :: eps

    integer :: i
    integer :: info
    integer, allocatable :: iwork(:)

    FLOAT, allocatable :: mat_tmp(:,:)

    PUSH_SUB(pcm_matrix)

!   Conforming the S_I matrix
    SAFE_ALLOCATE( s_mat_act(1:nts_act, 1:nts_act) )
    call s_i_matrix

!   Defining the matrix S_E=S_I/eps
    SAFE_ALLOCATE( Sigma(1:nts_act, 1:nts_act) )
    Sigma = s_mat_act/eps

!   Conforming the D_I matrix
    SAFE_ALLOCATE( d_mat_act(1:nts_act, 1:nts_act) )
    call d_i_matrix

!   Defining the matrix D_E=D_I 
    SAFE_ALLOCATE( Delta(1:nts_act, 1:nts_act) )
    Delta = d_mat_act

!   Start conforming the PCM matrix
    SAFE_ALLOCATE( pcm_mat(1:nts_act, 1:nts_act) )
    pcm_mat = M_ZERO

    pcm_mat = -d_mat_act

    do i=1, nts_act
       pcm_mat(i,i) = pcm_mat(i,i) + M_TWO*M_Pi  
    enddo

    SAFE_DEALLOCATE_A(d_mat_act) 
     
    SAFE_ALLOCATE( iwork(1:nts_act) )

!   Solving for X = S_I^-1*(2*Pi - D_I) 

    call dgesv(nts_act,nts_act,s_mat_act,nts_act,iwork,pcm_mat,nts_act,info)    

    SAFE_DEALLOCATE_A(iwork)

    SAFE_DEALLOCATE_A(s_mat_act)

    pcm_mat = -matmul( Sigma, pcm_mat ) ! Computing -S_E*S_I^-1*(2*Pi - D_I)
     
    do i=1, nts_act
       pcm_mat(i,i) = pcm_mat(i,i) + M_TWO*M_Pi
    enddo

    pcm_mat = pcm_mat - Delta
    
    SAFE_ALLOCATE( mat_tmp(1:nts_act,1:nts_act) )
    mat_tmp = M_ZERO

    SAFE_ALLOCATE( d_mat_act(1:nts_act,1:nts_act) )  
    call d_i_matrix

    mat_tmp = transpose(d_mat_act)        

    mat_tmp = matmul(Sigma, mat_tmp)

    mat_tmp = mat_tmp + M_TWO*M_Pi*Sigma

    SAFE_DEALLOCATE_A(d_mat_act)

    SAFE_ALLOCATE( s_mat_act(1:nts_act,1:nts_act) )
    call s_i_matrix

    mat_tmp = mat_tmp + M_TWO*M_Pi*s_mat_act - matmul(Delta, s_mat_act)

    SAFE_DEALLOCATE_A(s_mat_act)
    SAFE_DEALLOCATE_A(Sigma)
    SAFE_DEALLOCATE_A(Delta)

!   Solving for [(2*pi - D_E)*S_I + S_E*(2*Pi + D_I*)]*X = [(2*Pi - D_E) - S_E*S_I^-1*(2*Pi - D_I)]    

    SAFE_ALLOCATE( iwork(1:nts_act) )

    call dgesv(nts_act,nts_act,mat_tmp,nts_act,iwork,pcm_mat,nts_act,info)		  

    SAFE_DEALLOCATE_A(iwork)
    SAFE_DEALLOCATE_A(mat_tmp)

    pcm_mat = -pcm_mat
       
    POP_SUB(pcm_matrix)

  end subroutine pcm_matrix
!==============================================================
  subroutine s_i_matrix

    integer :: i
    integer :: j
    type(tess_pcm_t) :: tess_i
    type(tess_pcm_t) :: tess_j

    s_mat_act = M_ZERO 

    do i = 1, nts_act
       tess_i = cts_act(i)

     do j = i, nts_act
        tess_j = cts_act(j)

        s_mat_act(i,j) = s_mat_elem_I(tess_i,tess_j)

        if (i /= j) s_mat_act(j,i) = s_mat_act(i,j) !symmetric matrix

     enddo

    enddo

  end subroutine s_i_matrix
!==============================================================
  subroutine d_i_matrix

    integer :: i
    integer :: j
    TYPE(tess_pcm_t) :: tess_i
    TYPE(tess_pcm_t) :: tess_j

    d_mat_act = M_ZERO 

    do i = 1, nts_act
       tess_i = cts_act(i)

     do j = 1, nts_act !non-symmetric matrix
        tess_j = cts_act(j)

        d_mat_act(i,j) = d_mat_elem_I(tess_i,tess_j)

     enddo

    enddo

  end subroutine d_i_matrix
!==============================================================

  FLOAT function s_mat_elem_I( tessi, tessj )
    type(tess_pcm_t), intent(in) :: tessi
    type(tess_pcm_t), intent(in) :: tessj

    FLOAT, parameter :: M_SD_DIAG    = CNST(1.0694)
    FLOAT, parameter :: M_DIST_MIN   = CNST(0.1)

    FLOAT :: diff(3)
    FLOAT :: dist
    FLOAT :: s_diag
    FLOAT :: s_off_diag

!   GREEN FUNCTION IN VACUO G_I(r,r^\prime) = 1/|r-r^\prime|

    s_diag = M_ZERO
    s_off_diag = M_ZERO

    diff(1) = tessi%x - tessj%x    
    diff(2) = tessi%y - tessj%y
    diff(3) = tessi%z - tessj%z

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
!==============================================================

  FLOAT function d_mat_elem_I( tessi, tessj )
    type(tess_pcm_t), intent(in) :: tessi
    type(tess_pcm_t), intent(in) :: tessj


    FLOAT, parameter :: M_SD_DIAG    = CNST(1.0694)
    FLOAT, parameter :: M_DIST_MIN   = CNST(0.04)

    FLOAT :: diff(3)
    FLOAT :: dist
    FLOAT :: d_diag
    FLOAT :: d_off_diag

!!  GRADIENT OF THE GREEN FUNCTION IN VACUO GRAD[G_I(r,r^\prime)]

    d_diag = M_ZERO
    d_off_diag = M_ZERO

    diff(1) = tessi%x - tessj%x    
    diff(2) = tessi%y - tessj%y
    diff(3) = tessi%z - tessj%z

    dist = dot_product( diff, diff )
    dist = sqrt(dist)

    IF (dist == M_ZERO) THEN
!       Diagonal matrix elements  
        d_diag = M_SD_DIAG*sqrt(M_FOUR*M_Pi*tessi%area) 
        d_diag = -d_diag/(M_TWO*tessi%rsfe)
        d_mat_elem_I = d_diag
   
    else
!       off-diagonal matrix elements  
        if (dist > M_DIST_MIN) then
            d_off_diag = dot_product( diff, tessj%n(:) )	
            d_off_diag = d_off_diag*tessj%area/dist**3
        endif

        d_mat_elem_I = d_off_diag

    endif
   
    return
    end function d_mat_elem_I
!==============================================================

END MODULE pcm_m
