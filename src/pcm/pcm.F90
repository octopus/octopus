#include "global.h"

module pcm_m
  use datasets_m
  use global_m
  use geometry_m
  use io_m
  use messages_m
  use parser_m
  use profiling_m

  implicit none

  private

  public :: cts_act,    &
            nts_act,    &
            pcm_init,   &
            pcm_mat

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

  integer :: pcminfo_unit

  ! End of variable declaration.
  ! ----------------------------

contains

  !--------------------------------------------------------------------------------------------
  !> Initializes the PCM calculation: generate the molecular cavity and the PCM response matrix
  subroutine pcm_init(geo)
    type(geometry_t), intent(in) :: geo

    type(tess_pcm_t) :: dum2(1)
    integer :: ia
    integer :: itess,jtess
    integer :: nesf_act
    integer :: cav_unit, cav_unit_test
    integer :: pcmmat_unit
    integer, parameter :: mxts = 10000

    FLOAT   :: epsilon_static

    FLOAT   :: rcav_C
    FLOAT   :: rcav_O
    FLOAT   :: rcav_N
    FLOAT   :: rcav_S
    FLOAT   :: rcav_F

    logical :: run_pcm 

    PUSH_SUB(pcm_init)

    rcav_C = CNST(2.4)*P_Ang    ! 
    rcav_O = CNST(1.8)*P_Ang    !    
    rcav_N = CNST(1.9)*P_Ang    ! Agnstrom -> Bohr 
    rcav_S = CNST(2.0175)*P_Ang !
    rcav_F = CNST(1.682)*P_Ang  ! 

    ! -I need the RADII of the species within a block

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
    call parse_logical(datasets_check('Solvation'), .false., run_pcm)
    if (.not.(run_pcm)) then
        POP_SUB(pcm_init)
        return
    else
      call messages_experimental("polarizable continuum model")
    endif

    !%Variable SolventDielectricConstant
    !%Type float
    !%Default 1.0 (gas phase)
    !%Section Hamiltonian::PCM
    !%Description
    !% Static dielectric constant of the solvent (\epsilon_0).
    !%End
    call parse_float(datasets_check('SolventDielectricConstant'), CNST(1.0), epsilon_static)

    !%Variable PCMRadii
    !%Type float
    !%Section Hamiltonian::PCM
    !%Description
    !%
    !%
    !%End
!    if(parse_block(datasets_check('PCMRadii'), blk) == 0) then
!      call check_duplicated(done)

!      gf%n = parse_block_n(blk)

!      message(1) = "Reading " // trim(what) // " from " // trim(what) // " block"
!      call messages_info(1)

!      temp_SAFE_ALLOCATE(gf%atom(1:gf%n))

!      do ia = 1, gf%n
!        ncol = parse_block_cols(blk, ia - 1)
!        if((ncol .lt. space%dim + 1) .or. (ncol .gt. space%dim + 2)) then
!          write(message(1), '(3a,i2)') 'Error in block ', what, ' line #', ia
!          call messages_fatal(1)
!        end if
!        call parse_block_string (blk, ia - 1, 0, gf%atom(ia)%label)
!        call parse_block_float  (blk, ia - 1, 1, gf%atom(ia)%x(jdir))
!      end do

!      call parse_block_end(blk)
!    end if

    nesf_act = 0
    do ia = 1, geo%natoms
       if (geo%atom(ia)%label == 'H') cycle
       nesf_act = nesf_act + 1 !counting the number of species different from Hydrogen
    enddo

    SAFE_ALLOCATE(sfe_act(nesf_act))
    
    nesf_act = 0     
    do ia = 1, geo%natoms
       
       if (geo%atom(ia)%label == 'H') cycle
       nesf_act = nesf_act + 1
      
       ! This coordinates are already in atomic units (Bohr)
       sfe_act(nesf_act)%x = geo%atom(ia)%x(1)
       sfe_act(nesf_act)%y = geo%atom(ia)%x(2)
       sfe_act(nesf_act)%z = geo%atom(ia)%x(3)

       if (geo%atom(ia)%label == 'C') sfe_act(nesf_act)%r = rcav_C
       if (geo%atom(ia)%label == 'O') sfe_act(nesf_act)%r = rcav_O
       if (geo%atom(ia)%label == 'N') sfe_act(nesf_act)%r = rcav_N
       if (geo%atom(ia)%label == 'S') sfe_act(nesf_act)%r = rcav_S
       if (geo%atom(ia)%label == 'F') sfe_act(nesf_act)%r = rcav_F                

!       do ib = 1, n_species_solv
!          if (species_solv(ib)%symbol /= geo%atom(ia)%label) cycle
!          sfe_act(nesf_act)%r = species_solv(ib)%r_cav 
!       enddo

    enddo

    return

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

    cav_unit = io_open('pcm/cavity.xyz', action='write')
    cav_unit_test = io_open('pcm/cavity_mol.xyz', action='write')

    write (cav_unit,'(2X,I4)') nts_act
    write (cav_unit_test,'(2X,I4)') nts_act+geo%natoms

    write (cav_unit,'(2X)')
    write (cav_unit_test,'(2X)')

    do itess=1, nts_act

       write(cav_unit,'(2X,A2,3X,4f15.8,3X,4f15.8,3X,4f15.8)') 'H', cts_act(itess)%x, &
                                                                    cts_act(itess)%y, &
                                                                    cts_act(itess)%z

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

    call io_close(cav_unit)
    call io_close(cav_unit_test)

    call pcm_matrix(epsilon_static)

    pcmmat_unit = io_open('pcm/pcm_matrix.out', action='write')

     do jtess=1, nts_act
      do itess=1, nts_act
         write(pcmmat_unit,*) pcm_mat(itess,jtess)
      enddo
     enddo
	 
    call io_close(pcmmat_unit)

    if (nts_act.gt.mxts) then
        write(message(1),'(a,I5)') "WARNING: total number of tesserae = ", nts_act
        call messages_info(1)     
    endif

    message(1) = "Molecular cavity has been built"
    message(2) = "PCM response matrix has been evaluated"

    call messages_info(2)

    call io_close(pcminfo_unit)

    POP_SUB(pcm_init)
    return
  end subroutine pcm_init
!==============================================================

  subroutine pcm_matrix(eps)
    FLOAT, intent(in) :: eps

    integer :: i
    integer :: info
    integer, allocatable :: iwork(:)

    FLOAT, allocatable :: mat_tmp(:,:)

    PUSH_SUB(pcm_matrix)

!   Conforming the S_I matrix
    SAFE_ALLOCATE( s_mat_act(nts_act, nts_act) )
    call s_i_matrix

!   Defining the matrix S_E=S_I/eps
    SAFE_ALLOCATE( Sigma(nts_act, nts_act) )
    Sigma = s_mat_act/eps

!   Conforming the D_I matrix
    SAFE_ALLOCATE( d_mat_act(nts_act, nts_act) )
    call d_i_matrix

!   Defining the matrix D_E=D_I 
    SAFE_ALLOCATE( Delta(nts_act, nts_act) )
    Delta = d_mat_act

!   Start conforming the PCM matrix
    SAFE_ALLOCATE( pcm_mat(nts_act, nts_act) )
    pcm_mat = M_ZERO

    pcm_mat = -d_mat_act

    do i=1, nts_act
       pcm_mat(i,i) = pcm_mat(i,i) + M_TWO*M_Pi  
    enddo

    SAFE_DEALLOCATE_A(d_mat_act) 
     
    SAFE_ALLOCATE( iwork(nts_act) )

!   Solving for X = S_I^-1*(2*Pi - D_I) 

    call dgesv(nts_act,nts_act,s_mat_act,nts_act,iwork,pcm_mat,nts_act,info)    

    SAFE_DEALLOCATE_A(iwork)

    SAFE_DEALLOCATE_A(s_mat_act)

    pcm_mat = -matmul( Sigma, pcm_mat ) ! Computing -S_E*S_I^-1*(2*Pi - D_I)
     
    do i=1, nts_act
       pcm_mat(i,i) = pcm_mat(i,i) + M_TWO*M_Pi
    enddo

    pcm_mat = pcm_mat - Delta
    
    SAFE_ALLOCATE( mat_tmp(nts_act,nts_act) )
    mat_tmp = M_ZERO

    SAFE_ALLOCATE( d_mat_act(nts_act,nts_act) )  
    call d_i_matrix

    mat_tmp = transpose(d_mat_act)        

    mat_tmp = matmul(Sigma, mat_tmp)

    mat_tmp = mat_tmp + M_TWO*M_Pi*Sigma

    SAFE_DEALLOCATE_A(d_mat_act)

    SAFE_ALLOCATE( s_mat_act(nts_act,nts_act) )
    call s_i_matrix

    mat_tmp = mat_tmp + M_TWO*M_Pi*s_mat_act - matmul(Delta, s_mat_act)

    SAFE_DEALLOCATE_A(s_mat_act)
    SAFE_DEALLOCATE_A(Sigma)
    SAFE_DEALLOCATE_A(Delta)

!   Solving for [(2*pi - D_E)*S_I + S_E*(2*Pi + D_I*)]*X = [(2*Pi - D_E) - S_E*S_I^-1*(2*Pi - D_I)]    

    SAFE_ALLOCATE( iwork(nts_act) )

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
