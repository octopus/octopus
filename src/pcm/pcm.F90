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
    !% only for ground state calculations.
    !%End
    call parse_logical(datasets_check('Solvation'), .false., run_pcm)
    if (.not.(run_pcm)) then
        POP_SUB(pcm_init)
        return
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

    !subroutine pedra tesselates the cavity surface

    call pedra(0, 1, nesf_act, sfe_act, nts_act, dum2)

    SAFE_ALLOCATE(cts_act(1:nts_act))
	 
    call pedra(1, 1, nesf_act, sfe_act, nts_act, cts_act)

    cav_unit = io_open('pcm/cavity.out', action='write')

      write(cav_unit,'(i5)') nts_act

      do itess=1, nts_act
         write (cav_unit,'(4f15.8)') cts_act(itess)%x
      enddo
	
      do itess=1, nts_act
         write(cav_unit,'(4f15.8)') cts_act(itess)%y
      enddo
	
      do itess=1, nts_act
         write(cav_unit,'(4f15.8)') cts_act(itess)%z
      enddo

      do itess=1, nts_act
         write(cav_unit,'(4f15.8)') cts_act(itess)%area
      enddo

    call io_close(cav_unit)

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
    FLOAT :: eps ! FIXME: specify an intent here.

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

  subroutine pedra(i_count, n_tes, nesf, sfe, nts, cts)
    integer, intent(in) :: i_count, n_tes, nesf
    type(sphere_t), intent(inout) :: sfe(:)
    integer, intent(out) :: nts
    type(tess_pcm_t), intent(out) :: cts(:)

!    TYPE(sphere_t), ALLOCATABLE   :: sfe(:)

    FLOAT :: thev(24),fiv(24),fir,cv(122,3),th,fi,cth,sth
    FLOAT :: XCTST(240),YCTST(240),ZCTST(240),AST(240),nctst(3,240)
    FLOAT :: PTS(3,10),PP(3),PP1(3),CCC(3,10)
      
    integer :: idum(360), JVT1(6,60), isfet(240)
    integer :: i, ii, iii, j, k, nn, nsfe, its, n1,n2,n3,nv,i_tes

    FLOAT :: xen,yen,zen,ren,area,test,test2,rij,dnorm
    FLOAT :: xi,yi,zi,xj,yj,zj
    FLOAT :: vol,stot,prod

    character(LEN=4) :: add_sph='no  ' !this is not used. We built the Van der Waals surface.
!
!    FLOAT, PARAMETER :: MXTS=5000 
!    FLOAT :: DR=0.01  

     
!
!
!   Angoli che individuano i centri e i vertici di un poliedro
!   inscritto in una sfera di raggio unitario centrata nell origine
!
    DATA THEV/0.6523581398D+00,1.107148718D+00,1.382085796D+00, &
              &1.759506858D+00,2.034443936D+00,2.489234514D+00,   &
              &               0.3261790699D+00,0.5535743589D+00,   &
              &0.8559571251D+00,0.8559571251D+00,1.017221968D+00,   &
              &1.229116717D+00,1.229116717D+00,1.433327788D+00,   &
              &1.570796327D+00,1.570796327D+00,1.708264866D+00,   &
              &1.912475937D+00,1.912475937D+00,2.124370686D+00,   &
              &2.285635528D+00,2.285635528D+00,2.588018295D+00,   &
              &2.815413584D+00/
    DATA FIV/               0.6283185307D+00,0.0000000000D+00,   &
             0.6283185307D+00,0.0000000000D+00,0.6283185307D+00,   &
             0.0000000000D+00,0.6283185307D+00,0.0000000000D+00,   &
             0.2520539002D+00,1.004583161D+00,0.6283185307D+00,   &
             0.3293628477D+00,0.9272742138D+00,0.0000000000D+00,   &
             0.3141592654D+00,0.9424777961D+00,0.6283185307D+00,   &
             0.2989556830D+00,0.9576813784D+00,0.0000000000D+00,   &
             0.3762646305D+00,0.8803724309D+00,0.6283188307D+00,   &
             0.0000000000D+00/
    DATA FIR/1.256637061D+00/
!
!   Il vettore IDUM, ripreso nella matrice JVT1, indica quali sono
!   i vertici delle varie tessere (using less than 19 continuations)
!
    DATA (IDUM(III),III=1,280)/ &
       1, 6, 2, 32, 36, 37, 1, 2, 3, 33, 32, 38, 1, 3, 4, 34,&
       33, 39, 1, 4, 5, 35, 34, 40, 1, 5, 6, 36, 35, 41, 7, 2, 6, 51,&
       42, 37, 8, 3, 2, 47, 43, 38, 9, 4, 3, 48, 44, 39, 10, 5, 4,&
       49, 45, 40, 11, 6, 5, 50, 46, 41, 8, 2, 12, 62, 47, 52, 9,&
       3, 13, 63, 48, 53, 10, 4, 14, 64, 49, 54, 11, 5, 15, 65, 50,&
       55, 7, 6, 16, 66, 51, 56, 7, 12, 2, 42, 57, 52, 8, 13, 3,&
       43, 58, 53, 9, 14, 4, 44, 59, 54, 10, 15, 5, 45, 60, 55, 11,&
       16, 6, 46, 61, 56, 8, 12, 18, 68, 62, 77, 9, 13, 19, 69, 63,&
       78, 10, 14, 20, 70, 64, 79, 11, 15, 21, 71, 65, 80, 7, 16,&
       17, 67, 66, 81, 7, 17, 12, 57, 67, 72, 8, 18, 13, 58, 68, 73,&
       9, 19, 14, 59, 69, 74, 10, 20, 15, 60, 70, 75, 11, 21, 16,&
       61, 71, 76, 22, 12, 17, 87, 82, 72, 23, 13, 18, 88, 83, 73,&
       24, 14, 19, 89, 84, 74, 25, 15, 20, 90, 85, 75, 26, 16, 21,&
       91, 86, 76, 22, 18, 12, 82, 92, 77, 23, 19, 13, 83, 93, 78,&
       24, 20, 14, 84, 94, 79, 25, 21, 15, 85, 95, 80, 26, 17, 16,&
       86, 96, 81, 22, 17, 27, 102, 87, 97, 23, 18, 28, 103, 88, 98,&
       24, 19, 29, 104, 89, 99, 25, 20, 30, 105, 90, 100, 26, 21,&
       31, 106, 91, 101, 22, 28, 18, 92, 107, 98, 23, 29, 19, 93/
    DATA (IDUM(III),III=281,360)/&
       108, 99, 24, 30, 20, 94, 109, 100, 25, 31, 21, 95, 110, 101,&
       26, 27, 17, 96, 111, 97, 22, 27, 28, 107, 102, 112, 23, 28,&
       29, 108, 103, 113, 24, 29, 30, 109, 104, 114, 25, 30, 31,&
       110, 105, 115, 26, 31, 27, 111, 106, 116, 122, 28, 27, 117,&
       118, 112, 122, 29, 28, 118, 119, 113, 122, 30, 29, 119, 120,&
       114, 122, 31, 30, 120, 121, 115, 122, 27, 31, 121, 117, 116 /
!
!   It defines the cavity hosting the solute molecule and calculates vertices,
!   representative points and areas of tesserae with the
!   Gauss Bonnet Theorem.
!
!
!SC 19/8/03: Create New spheres, if necessary
! If added on 20/10/03 
! I SHOULD CHECK THIS WITH STEFANO!! 
!      IF (add_sph.eq.'add ') THEN
!       CALL new_sphere(0, nesf_in, sfe_in, nesf, sfe)
!       tmp_ALLOCATE ( sfe(nesf) )
!       CALL new_sphere(1, nesf_in, sfe_in, nesf, sfe)
!      ELSE
!       nesf=nesf_in
!       tmp_ALLOCATE ( sfe(nesf) )
!       sfe=sfe_in
!      ENDIF
!SC 19/8 End
!
    PUSH_SUB(pedra)

     if (i_count == 0) then
      if (n_tes == 1) then
          write(pcminfo_unit,'(2X,A32)') 'Number of tesserae / sphere = 60'
          write(pcminfo_unit,'(2X)') 
      else
          write(pcminfo_unit,'(2X,A33)') 'Number of tesserae / sphere = 240' 
          write(pcminfo_unit,'(2X)') 
      endif 
     endif

!      DR = DR*P_a_B ! A.U. -> Angstrom
!
!  PEDRA prevede che i dati geometrici siano espressi in ANGSTROM :
!  vengono trasformati, e solo alla fine i risultati tornano in bohr.
!
   90 CONTINUE !maybe it can be removed

      sfe(:)%X=sfe(:)%X*P_a_B  
      sfe(:)%Y=sfe(:)%Y*P_a_B  
      sfe(:)%Z=sfe(:)%Z*P_a_B  
      sfe(:)%R=sfe(:)%R*P_a_B  

!
!     ----- Partition of the cavity surface into tesserae -----
!
!      VOL=ZERO
!      STOT=ZERO

      VOL=M_ZERO
      STOT=M_ZERO

      jvt1=reshape(idum,(/6,60/))

! FIXME: much of the code below appears to come from some other (old) source.
! Where does it come from? Is it compatible with our GPL license?

!
!*****COORDINATES OF VERTICES OF TESSERAE IN A SPHERE WITH UNIT RADIUS.
!
!     Vengono memorizzati i vertici (nella matrice CV) e i centri (nei
!     vettori XC,YC,ZC) di 240 tessere (60 grandi divise in 4 piu
!     piccole) La matrice JVT1(i,j) indica quale e' il numero d'ordine
!     del vertice i-esimo della j-esima tessera grande. In ogni tessera
!     grande i 6 vertici sono cosi disposti:
!
!                                    1
!
!                                 4     5
!
!                              3     6     2
!
      CV(1,1)=0.0D+00
      CV(1,2)=0.0D+00
      CV(1,3)=1.0D+00
      CV(122,1)=0.0D+00
      CV(122,2)=0.0D+00
      CV(122,3)=-1.0D+00
      II=1
      DO 200 I=1,24
      TH=THEV(I)
      FI=FIV(I)
      CTH=COS(TH)
      STH=SIN(TH)
      DO 210 J=1,5
      FI=FI+FIR
      IF(J.EQ.1) FI=FIV(I)
      II=II+1
      CV(II,1)=STH*COS(FI)
      CV(II,2)=STH*SIN(FI)
      CV(II,3)=CTH
  210 CONTINUE
  200 CONTINUE
!
!     Controlla se ciascuna tessera e scoperta o va tagliata
!
      NN = 0
      DO 300 NSFE = 1, NESF
      XEN = sfe(NSFE)%x
      YEN = sfe(NSFE)%y
      ZEN = sfe(NSFE)%z
      REN = sfe(NSFE)%r
      XCTST(:) = M_ZERO
      YCTST(:) = M_ZERO
      ZCTST(:) = M_ZERO
      AST(:) = M_ZERO
!
      DO 310 ITS = 1, 60
!
!
      do i_tes=1,n_tes
      if (n_tes.eq.1) then
      N1 = JVT1(1,ITS)
      N2 = JVT1(2,ITS)
      N3 = JVT1(3,ITS)
      else
        if (i_tes.eq.1) then
          N1 = JVT1(1,ITS)
          N2 = JVT1(5,ITS)
          N3 = JVT1(4,ITS)
        elseif (i_tes.eq.2) then 
          N1 = JVT1(4,ITS)
          N2 = JVT1(6,ITS)
          N3 = JVT1(3,ITS)
        elseif (i_tes.eq.3)  then
          N1 = JVT1(4,ITS)
          N2 = JVT1(5,ITS)
          N3 = JVT1(6,ITS)
        elseif (i_tes.eq.4)  then
          N1 = JVT1(2,ITS)
          N2 = JVT1(6,ITS)
          N3 = JVT1(5,ITS)
         endif
      endif
      PTS(1,1)=CV(N1,1)*REN+XEN
      PTS(2,1)=CV(N1,3)*REN+YEN
      PTS(3,1)=CV(N1,2)*REN+ZEN
      PTS(1,2)=CV(N2,1)*REN+XEN
      PTS(2,2)=CV(N2,3)*REN+YEN
      PTS(3,2)=CV(N2,2)*REN+ZEN
      PTS(1,3)=CV(N3,1)*REN+XEN
      PTS(2,3)=CV(N3,3)*REN+YEN
      PTS(3,3)=CV(N3,2)*REN+ZEN
      PP(:) = M_ZERO
      PP1(:) = M_ZERO
      NV=3
!
!     Per ciascuna tessera, trova la porzione scoperta e ne
!     calcola l area con il teorema di Gauss-Bonnet; il punto
!     rappresentativo e definito come media dei vertici della
!     porzione scoperta di tessera e passato in PP (mentre in PP1
!     ci sono le coordinate del punto sulla normale interna).
!
!     I vertici di ciascuna tessera sono conservati in VERT(MXTS,10,3),
!     il numero di vertici di ciascuna tessera e in NVERT(MXTS), e i
!     centri dei cerchi di ciascun lato sono in CENTR(MXTS,10,3).
!     In INTSPH(numts,10) sono registrate le sfere a cui appartengono
!     i lati delle tessere.
!
      CALL SUBTESSERA(sfe,nsfe,nesf,NV,PTS,CCC,PP,PP1,AREA)
!
      IF(AREA.EQ.M_ZERO) cycle
      XCTST(n_tes*(ITS-1)+i_tes) = PP(1)
      YCTST(n_tes*(ITS-1)+i_tes) = PP(2)
      ZCTST(n_tes*(ITS-1)+i_tes) = PP(3)
      nctst(:,n_tes*(its-1)+i_tes)=pp1(:)
      AST(n_tes*(ITS-1)+i_tes) = AREA
      isfet(n_tes*(its-1)+i_tes) = nsfe
      enddo
 310  CONTINUE
!
!
!
      DO 320 ITS=1,60*n_tes
      IF(AST(ITS).EQ.0.0D+00) GO TO 320
      NN = NN + 1
!
!     check on the total number of tessera
!
!      IF(NN.GT.MXTS) THEN                        !! CHECK WITH STEFANO
!         WRITE(6,*) ' TOO MANY TESSERAE IN PEDRA'
!         write(message(1),'(a,I5)') "WARNING: total number of tesserae=", NN
!         call messages_info(1)     
!      END IF
!
      if (i_count.eq.1) then
       CTS(NN)%x = XCTST(ITS)
       CTS(NN)%y = YCTST(ITS)
       CTS(NN)%z = ZCTST(ITS)
       CTS(NN)%n(:)=nctst(:,its)
       CTS(NN)%area = AST(ITS)
       CTS(NN)%rsfe = sfe(isfet(ITS))%r
      endif
 320  CONTINUE
 300  CONTINUE
      NTS = NN
      if (i_count.eq.1) then
!
!
!     Verifica se due tessere sono troppo vicine
!
      TEST = 0.10D+00
      TEST2 = TEST*TEST
450   Continue
      DO 400 I = 1, NTS-1
      IF(cts(I)%area.EQ.M_ZERO) GO TO 400
      XI = CTS(I)%x
      YI = CTS(I)%y
      ZI = CTS(I)%z
      II = I + 1
      DO 410 J = II , NTS
      IF(cts(J)%area.EQ.M_ZERO) GO TO 410
      XJ = CTS(J)%x
      YJ = CTS(J)%y
      ZJ = CTS(J)%z
      RIJ = (XI-XJ)**2 + (YI-YJ)**2 + (ZI-ZJ)**2
      IF(RIJ.GT.TEST2) GO TO 410
!
!
      !WRITE(6,9010) I,J,SQRT(RIJ),TEST
      write(pcminfo_unit,9010) I,J,SQRT(RIJ),TEST
!SC 19/8
      XI=(XI*CTS(I)%area+XJ*CTS(J)%area)/(CTS(I)%area+CTS(J)%area)
      YI=(YI*CTS(I)%area+YJ*CTS(J)%area)/(CTS(I)%area+CTS(J)%area)
      ZI=(ZI*CTS(I)%area+ZJ*CTS(J)%area)/(CTS(I)%area+CTS(J)%area)
      CTS(I)%x=XI
      CTS(I)%y=YI
      CTS(I)%z=ZI
      CTS(I)%n=(CTS(I)%n*CTS(I)%area+CTS(J)%n*CTS(J)%area)
      DNORM=sqrt(dot_product(CTS(I)%n,CTS(I)%n))
      CTS(I)%n=CTS(I)%n/DNORM
      CTS(I)%rsfe=(CTS(I)%rsfe*CTS(I)%area+CTS(J)%rsfe*CTS(J)%area)/&
                  (CTS(I)%area+CTS(J)%area)
      CTS(I)%area=CTS(I)%area+CTS(J)%area
! Delete Tessera J
      Do K=J+1,NTS
       CTS(K-1)=CTS(K)
      Enddo
      NTS=NTS-1
      GoTo 450
      
      
 410  CONTINUE
 400  CONTINUE
!
!     Calcola il volume della cavita con la formula (t. di Gauss):
!                V=SOMMAsulleTESSERE{A r*n}/3
!     dove r e' la distanza del punto rappresentativo dall'origine,
!     n e' il versore normale alla tessera, A l'area della tessera,
!     e * indica il prodotto scalare.
!
      VOL = M_ZERO
      DO ITS = 1, NTS
!
!
!     Trova il prodotto scalare
!
         PROD = CTS(ITS)%x*cts(its)%n(1) + CTS(ITS)%y*cts(its)%n(2) + &
                CTS(ITS)%z*cts(its)%n(3)
         VOL = VOL + cts(ITS)%area * PROD / 3.0D+00
         stot = stot + cts(ITS)%area
      ENDDO
!
!     Stampa la geometria della cavita
!
!
         write(pcminfo_unit,9040) NTS,STOT,VOL
!
!
!     Trasform results in bohr
!     ========================
      cts(:)%area=cts(:)%area*(P_Ang)**2
      cts(:)%x=cts(:)%x*P_Ang
      cts(:)%y=cts(:)%y*P_Ang
      cts(:)%z=cts(:)%z*P_Ang
      cts(:)%rsfe=cts(:)%rsfe*P_Ang
      endif
!
      sfe(:)%x=sfe(:)%x*P_Ang 
      sfe(:)%y=sfe(:)%y*P_Ang
      sfe(:)%z=sfe(:)%z*P_Ang
      sfe(:)%r=sfe(:)%r*P_Ang

      POP_SUB(pedra)
      return
!
!
 9010 FORMAT(2X,'WARNING: the distance between tesserae ',I6, &
             ' and ',I6,' is ',F8.4,' A, less than ',F8.4,' A')

 9040 FORMAT(/2X,'Total number of tesserae = ',I8/ &
              2X,'Surface area = ',F14.8,'(A^2)',4X,'Cavity volume = ', &
                  F14.8,' (A^3)')

       END subroutine pedra
!==============================================================

      SUBROUTINE subtessera(sfe,ns,nesf,NV,PTS,CCC,PP,PP1,AREA)
!
      type(sphere_t), intent(in) :: sfe(:)

      ! FIXME: please put intents here.
      FLOAT  :: PTS(3,10),CCC(3,10),PP(3),PP1(3),area

      integer :: INTSPH(10),ns,nesf,nv,nsfe1,n,i,j,icop
      integer :: l,iv1,iv2,ii,icut,jj

      FLOAT  :: P1(3),P2(3),P3(3),P4(3),POINT(3),  &
                PSCR(3,10),CCCP(3,10),POINTL(3,10)

      integer :: IND(10),LTYP(10),INTSCR(10)
      integer, parameter :: MXTS=5000

      FLOAT  :: delr,delr2,rc,rc2,dnorm,dist,de2
      FLOAT, parameter :: tol= -1.d-10
!
! FIXME: please translate comments to English!
!
!     Coord. del centro che sottende l`arco tra i vertici
!     n e n+1 (per i primi tre vertici e sicuramente il centro della
!     sfera) e sfera alla cui intersezione con NS appartiene l arco (se
!     appartiene alla sfera originaria INTSPH(numts,N)=NS)
!
      AREA = 0.0D+00
      DO J=1, 3
        CCC(1,J) = sfe(NS)%x
        CCC(2,J) = sfe(NS)%y
        CCC(3,J) = sfe(NS)%z
      ENDDO
!
!     INTSPH viene riferito alla tessera -numts-, e in seguito riceve il
!     numero corretto.
!
      DO N = 1, 3
        INTSPH(N) = NS
      ENDDO
!
!     Loop sulle altre sfere
!
      DO 150 NSFE1=1,NESF
      IF(NSFE1.EQ.NS) GO TO 150
!
!     Memorizza i vertici e i centri che sottendono gli archi
!
      DO J =1, NV
        INTSCR(J) = INTSPH(J)
      DO I = 1,3
        PSCR(I,J) = PTS(I,J)
        CCCP(I,J) = CCC(I,J)
      ENDDO
      ENDDO
!
      ICOP = 0
      DO J =1, 10
        IND(J) = 0
        LTYP(J) = 0
      ENDDO
!
!     Loop sui vertici della tessera considerata
!
      DO 100 I=1,NV
        DELR2=(PTS(1,I)-sfe(NSFE1)%x)**2+(PTS(2,I)-sfe(NSFE1)%y)**2+&
        (PTS(3,I)-sfe(NSFE1)%z)**2
        DELR=SQRT(DELR2)
        IF(DELR.LT.sfe(NSFE1)%r) THEN
          IND(I) = 1
          ICOP = ICOP+1
        END IF
 100  CONTINUE
!     Se la tessera e completamente coperta, la trascura
      IF(ICOP.EQ.NV) RETURN
!                    ******
!
!     Controlla e classifica i lati della tessera: LTYP = 0 (coperto),
!     1 (tagliato con il II vertice coperto), 2 (tagliato con il I
!     vertice coperto), 3 (bitagliato), 4 (libero)
!     Loop sui lati
!
      DO L = 1, NV
        IV1 = L
        IV2 = L+1
        IF(L.EQ.NV) IV2 = 1
        IF(IND(IV1).EQ.1.AND.IND(IV2).EQ.1) THEN
          LTYP(L) = 0
        ELSE IF(IND(IV1).EQ.0.AND.IND(IV2).EQ.1) THEN
          LTYP(L) = 1
        ELSE IF(IND(IV1).EQ.1.AND.IND(IV2).EQ.0) THEN
          LTYP(L) = 2
        ELSE IF(IND(IV1).EQ.0.AND.IND(IV2).EQ.0) THEN
          LTYP(L) = 4
!
          RC2 = (CCC(1,L)-PTS(1,L))**2 + (CCC(2,L)-PTS(2,L))**2 + &
                (CCC(3,L)-PTS(3,L))**2
          RC = SQRT(RC2)
!
!     Su ogni lato si definiscono 11 punti equispaziati, che vengono
!     controllati
!
          DO II = 1, 11
          POINT(1) = PTS(1,IV1) + II * (PTS(1,IV2)-PTS(1,IV1)) / 11
          POINT(2) = PTS(2,IV1) + II * (PTS(2,IV2)-PTS(2,IV1)) / 11
          POINT(3) = PTS(3,IV1) + II * (PTS(3,IV2)-PTS(3,IV1)) / 11
          POINT(1) = POINT(1) - CCC(1,L)
          POINT(2) = POINT(2) - CCC(2,L)
          POINT(3) = POINT(3) - CCC(3,L)
          DNORM = SQRT(POINT(1)**2 + POINT(2)**2 + POINT(3)**2)
          POINT(1) = POINT(1) * RC / DNORM + CCC(1,L)
          POINT(2) = POINT(2) * RC / DNORM + CCC(2,L)
          POINT(3) = POINT(3) * RC / DNORM + CCC(3,L)
          DIST = SQRT( (POINT(1)-sfe(NSFE1)%x)**2 + &
          (POINT(2)-sfe(NSFE1)%y)**2 + (POINT(3)-sfe(NSFE1)%z)**2 )
          IF((DIST - sfe(NSFE1)%r) .LT. TOL) THEN
!         IF(DIST.LT.sfe(NSFE1)%r) then
            LTYP(L) = 3
            DO JJ = 1, 3
              POINTL(JJ,L) = POINT(JJ)
            ENDDO
            GO TO 160
          END IF
          ENDDO
        END IF
 160    CONTINUE
      ENDDO
!
!     Se la tessera e spezzata in due o piu tronconi, la trascura
!
      ICUT = 0
      DO L = 1, NV
        IF(LTYP(L).EQ.1.OR.LTYP(L).EQ.2) ICUT = ICUT + 1
        IF(LTYP(L).EQ.3) ICUT = ICUT + 2
      ENDDO
      ICUT = ICUT / 2
      IF(ICUT.GT.1) RETURN
!
!     Creazione dei nuovi vertici e lati della tessera
!     Loop sui lati
!
      N = 1
      DO 300 L = 1, NV
!     Se il lato L e coperto:
        IF(LTYP(L).EQ.0) GO TO 300
        IV1 = L
        IV2 = L+1
        IF(L.EQ.NV) IV2 = 1
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Se il lato L e tagliato (con il I vertice scoperto):
        IF(LTYP(L).EQ.1) THEN
        DO JJ = 1, 3
          PTS(JJ,N) = PSCR(JJ,IV1)
          CCC(JJ,N) = CCCP(JJ,IV1)
        ENDDO
        INTSPH(N) = INTSCR(IV1)
        N = N+1
!
!     Trova l intersezione tra i due vertici del lato L
!
!     P1 = coord. del primo vertice
!     P2 = coord. del secondo vertice
!     P3 = coord. del centro dell`arco sotteso
!     P4 = coord. dell intersezione
!
        DO JJ = 1, 3
          P1(JJ) = PSCR(JJ,IV1)
          P2(JJ) = PSCR(JJ,IV2)
          P3(JJ) = CCCP(JJ,IV1)
        ENDDO
        CALL INTER(sfe,P1,P2,P3,P4,NSFE1,0)
!     Aggiorna i vertici della tessera e il centro dell arco
        DO JJ = 1,3
          PTS(JJ,N) = P4(JJ)
        ENDDO
!
!     Il nuovo arco sara sotteso tra questo e il prossimo punto
!     di intersezione: il centro che lo sottende
!     sara il centro del cerchio di intersezione tra la sfera NS
!     e la sfera NSFE1.
!
        DE2 = (sfe(NSFE1)%x-sfe(NS)%x)**2+(sfe(NSFE1)%y-sfe(NS)%y)**2+ &
              (sfe(NSFE1)%z-sfe(NS)%z)**2
        CCC(1,N)=sfe(NS)%x+(sfe(NSFE1)%x-sfe(NS)%x)* &
                 (sfe(NS)%r**2-sfe(NSFE1)%r**2+DE2)/(2.0D+00*DE2)
        CCC(2,N)=sfe(NS)%y+(sfe(NSFE1)%y-sfe(NS)%y)* &
                 (sfe(NS)%r**2-sfe(NSFE1)%r**2+DE2)/(2.0D+00*DE2)
        CCC(3,N)=sfe(NS)%z+(sfe(NSFE1)%z-sfe(NS)%z)* &
                 (sfe(NS)%r**2-sfe(NSFE1)%r**2+DE2)/(2.0D+00*DE2)
        INTSPH(N) = NSFE1
        N = N+1
        END IF
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Se il lato L e tagliato (con il II vertice scoperto):
        IF(LTYP(L).EQ.2) THEN
!     Trova l intersezione tra i due vertici del lato L
!
!     P1 = coord. del primo vertice
!     P2 = coord. del secondo vertice
!     P3 = coord. del centro dell`arco sotteso
!     P4 = coord. dell intersezione
!
        DO JJ = 1, 3
          P1(JJ) = PSCR(JJ,IV1)
          P2(JJ) = PSCR(JJ,IV2)
          P3(JJ) = CCCP(JJ,IV1)
        ENDDO
        CALL INTER(sfe,P1,P2,P3,P4,NSFE1,1)
!     Aggiorna i vertici della tessera e il centro dell arco
        DO JJ = 1,3
          PTS(JJ,N) = P4(JJ)
          CCC(JJ,N) = CCCP(JJ,IV1)
        ENDDO
        INTSPH(N) = INTSCR(IV1)
        N = N+1
        END IF
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Se il lato e intersecato due volte:
        IF(LTYP(L).EQ.3) THEN
        DO JJ = 1, 3
          PTS(JJ,N) = PSCR(JJ,IV1)
          CCC(JJ,N) = CCCP(JJ,IV1)
        ENDDO
        INTSPH(N) = INTSCR(IV1)
        N = N+1
!
!     Trova l intersezione tra il primo vertice e un punto intermedio
!     coperto
!
!     P1 = coord. del primo vertice
!     P2 = coord. del secondo vertice
!     P3 = coord. del centro dell`arco sotteso
!     P4 = coord. dell intersezione
!
        DO JJ = 1, 3
          P1(JJ) = PSCR(JJ,IV1)
          P2(JJ) = POINTL(JJ,L)
          P3(JJ) = CCCP(JJ,IV1)
        ENDDO
        CALL INTER(sfe,P1,P2,P3,P4,NSFE1,0)
!     Aggiorna i vertici della tessera e il centro dell arco
        DO JJ = 1,3
          PTS(JJ,N) = P4(JJ)
        ENDDO
!
!     Il nuovo arco sara sotteso tra questo e il prossimo punto
!     di intersezione: il centro che lo sottende
!     sara il centro del cerchio di intersezione tra la sfera NS
!     e la sfera NSFE1.
!
        DE2 = (sfe(NSFE1)%x-sfe(NS)%x)**2+(sfe(NSFE1)%y-sfe(NS)%y)**2+ &
              (sfe(NSFE1)%z-sfe(NS)%z)**2
        CCC(1,N)=sfe(NS)%x+(sfe(NSFE1)%x-sfe(NS)%x)* &
                 (sfe(NS)%r**2-sfe(NSFE1)%r**2+DE2)/(2.0D+00*DE2)
        CCC(2,N)=sfe(NS)%y+(sfe(NSFE1)%y-sfe(NS)%y)* &
                 (sfe(NS)%r**2-sfe(NSFE1)%r**2+DE2)/(2.0D+00*DE2)
        CCC(3,N)=sfe(NS)%z+(sfe(NSFE1)%z-sfe(NS)%z)* &
                 (sfe(NS)%r**2-sfe(NSFE1)%r**2+DE2)/(2.0D+00*DE2)
        INTSPH(N) = NSFE1
        N = N+1
!
!     Trova l intersezione tra un punto intermedio coperto e il
!     secondo vertice
!
!     P1 = coord. del primo vertice
!     P2 = coord. del secondo vertice
!     P3 = coord. del centro dell`arco sotteso
!     P4 = coord. dell intersezione
!
        DO JJ = 1, 3
          P1(JJ) = POINTL(JJ,L)
          P2(JJ) = PSCR(JJ,IV2)
          P3(JJ) = CCCP(JJ,IV1)
        ENDDO
        CALL INTER(sfe,P1,P2,P3,P4,NSFE1,1)
!     Aggiorna il vertice e il centro dell arco
        DO JJ = 1,3
          PTS(JJ,N) = P4(JJ)
          CCC(JJ,N) = CCCP(JJ,IV1)
        ENDDO
        INTSPH(N) = INTSCR(IV1)
        N = N + 1
        END IF
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Se il lato e scoperto:
        IF(LTYP(L).EQ.4) THEN
        DO JJ = 1, 3
          PTS(JJ,N) = PSCR(JJ,IV1)
          CCC(JJ,N) = CCCP(JJ,IV1)
        ENDDO
        INTSPH(N) = INTSCR(IV1)
        N = N+1
        END IF
!
 300  CONTINUE
!
      NV = N - 1
!     Controlla che il numero di vertici creati non sia eccessivo
      IF(NV.GT.10) THEN
         WRITE(6,*) 'TROPPI VERTICI CREATI IN TESSERA: BYE BYE...'
         STOP
      END IF
 150  CONTINUE
!
!     Se la tessera non e stata scartata, a questo punto ne troviamo
!     l area e il punto rappresentativo
!
      CALL GAUBON(sfe,NV,NS,PTS,CCC,PP,PP1,AREA,INTSPH)
      RETURN
      END SUBROUTINE subtessera
!==============================================================

      SUBROUTINE INTER(sfe,P1,P2,P3,P4,NS,I)
!
      type(sphere_t), intent(in) :: sfe(:)

      FLOAT :: P1(:),P2(:),P3(:),P4(:)
      FLOAT :: r2,r,alpha,delta,dnorm,diff,diff2
      FLOAT, parameter :: TOL = 1.0D-08

      integer :: i,ns,m,jj
!
      FLOAT, parameter :: MXTS=5000
!
!
!     Trova il punto P4, sull`arco P1-P2 sotteso dal centro P3, che
!     si trova sulla superficie della sfera NS
!     P4 e definito come combinazioe lineare di P1 e P2, con
!     il parametro ALPHA ottimizzato per tentativi.
!
      R2 = (P1(1)-P3(1))**2+(P1(2)-P3(2))**2+(P1(3)-P3(3))**2
      R = SQRT(R2)
      ALPHA = 0.5D+00
      DELTA = 0.0D+00
      M = 1
  10  CONTINUE
      IF (M.GT.1000) THEN
         WRITE(6,*) 'TROPPE ITERAZIONI IN INTER! BYE BYE ....'
         STOP
      END IF
      ALPHA = ALPHA + DELTA
      DNORM = 0.0D+00
      DO JJ = 1,3
       P4(JJ)=P1(JJ)+ALPHA*(P2(JJ)-P1(JJ))-P3(JJ)
       DNORM = DNORM + P4(JJ)**2
      ENDDO
      DNORM = SQRT(DNORM)
      DO JJ = 1,3
       P4(JJ)= P4(JJ)*R/DNORM + P3(JJ)
      ENDDO
      DIFF2=(P4(1)-sfe(NS)%x)**2 + (P4(2)-sfe(NS)%y)**2 + (P4(3)-sfe(NS)%z)**2
      DIFF = SQRT(DIFF2) - sfe(NS)%r
!
      IF(ABS(DIFF).LT.TOL) RETURN
!                          ******
      IF(I.EQ.0) THEN
       IF(DIFF.GT.0.0D+00) DELTA =  1.0D+00/(2.0D+00**(M+1))
       IF(DIFF.LT.0.0D+00) DELTA = -1.0D+00/(2.0D+00**(M+1))
       M = M + 1
       GO TO 10
      END IF
      IF(I.EQ.1) THEN
       IF(DIFF.GT.0.0D+00) DELTA = -1.0D+00/(2.0D+00**(M+1))
       IF(DIFF.LT.0.0D+00) DELTA =  1.0D+00/(2.0D+00**(M+1))
       M = M + 1
       GO TO 10
      END IF
!          the code probably never reaches this return
      RETURN
      END subroutine inter
!==============================================================

      SUBROUTINE GAUBON(sfe,NV,NS,PTS,CCC,PP,PP1,AREA,INTSPH)
!
      FLOAT, parameter :: MXTS=5000
!
      type(sphere_t), intent(in) :: sfe(:)

      FLOAT :: PTS(3,10),CCC(3,10),PP(3),PP1(3),area

      integer :: INTSPH(:),nv,ns

      FLOAT ::  P1(3),P2(3),P3(3),U1(3),U2(3)
      FLOAT :: tpi,sum1,x1,y1,z1,x2,y2,z2,dnorm,dnorm1,dnorm2, &
                dnorm3,scal
      FLOAT :: cosphin,phin,costn,sum2,betan

      integer :: nsfe1,i,jj,n,n0,n1,n2
!
!     Sfrutta il teorema di Gauss-Bonnet per calcolare l area
!     della tessera con vertici PTS(3,NV). Consideriamo sempre
!     che il lato N della tessera e quello compreso tra i vertici
!     N e N+1 (oppure NV e 1). In CCC(3,NV) sono le posizioni dei
!     centri degli archi che sottendono i vari lati della tessera.
!     La formula di Gauss-Bonet per le sfere e':
!            Area=R^2[2pi+S(Phi(N)cosT(N))-S(Beta(N))]
!     dove Phi(N) e la lunghezza d arco (in radianti) del lato N,
!     T(N) e l angolo polare del lato N, Beta(N) l angolo esterno
!     relativo al vertice N.
!
!      TPI=2*PI
      TPI=2*M_Pi
!
!     Calcola la prima sommatoria
      SUM1 = M_ZERO
      DO 100 N = 1, NV
      X1 = PTS(1,N) - CCC(1,N)
      Y1 = PTS(2,N) - CCC(2,N)
      Z1 = PTS(3,N) - CCC(3,N)
      IF(N.LT.NV) THEN
        X2 = PTS(1,N+1) - CCC(1,N)
        Y2 = PTS(2,N+1) - CCC(2,N)
        Z2 = PTS(3,N+1) - CCC(3,N)
      ELSE
        X2 = PTS(1,1) - CCC(1,N)
        Y2 = PTS(2,1) - CCC(2,N)
        Z2 = PTS(3,1) - CCC(3,N)
      END IF
      DNORM1 = X1*X1 + Y1*Y1 + Z1*Z1
      DNORM2 = X2*X2 + Y2*Y2 + Z2*Z2
      SCAL = X1*X2 + Y1*Y2 + Z1*Z2
      COSPHIN = SCAL / (SQRT(DNORM1*DNORM2))
      IF(COSPHIN.GT.1.0D+00) COSPHIN = 1.0D+00
      IF(COSPHIN.LT.-1.0D+00) COSPHIN = -1.0D+00
      PHIN = ACOS(COSPHIN)
!
!     NSFE1 e la sfera con cui la sfera NS si interseca (eventualmente)
        NSFE1 = INTSPH(N)
        X1 = sfe(NSFE1)%x - sfe(NS)%x
        Y1 = sfe(NSFE1)%y - sfe(NS)%y
        Z1 = sfe(NSFE1)%z - sfe(NS)%z
      DNORM1 = SQRT(X1*X1 + Y1*Y1 + Z1*Z1)
      IF(DNORM1.EQ.M_ZERO) DNORM1 = 1.0D+00
        X2 = PTS(1,N) - sfe(NS)%x
        Y2 = PTS(2,N) - sfe(NS)%y
        Z2 = PTS(3,N) - sfe(NS)%z
      DNORM2 = SQRT(X2*X2 + Y2*Y2 + Z2*Z2)
      COSTN = (X1*X2+Y1*Y2+Z1*Z2)/(DNORM1*DNORM2)
      SUM1 = SUM1 + PHIN * COSTN
 100  CONTINUE
!
!     Calcola la seconda sommatoria: l angolo esterno Beta(N) e
!     definito usando i versori (u(N-1),u(N)) tangenti alla sfera
!     nel vertice N lungo le direzioni dei lati N-1 e N:
!                cos( Pi-Beta(N) )=u(N-1)*u(N)
!            u(N-1) = [V(N) x (V(N) x V(N-1))]/NORM
!            u(N) = [V(N) x (V(N) x V(N+1))]/NORM
!     dove V(I) e il vettore posizione del vertice I RISPETTO AL
!     CENTRO DELL ARCO CHE SI STA CONSIDERANDO.
!
      SUM2 = M_ZERO
!     Loop sui vertici
      DO 200 N = 1, NV
      DO JJ = 1, 3
      P1(JJ) = M_ZERO
      P2(JJ) = M_ZERO
      P3(JJ) = M_ZERO
      ENDDO
      N1 = N
      IF(N.GT.1) N0 = N - 1
      IF(N.EQ.1) N0 = NV
      IF(N.LT.NV) N2 = N + 1
      IF(N.EQ.NV) N2 = 1
!     Trova i vettori posizione rispetto ai centri corrispondenti
!     e i versori tangenti
!
!     Lato N0-N1:
      DO JJ = 1, 3
      P1(JJ) = PTS(JJ,N1) - CCC(JJ,N0)
      P2(JJ) = PTS(JJ,N0) - CCC(JJ,N0)
      ENDDO
!
      CALL VECP(P1,P2,P3,DNORM3)
      DO JJ = 1, 3
      P2(JJ) = P3(JJ)
      ENDDO
      CALL VECP(P1,P2,P3,DNORM3)
      DO JJ = 1, 3
      U1(JJ) = P3(JJ)/DNORM3
      ENDDO
!
!     Lato N1-N2:
      DO JJ = 1, 3
      P1(JJ) = PTS(JJ,N1) - CCC(JJ,N1)
      P2(JJ) = PTS(JJ,N2) - CCC(JJ,N1)
      ENDDO
!
      CALL VECP(P1,P2,P3,DNORM3)
      DO JJ = 1, 3
      P2(JJ) = P3(JJ)
      ENDDO
      CALL VECP(P1,P2,P3,DNORM3)
      DO JJ = 1, 3
      U2(JJ) = P3(JJ)/DNORM3
      ENDDO
!
      BETAN = ACOS(U1(1)*U2(1)+U1(2)*U2(2)+U1(3)*U2(3))
!      SUM2 = SUM2 + (PI - BETAN)
      SUM2 = SUM2 + (M_Pi - BETAN)
 200  CONTINUE
!     Calcola l area della tessera
        AREA = sfe(NS)%r*sfe(NS)%r*(TPI + SUM1 - SUM2)
!     Trova il punto rappresentativo (come media dei vertici)
      DO JJ = 1, 3
      PP(JJ) = M_ZERO
      ENDDO
      DO I = 1, NV
      PP(1) = PP(1) + (PTS(1,I)-sfe(NS)%x)
      PP(2) = PP(2) + (PTS(2,I)-sfe(NS)%y)
      PP(3) = PP(3) + (PTS(3,I)-sfe(NS)%z)
      ENDDO
      DNORM = M_ZERO
      DO JJ = 1, 3
      DNORM = DNORM + PP(JJ)*PP(JJ)
      ENDDO
      PP(1) = sfe(NS)%x + PP(1) * sfe(NS)%r / SQRT(DNORM)
      PP(2) = sfe(NS)%y + PP(2) * sfe(NS)%r / SQRT(DNORM)
      PP(3) = sfe(NS)%z + PP(3) * sfe(NS)%r / SQRT(DNORM)
!     Trova la normale (interna!) nel punto rappresentativo
      PP1(1) = (PP(1) - sfe(NS)%x) / sfe(NS)%r
      PP1(2) = (PP(2) - sfe(NS)%y) / sfe(NS)%r
      PP1(3) = (PP(3) - sfe(NS)%z) / sfe(NS)%r
!
!     A causa delle approssimazioni numeriche, l area di alcune piccole
!     tessere puo risultare negativa, e viene in questo caso trascurata
      IF(AREA.LT.M_ZERO)THEN
        WRITE(6,1000) NS,AREA
 1000   FORMAT(1X,'WARNING: THE AREA OF A TESSERA ON SPHERE ',I6, &
        ' IS NEGATIVE (',E10.3,' A^2 ), THUS DISCARDED')
        AREA = M_ZERO
      END IF
      RETURN
      END SUBROUTINE gaubon
!==============================================================

      SUBROUTINE VECP(P1,P2,P3,DNORM3)
!
      FLOAT :: P1(3),P2(3),P3(3)
      FLOAT :: dnorm3
!
!     Esegue il prodotto vettoriale P3 = P1 x P2
!
      P3(1) = P1(2)*P2(3) - P1(3)*P2(2)
      P3(2) = P1(3)*P2(1) - P1(1)*P2(3)
      P3(3) = P1(1)*P2(2) - P1(2)*P2(1)
      DNORM3 = SQRT(P3(1)*P3(1) + P3(2)*P3(2) + P3(3)*P3(3))

      RETURN
      END SUBROUTINE vecp

END MODULE pcm_m
