subroutine xc_gga(func, m, nspin, rho, rho_core, pot, energy)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: func, nspin
  real(r8), intent(IN) :: rho(m%np, nspin), rho_core(m%np)
  real(r8), intent(out) :: pot(m%np, nspin), energy

  real(r8) :: dedd(nspin), dedgd(3, nspin), e, dvol, den(3)
  real(r8), allocatable :: d(:,:), gd(:,:,:)
  integer :: i, is, in, ic, ind(3)

  sub_name = 'xc_gga'; call push_sub()

  allocate(d(0:m%np, nspin), gd(3, m%np, nspin))
  do is = 1, nspin
    d(0, is) = 0._r8
    d(1:m%np, is) = rho(1:m%np, is) + rho_core(1:m%np)/nspin
    call dmesh_derivatives(m, d(:, is), grad=gd(:,:, is))
  end do

  energy = 0._r8
  do i = 1, m%np
    select case(func)      
    case(X_FUNC_GGA_PBE)
      call xc_x_pbe(.false., nspin, d(i, :), gd(i,:,:), e, dedd, dedgd) 
      
    case(X_FUNC_GGA_PBER)
      call xc_x_pbe(.true.,  nspin, d(i, :), gd(i,:,:), e, dedd, dedgd) 

    case(X_FUNC_GGA_LB94)
      call xc_x_lb94(nspin, d(i, :), gd(i,:,:), e, dedd, dedgd) 
      
    case(C_FUNC_GGA_PBE)
      call xc_c_pbe(nspin, d(i, :), gd(i,:,:), e, dedd, dedgd) 
      
    end select
    energy = energy + sum(d(i, :)) * e * m%vol_pp

    do is = 1, nspin
      pot(i, is) = pot(i, is) + dedd(is)

      do in = -m%d%norder , m%d%norder
        ind(1) = m%Lxyz_inv(m%Lx(i)+in,m%Ly(i),m%Lz(i))
        ind(2) = m%Lxyz_inv(m%Lx(i),m%Ly(i)+in,m%Lz(i))
        ind(3) = m%Lxyz_inv(m%Lx(i),m%Ly(i),m%Lz(i)+in)

#ifndef BOUNDARIES_ZERO_DERIVATIVE
        den = 0.0_r8
#endif
        if(ind(1) > 0)den(1) = d(ind(1), is)
        if(ind(2) > 0)den(2) = d(ind(2), is)
        if(ind(3) > 0)den(3) = d(ind(3), is)
        
        do ic = 1, 3
          if(ind(ic) > 0) then
            pot(ind(ic), is) = pot(ind(ic), is) + &
                 dedgd(ic, is) * m%d%dgidfj(in)/ m%h
          end if
        end do
      end do
    end do
  end do

  if(func == X_FUNC_GGA_LB94) then ! we have to calculate the energy
    ! Levy-Perdew relation (PRA 32, 2010 (1985))
    energy = 0._r8
    do is = 1, nspin
      call dmesh_derivatives(m, pot(:, is), grad=gd(:, :, is))
      do i = 1, m%np
        energy = energy + d(i, is) * m%h * (  &
             m%Lx(i)*gd(i, 1, is) + m%Ly(i)*gd(i, 2, is) + m%Lz(i)*gd(i, 3, is))
      end do
    end do
    energy = - energy * m%vol_pp
  end if

  call pop_sub()

  return
end subroutine xc_gga


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Implements Perdew-Burke-Ernzerhof Generalized-Gradient-Approximation.       !
! Ref: J.P.Perdew, K.Burke & M.Ernzerhof, PRL 77, 3865 (1996)                 !
! Written by L.C.Balbas and J.M.Soler. December 1996. Version 0.5.            !
! ******** INPUT ******************************************************       !
! INTEGER IREL           : Relativistic-exchange switch (0=No, 1=Yes)         !
! INTEGER NSPIN          : Number of spin polarizations (1 or 2)              !
! REAL*8  DENS(NSPIN)    : Total electron density (if NSPIN=1) or             !
!                           spin electron density (if NSPIN=2)                !
! REAL*8  GDENS(3,NSPIN) : Total or spin density gradient                     !
! ******** OUTPUT *****************************************************       !
! REAL*8  EX             : Exchange energy density                            !
! REAL*8  EC             : Correlation energy density                         !
! REAL*8  DEXDD(NSPIN)   : Partial derivative                                 !
!                           d(DensTot*Ex)/dDens(nspin),                       !
!                           where DensTot = Sum_nspin( DENS(nspin) )          !
!                           For a constant density, this is the               !
!                          exchange potential                                 !
! REAL*8  DECDD(NSPIN)   : Partial derivative                                 !
!                           d(DensTot*Ec)/dDens(nspin),                       !
!                           where DensTot = Sum_nspin( DENS(nspin) )          !
!                          For a constant density, this is the                !
!                          correlation potential                              !
! REAL*8  DEXDGD(3,NSPIN): Partial derivative                                 !
!                           d(DensTot*Ex)/d(GradDens(i,nspin))                !
! REAL*8  DECDGD(3,NSPIN): Partial derivative                                 !
!                           d(DensTot*Ec)/d(GradDens(i,nspin))                !
! ********* UNITS ****************************************************        !
! Lengths in Bohr                                                             !
! Densities in electrons per Bohr**3                                          !
! Energies in Hartrees                                                        !
! Gradient vectors in cartesian coordinates                                   !
! ********* ROUTINES CALLED ******************************************        !
! EXCHNG, PW92C                                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine xc_x_pbe(IREL, NSPIN, DENS, GDENS, EX, DEXDD, DEXDGD)

  logical, intent(in)   :: irel
  integer, intent(in)   :: NSPIN
  real(r8), intent(IN)  :: DENS(NSPIN), GDENS(3, NSPIN)
  real(r8), intent(out) :: EX, DEXDD(NSPIN), DEXDGD(3, NSPIN)
  
  ! Internal variables
  integer :: IS, IX 
  real(r8) :: &
      D(2), DF1DD, DF1DGD, DFDD, DFDGD, DFXDD(2), DFXDGD(3,2), &
      DKFDD, DS(2), DSDD, DSDGD, DT, EXUNIF, F, F1, FX,        &
      GD(3,2), GDM(2), GDMS, GDMT, GDS, GDT(3), KFS, S, VXUNIF(2)

! Lower bounds of density and its gradient to avoid divisions by zero
! plus some numerical constants
  real(r8), parameter :: &
      DENMIN = 1.E-12_r8, GDMIN  = 1.E-12_r8,          &
      HALF=0.50_r8, THD=1.0_r8/3.0_r8, THRHLF=1.50_r8, &
      TWO=2.0_r8, BETA = 0.066725_r8,                  &
      MU = BETA * M_PI**2 / 3, KAPPA = 0.8040_r8

! Translate density and its gradient to new variables
  if (NSPIN .eq. 1) then
    D(1) = HALF * DENS(1)
    D(2) = D(1)
    DT = max( DENMIN, DENS(1) )
    do IX = 1,3
      GD(IX,1) = HALF * GDENS(IX,1)
      GD(IX,2) = GD(IX,1)
      GDT(IX) = GDENS(IX,1)
    end do
  else
    D(1) = DENS(1)
    D(2) = DENS(2)
    DT = max( DENMIN, DENS(1)+DENS(2) )
    do IX = 1,3
      GD(IX,1) = GDENS(IX,1)
      GD(IX,2) = GDENS(IX,2)
      GDT(IX) = GDENS(IX,1) + GDENS(IX,2)
    end do
  endif
  GDM(1) = sqrt( GD(1,1)**2 + GD(2,1)**2 + GD(3,1)**2 )
  GDM(2) = sqrt( GD(1,2)**2 + GD(2,2)**2 + GD(3,2)**2 )
  GDMT   = sqrt( GDT(1)**2  + GDT(2)**2  + GDT(3)**2  )
  GDMT = max( GDMIN, GDMT )
  
! Find exchange energy and potential
  FX = 0
  do IS = 1,2
    DS(IS)   = max( DENMIN, 2 * D(IS) )
    GDMS = max( GDMIN, 2 * GDM(IS) )
    KFS = (3 * M_PI**2 * DS(IS))**THD
    S = GDMS / (2 * KFS * DS(IS))
    F1 = 1 + MU * S**2 / KAPPA
    F = 1 + KAPPA - KAPPA / F1
    
    !       Note nspin=1 in call to exchng...
    
    call xc_x_lda( IREL, 1_i4, DS(is:is), VXUNIF(is:is), EXUNIF )
    FX = FX + DS(IS) * EXUNIF * F
    
    DKFDD = THD * KFS / DS(IS)
    DSDD = S * ( -(DKFDD/KFS) - 1/DS(IS) )
    DF1DD = 2 * (F1-1) * DSDD / S
    DFDD = KAPPA * DF1DD / F1**2
    DFXDD(IS) = VXUNIF(IS) * F + DS(IS) * EXUNIF * DFDD
    
    do IX = 1,3
      GDS = 2 * GD(IX,IS)
      DSDGD = (S / GDMS) * GDS / GDMS
      DF1DGD = 2 * MU * S * DSDGD / KAPPA
      DFDGD = KAPPA * DF1DGD / F1**2
      DFXDGD(IX,IS) = DS(IS) * EXUNIF * DFDGD
    end do
  end do
  FX = HALF * FX / DT

  ! Set output arguments
  EX = FX
  do IS = 1, NSPIN
    DEXDD(IS) = DFXDD(IS)
    do IX = 1,3
      DEXDGD(IX,IS) = DFXDGD(IX,IS)
    end do
  end do

  return
end subroutine xc_x_pbe

subroutine xc_c_pbe(NSPIN, DENS, GDENS, EC, DECDD, DECDGD)

  integer, intent(in)   :: NSPIN
  real(r8), intent(IN)  :: DENS(NSPIN), GDENS(3, NSPIN)
  real(r8), intent(out) :: EC, DECDD(NSPIN), DECDGD(3, NSPIN)

! Internal variables
  integer :: IS, IX
  real(r8) :: &
      A, D(2), DADD, DECUDD, DF1DD, DF2DD, DF3DD, DF4DD,  &
      DF3DGD, DF4DGD, DFCDD(2), DFCDGD(3,2), DHDD, DHDGD, &
      DKFDD, DKSDD, DPDD, DPDZ, DRSDD, DT, DTDD, DTDGD,   &
      DZDD(2), ECUNIF, F1, F2, F3, F4, FC, GAMMA, GD(3,2),&
      GDM(2), GDMT, GDT(3), H, KF, KS, PHI, RS, T,        &
      VCUNIF(2), ZETA

! Lower bounds of density and its gradient to avoid divisions by zero

  real(r8), parameter :: &
      DENMIN = 1.E-12, GDMIN  = 1.E-12,        &
      FOUTHD=4.0_r8/3.0_r8, HALF=0.50_r8,      &
      THD=1.0_r8/3.0_r8, TWOTHD=2.0_r8/3.0_r8, &
      BETA = 0.066725_r8
  
  GAMMA = (1 - log(2.0_r8)) / M_PI**2

! Translate density and its gradient to new variables
  if (NSPIN .eq. 1) then
    D(1) = HALF * DENS(1)
    D(2) = D(1)
    DT = max( DENMIN, DENS(1) )
    do IX = 1,3
      GD(IX,1) = HALF * GDENS(IX,1)
      GD(IX,2) = GD(IX,1)
      GDT(IX) = GDENS(IX,1)
    end do
  else
    D(1) = DENS(1)
    D(2) = DENS(2)
    DT = max( DENMIN, DENS(1)+DENS(2) )
    do IX = 1,3
      GD(IX,1) = GDENS(IX,1)
      GD(IX,2) = GDENS(IX,2)
      GDT(IX) = GDENS(IX,1) + GDENS(IX,2)
    end do
  endif
  GDM(1) = sqrt( GD(1,1)**2 + GD(2,1)**2 + GD(3,1)**2 )
  GDM(2) = sqrt( GD(1,2)**2 + GD(2,2)**2 + GD(3,2)**2 )
  GDMT   = sqrt( GDT(1)**2  + GDT(2)**2  + GDT(3)**2  )
  GDMT = max( GDMIN, GDMT )

! Find local correlation energy and potential
  call xc_c_pw92(2, D, VCUNIF, ECUNIF)

! Find total correlation energy
  RS = ( 3 / (4*M_PI*DT) )**THD
  KF = (3 * M_PI**2 * DT)**THD
  KS = sqrt( 4 * KF / M_PI )
  ZETA = ( D(1) - D(2) ) / DT
  ZETA = max( -1.0_r8+DENMIN, ZETA )
  ZETA = min(  1.0_r8-DENMIN, ZETA )
  PHI = HALF * ( (1+ZETA)**TWOTHD + (1-ZETA)**TWOTHD )
  T = GDMT / (2 * PHI * KS * DT)
  F1 = ECUNIF / GAMMA / PHI**3
  F2 = exp(-F1)
  A = BETA / GAMMA / (F2-1)
  F3 = T**2 + A * T**4
  F4 = BETA/GAMMA * F3 / (1 + A*F3)
  H = GAMMA * PHI**3 * log( 1 + F4 )
  FC = ECUNIF + H


! Find correlation energy derivatives
  DRSDD = - (THD * RS / DT)
  DKFDD =   THD * KF / DT
  DKSDD = HALF * KS * DKFDD / KF
  DZDD(1) =   1 / DT - ZETA / DT
  DZDD(2) = - (1 / DT) - ZETA / DT
  DPDZ = HALF * TWOTHD * ( 1/(1+ZETA)**THD - 1/(1-ZETA)**THD )
  do IS = 1,2
     DECUDD = ( VCUNIF(IS) - ECUNIF ) / DT
     DPDD = DPDZ * DZDD(IS)
     DTDD = (- T) * ( DPDD/PHI + DKSDD/KS + 1/DT )
     DF1DD = F1 * ( DECUDD/ECUNIF - 3*DPDD/PHI )
     DF2DD = (- F2) * DF1DD
     DADD = (- A) * DF2DD / (F2-1)
     DF3DD = (2*T + 4*A*T**3) * DTDD + DADD * T**4
     DF4DD = F4 * ( DF3DD/F3 - (DADD*F3+A*DF3DD)/(1+A*F3) )
     DHDD = 3 * H * DPDD / PHI
     DHDD = DHDD + GAMMA * PHI**3 * DF4DD / (1+F4)
     DFCDD(IS) = VCUNIF(IS) + H + DT * DHDD

     do IX = 1,3
        DTDGD = (T / GDMT) * GDT(IX) / GDMT
        DF3DGD = DTDGD * ( 2 * T + 4 * A * T**3 )
        DF4DGD = F4 * DF3DGD * ( 1/F3 - A/(1+A*F3) ) 
        DHDGD = GAMMA * PHI**3 * DF4DGD / (1+F4)
        DFCDGD(IX,IS) = DT * DHDGD
      end do
    end do

! Set output arguments
  EC = FC
  do IS = 1,NSPIN
    DECDD(IS) = DFCDD(IS)
    do IX = 1,3
      DECDGD(IX,IS) = DFCDGD(IX,IS)
    end do
  end do
  
  return
end subroutine xc_c_pbe

subroutine xc_x_lb94(nspin, dens, gdens, ex, dexdd, dexdgd)
  integer, intent(in)  :: nspin
  real(r8), intent(IN) :: dens(nspin), gdens(3, nspin)
  real(r8), intent(out):: ex, dexdd(nspin), dexdgd(3, nspin)
  
  ! Internal variables
  integer(i4) :: is 
  real(r8)    :: d(nspin), gd(3, nspin), gdm, x, f

! Lower bounds of density and its gradient to avoid divisions by zero
! plus some numerical constants
  real(r8), parameter :: &
      DENMIN = 1.E-20_r8,    GDMIN  = 1.E-20_r8,          &
      HALF   = 0.5_r8,      THRD   = 1._r8/3._r8,        &
      FTHRD  = 4._r8/3._r8, BETA   = 0.05_r8

! first we add the LDA potential
  call xc_x_lda( .false., nspin, dens, dexdd, ex)

! Translate density and its gradient to new variables
  if (nspin .eq. 1) then
    d(1)     = dens(1) * HALF
    gd(:, 1) = gdens(:, 1) * HALF
  else
    d  = dens
    gd = gdens
  endif

  do is = 1, nspin
    gdm   = sqrt( gd(1,is)**2 + gd(2,is)**2 + gd(3,is)**2 )
    
    if(d(is) >= DENMIN .and. gdm >=GDMIN) then
      x = gdm / d(is)**FTHRD
      f = x**2/(1._r8 + 3._r8*BETA*x*oct_asinh(x))
      dexdd(is) = dexdd(is) - BETA * d(is)**THRD * f
    end if
  end do

  ex    = 0._r8 ! this is to be calculated afterwards
  dexdgd = 0._r8

  return
end subroutine xc_x_lb94
  
