#include "config_F90.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module includes a set of subroutines needed by (and only by) init_kb   !
! it should be rewritten (implicit none does not work!)                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  module kb

  use global

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine to generate a smooth local pseudopotential,                      !
! and its charge density.                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine vps_loc(zval,nrval,a,b,rofi,drdi,s,rcmax,                 &
                          rclocal,vlocal,nchloc,chlocal)
!
! subroutine to generate a smooth local pseudopotential and its charge density.
!

       implicit none
!
! Arguments
!
       integer nrval, nchloc
       real(r8) a, b, rclocal, rcmax, zval, rofi(nrval),vlocal(nrval),        &
                chlocal(nrval),drdi(nrval),s(nrval) 
!
! Local variables
!
       integer ir, nrc
       real(r8), parameter :: exp_range=120.d0 
       real(r8) pi, qtot, vconst
       real(r8) dev1, dev2, dev3, var1, var2, var3, v1, v2, v3, v4,         &
                dm11, dm12,dm13,dm21,dm22,dm23,dm31,dm32,dm33,              &
                g0, g1, g2, g3, g4, g5
!
       pi = 4.0d0*datan(1.0d0)
!
! We fit either an analytic function for the potential 
!
        write(6,'(/,2a)')'vps_local:',                                        &
                         ' Fitted the Vlocal to v4*exp(v1*r^2+v2*r^3+v3*r^4)'
!
         nrc=nint(dlog(rcmax/b+1.0d0)/a) + 1 
!
         vconst = vlocal(nchloc+1)
!
! local potential is Vloc(r)=v4*exp(v1*r^2+v2*r^4+v3*r^6) inside rcmax and equal
! to the all electron outside
!                             fit upto third derivative at rcmax
!
         dev1 = ( - vlocal(nrc+2)/12.0d0                        &
                  + 2.0d0/3.0d0*vlocal(nrc+1)                   &
                  - 2.0d0/3.0d0*vlocal(nrc-1)                   &
                  + vlocal(nrc-2)/12.0d0      )
         dev2 = ( - vlocal(nrc+2)/12.0d0                        &
                  + 4.0d0/3.0d0*vlocal(nrc+1)                   &
                  - 5.0d0/2.0d0*vlocal(nrc)                     &
                  + 4.0d0/3.0d0*vlocal(nrc-1)                   &
                  - vlocal(nrc-2)/12.0d0      )
         dev3 = ( - vlocal(nrc+3)/8.0d0                         &
                  + vlocal(nrc+2)                               &
                  - 13.0d0/8.0d0*vlocal(nrc+1)                  &
                  + 13.0d0/8.0d0*vlocal(nrc-1)                  &
                  - vlocal(nrc-2)                               &
                  + vlocal(nrc-3)/8.0d0       )
!
         dev3=(dev3-3.0d0*a*dev2+2.0d0*a**2*dev1)/(drdi(nrc)**3)
         dev2=(dev2 - a*dev1)/(drdi(nrc)**2)
         dev1=dev1/drdi(nrc)
!
         var1 = dev1/vlocal(nrc)
         var2 = dev2/vlocal(nrc) - var1**2
         var3 = dev3/vlocal(nrc) - var1**3 - 2.0d0*var2*var1
!
         v3 = (var3*rofi(nrc)/3.0d0-var2+var1/rofi(nrc))              &
              /(16.0d0*rofi(nrc)**4)
         v2 = (var2-var1/rofi(nrc)-24.0d0*rofi(nrc)**4*v3)            &
              /(8.0d0*rofi(nrc)**2)
         v1 = (var1-4.0d0*v2*rofi(nrc)**3-6.0d0*v3*rofi(nrc)**5)      &
              /(2.0d0*rofi(nrc))
         v4=vlocal(nrc)/(exp((v1+v2*rofi(nrc)**2+v3*rofi(nrc)**4)*rofi(nrc)**2))
!
         vlocal(1:nrc) = v4*exp((v1+v2*rofi(1:nrc)**2+v3*rofi(1:nrc)**4)*     &
                                 rofi(1:nrc)**2)
!
! Now we get the local pseudopotential charge
!
          qtot=0.d0
          chlocal = 0.0d0
          do ir=1,nchloc
            if(ir > nrc) then
               dev1 = ( - vlocal(ir+2)/12.0d0                        &
                        + 2.0d0/3.0d0*vlocal(ir+1)                   &
                        - 2.0d0/3.0d0*vlocal(ir-1)                   &
                        + vlocal(ir-2)/12.0d0      )
               dev2 = ( - vlocal(ir+2)/12.0d0                        &
                        + 4.0d0/3.0d0*vlocal(ir+1)                   &
                        - 5.0d0/2.0d0*vlocal(ir)                     &
                        + 4.0d0/3.0d0*vlocal(ir-1)                   &
                        - vlocal(ir-2)/12.0d0      )
               dev2 = (dev2 - a*dev1)/(drdi(ir)**2)
               dev1 = dev1/drdi(ir)
               g3 = (2.0d0*dev1/rofi(ir) + dev2)
            else
              g0=v4*exp((v1+v2*rofi(ir)**2+v3*rofi(ir)**4)*rofi(ir)**2)
              g1=(2.0d0*v1+4.0d0*v2*rofi(ir)**2+6.0d0*v3*rofi(ir)**4)
              g2=(2.0d0*v1+12.0d0*v2*rofi(ir)**2+30.0d0*v3*rofi(ir)**4)
              g3=(g2+g1*g1*rofi(ir)**2+2.0d0*g1)*g0
            endif
            chlocal(ir)= (-g3)/ (8.0d0*pi)
            !write(8,*)rofi(ir),chlocal(ir)
         enddo
            dev1 = ( - chlocal(nrc+2)/12.0d0                        &
                     + 2.0d0/3.0d0*chlocal(nrc+1)                   &
                     - 2.0d0/3.0d0*chlocal(nrc-1)                   &
                     + chlocal(nrc-2)/12.0d0      )
            dev1=dev1/drdi(nrc)
!
! now we make the density going to zero at nchloc (also the 1st and 2nd 
! derivatives) with a function g1*exp(-g0x)*(1+g2*x+g3*x^2+g4*x^3)
!
            g3 = 3.5d0/rofi(nchloc)**2 
            g2 = -rofi(nchloc)*g3
            g4 = -g3/(3.50d0*rofi(nchloc))
! 
            g0 = (- dev1/chlocal(nrc) +                                     &
                 (g2+2.0d0*g3*rofi(nrc)+3.0d0*g4*rofi(nrc)**3)/           &
                 (1.0d0+g2*rofi(nrc)+g3*rofi(nrc)**2+g4*rofi(nrc)**4) ) 
            g1 = chlocal(nrc)/( exp(-g0*rofi(nrc)) *                      &
                 (1.0d0+g2*rofi(nrc)+g3*rofi(nrc)**2+g4*rofi(nrc)**4) ) 
!
            qtot = 0.0d0
            do ir =1, nrc
              qtot = qtot  - (chlocal(ir)*4.0d0*pi)*rofi(ir)**2*drdi(ir)
            enddo
            do ir = nrc+1,nchloc
             chlocal(ir)=g1*(1.0d0+g2*rofi(ir)+g3*rofi(ir)**2+g4*rofi(ir)**4)* &
                             exp(-g0*rofi(ir))
             qtot = qtot  - (chlocal(ir)*4.0d0*pi)*rofi(ir)**2*drdi(ir)
            enddo
!
          chlocal(:) = zval*chlocal(:)/qtot
!
          qtot = chlocal(1)
          chlocal(1)=chlocal(2)-(chlocal(3)-chlocal(2))*      &
                     rofi(2)/(rofi(3)-rofi(2))
          chlocal(2:nrval) = 4.0d0*pi*rofi(2:nrval)**2*chlocal(2:nrval)
!cang          call vhrtre(chlocal,vlocal,rofi,drdi,s,nrval,a)
          chlocal(2:nrval) = chlocal(2:nrval) / (4.0d0*pi*rofi(2:nrval)**2)
          chlocal(1) = qtot 
!          


  return


  end subroutine vps_loc






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   VHRTRE CONSTRUCTS THE ELECTROSTATIC POTENTIAL DUE TO A SUPPLIED           !
!   ELECTRON DENSITY.  THE NUMEROV METHOD IS USED TO INTEGRATE                !
!   POISSONS EQN.                                                            !
!                                                                             !
!   DESCRIPTION OF ARGUMENTS:                                                 !
!      RHO....4*PI*R**2 * THE ELECTRON DENSITY FOR WHICH WE CALCULATING       !
!             THE ELECTROSTATIC POTENTIAL                                     !
!      V......THE ELECTROSTATIC POTENTIAL DUE TO THE ELECTRON DENSITY         !
!             RHO.  THE CONSTANTS OF INTEGRATION ARE FIXED SO THAT THE        !
!             POTENTIAL TENDS TO A CONSTANT AT THE ORIGIN AND TO              !
!             2*Q/R AT R=R(NR), WHERE Q IS THE INTEGRATED CHARGE              !
!             CONTAINED IN RHO(R)                                             !
!      R......THE RADIAL MESH R(I) = B*(EXP(A(I-1))-1)                        !
!      NR.....THE NUMBER OF RADIAL MESH POINTS                                !
!      DRDI...DR(I)/DI                                                        !
!     SRDRDI.SQRT(DR/DI)                                                      !
!      A......THE PARAMETER APPEARING IN R(I) = B*(EXP(A(I-1))-1)             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine vhrtre(rho, v, r, drdi, srdrdi, nr, a)
  real(r8), intent(IN), dimension(*) :: rho, r, drdi, srdrdi
  real(r8), intent(in) :: a
  real(r8), intent(out), dimension(*) :: v
  integer, intent(in) :: nr

  real(r8) :: x,y,q,a2by4,ybyq,qbyy,qpartc,v0,qt,dz,t,beta,dv
  integer :: nrm1,nrm2,ir

  NRM1 = NR - 1
  NRM2 = NR - 2
  A2BY4 = A*A/4.0_r8
  YBYQ = 1.0_r8 - A*A/48.0_r8
  QBYY = 1.0_r8/YBYQ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  SIMPSONS RULE IS USED TO PERFORM TWO INTEGRALS OVER THE ELECTRON          ! 
!  DENSITY.  THE TOTAL CHARGE QT IS USED TO FIX THE POTENTIAL AT R=R(NR)      !
!  AND V0 (THE INTEGRAL OF THE ELECTRON DENSITY DIVIDED BY R) FIXES           !
!  THE ELECTROSTATIC POTENTIAL AT THE ORIGIN                                  !
!  Modified to be consistent with pseudopotential generation (use the         !
!  trapeziodal rule for integration). A. Rubio. Jan. 2000                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  V0 = 0.0_r8
  QT = 0.0_r8
  do IR = 2, NRM1, 2
    DZ = DRDI(IR)*RHO(IR)
    QT = QT + DZ
    V0 = V0 + DZ/R(IR)
  end do

  do IR=3, NRM2, 2
    DZ = DRDI(IR)*RHO(IR)
    QT = QT + DZ
    V0 = V0 + DZ/R(IR)
  end do
  DZ = DRDI(NR)*RHO(NR)

  QT = (QT + QT + DZ)*0.50d0
  V0 = (V0 + V0 + DZ/R(NR))
  V(1) = 2.0_r8*V0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  THE ELECTROSTATIC POTENTIAL AT R=0 IS SET EQUAL TO                         !
!                       THE AVERAGE VALUE OF RHO(R)/R                         !
!  BEGIN CONSTRUCTION OF THE POTENTIAL AT FINITE                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IR   = 2
  T    = SRDRDI(IR)/R(IR)
  BETA = DRDI(IR)*T*RHO(IR)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  THE NEXT 4 STATEMENTS INDICATE THAT WE FIRST FIND THE PARTICULAR           !
!  SOLUTION TO THE INHOMOGENEOUS EQN. FOR WHICH Q(2)=0, WE THEN               !
!  ADD TO THIS PARTICULAR SOLUTION A SOLUTION OF THE HOMOGENEOUS EQN.         !
!  (A CONSTANT IN V OR A Q PROPORTIONAL TO R)                                 !
!  WHICH WHEN DIVIDED BY R IN GOING FROM Q TO V GIVES                         !
!  THE POTENTIAL THE DESIRED COULOMB TAIL OUTSIDE THE ELECTRON DENSITY.       !
!  THE SIGNIFICANCE OF THE SOLUTION VANISHING AT THE SECOND RADIAL            !
!  MESH POINT IS THAT, SINCE ALL REGULAR SOLUTIONS OF THE EQUATION            !
!  FOR Q=R*V VANISH AT THE ORIGIN, THE KNOWLEDGE OF THE SOLUTION              !
!  VALUE AT THE SECOND MESH POINT PROVIDES THE TWO SOLUTION VALUES            !
!  REQUIRED TO START THE NUMEROV PROCEDURE.                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  X = 0.0_r8
  Y = 0.0_r8
  Q = (Y - BETA/12.0_r8)*QBYY
  V(IR) = 2.0_r8*T*Q

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  BEGINNING OF THE NUMEROV ALGORITHM                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

3 X  = X + A2BY4*Q - BETA
  Y  = Y + X
  IR = IR + 1
  T  = SRDRDI(IR)/R(IR)
  BETA = T*DRDI(IR)*RHO(IR)
  Q = (Y-BETA/12.0_r8)*QBYY
  V(IR) = 2.0_r8*T*Q
  IF(IR.LT.NR) GO TO 3


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  END OF THE NUMEROV ALGORITHM                                               !
!                                                                             !
!  WE HAVE NOW FOUND A PARTICULAR SOLUTION TO THE INHOMOGENEOUS EQN.          !
!  FOR WHICH Q(R) AT THE SECOND RADIAL MESH POINT EQUALS ZERO.                !
!  NOTE THAT ALL REGULAR SOLUTIONS TO THE EQUATION FOR Q=R*V                  !
!  VANISH AT THE ORIGIN.                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  QPARTC = R(NR)*V(NR)/2.0_r8
  DZ = QT - QPARTC
  DV = 2.0_r8*DZ/R(NR)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  THE LOOP FOLLOWING ADDS THE CONSTANT SOLUTION OF THE HOMOGENEOUSi          !
!  EQN TO THE PARTICULAR SOLUTION OF THE INHOMOGENEOUS EQN.                   !
!  NOTE THAT V(1) IS CONSTRUCTED INDEPENDENTLY                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do IR = 2, NR
    V(IR) = V(IR) + DV
  end do

  return
  end subroutine vhrtre






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Next comes a set of subroutines to solve the readial schrodinger equation   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






  subroutine egofv(h,s,n,e,g,y,l,z,a,b,rmax,nprin,nnode,dr)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  egofv determines the eigenenergy and wavefunction corresponding            !
!  to a particular l, principal quantum number and boundary condition.        !
!                                                                             !
!  two fundamental techniques are used to locate the solution:                ! 
!      1) node counting and bisection                                         !
!       2) variational estimate based on a slope discontinuity in psi         !
!  the arguments are defined as follows:                                      !
!       h,s: g = (h-e*s)*g                                                   !
!       nr: maximum allowed number of radial points                           !
!       e: e(i) is the i-th energy found                                      !
!       ne: number of energies found                                          !
!       l: the angular momentum                                               !
!       ncor: the number of lower-energy state                                !
!                                                                             ! 
!  the individual energies are resolved by performing a fixed number          !
!  of bisections after a given eigenvalue has been isolated                   ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  implicit real(r8) (a-h,o-z)
  integer :: i,n,l,nprin,nnode,ncor,n1,n2,niter,nt

  real(r8),dimension(*) ::h,s,g,y

  data tol   /1.d-5/
  ncor=nprin-l-1
  n1=nnode
  n2=nnode-1
  e1=e
  e2=e


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  the labels 1 and 2 refer to the bisection process, defining the            !
!  range in which the desired solution is located.  the initial               !
!  settings of n1, n2, e1 and e2 are not consistent with the bisection        !
!  algorithm; they are set to consistent values when the desired              !
!  energy interval has been located.                                          !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  del = 5.d-1
  de  = 0.d0
  niter = 0
1 niter = niter + 1
  if(niter.gt.40) go to 3
  et = e + de
! the following line is the fundamental "bisection"
  e = 0.5*(e1+e2)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  the following concatenation of logical ors ensures that node               !
!  counting is used unless the previous integration of the radial             !
!  eq produced both the correct number of nodes and a sensible                !
!  prediction for the energy.                                                 !
!                                                                             !
!     sensible means that et must be greater than e1 and less than e2         !
!     correct number of nodes means that nt = nnode or nnode-1.               !
!                                                                             !
!     leaving e set to its bisection value, and transfering to                !
!     the call to yofe means that we are performing bisection,                !
!     whereas setting e to et is use of the variational estimate.             !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  if(et.le.e1 .or. et.ge.e2 .or.                                              &
     nt.lt.nnode-1 .or. nt.gt.nnode) go to 2
  e=et
  if(dabs(de).lt.tol) go to 6
2 call yofe(e,de,dr,rmax,h,s,y,n,l,ncor,nt,z,a,b)
!     write(6,101) l,dr,n1,nt,nnode,n2,e1,e,e2,de
!101  format('  l     dr     n1  nt   n  n2       e1           e',
!    1       '          e2          de'/i3,d10.3,4i4,4f12.5)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  yofe integrates the schro eq.; now the bisection logic                     !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  if(nt.ge.nnode) go to 5
!  too few nodes; set e1 and n1
  e1=e
  n1=nt


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  at this point, we have just set the bottom of the bisection range;         !
!  if the top is also set, we procede.  if the top of the range has not       !
!  been set, it means that we have yet to find an e greater than the          !
!  desired energy.  the upper end of the range is extended.                   !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  if(n2.ge.nnode) go to 1
  del=del*2.d0
  e2=e1+del
  go to 1
!  too many nodes; set e2 and n2
5 e2=e
  n2=nt


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  at this point, we have just set the top of the bisection range;            !
!  if the top is also set, we procede.  if the top of the range has           !
!  not been set, it means that we have yet to find an e less than the         !
!  desired energy.  the lower end of the range is extended.                   !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  if(n1.lt.nnode) go to 1
  del=del*2.d0
  e1=e2-del
  go to 1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  the numerov method uses a transformation of the radial wave fcn.           !
!  that we call "y".  having located the eigenenergy, we transform            !
!  y to "g", from which the density is easily constructed.                    !
!  finally, the call to "nrmlzg" normalizes g to one electron.                !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


6 g(1) = 0.d0
  do 7 i=2,n
  t=h(i)-e*s(i)
  g(i)=y(i)/(1.d0-t/12.d0)
7 continue
  call nrmlzg(g,s,n)
  return
3 write(6,4) z,l,nnode,e,de
4 format(' egofv: too many iterations; execution stopping'/                   &
            ' z=',f3.0,'  l=',i2,'  nnode=',i2,'  e=',f12.5,                  &
            '  de=',f12.5)
  stop 8


  end subroutine egofv








  subroutine yofe(e,de,dr,rmax,h,s,y,nmax,l,ncor,nnode,z,a,b)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   yofe integrates the radial schrodinger eqn using the numerov              !
!   method.                                                                   !
!                                                                             !
!   the arguments are defined as follows:                                     !
!       e is the old energy(overwritten) by the new energy                    !
!      de is the e change predicted to elim the kink in psi                   !
!       dr is the log deriv (the boundary condition)                          !
!       gpp = (h-es)g (all diagonal in i (radius) )                            !
!       y is the numerov independent variable y = g - gpp/12                   !
!       n is the number of radial mesh points                                 !
!       l is the angular momentum                                             !
!       ncor is the number of states of lower energy                          !
!       nnode is 1 + the number of interior nodes in psi                      !
!       z is the atomic number                                                !
!       a and b specify the radial mesh r(i)=(exp(a*(i-1))-1)*b               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  implicit real(r8) (a-h,o-z)
  dimension h(*),s(*),y(*)
  integer :: nmax,l,ncor,nnode,n,knk,nndin,i

  zdr = z*a*b
  n=nmax
8 if( h(n)-e*s(n) .lt. 1.d0 ) go to 9
  y(n)=0.d0
  n=n-1
  go to 8
9 continue
  call bcorgn(e,h,s,l,zdr,y2)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  bcorgn computes y2, which embodies the boundary condition                  !
!  satisfied by the radial wave function at the origin                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  knk=n
  call numout(e,h,s,y,ncor,knk,nnode,y2,g,gsg,x)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  the outward integration is now complete                                    !
!                                                                             !
!     we first decide if the kinetic energy is sufficiently non               !
!     negative to permit use of the numerov eq at rmax.  if                   !
!     it is not, then zero-value boundary condition is used                   !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  yn=0.d0
  if(n.lt.nmax .or. dabs(dr).gt.1.d3) go to 7
  call bcrmax(e,dr,rmax,h,s,n,yn,a,b)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  bcrmax computes yn, which embodies the boundary condition                  !
!  satisfied by the radial wave function at rmax                              !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


7 call numin(e,h,s,y,n,nndin,yn,gin,gsgin,xin,knk)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  numin performs the inward integration by the numerov method                !
!                                                                             !
!  the energy increment is now evaluated from the kink in psi                 !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ratio = g/gin
  xin=xin*ratio
  gsg=gsg+gsgin*ratio*ratio
  t=h(knk)-e*s(knk)
  de=g*(x+xin+t*g)/gsg
  nnode=nnode+nndin
  if(de.lt.0.d0) nnode=nnode+1
  do 6 i=knk,n
  y(i) = y(i)*ratio
6 continue
  return


  end subroutine yofe




  subroutine nrmlzg(g,s,n)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   nrmlzg normalizes the supplied radial wave function                       !
!                                                                             !
!   the arguments are defined as follows:                                     !
!       g is the radial wave function appropriate to the numerov              !
!             representation of the radial schrodinger equation               !
!             that is, the radial fcn r(r) = (drdi)**1/2 g(i) / r(i)          !
!       gpp = (h-es)g (all diagonal in i (radius) )                           !
!       n1 is the number of radial mesh points corresponding to               !
!             the portion of the radial mesh on which the norm                !
!             is defined                                                      !
!       n2 is the number of radial mesh points corresponding to               !
!             the portion of the radial mesh on which the wave                !
!             function is defined.  for the intended use of this              !
!             routine, n1 = nrval and n2 = nrcor                              !
!       a and b are the radial mesh parameters                                !
!             r(i) = ( exp(a*(i-1)) - 1 ) * b                                 !
!             (dr/di = a*b at the origin)                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  implicit real(r8) (a-h,o-z)
  real(r8) s(*),g(*),norm
  integer :: n,nm1,nm2,i


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  determine the norm of g(i) using simpsons rule                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  if (mod(n,2).ne.1) write(6,*) ' nrmlzg: n should be odd. n =',n
  norm = 0.d0
  nm1 = n - 1
  do 2 i = 2,nm1,2
  norm=norm + g(i)*s(i)*g(i)
2 continue
  norm = norm * 2.d0
  nm2  = n - 2
  do 3 i = 3,nm2,2
  norm=norm + g(i)*s(i)*g(i)
3 continue
  norm = norm * 2.d0
  nm2  = n - 2
  do 4 i = 1,n,nm1
  norm=norm + g(i)*s(i)*g(i)
4 continue
  norm = norm/3.d0
  srnrm = dsqrt(norm)
  do 5 i=1,n
  g(i) = g(i)/srnrm
5 continue
  return


  end subroutine nrmlzg






  subroutine bcorgn(e,h,s,l,zdr,y2)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   yofe integrates the radial schrodinger eqn using the numerov              !
!   method.                                                                   !
!   the arguments are defined as follows:                                     !
!       e is the old energy(overwritten) by the new energy                    !
!       de is the e change predicted to elim the kink in psi                  !
!       dr is the log deriv (the boundary condition)                          !
!       gpp = (h-es)g (all diagonal in i (radius) )                           !
!       y is the numerov independent variable y = g - gpp/12                  !
!       n is the number of radial mesh points                                 !
!       l is the angular momentum                                             !
!       nnode is 1 + the number of interior nodes in psi                      !
!       z is the atomic number                                                !
!       a and b specify the radial mesh r(i)=(exp(a*(i-1))-1)*b               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  implicit real(r8) (a-h,o-z)
  real(r8) h(*),s(*)
  integer :: l


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   the quantity called d(i) in the program is actually the inverse           !
!   of the diagonal of the tri-diagonal numerov matrix                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  t2=h(2)-e*s(2)
  d2=-((24.d0+10.d0*t2)/(12.d0-t2))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  the following section deals with the fact that the independent             !
!  variable "y" in the numerov equation is not zero at the origin             !
!  for l less than 2                                                          !
!  the l=0 solution g vanishes, but the first and second                      !
!  derivatives are finite, making the numerov variable y finite               !
!  the l=1 solution g vanishes, and gprime also vanishes, but                 !
!  the second derivative gpp is finite making y finite.  for l > 1,           !
!  g and its first two derivatives vanish, making y zero.                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  if(l.ge.2) goto 3
  if(l.gt.0) goto 1
  c0=zdr/6.d0
  c0=c0/(1.d0-0.75*zdr)
  go to 2
1 c0=1.d0/12.d0
  c0=(-c0)*8.d0/3.d0
2 c1=c0*(12.d0+13.d0*t2)/(12.d0-t2)
  t3=h(3)-e*s(3)
  c2=(-5.d-1)*c0*(24.d0-t3)/(12.d0-t3)
  d2=(d2-c1)/(1.d0-c2)
3 y2=(-1.d0)/d2
  return

  end subroutine bcorgn



  subroutine bcrmax(e,dr,rmax,h,s,n,yn,a,b)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 22.7.85                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  implicit real(r8) (a-h,o-z)
  real(r8)  h(*),s(*),                                                          &
   e,dr,rmax,yn,a,b,tnm1,tn,tnp1,beta,dg,c1,c2,c3,dn
  integer :: n

!     write(6,*) 'bcrmax:',dr
  tnm1=h(n-1)-e*s(n-1)
  tn  =h(n  )-e*s(n  )
  tnp1=h(n+1)-e*s(n+1)
  beta=1.d0+b/rmax
  dg=a*beta*(dr+1.d0-5.d-1/beta)


  c2=24.d0*dg/(12.d0-tn)
  dn=-((24.d0+10.d0*tn)/(12.d0-tn))

  c1= (1.d0-tnm1/6.d0)/(1.d0-tnm1/12.d0)
  c3=-((1.d0-tnp1/6.d0)/(1.d0-tnp1/12.d0))
  yn=-((1.d0-c1/c3)/(dn-c2/c3))


  return

  end subroutine bcrmax



  subroutine numin(e,h,s,y,n,nnode,yn,g,gsg,x,knk)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   yofe integrates the radial schrodinger eqn using the numerov              !
!   method.                                                                   !
!                                                                             !
!   the arguments are defined as follows:                                     !
!       e is the old energy(overwritten) by the new energy                    !
!       de is the e change predicted to elim the kink in psi                  !
!       dr is the log deriv (the boundary condition)                          !
!       gpp = (h-es)g (all diagonal in i (radius) )                            !
!       y is the numerov independent variable y = g - gpp/12                   !
!       n is the number of radial mesh points                                 !
!       l is the angular momentum                                             !
!       nnode is 1 + the number of interior nodes in psi                      !
!       z is the atomic number                                                !
!       a and b specify the radial mesh r(i)=(exp(a*(i-1))-1)*b               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  implicit real(r8) (a-h,o-z)
  real(r8) h(n),s(n),y(n)
  integer :: i,n,nnode,knk

  y(n)=yn
  t=h(n)-e*s(n)
  g=y(n)/(1.d0-t/12.d0)
  gsg=g*s(n)*g
  i=n-1
  y(i)=1.d0
  t=h(i)-e*s(i)
  g=y(i)/(1.d0-t/12.d0)
  gsg=gsg+g*s(i)*g
  x=y(i)-y(n)
  nnode=0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  begin the inward integrationby the numerov method                          !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


1 x=x+t*g
  i=i-1
  y(i)=y(i+1)+x
  if( y(i)*y(i+1) .lt. 0.d0) nnode=nnode+1
  t=h(i)-e*s(i)
  g=y(i)/(1.d0-t/12.d0)
  gsg=gsg+g*s(i)*g
  if(i.gt.knk) go to 1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  the last statement defines the kink radius as the point where              !
!  psi first turns downward.  this usually means at the outermost             !
!  maximum                                                                    !
!                                                                             !
!  the inward integration is now complete                                     !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  knk=i
  return

  end subroutine numin




  subroutine numout(e,h,s,y,ncor,knk,nnode,y2,g,gsg,x)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   yofe integrates the radial schrodinger eqn using the numerov              !
!   method.                                                                   !
!                                                                             ! 
!   the arguments are defined as follows:                                     !
!       e is the old energy(overwritten) by the new energy                    !
!       de is the e change predicted to elim the kink in psi                  !
!       dr is the log deriv (the boundary condition)                          !
!       gpp = (h-es)g (all diagonal in i (radius) )                            !
!       y is the numerov independent variable y = g - g/12                   !
!       n is the number of radial mesh points                                 !
!       l is the angular momentum                                             !
!       nnode is 1 + the number of interior nodes in psi                      !
!       z is the atomic number                                                !
!       a and b specify the radial mesh r(i)=(exp(a*(i-1))-1)*b               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  implicit real(r8) (a-h,o-z)
  real(r8) h(knk),s(knk),y(knk)
  integer :: ncor,nnode,knk,i,nm4

  y(1)=0.d0
  y(2)=y2
  t=h(2)-e*s(2)
  g=y(2)/(1.d0-t/12.d0)
  gsg=g*s(2)*g
  y(3)=1.d0
  t=h(3)-e*s(3)
  g=y(3)/(1.d0-t/12.d0)
  gsg=gsg+g*s(3)*g
  x=y(3)-y(2)
  i=3
  nnode=0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  begin the outward integrationby the numerov method                         !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  nm4=knk-4
1 xl=x
  x=x+t*g
  i=i+1
  y(i)=y(i-1)+x
!     write(6,300) i,y(i),x,t,h(i),s(i)
!300  format(i5,5d14.5)
  if( y(i)*y(i-1) .lt. 0.d0) nnode=nnode+1
  t=h(i)-e*s(i)
  g=y(i)/(1.d0-t/12.d0)
  gsg=gsg+g*s(i)*g
  if(i.eq.nm4) go to 2
  if(nnode.lt.ncor) go to 1
  if(xl*x.gt.0.d0) go to 1
2 knk=i


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!  the outward integration is now complete                                    !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  return

  end subroutine numout





  end module kb

