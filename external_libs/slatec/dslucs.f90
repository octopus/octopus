!DECK DSLUCS                                                            
!***BEGIN PROLOGUE  DSLUCS                                              
!***PURPOSE  Incomplete LU BiConjugate Gradient Squared Ax=b Solver.    
!            Routine to solve a linear system  Ax = b  using the        
!            BiConjugate Gradient Squared method with Incomplete LU     
!            decomposition preconditioning.                             
!***LIBRARY   SLATEC (SLAP)                                             
!***CATEGORY  D2A4, D2B4                                                
!***TYPE      DOUBLE PRECISION (SSLUCS-S, DSLUCS-D)                     
!***KEYWORDS  ITERATIVE INCOMPLETE LU PRECONDITION,                     
!             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE                 
!***AUTHOR  Greenbaum, Anne, (Courant Institute)                        
!           Seager, Mark K., (LLNL)                                     
!             Lawrence Livermore National Laboratory                    
!             PO BOX 808, L-60                                          
!             Livermore, CA 94550 (510) 423-3141                        
!             seager@llnl.gov                                           
!***DESCRIPTION                                                         
!                                                                       
! *Usage:                                                               
!     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX            
!     INTEGER ITER, IERR, IUNIT, LENW, IWORK(NL+NU+4*N+2), LENIW        
!     DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, RWORK(NL+NU+8*N)  
!                                                                       
!     CALL DSLUCS(N, B, X, NELT, IA, JA, A, ISYM, ITOL, TOL,            
!    $     ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW)    
!                                                                       
! *Arguments:                                                           
! N      :IN       Integer.                                             
!         Order of the Matrix.                                          
! B      :IN       Double Precision B(N).                               
!         Right-hand side vector.                                       
! X      :INOUT    Double Precision X(N).                               
!         On input X is your initial guess for solution vector.         
!         On output X is the final approximate solution.                
! NELT   :IN       Integer.                                             
!         Number of Non-Zeros stored in A.                              
! IA     :INOUT    Integer IA(NELT).                                    
! JA     :INOUT    Integer JA(NELT).                                    
! A      :INOUT    Double Precision A(NELT).                            
!         These arrays should hold the matrix A in either the SLAP      
!         Triad format or the SLAP Column format.  See "Description",   
!         below.  If the SLAP Triad format is chosen it is changed      
!         internally to the SLAP Column format.                         
! ISYM   :IN       Integer.                                             
!         Flag to indicate symmetric storage format.                    
!         If ISYM=0, all non-zero entries of the matrix are stored.     
!         If ISYM=1, the matrix is symmetric, and only the upper        
!         or lower triangle of the matrix is stored.                    
! ITOL   :IN       Integer.                                             
!         Flag to indicate type of convergence criterion.               
!         If ITOL=1, iteration stops when the 2-norm of the residual    
!         divided by the 2-norm of the right-hand side is less than TOL.
!         This routine must calculate the residual from R = A*X - B.    
!         This is unnatural and hence expensive for this type of iter-  
!         ative method.  ITOL=2 is *STRONGLY* recommended.              
!         If ITOL=2, iteration stops when the 2-norm of M-inv times the 
!         residual divided by the 2-norm of M-inv times the right hand  
!         side is less than TOL, where M-inv time a vector is the pre-  
!         conditioning step.  This is the *NATURAL* stopping for this   
!         iterative method and is *STRONGLY* recommended.               
! TOL    :INOUT    Double Precision.                                    
!         Convergence criterion, as described above.  (Reset if IERR=4.)
! ITMAX  :IN       Integer.                                             
!         Maximum number of iterations.                                 
! ITER   :OUT      Integer.                                             
!         Number of iterations required to reach convergence, or        
!         ITMAX+1 if convergence criterion could not be achieved in     
!         ITMAX iterations.                                             
! ERR    :OUT      Double Precision.                                    
!         Error estimate of error in final approximate solution, as     
!         defined by ITOL.                                              
! IERR   :OUT      Integer.                                             
!         Return error flag.                                            
!           IERR = 0 => All went well.                                  
!           IERR = 1 => Insufficient space allocated for WORK or IWORK. 
!           IERR = 2 => Method failed to converge in ITMAX steps.       
!           IERR = 3 => Error in user input.                            
!                       Check input values of N, ITOL.                  
!           IERR = 4 => User error tolerance set too tight.             
!                       Reset to 500*D1MACH(3).  Iteration proceeded.   
!           IERR = 5 => Breakdown of the method detected.               
!                       (r0,r) approximately 0.                         
!           IERR = 6 => Stagnation of the method detected.              
!                       (r0,v) approximately 0.                         
!           IERR = 7 => Incomplete factorization broke down and was     
!                       fudged.  Resulting preconditioning may be less  
!                       than the best.                                  
! IUNIT  :IN       Integer.                                             
!         Unit number on which to write the error at each iteration,    
!         if this is desired for monitoring convergence.  If unit       
!         number is 0, no writing will occur.                           
! RWORK  :WORK     Double Precision RWORK(LENW).                        
!         Double Precision array used for workspace.  NL is the number  
!         of non-zeros in the lower triangle of the matrix (including   
!         the diagonal).  NU is the number of non-zeros in the upper    
!         triangle of the matrix (including the diagonal).              
! LENW   :IN       Integer.                                             
!         Length of the double precision workspace, RWORK.              
!         LENW >= NL+NU+8*N.                                            
! IWORK  :WORK     Integer IWORK(LENIW).                                
!         Integer array used for workspace.  NL is the number of non-   
!         zeros in the lower triangle of the matrix (including the      
!         diagonal).  NU is the number of non-zeros in the upper        
!         triangle of the matrix (including the diagonal).              
!         Upon return the following locations of IWORK hold information 
!         which may be of use to the user:                              
!         IWORK(9)  Amount of Integer workspace actually used.          
!         IWORK(10) Amount of Double Precision workspace actually used. 
! LENIW  :IN       Integer.                                             
!         Length of the integer workspace, IWORK.                       
!         LENIW >= NL+NU+4*N+12.                                        
!                                                                       
! *Description:                                                         
!       This routine is simply a  driver for the DCGSN  routine.  It    
!       calls the DSILUS routine to set  up the  preconditioning and    
!       then  calls DCGSN with  the appropriate   MATVEC, MTTVEC and    
!       MSOLVE, MTSOLV routines.                                        
!                                                                       
!       The Sparse Linear Algebra Package (SLAP) utilizes two matrix    
!       data structures: 1) the  SLAP Triad  format or  2)  the SLAP    
!       Column format.  The user can hand this routine either of the    
!       of these data structures and SLAP  will figure out  which on    
!       is being used and act accordingly.                              
!                                                                       
!       =================== S L A P Triad format ===================    
!                                                                       
!       This routine requires that the  matrix A be   stored in  the    
!       SLAP  Triad format.  In  this format only the non-zeros  are    
!       stored.  They may appear in  *ANY* order.  The user supplies    
!       three arrays of  length NELT, where  NELT is  the number  of    
!       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For    
!       each non-zero the user puts the row and column index of that    
!       matrix element  in the IA and  JA arrays.  The  value of the    
!       non-zero   matrix  element is  placed  in  the corresponding    
!       location of the A array.   This is  an  extremely  easy data    
!       structure to generate.  On  the  other hand it   is  not too    
!       efficient on vector computers for  the iterative solution of    
!       linear systems.  Hence,   SLAP changes   this  input    data    
!       structure to the SLAP Column format  for  the iteration (but    
!       does not change it back).                                       
!                                                                       
!       Here is an example of the  SLAP Triad   storage format for a    
!       5x5 Matrix.  Recall that the entries may appear in any order.   
!                                                                       
!           5x5 Matrix      SLAP Triad format for 5x5 matrix on left.   
!                              1  2  3  4  5  6  7  8  9 10 11          
!       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21          
!       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2          
!       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1          
!       | 0  0  0 44  0|                                                
!       |51  0 53  0 55|                                                
!                                                                       
!       =================== S L A P Column format ==================    
!                                                                       
!       This routine  requires that  the matrix A  be stored in  the    
!       SLAP Column format.  In this format the non-zeros are stored    
!       counting down columns (except for  the diagonal entry, which    
!       must appear first in each  "column")  and are stored  in the    
!       double precision array A.   In other words,  for each column    
!       in the matrix put the diagonal entry in  A.  Then put in the    
!       other non-zero  elements going down  the column (except  the    
!       diagonal) in order.   The  IA array holds the  row index for    
!       each non-zero.  The JA array holds the offsets  into the IA,    
!       A arrays  for  the  beginning  of each   column.   That  is,    
!       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the    
!       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1),    
!       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column.    
!       Note that we always have  JA(N+1) = NELT+1,  where N is  the    
!       number of columns in  the matrix and NELT  is the number  of    
!       non-zeros in the matrix.                                        
!                                                                       
!       Here is an example of the  SLAP Column  storage format for a    
!       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a    
!       column):                                                        
!                                                                       
!           5x5 Matrix      SLAP Column format for 5x5 matrix on left.  
!                              1  2  3    4  5    6  7    8    9 10 11  
!       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35  
!       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3  
!       | 0  0 33  0 35|  JA:  1  4  6    8  9   12                     
!       | 0  0  0 44  0|                                                
!       |51  0 53  0 55|                                                
!                                                                       
! *Side Effects:                                                        
!       The SLAP Triad format (IA, JA,  A) is modified internally to    
!       be the SLAP Column format.  See above.                          
!                                                                       
! *Cautions:                                                            
!     This routine will attempt to write to the Fortran logical output  
!     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that  
!     this logical unit is attached to a file or terminal before calling
!     this routine with a non-zero value for IUNIT.  This routine does  
!     not check for the validity of a non-zero IUNIT unit number.       
!                                                                       
!***SEE ALSO  DCGS, DSDCGS                                              
!***REFERENCES  1. P. Sonneveld, CGS, a fast Lanczos-type solver        
!                  for nonsymmetric linear systems, Delft University    
!                  of Technology Report 84-16, Department of Mathe-     
!                  matics and Informatics, Delft, The Netherlands.      
!               2. E. F. Kaasschieter, The solution of non-symmetric    
!                  linear systems by biconjugate gradients or conjugate 
      SUBROUTINE DSLUCS (N, B, X, NELT, IA, JA, A, ISYM, ITOL, TOL,     &
     &   ITMAX, ITER, ITERMIN, ERR, IERR, IUNIT, RWORK, LENW, IWORK,    &
     &   LENIW)                                                         
!                  gradients squared,  Delft University of Technology   
!                  Report 86-21, Department of Mathematics and Informa- 
!                  tics, Delft, The Netherlands.                        
!***ROUTINES CALLED  DCGS, DCHKW, DS2Y, DSILUS, DSLUI, DSMV             
!***REVISION HISTORY  (YYMMDD)                                          
!   890404  DATE WRITTEN                                                
!   890404  Previous REVISION DATE                                      
!   890915  Made changes requested at July 1989 CML Meeting.  (MKS)     
!   890921  Removed TeX from comments.  (FNF)                           
!   890922  Numerous changes to prologue to make closer to SLATEC       
!           standard.  (FNF)                                            
!   890929  Numerous changes to reduce SP/DP differences.  (FNF)        
!   910411  Prologue converted to Version 4.0 format.  (BAB)            
!   920511  Added complete declaration section.  (WRB)                  
!   920929  Corrected format of references.  (FNF)                      
!   921113  Corrected C***CATEGORY line.  (FNF)                         
!***END PROLOGUE  DSLUCS                                                
!     .. Parameters ..                                                  
      INTEGER LOCRB, LOCIB 
      PARAMETER (LOCRB=1, LOCIB=11) 
!     .. Scalar Arguments ..                                            
      DOUBLE PRECISION ERR, TOL 
      INTEGER IERR, ISYM, ITER, ITERMIN, ITMAX, ITOL, IUNIT, LENIW, LENW, N, NELT 
!     .. Array Arguments ..                                             
      DOUBLE PRECISION A(NELT), B(N), RWORK(LENW), X(N) 
      INTEGER IA(NELT), IWORK(LENIW), JA(NELT) 
!     .. Local Scalars ..                                               
      INTEGER ICOL, J, JBGN, JEND, LOCDIN, LOCIL, LOCIU, LOCIW, LOCJL,  &
     &        LOCJU, LOCL, LOCNC, LOCNR, LOCP, LOCQ, LOCR, LOCR0, LOCU, &
     &        LOCUU, LOCV1, LOCV2, LOCW, NL, NU                         
!     .. External Subroutines ..                                        
      EXTERNAL DCGS, DCHKW, DS2Y, DSILUS, DSLUI, DSMV 
!***FIRST EXECUTABLE STATEMENT  DSLUCS                                  
!                                                                       
      IERR = 0 
      IF( N.LT.1 .OR. NELT.LT.1 ) THEN 
         IERR = 3 
         RETURN 
      ENDIF 
      write(6,*)' ierr ',ierr,' in dslucs ' 
      write(6,*)n,nelt,isym,itol,lenw,leniw,iunit 
!                                                                       
!         Change the SLAP input matrix IA, JA, A to SLAP-Column format. 
      CALL DS2Y( N, NELT, IA, JA, A, ISYM ) 
!                                                                       
!         Count number of Non-Zero elements preconditioner ILU matrix.  
!         Then set up the work arrays.                                  
      NL = 0 
      NU = 0 
      DO 20 ICOL = 1, N 
!         Don't count diagonal.                                         
         JBGN = JA(ICOL)+1 
         JEND = JA(ICOL+1)-1 
         IF( JBGN.LE.JEND ) THEN 
!VD$ NOVECTOR                                                           
            DO 10 J = JBGN, JEND 
               IF( IA(J).GT.ICOL ) THEN 
                  NL = NL + 1 
                  IF( ISYM.NE.0 ) NU = NU + 1 
               ELSE 
                  NU = NU + 1 
               ENDIF 
   10       CONTINUE 
         ENDIF 
   20 END DO 
!                                                                       
      LOCIL = LOCIB 
      LOCJL = LOCIL + N+1 
      LOCIU = LOCJL + NL 
      LOCJU = LOCIU + NU 
      LOCNR = LOCJU + N+1 
      LOCNC = LOCNR + N 
      LOCIW = LOCNC + N 
!                                                                       
      LOCL   = LOCRB 
      LOCDIN = LOCL + NL 
      LOCUU  = LOCDIN + N 
      LOCR   = LOCUU + NU 
      LOCR0  = LOCR + N 
      LOCP   = LOCR0 + N 
      LOCQ   = LOCP + N 
      LOCU   = LOCQ + N 
      LOCV1  = LOCU + N 
      LOCV2  = LOCV1 + N 
      LOCW   = LOCV2 + N 
!                                                                       
!         Check the workspace allocations.                              
      CALL DCHKW( 'DSLUCS', LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR ) 
      write(6,*)' ierr ',ierr 
      IF( IERR.NE.0 ) RETURN 
!                                                                       
      IWORK(1) = LOCIL 
      IWORK(2) = LOCJL 
      IWORK(3) = LOCIU 
      IWORK(4) = LOCJU 
      IWORK(5) = LOCL 
      IWORK(6) = LOCDIN 
      IWORK(7) = LOCUU 
      IWORK(9) = LOCIW 
      IWORK(10) = LOCW 
!                                                                       
!         Compute the Incomplete LU decomposition.                      
      CALL DSILUS( N, NELT, IA, JA, A, ISYM, NL, IWORK(LOCIL),          &
     &     IWORK(LOCJL), RWORK(LOCL), RWORK(LOCDIN), NU, IWORK(LOCIU),  &
     &     IWORK(LOCJU), RWORK(LOCUU), IWORK(LOCNR), IWORK(LOCNC) )     
!                                                                       
!         Perform the incomplete LU preconditioned                      
!         BiConjugate Gradient Squared algorithm.                       
      CALL DCGS(N, B, X, NELT, IA, JA, A, ISYM, DSMV,                   &
     &     DSLUI, ITOL, TOL, ITMAX, ITER, ITERMIN,                      &
     &     ERR, IERR, IUNIT,                                            &
     &     RWORK(LOCR), RWORK(LOCR0), RWORK(LOCP),                      &
     &     RWORK(LOCQ), RWORK(LOCU), RWORK(LOCV1),                      &
     &     RWORK(LOCV2), RWORK, IWORK )                                 
      write(6,*)' dslucs iter,itermin,ierr ',iter,itermin,ierr 
      RETURN 
!------------- LAST LINE OF DSLUCS FOLLOWS ---------------------------- 
      END                                           
