C SL_GMINR.F - Gateway function for SLICOT minimal realization routines
C              TG01JD.F, TB01PX.F and TG01GD.F.
C
C MATLAB call:
C   For meth = 1
C   [Ar,Er,Br,Cr,Dr,redinfo,tau] = ...
C                 sl_gminr(meth,A,E,B,C,D,tol,job,systype)
C   For meth = 2
C   [Ar,Er,Br,Cr,Dr,redinfo] = sl_gminr(meth,A,E,B,C,D,tol,job)
C   For meth = 3
C   [Ar,Er,Br,Cr,Dr,redinfo,tau] = ...
C                 sl_gminr(meth,A,E,B,C,D,tol,job,systype)
C
C Purpose:
C   To compute a reduced (controllable, observable, or minimal)
C   state-space realization (Ar-lambda*Er,Br,Cr) for a given descriptor
C   representation (A-lambda*E,B,C). For meth = 1, the pencil
C   Ar-lambda*Er is in upper block Hessenberg form, whose structure
C   is defined by redinfo and tau. If E = [] or E = I minimal
C   realization for standard systems is performed.
C
C Input parameters:
C   meth    - type of order reduction to be performed
C               = 1 : irreducible descriptor realization
C               = 2 : non-dynamic modes reduction
C               = 3 : minimal descriptor realization
C   A,E,
C   B,C,D   - descriptor system matrices
C   tol     - (optional) tolerance parameter for controllability
C             and observability tests
C   job     - (optional) option parameter
C                for meth = 1,3, reduction/scaling option of decimal
C                form IJ
C                  J = 0 - remove both uncontrollable and
C                          unobservable eigenvalues (default)
C                  J = 1 - remove all uncontrollable eigenvalues
C                  J = 2 - remove all unobservable eigenvalues
C                  I = 0 - no scaling (default)
C                  I = 1 - perform scaling
C                for meth = 2, system output form option
C                  = 0 - descriptor system form
C                  = 1 - standard system form with leading block
C                        of E identity
C   systype - (optional) system type
C                = 0 - general improper system (default)
C                = 1 - proper system
C                = 2 - polynomial system
C
C  Output parameters:
C  Ar, Er,
C  Br, Cr,
C  Dr       - matrices of the reduced descriptor system
C  redinfo  - reduction information
C             redinfo(k) > 0 - reduction achieved in Phase k
C                              k = 1..4 (see TG01JD)
C                              k = 5 (see TG01GD)
C             redinfo(6) - number of nonzero subdiagonals of Ar
C             redinfo(7) - number of nonzero subdiagonals of Er
C             redinfo(8) >= 0 - rank of E
C             redinfo(8)  < 0 - rank of E not computed
C  tau      - for meth = 1 or 3: vector containing the dimension of full
C                        row/column rank matrices in the resulting
C                        staircase form (see TG01JD)
C
C Contributor
C   A. Varga, May 1999.
C
C Revisions:
C   A. Varga,  13-01-2016.
C   V. Sima, Dec. 2016, Feb. 2017, Apr. 2017.
C
C **********************************************************************
C
#include "fintrf.h"
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C .. Parameters ..
      IMPLICIT          NONE
      DOUBLE PRECISION  ZERO
      PARAMETER         ( ZERO = 0.0D0 )
C
C .. Mex-file interface parameters ..
      MWPOINTER         PLHS(*), PRHS(*)
      INTEGER*4         NLHS, NRHS
C
C .. Mex-file integer functions ..
      MWPOINTER         mxCreateDoubleMatrix, mxGetPr
      MWSIZE            mxGetM, mxGetN
      INTEGER*4         mxIsNumeric, mxIsComplex
C
C .. Scalar parameters used by SLICOT subroutines ..
      CHARACTER         EQUIL, JOB, JOBS, SYSTYP
      INTEGER           INFO, INFOKD, LDA, LDB, LDC, LDD, LDE, LDWORK,
     $                  M, N, NR, P, RANKE
      MWSIZE            MW, NRW, NW, PW
      DOUBLE PRECISION  TOL
C
C ..  Fixed dimension arrays used by SLICOT subroutines ..
      INTEGER           INFRED(8)
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      DOUBLE PRECISION, ALLOCATABLE:: A(:,:), B(:,:), C(:,:), D(:,:),
     $                                DWORK(:), E(:,:)
      INTEGER, ALLOCATABLE:: IWORK(:), RTAU(:)
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      LOGICAL           SSTYPE
      INTEGER           I, ISCALE, ISTYPE, JOBN, L, L1, L2, LE, M4,
     $                  METH, NBLCK, N1, N3, NN, P4
      MWSIZE            L1W, L2W, LW, M4W, N1W, N3W, NBLCKW, NBLW, P4W
      MWSIZE            ZEROW, ONEW, EIGHTW
      DOUBLE PRECISION  DUM(8)
C
C     .. matlab integer*4 flags
C
      INTEGER*4         SZERO
      PARAMETER         ( SZERO = 0 )
C
C .. External functions ..
      LOGICAL           CHCKID
      EXTERNAL          CHCKID
C
C .. External Subroutines ..
      EXTERNAL          DLACPY, DLASET, TB01PX, TG01GD, TG01JD
C
C .. Intrinsic functions ..
      INTRINSIC         INT, MAX, MIN, MOD
C
      ZEROW  = 0
      ONEW   = 1
      EIGHTW = 8
C
*
* Debuging by writing in a file
*
*      OPEN (1, FILE = 'test.txt', STATUS = 'REPLACE')
*
      IF( NRHS.EQ.0 ) RETURN
C
C Check dimensions of input parameters and read/set scalar parameters
C
C  meth
C
      IF( mxGetM( PRHS(1) ).NE.1 .OR. mxGetN( PRHS(1) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'METH must be a scalar' )
      END IF
      IF( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'METH must be an integer scalar' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), DUM, ONEW )
      METH = DUM(1)
      IF( METH.LE.0 .OR. METH.GT.3 ) THEN
         CALL mexErrMsgTxt
     $        ( 'METH has 1 ... 3 the only admissible values' )
      END IF
C
C Check for proper number of arguments
C
      IF( NRHS.LT.6 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SL_GMINR requires at least 6 input arguments' )
      END IF
      IF( METH.EQ.1 .OR. METH.EQ.3) THEN
         IF( NLHS.GT.7 ) THEN
            CALL mexErrMsgTxt
     $           ( 'SL_GMINR has at most 7 output arguments' )
         END IF
      ELSE IF( METH.EQ.2 ) THEN
         IF( NLHS.GT.6 ) THEN
            CALL mexErrMsgTxt
     $           ( 'SL_GMINR has at most 6 output arguments' )
         END IF
      END IF
C
C   A(NxN), E(NxN), B(NxM), C(PxN), D(PxM)
C
      LW  = mxGetM( PRHS(2) )
      L   = LW
      NW  = mxGetN( PRHS(2) )
      N   = NW
      L1W = mxGetM( PRHS(3) )
      L1  = L1W
      N1W = mxGetN( PRHS(3) )
      N1  = N1W
      L2W = mxGetM( PRHS(4) )
      L2  = L2W
      MW  = mxGetN( PRHS(4) )
      M   = MW
      PW  = mxGetM( PRHS(5) )
      P   = PW
      N3W = mxGetN( PRHS(5) )
      N3  = N3W
      P4W = mxGetM( PRHS(6) )
      P4  = P4W
      M4W = mxGetN( PRHS(6) )
      M4  = M4W
C
      IF( L.NE.N .OR. L1.NE.N1 ) THEN
         CALL mexErrMsgTxt
     $        ( 'A and E must be square matrices' )
      END IF
      SSTYPE = L1.EQ.0
      IF( .NOT.SSTYPE .AND. ( (N1.NE.N) .OR. (L1.NE.L) ) ) THEN
         CALL mexErrMsgTxt
     $        ( 'A and E must have the same dimensions' )
      END IF
      IF( M.GT.0 .AND. L2.NE.L ) THEN
         CALL mexErrMsgTxt
     $        ( 'B must have the same row dimension as A ' )
      END IF
      IF( P.GT.0 .AND. N3.NE.N ) THEN
         CALL mexErrMsgTxt
     $        ( 'C must have the same column dimension as A ' )
      END IF
      IF( MIN( P4, M4 ).GT.0 .AND. N.GT.0 .AND.
     $    P4.NE.P .AND. M4.NE.M ) THEN
         CALL mexErrMsgTxt
     $        ( 'D must have compatible dimensions with B and C' )
      END IF
      IF( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a real matrix' )
      END IF
      IF( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'E must be a real matrix' )
      END IF
      IF( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(4) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'B must be a real matrix' )
      END IF
      IF( mxIsNumeric( PRHS(5) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(5) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'C must be a real matrix' )
      END IF
      IF( mxIsNumeric( PRHS(6) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(6) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'D must be a real matrix' )
      END IF
      M = MAX( M, M4 )
      MW = M
      P = MAX( P, P4 )
      PW = P
C
C   tol
C
      TOL = ZERO
      IF( NRHS.GE.7 ) THEN
         IF( mxGetM( PRHS(7) )*mxGetN( PRHS(7) ).GT.1 ) THEN
            CALL mexErrMsgTxt( 'TOL must be scalar' )
         END IF
         IF( mxIsNumeric( PRHS(7) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(7) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'TOL must be a real scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(7) ), TOL, ONEW )
      END IF
C
C     job
C
      JOBN = 0
      ISCALE = 0
      IF( NRHS.GE.8 ) THEN
         IF( mxGetM( PRHS(8) )*mxGetN( PRHS(8) ).GT.1 ) THEN
            CALL mexErrMsgTxt( 'JOB must be scalar' )
         END IF
         IF( mxIsNumeric( PRHS(8) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(8) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'JOB must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(8) ), DUM, ONEW )
         JOBN = DUM(1)
         ISCALE = JOBN/10
         JOBN = MOD( JOBN, 10 )
         IF( METH.EQ.2 ) THEN
            IF( JOBN.NE.0 .AND. JOBN.NE.1 ) THEN
               CALL mexErrMsgTxt
     $         ('JOB has 0 and 1 the only admissible values')
            END IF
         ELSE
            IF( .NOT.( JOBN.GE.0 .AND. JOBN.LE.2 ) .OR.
     $          .NOT.( ISCALE.EQ.0 .OR. ISCALE.EQ.1 ) ) THEN
                CALL mexErrMsgTxt
     $          ('JOB has 0..2 and 10..12 the only admissible values')
            END IF
         END IF
      END IF
C
C     systype
C
      ISTYPE = 0
      IF( NRHS.GE.9 ) THEN
         IF( METH.EQ.2 ) THEN
            CALL mexErrMsgTxt
     $           ( 'SL_GMINR has at most 8 input arguments' )
         END IF
         IF( mxGetM( PRHS(9) )*mxGetN( PRHS(9) ).GT.1 ) THEN
            CALL mexErrMsgTxt( 'SYSTYPE must be scalar' )
         END IF
         IF( mxIsNumeric( PRHS(9) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(9) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'SYSTYPE must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(9) ), DUM, ONEW )
         ISTYPE = DUM(1)
         IF( ISTYPE.LT.0 .OR. ISTYPE.GT.2 ) THEN
            CALL mexErrMsgTxt
     $      ( 'SYSTYPE has 0 ... 2 the only admissible values' )
         END IF
      END IF
C
      IF( METH.EQ.2 ) THEN
         IF( JOBN.EQ.0 ) THEN
            JOBS = 'D'
         ELSE
            JOBS = 'S'
         END IF
      ELSE
         JOBS = 'D'
      END IF
C
C  Set Problem Parmeters
C
      LDA = MAX( 1, N )
      IF( SSTYPE ) THEN
         LDE = 1
         LE  = 1
      ELSE
         LDE = LDA
         LE  = LDE
      END IF
      LDB = LDA
      LDC = MAX( 1, M, P )
      LDD = LDC
C
C Determine the lengths of working arrays
C
      LDWORK = 0
      IF( METH.EQ.1 .OR. METH.EQ.3 ) THEN
         IF( SSTYPE ) THEN
            LDWORK = MAX( 1, N + MAX( N, 3*M, 3*P ) + N*( N + M + P ) )
         ELSE
            LDWORK = MAX( 1, 8*N, N*( 2*N + M + P ) +
     $                            MAX( N, 2*M, 2*P ) )
         END IF
      END IF
      IF( METH.GE.2 .AND. .NOT.SSTYPE ) THEN
         CALL TG01GD( JOBS, N, N, M, P, A, LDA, E, LDE, B, LDB, C, LDC,
     $                D, LDD, NR, NR, RANKE, INFRED, TOL, IWORK, DUM,
     $                -1, INFOKD )
         LDWORK = MAX( LDWORK, INT( DUM(1) ), N*( 2*N + M + P ) + 
     $                 MAX( N+P, MIN( L, N ) + MAX( 3*N-1, M, L ) ) )
      END IF
      LDWORK = MAX( LDWORK, MAX( N, M )*P )
C
C Allocate variable dimension local arrays
C !Fortran 90/95
C
      ALLOCATE ( A(LDA,MAX( 1, N )), B(LDB,LDC),
     $           C(LDC,MAX( 1, N )), E(LDE,MAX( 1, N )),
     $           D(LDD,LDC), DWORK(LDWORK),
     $           IWORK(MAX( 1, 2*N+MAX( M, P ) )) )
C
C Copy inputs from MATLAB workspace to locally allocated arrays
C
      CALL mxCopyPtrToReal8( mxGetPR( PRHS(2) ), A, NW*NW )
      IF( .NOT. SSTYPE )
     $   CALL mxCopyPtrToReal8( mxGetPR( PRHS(3) ), E, NW*NW )
      CALL mxCopyPtrToReal8( mxGetPR( PRHS(4) ), B, NW*MW )
      CALL mxCopyPtrToReal8( mxGetPR( PRHS(5) ), DWORK, PW*NW )
      CALL DLACPY( 'F', P, N, DWORK, MAX( 1, P ), C, LDC )
      IF( MIN( P4, M4 ).EQ.0 ) THEN
         CALL DLASET( 'F', P, M, ZERO, ZERO, D, LDD )
      ELSE
         CALL mxCopyPtrToReal8( mxGetPR( PRHS(6) ), DWORK, PW*MW )
         CALL DLACPY( 'F', P, M, DWORK, MAX( 1, P ), D, LDD )
      END IF
C
C  Set options.
C
      IF( .NOT.SSTYPE ) SSTYPE = CHCKID( N, E, LDE )
      IF( METH.NE.2 ) THEN
         IF( JOBN.EQ.0 ) THEN
            JOB = 'I'
            IF( SSTYPE ) JOB = 'M'
         ELSE IF( JOBN.EQ.1 ) THEN
            JOB = 'C'
         ELSE
            JOB = 'O'
         END IF
      END IF
      IF( ISTYPE.EQ.0 ) THEN
         SYSTYP = 'R'
      ELSE IF( ISTYPE.EQ.1 ) THEN
         SYSTYP = 'S'
      ELSE
         SYSTYP = 'P'
      END IF
      IF( ISCALE.EQ.0 ) THEN
         EQUIL = 'N'
      ELSE
         EQUIL = 'S'
      END IF
C
C Do the actual computations
C
      DO 10 I = 1, 5
         DUM(I) = -1
   10 CONTINUE
      DUM(6) = MAX( 0, N-1 )
      IF( SSTYPE) THEN
         DUM(7) = 0
         DUM(8) = N
      ELSE
         DUM(7) = MAX( 0, N-1 )
         DUM(8) = -1
      END IF
      NR = N
      NN = N
      INFO = 0
      INFOKD = 0
      IF( METH.EQ.1 .OR. METH.EQ.3 ) THEN
C
C        Irreducible realization
C
         IF( SSTYPE ) THEN
            CALL TB01PX( JOB, EQUIL, N, M, P, A, LDA, B, LDB, C, LDC,
     $                   NR, INFRED, TOL, IWORK, DWORK, LDWORK, INFO )
            DUM(1) = INFRED(1)
            DUM(3) = INFRED(2)
            DUM(6) = INFRED(3)
            DUM(8) = N
         ELSE
            CALL TG01JD( JOB, SYSTYP, EQUIL, N, M, P, A, LDA, E, LDE,
     $                   B, LDB, C, LDC, NR, INFRED, TOL, IWORK, DWORK,
     $                   LDWORK, INFO )
            DUM(1) = INFRED(1)
            DUM(2) = INFRED(2)
            DUM(3) = INFRED(3)
            DUM(4) = INFRED(4)
            DUM(6) = INFRED(5)
            DUM(7) = INFRED(6)
            IF( INFO.EQ.0 .AND. METH.EQ.3 .AND. NLHS.GE.7 ) THEN
               ALLOCATE ( RTAU(INFRED(7)) )
               DO 20 I = 1, INFRED(7)
                  RTAU(I) = IWORK(I)
   20          CONTINUE
            END IF
         END IF
         NN = NR
      END IF
      IF( INFO.EQ.0 .AND. METH.GE.2 .AND. .NOT.SSTYPE ) THEN
C
C        Elliminate non-dynamic modes.
C
         CALL TG01GD( JOBS, NN, NN, M, P, A, LDA, E, LDE, B, LDB,
     $                C, LDC, D, LDD, NR, NR, RANKE, INFRED, TOL,
     $                IWORK, DWORK, LDWORK, INFOKD )
         DUM(5) = INFRED(1)
         IF( INFRED(1).GE.0 ) THEN
            DUM(6) = MAX( 0, NR-1 )
            DUM(7) = 0
         END IF
         DUM(8) = RANKE
      END IF
C
C Copy output to MATLAB workspace
C
      IF( INFO.EQ.0 .AND. INFOKD.EQ.0 ) THEN
         IF( NLHS.GE.1 ) THEN
            IF( NR.LT.N )
     $         CALL DLACPY( 'F', NR, NR, A, LDA, A, MAX( 1, NR ) )
            NRW = NR
            PLHS(1) = mxCreateDoubleMatrix( NRW, NRW, SZERO )
            CALL mxCopyReal8ToPtr( A, mxGetPR( PLHS(1) ), NRW*NRW )
         END IF
         IF( NLHS.GE.2 ) THEN
            IF( L1.GT.0 ) THEN
               IF( NR.LT.N )
     $            CALL DLACPY( 'F', NR, NR, E, LDE, E, MAX( 1, NR ) )
               NRW = NR
               PLHS(2) = mxCreateDoubleMatrix( NRW, NRW, SZERO )
               CALL mxCopyReal8ToPtr( E, mxGetPR( PLHS(2) ), NRW*NRW )
            ELSE
               PLHS(2) = mxCreateDoubleMatrix( ZEROW, ZEROW, SZERO )
            END IF
         END IF
         IF( NLHS.GE.3 ) THEN
            IF( NR.LT.N )
     $         CALL DLACPY( 'F', NR, M, B, LDB, B, MAX( 1, NR ) )
            NRW = NR
            MW  = M
            PLHS(3) = mxCreateDoubleMatrix( NRW, MW, SZERO )
            CALL mxCopyReal8ToPtr( B, mxGetPR( PLHS(3) ), NRW*MW )
         END IF
         IF( NLHS.GE.4 ) THEN
            IF( M.GT.P )
     $         CALL DLACPY( 'F', P, NR, C, LDC, C, MAX( 1, P ) )
            NRW = NR
            PW  = P
            PLHS(4) = mxCreateDoubleMatrix( PW, NRW, SZERO )
            CALL mxCopyReal8ToPtr( C, mxGetPR( PLHS(4) ), PW*NRW )
         END IF
         IF( NLHS.GE.5 ) THEN
            IF( M.GT.P )
     $         CALL DLACPY( 'F', P, M, D, LDD, D, MAX( 1, P ) )
            MW  = M
            PW  = P
            PLHS(5) = mxCreateDoubleMatrix( PW, MW, SZERO )
            CALL mxCopyReal8ToPtr( D, mxGetPR( PLHS(5) ), PW*MW )
         END IF
         IF( NLHS.GE.6 ) THEN
            PLHS(6) = mxCreateDoubleMatrix( ONEW, EIGHTW, SZERO )
            CALL mxCopyReal8ToPtr( DUM, mxGetPR( PLHS(6) ), EIGHTW )
         END IF
         IF( METH.NE.2) THEN
            IF( NLHS.GE.7 ) THEN
               IF( SSTYPE ) THEN
                  NBLCK = INFRED(4)
               ELSE
                  NBLCK = INFRED(7)
               END IF
               IF( METH.EQ.1 .OR. SSTYPE ) THEN
                  DO 30 I = 1, NBLCK
                     DWORK(I) = IWORK(I)
   30             CONTINUE
               ELSE
                  DO 40 I = 1, NBLCK
                     DWORK(I) = RTAU(I)
   40             CONTINUE
                  DEALLOCATE ( RTAU )
               END IF
               NBLCKW = NBLCK
               NBLW   = MIN( 1, NBLCK )
               PLHS(7) = mxCreateDoubleMatrix( NBLW, NBLCKW, SZERO )
               CALL mxCopyReal8ToPtr( DWORK, mxGetPR( PLHS(7) ),
     $                                NBLCKW )
            END IF
         END IF
      END IF
C
C Deallocate local arrays
C !Fortran 90/95
C
      DEALLOCATE (A, B, C, D, DWORK, E, IWORK )
*      close(1)
C
C Error handling
C
      IF( INFO.NE.0 .OR. INFOKD.NE.0 ) THEN
         IF( INFO.NE.0 ) THEN
            IF( SSTYPE ) THEN
               write( TEXT,
     $                '( " INFO =", I4, " ON EXIT FROM TB01PX" )' ) INFO
            ELSE
               write( TEXT,
     $                '( " INFO =", I4, " ON EXIT FROM TGO1JD" )' ) INFO
            END IF
         ELSE IF( INFOKD.NE.0 ) THEN
               write( TEXT,
     $              '( " INFO =", I4, " ON EXIT FROM TGO1KD" )' ) INFOKD
         END IF
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of SL_GMINR ***
      END
