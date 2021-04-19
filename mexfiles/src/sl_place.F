C SL_PLACE.F - Gateway function for SLICOT eigenvalue assignment
C              routine SB01BD.F.
C
C MATLAB call:
C   [F, Z, Acl,info] = sl_place(meth,A,B,poles,tol,discr,alpha)
C
C Purpose:
C   To compute a state feedback F which assigns the closed-loop
C   system poles, i.e. the eigenvalues of A+B*F to desired
C   locations. Partial pole placement is possible by keeping
C   unaltered eigenvalues having an acceptable stability degree.
C
C Input parameters:
C   meth    - eigenvalue assignement method
C               = 1 : Schur method
C   A,B     - state-space system matrices
C   poles   - desired eigenvalues for A+B*F
C   tol     - (optional) controllability tolerance
C             default: tol = epsilon_machine*max(norm(A),norm(B))
C   discr   - (optional) type of system
C                = 0 : continuous-time (default)
C                = 1 : discrete-time
C   alpha   - (optional) acceptable stability degree
C             default: -1.e7  for continuous-time
C                        0.0  for discrete-time
C
C Output parameters:
C   F       - feedback matrix assigning the closed loop eiegenvalues
C   Z       - orthonal matrix such that Z'*(A+B*F)*Z is a real
C             Schur matrix
C   Acl     - closed-loop state matrix Z'*(A+B*F)*Z
C   info    - information vector assembled as
C             [nfp, nap, nup, iwarn, ierr], where
C             nfp - number of fixed eigenvalues of A
C             nap - number of assigned eigenvalues from provided poles
C             nup - number of detected uncontrollable eigenvalues of A
C             iwarn - warning flag from SB01BD
C             ierr  - error flag from SB01BD
C
C Contributor
C   A. Varga, May 1999.
C
C Revisions:
C   A. Varga, March 2002.
C   V. Sima, Dec. 2016.
C
C **********************************************************************
C
#include "fintrf.h"
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         ( ZERO = 0.0D0 )
C
C .. Mex-file interface parameters ..
      MWPOINTER         PLHS(*), PRHS(*)
      INTEGER*4         NLHS, NRHS
C
C .. Mex-file integer functions ..
      MWPOINTER         mxCreateDoubleMatrix, mxGetPr, mxGetPi
      MWSIZE            mxGetM, mxGetN
      INTEGER*4         mxIsNumeric, mxIsComplex
C
C .. Scalar parameters used by SLICOT subroutines ..
      CHARACTER         DICO
      INTEGER           INFO, IWARN, LDA, LDB, LDF, LDWORK, LDZ, M, N,
     $                  NAP, NFP, NP, NUP
      MWSIZE            MW, NW, NPW
      DOUBLE PRECISION  ALPHA, TOL
C
C .. Allocatable local arrays ..
C !Fortran 90/95
      DOUBLE PRECISION, ALLOCATABLE:: A(:,:), B(:,:), DWORK(:), F(:,:),
     $                                WI(:), WR(:), Z(:,:)
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      LOGICAL           DISCR, RFLAG
      INTEGER           I, IP, METH, N1, N2
      MWPOINTER         IPW
      MWSIZE            N1W, N2W
      MWSIZE            ONEW, FIVEW
      DOUBLE PRECISION  DUM
C
C     .. matlab integer*4 flags
C
      INTEGER*4         SZERO
      PARAMETER         ( SZERO = 0 )
C
C .. External functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH
C
C .. External Subroutines ..
      EXTERNAL          SB01BD
C
C .. Intrinsic functions ..
      INTRINSIC         MAX
*
* Debuging by writing in a file
*
*      OPEN (1, FILE = 'test.txt', STATUS = 'REPLACE')
*
      ONEW  = 1
      FIVEW = 5
      IF( NRHS.EQ.0 ) RETURN
C
C Check for proper number of arguments
C
      IF( NRHS.LT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SL_PLACE requires at least 4 input arguments' )
      ELSE IF( NLHS.GT.4 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SL_PLACE has at most 4 output argument' )
      END IF
C
C Check dimensions of input parameters and read/set scalar parameters
C
C  meth
C
      METH = 1
      IF( mxGetM( PRHS(1) ).NE.1 .OR. mxGetN( PRHS(1) ).NE.1 ) THEN
         CALL mexErrMsgTxt( 'METH must be a scalar' )
      END IF
      IF( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'METH must be an integer scalar 0 or 1' )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), DUM, ONEW )
      METH = DUM
      IF( METH.LE.0 .OR. METH.GT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 'METH has 1 ... 1 the only admissible values' )
      END IF
C
C   A(NxN), B(NxM)
C
      NW  = mxGetM( PRHS(2) )
      N   = NW
      MW  = mxGetN( PRHS(3) )
      M   = MW
      NPW = mxGetM( PRHS(4) )*mxGetN( PRHS(4) )
      NP  = NPW
      N1W = mxGetN( PRHS(2) )
      N1  = N1W
      N2W = mxGetM( PRHS(3) )
      N2  = N2W
C
      IF( (N1.NE.N) ) THEN
         CALL mexErrMsgTxt
     $        ( 'A must be a square matrix' )
      END IF
      IF( N2.NE.N ) THEN
         CALL mexErrMsgTxt
     $        ( 'B must have the same row dimension as A' )
      END IF
      IF( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a real matrix' )
      END IF
      IF( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'B must be a real matrix' )
      END IF
      IF( mxIsNumeric( PRHS(4) ).EQ.0 ) THEN
         CALL mexErrMsgTxt( 'POLES must be a numeric matrix' )
      END IF
C
C   tol
C
      TOL = ZERO
      IF( NRHS.GT.4 ) THEN
        IF( mxGetM( PRHS(5) ).NE.1 .OR. mxGetN( PRHS(5) ).NE.1 ) THEN
           CALL mexErrMsgTxt( 'TOL must be a scalar' )
        END IF
      IF( mxIsNumeric( PRHS(5) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(5) ).EQ.1 ) THEN
           CALL mexErrMsgTxt( 'TOL must be a real scalar' )
        END IF
        CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), TOL, ONEW )
      END IF
C
C   discr
C
      DICO = 'C'
      IF( NRHS.GT.5 ) THEN
         IF( mxGetM( PRHS(6) ).NE.1 .OR. mxGetN( PRHS(6) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'DISCR must be a scalar' )
         END IF
         IF( mxIsNumeric( PRHS(6) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(6) ).EQ.1 ) THEN
            CALL mexErrMsgTxt
     $           ( 'DISCR must be an integer scalar 0 or 1' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(6) ), DUM, ONEW )
         IF( DUM.NE.ZERO ) DICO = 'D'
      END IF
      DISCR = LSAME( DICO, 'D' )
C
C   alpha
C
      IF( NRHS.GT.6 ) THEN
         IF( mxGetM( PRHS(7) ).NE.1 .OR. mxGetN( PRHS(7) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'ALPHA must be a scalar' )
         END IF
         IF( mxIsNumeric( PRHS(7) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(7) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'ALPHA must be a real scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr(PRHS(7) ), ALPHA, ONEW )
      ELSE
         ALPHA = -1.D7
         IF( DISCR ) ALPHA = ZERO
      END IF
C
C Set the length of working array.
C
      LDWORK = MAX( 1, MAX( 5*M, 5*N, 2*N+4*M ) )
C
      LDA = MAX( 1, N )
      LDB = MAX( 1, N )
      LDF = MAX( 1, M )
      LDZ = LDA
C
C Allocate variable dimension local arrays
C !Fortran 90/95
C
      ALLOCATE ( A(LDA,N), B(LDB,M), DWORK(LDWORK), F(LDF,N), WI(NP),
     $           WR(NP), Z(LDZ,N) )
C
C Copy inputs from MATLAB workspace to locally allocated arrays
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), A, NW*NW )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), B, NW*MW )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), WR, NPW )
      IPW = mxGetPi( PRHS(4) )
      IP  = IPW
      IF( IP.EQ.0 ) THEN
         DO 10 I = 1, NP
            WI(I) = ZERO
   10    CONTINUE
      ELSE
         CALL mxCopyPtrToReal8( IPW, WI, NPW )
      END IF
      RFLAG = .TRUE.
      INFO = 0
      DO I = 1, NP
         IF( RFLAG .AND. WI(I).NE.ZERO ) THEN
            RFLAG = .FALSE.
            IF( I.LT.NP ) THEN
               IF( WR(I).NE.WR(I+1) .OR. (WI(I)+WI(I+1)).NE.ZERO )
     $            INFO = 1
            ELSE
               INFO = 1
            END IF
         ELSE
            RFLAG = .TRUE.
         END IF
         IF( INFO.GT.0 ) THEN
            CALL mexErrMsgTxt(
     $      'POLES must be a self-conjugated vector of complex numbers')
         END IF
      END DO
C
C Do the actual computations
C
C
C Schur method
C
      IF( METH.EQ.1 )
     $   CALL SB01BD( DICO, N, M, NP, ALPHA, A, LDA, B, LDB, WR, WI,
     $                NFP, NAP, NUP, F, LDF, Z, LDZ, TOL, DWORK,
     $                LDWORK, IWARN, INFO )
*      close(1)
C
C Copy output to MATLAB workspace
C
      IF( INFO.EQ.0 .OR. (METH.EQ.1 .AND. INFO.EQ.3) ) THEN
         MW = M
         NW = N
         IF( NLHS.GE.1 ) THEN
            PLHS(1) = mxCreateDoubleMatrix( MW, NW, SZERO )
            CALL mxCopyReal8ToPtr( F, mxGetPr( PLHS(1) ), MW*NW )
      END IF
         IF( NLHS.GE.2 ) THEN
            PLHS(2) = mxCreateDoubleMatrix( NW, NW, SZERO )
            CALL mxCopyReal8ToPtr( Z, mxGetPr( PLHS(2) ), NW*NW )
      END IF
         IF( NLHS.GE.3 ) THEN
            PLHS(3) = mxCreateDoubleMatrix( NW, NW, SZERO )
            CALL mxCopyReal8ToPtr( A, mxGetPr( PLHS(3) ), NW*NW )
      END IF
         IF( NLHS.GE.4 ) THEN
            PLHS(4)  = mxCreateDoubleMatrix( FIVEW, ONEW, SZERO )
            DWORK(1) = NFP
            DWORK(2) = NAP
            DWORK(3) = NUP
            DWORK(4) = IWARN
            DWORK(5) = INFO
            CALL mxCopyReal8ToPtr( DWORK, mxGetPr( PLHS(4) ), FIVEW )
         END IF
C
      END IF
C
C Deallocate local arrays
C !Fortran 90/95
C
      DEALLOCATE ( A, B, F, Z, WR, WI, DWORK )
C
C Error handling
C
      IF( INFO.NE.0 ) THEN
          IF( METH.EQ.1 .AND. INFO.NE.3 ) THEN
      write( TEXT,'( " INFO =", I4, " ON EXIT FROM SB01BD" )' )
     $              INFO
             CALL mexErrMsgTxt( TEXT )
          END IF
      END IF
C
      RETURN
C *** Last line of SL_PLACE ***
      END
