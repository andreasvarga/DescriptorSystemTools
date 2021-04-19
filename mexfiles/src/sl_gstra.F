C SL_GSTRA.F - Gateway function for SLICOT system transformation routines
C              TG01AD, TG01BD, TG01CD, TG01DD, TG01ED, TG01FD, TG01HD,
C              TG01HU, TG01ID, TB01ID, TB01UD, TB01UX, TB01UY, TB01WX.
C
C MATLAB call:
C   For meth = 1
C   [At,Et,Bt,Ct,Q,Z] = sl_gstra(meth,A,E,B,C,job,thresh)
C   For meth = 2
C   [At,Et,Bt,Q] = sl_gstra(meth,A,E,B)
C   For meth = 3
C   [At,Et,Ct,Z] = sl_gstra(meth,A,E,C)
C   For meth = 4, 5, 6, 7
C   [At,Et,Bt,Ct,ranks,Q,Z] = sl_gstra(meth,A,E,B,C,tol)
C   For meth = 8, 9
C   [At,Et,Bt,Ct,dims,tau,Q,Z] = sl_gstra(meth,A,E,B,C,job,tol)
C   For meth = 10
C   [At,Et,Bt,Ct,ncont,tau,Q,Z] = sl_gstra(meth,A,E,B,C,m1,tol)
C   For meth = 11
C   [At,Et,Bt,Ct,Q,Z] = sl_gstra(meth,A,E,B,C)
C
C Purpose:
C   To compute several similarity transformations on a
C   a given descriptor system (A-lambda*E,B,C).
C
C Input parameters:
C   meth    - type of system transformation to be performed
C               = 1 : descriptor system scaling
C               = 2 : reduction to QR-coordinate form
C               = 3 : reduction to RQ-coordinate form
C               = 4 : reduction to SVD-coordinate form
C               = 5 : reduction to more detailed SVD-coordinate form
C               = 6 : reduction to SVD-like coordinate form
C               = 7 : reduction to more detailed SVD-like
C                     coordinate form
C               = 8 : reduction to controllability forms
C               = 9 : reduction to observability forms
C               = 10: reduction to a special controllability form
C               = 11: reduction to generalized Hessenberg form (*)
C               = 12: reduction to generalized Schur form (*)
C   A,B,
C   C,E     - descriptor system matrices.
C   job     - (optional) option parameter
C             for meth = 1, scaling option
C                = 0 - scale using all matrices A, E, B, C (default)
C                = 1 - scale using only A, E, B
C                = 2 - scale using only A, E, C
C                = 3 - scale using only A, E;
C             for meth = 8, controllability form option
C                = 0 - complete controllability form (default)
C                = 1 - finite controllability form
C                = 2 - infinite controllability form;
C             for meth = 9, observability form option
C                = 0 - complete observability form (default)
C                = 1 - finite observability form
C                = 2 - infinite observability form.
C   m1      - (optional) number of columns of matrix B1 in the
C             partition of B as B = [B1 B2].
C   tol     - (optional) tolerance parameter
C             for meth = 4, 5, 6, 7 tolerance for determining the
C                rank of E and of submatrix A22 (see TG01FD routine);
C                if meth = 4, 5 1/tol is a lower bound for the
C                reciprocal condition numbers of the triangular
C                matrices from the QR decomposition of E and A22;
C                if meth = 6, 7 1/tol is a lower bound for the
C                reciprocal condition numbers of the diagonal
C                matrices from the SVD decomposition of E and A22;
C             for meth = 8, 10 tolerance for controllability tests;
C             for meth = 9 tolerance for observability tests.
C   thresh  - (optional) threshold value for scaling.
C
C Output parameters:
C   At, Et,
C   Bt, Ct  - matrices of the transformed descriptor system
C               (Q'*A*Z-lambda Q'*E*Z, Q'*B, C*Z).
C   Q,Z     - orthogonal transformation matrices.
C   ranks   - for meth = 4: rank of matrix E;
C             for meth = 5: vector [ranke ranka22], where ranke
C                           is the rank of E and ranka22 is the rank
C                           of submatrix A22 (see TG01FD routine);
C             for meth = 6: rank of matrix E;
C             for meth = 7: vector [ranke ranka22], where ranke
C                           is the rank of E and ranka22 is the rank
C                           of submatrix A22 (see TG01FD routine).
C   dims    - for meth = 8: vector [ncont,niucon], where :
C                           if job = 0, ncont is the order of
C                           controllable subsystem and niucon is the
C                           number of uncontrollable infinite
C                           eigenvalues (see TG01HD routine);
C                           if job = 1, ncont is the order of finite
C                           controllable subsystem and niucon = 0;
C                           if job = 2, ncont is the order of a
C                           susbsystem having controllable infinite
C                           and non-zero eigenvalues and niucon = 0.
C           - for meth = 9: vector [nobs,niuobs], where:
C                           if job = 0, nobs is the order of
C                           observable subsystem and niuobs is the
C                           number of unobservable infinite
C                           eigenvalues (see TG01ID routine);
C                           if job = 1, nobs is the order of finite
C                           observable subsystem and niuobs = 0;
C                           if job = 2, nobs is the order of a
C                           susbsystem having observable infinite
C                           and non-zero eigenvalues and niuobs = 0.
C   ncont   - for meth =10: the order of the controllable subsystem.
C   tau     - for meth = 8: vector containing the dimension of full
C                           row rank matrices in the resulting
C                           staircase form (see TB01UD and TG01HD);
C           - for meth = 9: vector containing the dimension of full
C                           column rank matrices in the resulting
C                           staircase form (see TG01ID);
C           - for meth =10: vector containing the dimension of full
C                           row rank matrices in the resulting
C                           staircase form (see TB01UY).
C
C Contributor
C   A. Varga, April 1999.
C
C Revisions:
C   A. Varga, March 2002-2004, 13-01-2016.
C   V. Sima, Dec. 2016, Feb. 2018.
C
C **********************************************************************
C
#include "fintrf.h"
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C .. Parameters ..
      IMPLICIT NONE
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         ( ONE = 1.0D0, ZERO = 0.0D0 )
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
      CHARACTER         COMPQ, COMPZ, JOBA, JOBCO, JOBS
      INTEGER           INFO, L, LDA, LDB, LDC, LDE, LDQ, LDWORK, LDZ,
     $                  M, M1, M2, N, NCONT, NIUCON, NRBLCK, P, RANKE,
     $                  RNKA22
      DOUBLE PRECISION  MAXRED, TOL
C
C .. Allocatable local arrays ..
C !Fortran 90/95
      DOUBLE PRECISION, ALLOCATABLE:: A(:,:), B(:,:), C(:,:), DWORK(:),
     $                                E(:,:), Q(:,:), Z(:,:)
      INTEGER, ALLOCATABLE::          IWORK(:)
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      LOGICAL           SSTYPE
      INTEGER           ILEN, I, JOBN, L1, L2, LB, LC, LQ, LZ, METH, N1,
     $                  N3, NR, NRSLIM
      MWSIZE            IW, L1W, L2W, LIWORK, LW, MW, N1W, N3W, NRBLCKW,
     $                  NRW, NW, ONEW, PW
      DOUBLE PRECISION  DUM(2)
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
      EXTERNAL          DCOPY, DLACPY, DLASET, TB01ID, TB01UD, TB01UX,
     $                  TB01UY, TB01WX, TG01AD, TG01BD, TG01CD, TG01DD,
     $                  TG01ED, TG01FD, TG01HD, TG01HU, TG01ID
C
C .. Intrinsic functions ..
      INTRINSIC         MAX, MIN
C
      ONEW = 1
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
      IF( METH.LE.0 .OR. METH.GT.11 ) THEN
         CALL mexErrMsgTxt
     $        ( 'METH has 1 ... 11 the only admissible values' )
      END IF
C
C Check for proper number of arguments
C
      IF( METH.EQ.2 .OR. METH.EQ.3 ) THEN
         IF( NRHS.LT.4 ) THEN
            CALL mexErrMsgTxt
     $           ( 'SL_GSTRA requires at least 4 input arguments' )
         ELSE IF( NLHS.GT.4 ) THEN
            CALL mexErrMsgTxt
     $           ( 'SL_GSTRA has at most 4 output arguments' )
         END IF
      ELSE IF( METH.EQ.8 .OR. METH.EQ.9 .OR. METH.EQ.10 ) THEN
         IF( NRHS.LT.5 ) THEN
            CALL mexErrMsgTxt
     $           ( 'SL_GSTRA requires at least 5 input arguments' )
         ELSE IF( NLHS.GT.8 ) THEN
            CALL mexErrMsgTxt
     $           ( 'SL_GSTRA has at most 8 output arguments' )
         END IF
      ELSE IF( METH.EQ.11 ) THEN
         IF( NRHS.NE.5 ) THEN
            CALL mexErrMsgTxt
     $           ( 'SL_GSTRA requires 5 input arguments' )
         ELSE IF( NLHS.GT.6 ) THEN
            CALL mexErrMsgTxt
     $           ( 'SL_GSTRA has at most 6 output arguments' )
         END IF
      ELSE
         IF( NRHS.LT.5 ) THEN
            CALL mexErrMsgTxt
     $           ( 'SL_GSTRA requires at least 5 input arguments' )
         ELSE IF( METH.EQ.1 ) THEN
            IF( NLHS.GT.6 ) THEN
               CALL mexErrMsgTxt
     $              ( 'SL_GSTRA has at most 6 output arguments' )
            END IF
         ELSE IF( NLHS.GT.7 ) THEN
            CALL mexErrMsgTxt
     $           ( 'SL_GSTRA has at most 7 output arguments' )
         END IF
      END IF
C
C   A(LxN), E(LxN), B(LxM), C(PxN)
C
      LW  = mxGetM( PRHS(2) )
      L   = LW
      NW  = mxGetN( PRHS(2) )
      N   = NW
      L1W = mxGetM( PRHS(3) )
      L1  = L1W
      N1W = mxGetN( PRHS(3) )
      N1  = N1W
      IF( METH.EQ.2 ) THEN
         MW = mxGetN( PRHS(4) )
         M  = MW
         P  = 0
         PW = P
         L2W = mxGetM( PRHS(4) )
         L2 = L2W
         N3 = 0
         N3W = N3
      ELSE IF( METH.EQ.3 ) THEN
         M  = 0
         MW = M
         PW = mxGetM( PRHS(4) )
         P  = PW
         L2 = 0
         L2W = L2
         N3W = mxGetN( PRHS(4) )
         N3 = N3W
      ELSE
         MW  = mxGetN( PRHS(4) )
         M   = MW
         PW  = mxGetM( PRHS(5) )
         P   = PW
         L2W = mxGetM( PRHS(4) )
         L2  = L2W
         N3W = mxGetN( PRHS(5) )
         N3  = N3W
      END IF
C
      SSTYPE = L1.EQ.0 .AND. N1.EQ.0
      IF( .NOT.SSTYPE .AND. ( (N1.NE.N) .OR. (L1.NE.L) ) ) THEN
         CALL mexErrMsgTxt
     $        ( 'A and E must have the same dimensions' )
      END IF
      IF( (METH.GE.8 .AND. METH.LE.11) .AND. L.NE.N ) THEN
         CALL mexErrMsgTxt
     $        ( 'A and E must be square matrices' )
      END IF
      IF( M.GT.0 .AND. L2.NE.L ) THEN
         CALL mexErrMsgTxt
     $        ( 'B must have the same row dimension as A and E' )
      END IF
      IF( P.GT.0 .AND. N3.NE.N ) THEN
         CALL mexErrMsgTxt
     $        ( 'C must have the same column dimension as A and E' )
      END IF
      IF( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a real matrix' )
      END IF
      IF( .NOT.SSTYPE .AND. ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(3) ).EQ.1 ) ) THEN
         CALL mexErrMsgTxt( 'E must be a real matrix' )
      END IF
      IF( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(4) ).EQ.1 ) THEN
         IF( METH.EQ.3 ) THEN
            CALL mexErrMsgTxt( 'C must be a real matrix' )
         ELSE
            CALL mexErrMsgTxt( 'B must be a real matrix' )
         END IF
      END IF
      IF( METH.NE.2 .AND. METH.NE.3 ) THEN
         IF( mxIsNumeric( PRHS(5) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(5) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'C must be a real matrix' )
        END IF
      END IF
C
C   tol, thresh
C
      TOL = ZERO
      IF( METH.EQ.1 .OR. (METH.GE.4 .AND. METH.LE.10) ) THEN
         NRSLIM = 6
         IF( METH.EQ.1 .OR. METH.GT.7 ) NRSLIM = 7
         IF( NRHS.GE.NRSLIM ) THEN
            IF( mxGetM( PRHS(NRSLIM) )*mxGetN( PRHS(NRSLIM) ).GT.1 )
     $         THEN
               IF( METH.NE.1) THEN
                  CALL mexErrMsgTxt( 'TOL must be scalar' )
               ELSE
                  CALL mexErrMsgTxt( 'THRESH must be scalar' )
               END IF
            END IF
            IF( mxIsNumeric( PRHS(NRSLIM) ).EQ.0 .OR.
     $          mxIsComplex( PRHS(NRSLIM) ).EQ.1 ) THEN
               IF( METH.NE.1) THEN
                  CALL mexErrMsgTxt( 'TOL must be a real scalar' )
               ELSE
                  CALL mexErrMsgTxt( 'THRESH must be a real scalar' )
               END IF
            END IF
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(NRSLIM) ), TOL, ONEW )
         END IF
      END IF
C
C     job
C
      JOBN = 0
      IF( METH.EQ.1 .OR. METH.EQ.8 .OR. METH.EQ.9 .OR. METH.EQ.10 ) THEN
         IF( NRHS .GT. 5 ) THEN
            ILEN = mxGetM( PRHS(6) )*mxGetN( PRHS(6) )
            IF( ILEN.GT.1 ) THEN
               IF( METH.EQ.10 ) THEN
                  CALL mexErrMsgTxt( 'M1 must be scalar' )
               ELSE
                  CALL mexErrMsgTxt( 'JOB must be scalar' )
               END IF
            END IF
            IF( mxIsNumeric( PRHS(6) ).EQ.0 .OR.
     $          mxIsComplex( PRHS(6) ).EQ.1 ) THEN
               IF( METH.EQ.10 ) THEN
                  CALL mexErrMsgTxt( 'M1 must be an integer scalar' )
               ELSE
                  CALL mexErrMsgTxt( 'JOB must be an integer scalar' )
               END IF
            END IF
            IF( ILEN.EQ.0 ) THEN
               JOBN = 0
            ELSE
               CALL mxCopyPtrToReal8( mxGetPr( PRHS(6) ), DUM, ONEW )
               JOBN = DUM(1)
            END IF
            IF( METH.EQ.1 ) THEN
               IF( JOBN.LT.0 .OR. JOBN.GT.3 ) THEN
                  CALL mexErrMsgTxt
     $                 ( 'JOB has 0 ... 3 the only admissible values' )
               END IF
            ELSE IF( METH.EQ.10 ) THEN
               IF( JOBN.LT.0 ) THEN
                  CALL mexErrMsgTxt
     $                 ( 'M1 must be non-negative' )
               END IF
            ELSE
               IF( JOBN.LT.0 .OR. JOBN.GT.2 ) THEN
                  CALL mexErrMsgTxt
     $                 ( 'JOB has 0 ... 2 the only admissible values' )
               END IF
            END IF
         END IF
      END IF
C
C Determine the lengths of working arrays
C
      LDA = MAX( 1, L )
      LDE = MAX( 1, L )
      LDB = MAX( 1, L )
      LDC = MAX( 1, P )
      IF( METH.EQ.9 ) LDC = MAX( LDC, M )
      LB = MAX( 1, M )
      LC = MAX( 1, N )
      IF( METH.EQ.1 ) THEN
         IF( NLHS.LT.5 ) THEN
            LQ = L
            LDQ = 1
         ELSE
            LQ = L
            LDQ = MAX( 1, L )
         END IF
         IF( NLHS.LT.6 ) THEN
            LZ = N
            LDZ = 1
         ELSE
            LZ = N
            LDZ = MAX( 1, N )
         END IF
      ELSE IF( METH.EQ.2 ) THEN
         IF( NLHS.LT.4 ) THEN
            LQ = 1
            COMPQ = 'N'
            LDQ = 1
         ELSE
            LQ = L
            COMPQ = 'I'
            LDQ = MAX( 1, L )
         END IF
         LC = 1
         LDZ = 1
         LZ = 1
      ELSE IF( METH.EQ.3 ) THEN
         IF( NLHS.LT.4 ) THEN
            LZ = 1
            COMPZ = 'N'
            LDZ = 1
         ELSE
            LZ = N
            COMPZ = 'I'
            LDZ = MAX( 1, N )
         END IF
         LB = 1
         LDQ = 1
         LQ = 1
      ELSE IF( METH.EQ.4 .OR. METH.EQ.5 ) THEN
         LQ = L
         LDQ = MAX( 1, L )
         LZ = N
         LDZ = MAX( 1, N )
      ELSE IF( METH.EQ.6 .OR. METH.EQ.7 ) THEN
         IF( NLHS.LT.6 ) THEN
            LQ = 1
            COMPQ = 'N'
            LDQ = 1
         ELSE
            LQ = L
            COMPQ = 'I'
            LDQ = MAX( 1, L )
         END IF
         IF( NLHS.LT.7 ) THEN
            LZ = 1
            COMPZ = 'N'
            LDZ = 1
         ELSE
            LZ = N
            COMPZ = 'I'
            LDZ = MAX( 1, N )
         END IF
      ELSE IF( METH.GE.8 .AND. METH.LE.10 ) THEN
         IF( NLHS.LT.7 ) THEN
            LQ = 1
            COMPQ = 'N'
            LDQ = 1
         ELSE
            LQ = L
            COMPQ = 'I'
            LDQ = MAX( 1, L )
         END IF
         IF( NLHS.LT.8 ) THEN
            LZ = 1
            COMPZ = 'N'
            LDZ = 1
         ELSE
            LZ = N
            COMPZ = 'I'
            LDZ = MAX( 1, N )
         END IF
         IF( METH.EQ.9 ) LB = MAX( LB, P )
      ELSE IF( METH.EQ.11 ) THEN
         IF( NLHS.LT.5 ) THEN
            LQ = 1
            COMPQ = 'N'
            LDQ = 1
         ELSE
            LQ = L
            COMPQ = 'I'
            LDQ = MAX( 1, L )
         END IF
         IF( NLHS.LT.6 ) THEN
            LZ = 1
            COMPZ = 'N'
            LDZ = 1
         ELSE
            LZ = N
            COMPZ = 'I'
            LDZ = MAX( 1, N )
         END IF
      END IF
C
      ALLOCATE ( E(LDE,MAX( 1, N )) )
      IF( SSTYPE ) THEN
         CALL DLASET( 'F', L, N, ZERO, ONE, E, LDE )
      ELSE
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), E, LW*NW )
         SSTYPE = L.EQ.N .AND. CHCKID( MIN( L, N ), E, LDE )
      END IF
C
      INFO   = 0
      LIWORK = 0
      IF( SSTYPE ) THEN
         LDWORK = 0
      ELSE
         IF( METH.EQ.1 ) THEN
            LDWORK = 3*( L + N )
         ELSE IF( METH.EQ.2 ) THEN
            LDWORK = MAX( 1, MIN( L, N ) + MAX( L, N, M ) )
         ELSE IF( METH.EQ.3 ) THEN
            LDWORK = MAX( 1, MIN( L, N ) + MAX( L, N, P ) )
         ELSE IF( METH.EQ.4 .OR. METH.EQ.5 ) THEN
            LDWORK = MAX( 1, MIN( L, N ) + MAX( M, P, 5*MIN( L, N ),
     $                     3*MIN( L, N ) + MAX( L, N ) ) )
         ELSE IF( METH.EQ.6 .OR. METH.EQ.7 ) THEN
            IF( METH.EQ.6 ) THEN
               JOBA = 'N'
            ELSE
               JOBA = 'R'
            END IF
            CALL TG01FD( COMPQ, COMPZ, JOBA, L, N, M, P, A, LDA, E, LDE,
     $                   B, LDB, C, LDC, Q, LDQ, Z, LDZ, RANKE, RNKA22,
     $                   TOL, IWORK, DUM, -1, INFO )
            LDWORK = DUM(1)
            LIWORK = N
         END IF
      END IF
      IF( METH.EQ.8 ) THEN
         IF( SSTYPE .OR. N.EQ.1 ) THEN
            LDWORK = N + MAX( 1, N, 3*M, P )
         ELSE
            LDWORK = MAX( N, 2*M )
         END IF
         LIWORK = N + M
      ELSE IF( METH.EQ.9 ) THEN
         IF( SSTYPE ) THEN
            LDWORK = N + MAX( 1, N, 3*P, M )
         ELSE
            LDWORK = MAX( N, 2*P )
         END IF
         LIWORK = N + P
      ELSE IF( METH.EQ.10 ) THEN
         M1 = MIN( JOBN, M )
         M2 = M - M1
         IF( SSTYPE ) THEN
            CALL TB01UY( COMPQ, N, M1, M2, P, A, LDA, B, LDB, C, LDC,
     $                   NCONT, NRBLCK, IWORK, Q, LDQ, DWORK, TOL,
     $                   IWORK, DUM, -1, INFO )
            LIWORK = 2*N + MAX( M1, M2 )
            LDWORK = N + DUM(1)
         ELSE
            CALL TG01HU( COMPQ, COMPZ, N, N, M1, M2, P, N,
     $                   MAX( 0, N-1 ), A, LDA,  E, LDE, B, LDB, C, LDC,
     $                   Q, LDQ, Z, LDZ, NCONT, NRBLCK, IWORK, TOL,
     $                   IWORK, DUM, -1, INFO )
            LIWORK = 2*N + M
            LDWORK = DUM(1)
         END IF
      ELSE IF( METH.EQ.11 ) THEN
         IF( SSTYPE ) THEN
            CALL TB01WX( COMPQ, N, M, P, A, LDA, B, LDB, C, LDC, Q, LDQ,
     $                   DUM, -1, INFO )
            LDWORK = DUM(1)
         ELSE
            LDWORK = N + MAX( N, M )
         END IF
      END IF
C
C Allocate variable dimension local arrays
C !Fortran 90/95
C
      ALLOCATE ( A(LDA,MAX( 1, N )), B(LDB,LB), C(LDC,LC),
     $           Q(LDQ,LQ), Z(LDZ,LZ), DWORK(LDWORK), IWORK(LIWORK) )
C
C Copy inputs from MATLAB workspace to locally allocated arrays
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), A, LW*NW )
      IF( METH.EQ.2 ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), B, LW*MW )
      ELSE IF( METH.EQ.3 ) THEN
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), C, PW*NW )
      ELSE
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), B, LW*MW )
         IF( METH.EQ.9 .AND. M.GT.P ) THEN
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), DWORK, PW*NW )
            CALL DLACPY( 'F', P, N, DWORK, MAX( 1, P ), C, LDC )
         ELSE
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), C, PW*NW )
         END IF
      END IF
C
C Do the actual computations
C
      IF( METH.EQ.1 ) THEN
C
C        Scaling
C
         IF( JOBN.EQ.0 ) THEN
            JOBS = 'A'
         ELSE IF( JOBN.EQ.1 ) THEN
            JOBS = 'B'
         ELSE IF( JOBN.EQ.2 ) THEN
            JOBS = 'C'
         ELSE
            JOBS = 'N'
         END IF
         IF( SSTYPE ) THEN
            MAXRED = ZERO
            CALL TB01ID( JOBS, N, M, P, MAXRED, A, LDA, B, LDB, C, LDC,
     $                   Z, INFO )
         ELSE
            CALL TG01AD( JOBS, L, N, M, P, TOL, A, LDA, E, LDE, B, LDB,
     $                   C, LDC, Q, Z, DWORK, INFO )
         END IF
      END IF
C
      IF( SSTYPE) THEN
         RANKE = N
         RNKA22 = 0
      ELSE
C
C        QR-coordinate form
C
         IF( METH.EQ.2 )
     $      CALL TG01CD( COMPQ, L, N, M, A, LDA, E, LDE, B, LDB, Q, LDQ,
     $                   DWORK, LDWORK, INFO )
C
C        RQ-coordinate form
C
         IF( METH.EQ.3 )
     $      CALL TG01DD( COMPZ, L, N, P, A, LDA, E, LDE, C, LDC, Z, LDZ,
     $                   DWORK, LDWORK, INFO )
C
C        SVD-coordinate form
C
         IF( METH.EQ.4 )
     $      CALL TG01ED( 'N', L, N, M, P, A, LDA, E, LDE,
     $                   B, LDB, C, LDC, Q, LDQ, Z, LDZ,
     $                   RANKE, RNKA22, TOL, DWORK, LDWORK, INFO )
C
C        SVD-coordinate form with A22 further reduced.
C
         IF( METH.EQ.5 )
     $      CALL TG01ED( 'R', L, N, M, P, A, LDA, E, LDE,
     $                   B, LDB, C, LDC, Q, LDQ, Z, LDZ,
     $                   RANKE, RNKA22, TOL, DWORK, LDWORK, INFO )
C
C        SVD-like coordinate form
C
         IF( METH.EQ.6 )
     $      CALL TG01FD( COMPQ, COMPZ, 'N', L, N, M, P, A, LDA, E, LDE,
     $                   B, LDB, C, LDC, Q, LDQ, Z, LDZ, RANKE, RNKA22,
     $                   TOL, IWORK, DWORK, LDWORK, INFO )
C
C        SVD-like coordinate form with A22 further reduced.
C
         IF( METH.EQ.7 )
     $      CALL TG01FD( COMPQ, COMPZ, 'R', L, N, M, P, A, LDA, E, LDE,
     $                   B, LDB, C, LDC, Q, LDQ, Z, LDZ, RANKE, RNKA22,
     $                   TOL, IWORK, DWORK, LDWORK, INFO )
      END IF
C
      IF( METH.EQ.8 ) THEN
C
C        Controllability forms.
C
         IF( JOBN.EQ.0 ) THEN
            JOBCO = 'C'
         ELSE IF( JOBN.EQ.1 ) THEN
            JOBCO = 'F'
         ELSE
            JOBCO = 'I'
         END IF
         IF( SSTYPE ) THEN
            CALL TB01UD( COMPQ, N, M, P, A, LDA, B, LDB, C, LDC, NCONT,
     $                   NRBLCK, IWORK, Q, LDQ, DWORK, TOL, IWORK(N+1),
     $                   DWORK(N+1), LDWORK-N+1, INFO )
            NIUCON = 0
         ELSE
            CALL TG01HD( JOBCO, COMPQ, COMPZ, N, M, P, A, LDA, E, LDE,
     $                   B, LDB, C, LDC, Q, LDQ, Z, LDZ, NCONT, NIUCON,
     $                   NRBLCK, IWORK, TOL, IWORK(N+1), DWORK, INFO )
         END IF
      END IF
C
      IF( METH.EQ.9 ) THEN
C
C        Observability forms.
C
         IF( JOBN.EQ.0 ) THEN
            JOBCO = 'O'
         ELSE IF( JOBN.EQ.1 ) THEN
            JOBCO = 'F'
         ELSE
            JOBCO = 'I'
         END IF
         IF( SSTYPE ) THEN
            CALL TB01UX( COMPQ, N, M, P, A, LDA, B, LDB, C, LDC,
     $                   Q, LDQ, NCONT, NRBLCK, IWORK, TOL, IWORK(N+1),
     $                   DWORK, INFO )
            NIUCON = 0
         ELSE
            CALL TG01ID( JOBCO, COMPQ, COMPZ, N, M, P, A, LDA, E, LDE,
     $                   B, LDB, C, LDC, Q, LDQ, Z, LDZ, NCONT, NIUCON,
     $                   NRBLCK, IWORK, TOL, IWORK(N+1), DWORK, INFO )
         END IF
      END IF
C
      IF( METH.EQ.10 ) THEN
C
C        Special controllability forms.
C
         M1 = MIN( JOBN, M )
         M2 = M - M1
         IF( SSTYPE ) THEN
            CALL TB01UY( COMPQ, N, M1, M2, P, A, LDA, B, LDB, C, LDC,
     $                   NCONT, NRBLCK, IWORK, Q, LDQ, DWORK, TOL,
     $                   IWORK(2*N+1), DWORK(N+1), LDWORK-N+1, INFO )
         ELSE
            CALL TG01HU( COMPQ, COMPZ, N, N, M1, M2, P, N,
     $                   MAX( 0, N-1 ), A, LDA,  E, LDE, B, LDB, C, LDC,
     $                   Q, LDQ, Z, LDZ, NCONT, NRBLCK, IWORK, TOL,
     $                   IWORK(2*N+1), DWORK, LDWORK, INFO )
         END IF
      END IF
C
      IF( METH.EQ.11 ) THEN
C
C        Hessenberg forms.
C
         IF( SSTYPE ) THEN
            CALL TB01WX( COMPQ, N, M, P, A, LDA, B, LDB, C, LDC, Q, LDQ,
     $                   DWORK, LDWORK, INFO )
         ELSE
            CALL TG01BD( 'G', COMPQ, COMPZ, N, M, P, 1, N, A, LDA,
     $                   E, LDE, B, LDB, C, LDC, Q, LDQ, Z, LDZ, DWORK,
     $                   LDWORK, INFO )
         END IF
      END IF
C
*      close(1)
C
C Copy output to MATLAB workspace
C
      IF( INFO.EQ.0 ) THEN
         LW = L
         MW = M
         NW = N
         PW = P
         IF( NLHS.GE.1 ) THEN
            PLHS(1) = mxCreateDoubleMatrix( LW, NW, SZERO )
            CALL mxCopyReal8ToPtr( A, mxGetPr( PLHS(1) ), LW*NW )
         END IF
         IF( NLHS.GE.2 ) THEN
            PLHS(2) = mxCreateDoubleMatrix( LW, NW, SZERO )
            CALL mxCopyReal8ToPtr( E, mxGetPr( PLHS(2) ), LW*NW )
         END IF
         IF( METH.EQ.2 ) THEN
            IF( NLHS.GE.3 ) THEN
               PLHS(3) = mxCreateDoubleMatrix( LW, MW, SZERO )
               CALL mxCopyReal8ToPtr( B, mxGetPr( PLHS(3) ), LW*MW )
            END IF
            IF( NLHS.GE.4 ) THEN
               PLHS(4) = mxCreateDoubleMatrix( LW, LW, SZERO )
               CALL mxCopyReal8ToPtr( Q, mxGetPr( PLHS(4) ), LW*LW )
            END IF
         ELSE IF( METH.EQ.3 ) THEN
            IF( NLHS.GE.3 ) THEN
               PLHS(3) = mxCreateDoubleMatrix( PW, NW, SZERO )
               CALL mxCopyReal8ToPtr(C, mxGetPr( PLHS(3) ), PW*NW )
            END IF
            IF( NLHS.GE.4 ) THEN
               PLHS(4) = mxCreateDoubleMatrix( NW, NW, SZERO )
               CALL mxCopyReal8ToPtr( Z, mxGetPr( PLHS(4) ), NW*NW )
            END IF
         ELSE IF( METH.EQ.1 .OR. METH.GE.4 ) THEN
            IF( NLHS.GE.3 ) THEN
               PLHS(3) = mxCreateDoubleMatrix( LW, MW, SZERO )
               CALL mxCopyReal8ToPtr( B, mxGetPr( PLHS(3) ), LW*MW )
            END IF
            IF( NLHS.GE.4 ) THEN
               PLHS(4) = mxCreateDoubleMatrix( PW, NW, SZERO )
               CALL mxCopyReal8ToPtr( C(1:P,1:N), mxGetPr( PLHS(4) ),
     $                                PW*NW )
            END IF
            IF( METH.EQ.1 ) THEN
               IF( NLHS.GE.5 ) THEN
                  IF( SSTYPE ) THEN
                     DO I = 1, L
                        DWORK(I) = ONE / Z(I,1)
                     END DO
                  ELSE
                     CALL DCOPY( L, Q, 1, DWORK, 1 )
                  END IF
                  CALL DLASET( 'F', L, L, ZERO, ZERO, Q, LDQ )
                  CALL DCOPY( L, DWORK, 1, Q, LDQ+1 )
                  PLHS(5) = mxCreateDoubleMatrix( LW, LW, SZERO )
                  CALL mxCopyReal8ToPtr( Q, mxGetPr( PLHS(5) ), LW*LW )
               END IF
               IF( NLHS.GE.6 ) THEN
                  CALL DCOPY( N, Z, 1, DWORK, 1 )
                  CALL DLASET( 'F', N, N, ZERO, ZERO, Z, LDZ )
                  CALL DCOPY( N, DWORK, 1, Z, LDZ+1 )
                  PLHS(6) = mxCreateDoubleMatrix( NW, NW, SZERO )
                  CALL mxCopyReal8ToPtr( Z, mxGetPr( PLHS(6) ), NW*NW )
               END IF
            ELSE IF( METH.GE.4 .AND. METH.LE.7 ) THEN
               IF( NLHS.GE.5 ) THEN
                  NR = 1
                  DUM(1) = RANKE
                  IF( METH.EQ.5 .OR. METH.EQ.7 ) THEN
                     NR = 2
                     DUM(2) = RNKA22
                  END IF
                  NRW = NR
                  PLHS(5) = mxCreateDoubleMatrix( NRW, ONEW, SZERO )
                  CALL mxCopyReal8ToPtr( DUM, mxGetPr( PLHS(5) ), NRW )
               END IF
               IF( NLHS.GE.6 ) THEN
                  PLHS(6) = mxCreateDoubleMatrix( LW, LW, SZERO )
                  IF( SSTYPE )
     $               CALL DLASET( 'F', L, L, ZERO, ONE, Q, LDQ )
                  CALL mxCopyReal8ToPtr( Q, mxGetPr( PLHS(6) ), LW*LW )
               END IF
               IF( NLHS.GE.7 ) THEN
                  PLHS(7) = mxCreateDoubleMatrix( NW, NW, SZERO )
                  IF( SSTYPE )
     $               CALL DLASET( 'F', N, N, ZERO, ONE, Z, LDZ )
                  CALL mxCopyReal8ToPtr( Z, mxGetPr( PLHS(7) ), NW*NW )
               END IF
            ELSE IF( METH.EQ.11 ) THEN
               IF( NLHS.GE.5 ) THEN
                  PLHS(5) = mxCreateDoubleMatrix( NW, NW, SZERO )
                  CALL mxCopyReal8ToPtr( Q, mxGetPr( PLHS(5) ), NW*NW )
               END IF
               IF( NLHS.GE.6 ) THEN
                  PLHS(6) = mxCreateDoubleMatrix( NW, NW, SZERO )
                  IF( SSTYPE ) THEN
                     CALL mxCopyReal8ToPtr( Q, mxGetPr( PLHS(6) ),
     $                                      NW*NW )
                  ELSE
                     CALL mxCopyReal8ToPtr( Z, mxGetPr( PLHS(6) ),
     $                                      NW*NW )
                  END IF
               END IF
            ELSE
               IF( NLHS.GE.5 ) THEN
                  DUM(1) = NCONT
                  I = 2
                  IF( METH.EQ.10 ) THEN
                     I = 1
                  ELSE
                     I = 2
                     DUM(2) = NIUCON
                  END IF
                  IW = I
                  PLHS(5) = mxCreateDoubleMatrix( IW, ONEW, SZERO )
                  CALL mxCopyReal8ToPtr( DUM, mxGetPr( PLHS(5) ), IW )
               END IF
               IF( NLHS.GE.6 ) THEN
                  DO 10 I = 1, NRBLCK
                     DWORK(I) = IWORK(I)
   10             CONTINUE
                  NRBLCKW = NRBLCK
                  PLHS(6) = mxCreateDoubleMatrix( NRBLCKW, ONEW, SZERO )
                  CALL mxCopyReal8ToPtr
     $                 ( DWORK, mxGetPr( PLHS(6) ), NRBLCKW )
               END IF
               IF( NLHS.GE.7 ) THEN
                  PLHS(7) = mxCreateDoubleMatrix( LW, LW, SZERO )
                  CALL mxCopyReal8ToPtr( Q, mxGetPr( PLHS(7) ), LW*LW )
               END IF
               IF( NLHS.GE.8 ) THEN
                  PLHS(8) = mxCreateDoubleMatrix( NW, NW, SZERO )
                  IF( SSTYPE ) THEN
                     CALL mxCopyReal8ToPtr( Q, mxGetPr( PLHS(8) ),
     $                                      NW*NW )
                  ELSE
                     CALL mxCopyReal8ToPtr( Z, mxGetPr( PLHS(8) ),
     $                                      NW*NW )
                  END IF
               END IF
            END IF
         END IF
C
      END IF
C Deallocate local arrays
C !Fortran 90/95
C
      DEALLOCATE ( A, B, C, E, Q, Z, DWORK, IWORK )
C
C Error handling
C
      IF( INFO.NE.0 ) THEN
         IF( METH.EQ.1 ) THEN
           write( TEXT,'( " INFO =", I4, " ON EXIT FROM TGO1AD" )' )
     $            INFO
         ELSE IF( METH.EQ.2 ) THEN
           write( TEXT,'( " INFO =", I4, " ON EXIT FROM TGO1CD" )' )
     $            INFO
         ELSE IF( METH.EQ.3 ) THEN
           write( TEXT,'( " INFO =", I4, " ON EXIT FROM TGO1DD" )' )
     $            INFO
         ELSE IF( METH.EQ.4 .OR. METH.EQ.5 ) THEN
           write( TEXT,'( " INFO =", I4, " ON EXIT FROM TGO1ED" )' )
     $            INFO
         ELSE IF( METH.EQ.6 .OR. METH.EQ.7 ) THEN
           write( TEXT,'( " INFO =", I4, " ON EXIT FROM TGO1FD" )' )
     $            INFO
         ELSE IF( METH.EQ.8 ) THEN
           write( TEXT,'( " INFO =", I4, " ON EXIT FROM TGO1HD" )' )
     $            INFO
         ELSE IF( METH.EQ.9 ) THEN
           write( TEXT,'( " INFO =", I4, " ON EXIT FROM TGO1ID" )' )
     $            INFO
         ELSE IF( METH.EQ.10 ) THEN
           write( TEXT,'( " INFO =", I4, " ON EXIT FROM TBO1UY" )' )
     $            INFO
         ELSE IF( METH.EQ.11 ) THEN
           write( TEXT,'( " INFO =", I4, " ON EXIT FROM TBO1BD" )' )
     $            INFO
         END IF
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of SL_GSTRA ***
      END
