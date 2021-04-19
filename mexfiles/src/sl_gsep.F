C SL_GSEP.F - Gateway function for SLICOT system transformation routines
C             TG01LD, TG01NX, TG01QD
C             TB01KX, TB01LD, TB01WD, TB01WX.
C
C MATLAB call:
C   [At,Et,Bt,Ct,dims,tau,Q,Z] = sl_gsep(meth,A,E,B,C,job,tol,disc,sdeg)
C
C Purpose:
C   To compute several similarity transformations on a
C   a given descriptor system (A-lambda*E,B,C) to separate the
C   spectrum of the pair (A,E).
C
C Input parameters:
C   meth    - type of system transformation to be performed
C               = 1 : finite-infinite separation
C               = 2 : reduction to Hessenberg form with
C                     finite-infinite separation
C               = 3 : reduction to generalized Schur form with
C                     finite-infinite separation
C               = 4 : reduction to block diagonal generalized Schur
C                     form with finite-infinite separation
C               = 5:  reduction to generalized Schur form with
C                     ordered finite-infinite separation
C               = 6 : reduction to block diagonal generalized Schur
C                     form with ordered finite-infinite separation
C   A,B,
C   C,E     - descriptor system matrices
C   job     - (optional) option parameter
C              =  0 - for meth = 1, 2, 3, or 4, the leading block
C                     contains the infinite spectrum and trailing block
C                     contains the finite spectrum (default);
C                     for meth = 5 or 6, the leading block contains the
C                     infinite spectrum, the midle block contains the
C                     unstable finite spectrum and the trailing block
C                     contains the stable finite spectrum (default);
C              =  1 - for meth = 1, 2, 3, or 4, the leading block
C                     contains the finite spectrum and trailing block
C                     contains the infinite spectrum;
C                     for meth = 5 or 6, the leading block contains the
C                     stable finite spectrum, the midle block contains
C                     the unstable finite spectrum, and the trailing
C                     block contains the infinite spectrum;
C              =  2 - for meth = 1, 2, 3, or 4, the leading block
C                     contains the infinite spectrum and trailing block
C                     contains the finite spectrum;
C                     for meth = 5 or 6, the leading block contains the
C                     infinite spectrum, the midle block contains the
C                     stable finite spectrum, and the trailing block
C                     contains the unstable finite spectrum;
C              =  3 - for meth = 1, 2, 3, or 4, the leading block
C                     contains the finite spectrum and trailing block
C                     contains the infinite spectrum;
C                     for meth = 5 or 6, the leading block contains the
C                     unstable finite spectrum, the midel block contains
C                     the stable finite spectrum and the trailing block
C                     contains the infinite spectrum;
C   tol     - (optional) tolerance for determining the rank of submatrices
C             in the system matrices A and E; generally, 1/tol is a
C             lower bound for the reciprocal condition numbers used for
C             the triangular matrices in the rank revealing QR decompositions.
C   disc    - (optional) continuous/disrete flag
C             = 0 - continuous-time system (default)
C             = 1 - discrete-time system
C   sdeg    - (optional) desired stability degree for the 'stable' part
C              default:  -sqrt(eps) for continuous-time
C                       1-sqrt(eps) for discrete-time
C
C Output parameters:
C   At, Et,
C   Bt, Ct  - matrices of the transformed descriptor system
C               (Q*A*Z-lambda Q*E*Z, Q*B, C*Z)
C   dims    - vector [n1,n2,n3,rankE], where :
C             n1 is the order of the leading diagonal block,
C             n2 is the order of the midle diagonal block,
C             n3 is the order of the trailing diagonal block, and
C             rankE is the rank of E;
C   tau     - vector containing the dimensions of diagonal blocks
C             of the subpencil containing infinite eigenvalues
C   Q,Z     - transformation matrices (orthogonal for meth <= 3)
C
C Contributor
C   A. Varga, April 2002.
C
C Revisions:
C   A. Varga,  13-01-2016.
C   V. Sima, Dec. 2016.
C
C **********************************************************************
C
#include "fintrf.h"
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C .. Parameters ..
      IMPLICIT NONE
      DOUBLE PRECISION  ONE, ZERO, THREE
      PARAMETER         ( ONE = 1.0D0, ZERO = 0.0D0, THREE = 3.0D0 )
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
      CHARACTER         COMPQ, COMPZ, DICO, JOB, JOBA, STDOM
      INTEGER           INFO, LDA, LDB, LDC, LDE, LDQ, LDWORK, LDZ, M,
     $                  N, NDIM, P
      MWSIZE            MW, NW, PW
      DOUBLE PRECISION  TOL
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
      INTEGER           I, JOBN, L1, L2, L3, LB, LC, LIWORK, LQ, LZ,
     $                  METH, N1, N2, N3, ND, NF, NIBLCK
      MWSIZE            L1W, L2W, L3W, N2W, N3W, NIBLCKW
      MWSIZE            ONEW, FOURW
      DOUBLE PRECISION  ALPHA, DUM(4)
C
C     .. matlab integer*4 flags
C
      INTEGER*4         SZERO
      PARAMETER         ( SZERO = 0 )
C
C .. External functions ..
      LOGICAL           CHCKID
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          CHCKID, DLAMCH
C
C .. External Subroutines ..
      EXTERNAL          DLASET, DSWAP, TB01KX, TB01LD, TB01WD, TB01WX,
     $                  TG01LD, TG01NX, TG01QD
C
C .. Intrinsic functions ..
      INTRINSIC         INT, MAX, SQRT
C
      ONEW  = 1
      FOURW = 4
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
      IF( METH.LE.0 .OR. METH.GT.6 ) THEN
         CALL mexErrMsgTxt
     $        ( 'METH has 1 ... 6 the only admissible values' )
      END IF
C
C Check for proper number of arguments
C
      IF( NRHS.LT.5 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SL_GSEP requires at least 5 input arguments' )
      ELSE IF( NRHS.GT.9 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SL_GSEP has at most 9 input argument' )
      ELSE IF( NLHS.GT.8 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SL_GSEP has at most 8 output argument' )
      END IF
C
C   A(LxN), E(LxN), B(LxM), C(PxN)
C
      L1W = mxGetM( PRHS(2) )
      L1  = L1W
      NW  = mxGetN( PRHS(2) )
      N   = NW
      L2W = mxGetM( PRHS(3) )
      L2  = L2W
      N2W = mxGetN( PRHS(3) )
      N2  = N2W
      L3W = mxGetM( PRHS(4) )
      L3  = L3W
      MW  = mxGetN( PRHS(4) )
      M   = MW
      PW  = mxGetM( PRHS(5) )
      P   = PW
      N3W = mxGetN( PRHS(5) )
      N3  = N3W
C
      SSTYPE = L2.EQ.0 .AND. N2.EQ.0
      IF( .NOT.SSTYPE .AND. ( (N.NE.N2) .OR. (L1.NE.L2) ) ) THEN
         CALL mexErrMsgTxt
     $        ( 'A and E must have the same dimensions' )
      END IF
      IF( L1.NE.N ) THEN
         CALL mexErrMsgTxt
     $        ( 'A and E must be square matrices' )
      END IF
      IF( M.GT.0 .AND. N.NE.L3 ) THEN
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
         CALL mexErrMsgTxt( 'B must be a real matrix' )
      END IF
      IF( mxIsNumeric( PRHS(5) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(5) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'C must be a real matrix' )
      END IF
C
C     job
C
      JOBN = 0
      IF( NRHS .GT. 5 ) THEN
         IF( mxGetM( PRHS(6) )*mxGetN( PRHS(6) ).GT.1 ) THEN
            CALL mexErrMsgTxt( 'JOB must be scalar' )
         END IF
         IF( mxIsNumeric( PRHS(6) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(6) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'JOB must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(6) ), DUM, ONEW )
         JOBN = DUM(1)
         IF( DUM(1).LT.ZERO .OR. DUM(1).GT.THREE ) THEN
            CALL mexErrMsgTxt
     $           ( 'JOB has 0 ... 3 the only admissible values' )
         END IF
      END IF
C
C   tol
C
      TOL = ZERO
      IF( NRHS.GT.6 ) THEN
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
C   disc
C
      DICO = 'C'
      DUM(1) = ZERO
      IF( NRHS.GT.7 ) THEN
         IF( mxGetM( PRHS(8) )*mxGetN( PRHS(8) ).GT.1 ) THEN
            CALL mexErrMsgTxt( 'DISC must be scalar' )
         END IF
         IF( mxIsNumeric( PRHS(8) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(8) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'DISC must be a real scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(8) ), DUM, ONEW )
         IF( DUM(1).NE.ZERO .AND. DUM(1).NE.ONE ) THEN
            CALL mexErrMsgTxt( 'DISC has 0 or 1 the only valid values' )
         END IF
         IF( DUM(1).NE.ZERO ) DICO = 'D'
      END IF
C
C   alpha
C
      IF( NRHS.GT.8 ) THEN
         IF( mxGetM( PRHS(9) ).NE.1 .OR. mxGetN( PRHS(9) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'ALPHA must be a scalar' )
         END IF
         IF( mxIsNumeric( PRHS(9) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(9) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'ALPHA must be a real scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(9) ), ALPHA, ONEW )
      ELSE
         ALPHA = -SQRT( DLAMCH( 'Epsilon' ) )
         IF( DUM(1).NE.ZERO ) ALPHA = ONE + ALPHA
      END IF
C
C Determine the lengths of working arrays
C
      LDA = MAX( 1, N )
      LDE = MAX( 1, N )
      LDB = MAX( 1, N )
      LDC = MAX( 1, M, P )
      LB  = MAX( 1, M, P )
      LC  = MAX( 1, N )
      IF( NLHS.LT.7 .AND. METH.LE.2 ) THEN
         LQ = 1
         COMPQ = 'N'
         LDQ = 1
      ELSE
         LQ = N
         COMPQ = 'I'
         LDQ = MAX( 1, N )
      END IF
      IF( NLHS.LT.8 .AND. METH.LE.2 ) THEN
         LZ = 1
         COMPZ = 'N'
         LDZ = 1
      ELSE
         LZ = N3
         IF( N3.EQ.0 ) THEN
            COMPZ = 'N'
         ELSE
            COMPZ = 'I'
         END IF
         LDZ = MAX( 1, N )
      END IF
C
      IF( JOBN.EQ.0 .OR. JOBN.EQ.2 ) THEN
         JOB = 'I'
      ELSE
         JOB = 'F'
      END IF
      IF( NLHS.GE.6 )
     $   LIWORK = N
      IF( METH.LE.2 ) THEN
         IF( SSTYPE ) THEN
            IF( METH.EQ.2 ) THEN
               CALL TB01WX( COMPQ, N, M, P, A, LDA, B, LDB, C, LDC,
     $                      Q, LDQ, DUM, -1, INFO )
               LDWORK = DUM(1)
            ELSE
               LDWORK = 1
            END IF
         ELSE
            IF( METH.EQ.1 ) THEN
               JOBA = 'N'
            ELSE
               JOBA = 'H'
            END IF
            CALL TG01LD( JOB, JOBA, COMPQ, COMPZ, N, M, P, A, LDA,
     $                   E, LDE, B, LDB, C, LDC, Q, LDQ, Z, LDZ, NF, ND,
     $                   NIBLCK, IWORK, TOL, IWORK, DUM, -1, INFO )
            LIWORK = 2*N
            LDWORK = DUM(1)
         END IF
      ELSE
         IF( SSTYPE ) THEN
            LDWORK = MAX( 1, 5*N )
         ELSE
            IF( METH.EQ.3 .OR. METH.EQ.4 ) THEN
               STDOM = 'N'
            ELSE
               IF( JOBN.EQ.0 .OR. JOBN.EQ.3 ) THEN
                  STDOM = 'U'
               ELSE
                  STDOM = 'S'
               END IF
            END IF
            CALL TG01QD( DICO, STDOM, JOB, N, M, P, ALPHA, A, LDA,
     $                   E, LDE, B, LDB, C, LDC, N1, N2, N3, ND,
     $                   NIBLCK, IWORK, Q, LDQ, Z, LDZ, DUM, DUM,
     $                   DUM, TOL, IWORK, DUM, -1, INFO )
            LDWORK = INT( DUM(1) ) + 3*N
            IF( METH.EQ.3 .OR. METH.EQ.5 ) THEN
               LIWORK = 2*N
            ELSE
               LIWORK = 2*N + 6
            END IF
         END IF
      END IF
C
C Allocate variable dimension local arrays
C !Fortran 90/95
C
      ALLOCATE ( A(LDA,MAX( 1, N )), B(LDB,LB), C(LDC,LC),
     $           E(LDE,MAX( 1, N )), Q(LDQ,LQ), Z(LDZ,LZ),
     $           DWORK(LDWORK), IWORK(LIWORK) )
C
C Copy inputs from MATLAB workspace to locally allocated arrays
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), A, NW*NW )
      IF( SSTYPE ) THEN
         CALL DLASET( 'F', N, N, ZERO, ONE, E, LDE )
      ELSE
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), E, NW*NW )
         SSTYPE = CHCKID( N, E, LDE )
      END IF
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ),
     $            B(1:MAX( 1, N ),1:MAX( 1, M)), NW*MW )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ),
     $            C(1:MAX( 1, P),1:MAX( 1, N )), PW*NW )
C
C     Do actual computations.
C
      NDIM = 0
      N1 = 0
      N2 = 0
      N3 = 0
      ND = 0
      NIBLCK = 0
      INFO = 0
      IF( METH.EQ.1 .OR. METH.EQ.2 ) THEN
         IF( SSTYPE ) THEN
            IF( METH.EQ.2 ) THEN
C
C              Reduction to Hessenberg form.
C
               CALL TB01WX( COMPQ, N, M, P, A, LDA, B, LDB, C, LDC,
     $                      Q, LDQ, DWORK, LDWORK, INFO )
            ELSE
               IF( NLHS.GE.7 )
     $            CALL DLASET( 'F', N, N, ZERO, ONE, Q, LDQ )
            END IF
            IF( JOBN.EQ.0 .OR. JOBN.EQ.2 ) THEN
               N3 = N
            ELSE
               N1 = N
            END IF
         ELSE
C
C           Finite-infinite separation in Hessenberg form
C
            CALL TG01LD( JOB, JOBA, COMPQ, COMPZ, N, M, P, A, LDA,
     $                   E, LDE, B, LDB, C, LDC, Q, LDQ, Z, LDZ, NF, ND,
     $                   NIBLCK, IWORK, TOL, IWORK(N+1), DWORK, LDWORK,
     $                   INFO )
           IF( JOBN.EQ.0 .OR. JOBN.EQ.2 ) THEN
               N1 = N - NF
               N3 = NF
            ELSE
               N1 = NF
               N3 = N - NF
            END IF
C
         END IF
      ELSE IF( METH.EQ.3 .OR. METH.EQ.5 ) THEN
         IF( SSTYPE ) THEN
            IF( METH.EQ.3 ) THEN
C
C              Reduction to real Schur form.
C
               CALL TB01WD( N, M, P, A, LDA, B, LDB, C, LDC, Q, LDQ,
     $                      DWORK(1), DWORK(N+1), DWORK(2*N+1),
     $                      LDWORK-2*N, INFO )
               IF( JOBN.EQ.0 .OR. JOBN.EQ.2 ) THEN
                  N3 = N
               ELSE
                  N1 = N
               END IF
            ELSE
C
C              Reduction to ordered real Schur form.
C
               IF( JOBN.EQ.0 .OR. JOBN.EQ.3 ) THEN
                  STDOM = 'U'
               ELSE
                  STDOM = 'S'
               END IF
               CALL TB01LD( DICO, STDOM, 'General', N, M, P, ALPHA,
     $                      A, LDA, B, LDB, C, LDC, NDIM, Q, LDQ,
     $                      DWORK(1), DWORK(N+1), DWORK(2*N+1),
     $                      LDWORK-2*N, INFO )
               IF( JOBN.EQ.0 .OR. JOBN.EQ.2 ) THEN
                  N2 = NDIM
                  N3 = N - NDIM
               ELSE
                  N1 = NDIM
                  N2 = N - NDIM
               END IF
            END IF
         ELSE
C
C           Finite-infinite separation in ordered generalized
C           real Schur form.
C
            CALL TG01QD( DICO, STDOM, JOB, N, M, P, ALPHA, A, LDA,
     $                   E, LDE, B, LDB, C, LDC, N1, N2, N3, ND,
     $                   NIBLCK, IWORK, Q, LDQ, Z, LDZ, DWORK(1),
     $                   DWORK(N+1), DWORK(2*N+1), TOL, IWORK(N+1),
     $                   DWORK(3*N+1), LDWORK-3*N, INFO )
         END IF
      ELSE IF( METH.EQ.4 .OR. METH.EQ.6 ) THEN
         IF( SSTYPE ) THEN
            IF( METH.EQ.4 ) THEN
C
C              Reduction to real Schur form.
C
               CALL TB01WD( N, M, P, A, LDA, B, LDB, C, LDC, Q, LDQ,
     $                      DWORK(1), DWORK(N+1), DWORK(2*N+1),
     $                      LDWORK-2*N, INFO )
               IF( JOBN.EQ.0 .OR. JOBN.EQ.2 ) THEN
                  N3 = N
               ELSE
                  N1 = N
               END IF
            ELSE
C
C              Stable-unstable separation in block diagonal
C              real Schur form form.
C
               IF( JOBN.EQ.0 .OR. JOBN.EQ.3 ) THEN
                  STDOM = 'Unstable'
               ELSE
                  STDOM = 'Stable'
               END IF
               CALL TB01LD( DICO, STDOM, 'General', N, M, P, ALPHA,
     $                      A, LDA, B, LDB, C, LDC, NDIM, Z, LDZ,
     $                      DWORK(1), DWORK(N+1), DWORK(2*N+1),
     $                      LDWORK-2*N, INFO )
               CALL TB01KX( N, M, P, NDIM, A, LDA, B, LDB, C, LDC,
     $                      Z, LDZ, Q, LDQ, INFO )
               IF( JOBN.EQ.0 .OR. JOBN.EQ.2 ) THEN
                  N2 = NDIM
                  N3 = N - NDIM
               ELSE
                  N1 = NDIM
                  N2 = N - NDIM
               END IF
            END IF
         ELSE
C
C           Finite-infinite separation in ordered generalized
C           real Schur form.
C
            CALL TG01QD( DICO, STDOM, JOB, N, M, P, ALPHA, A, LDA,
     $                   E, LDE, B, LDB, C, LDC, N1, N2, N3, ND,
     $                   NIBLCK, IWORK, Q, LDQ, Z, LDZ, DWORK(1),
     $                   DWORK(N+1), DWORK(2*N+1), TOL, IWORK(N+1),
     $                   DWORK(3*N+1), LDWORK-3*N, INFO )
            IF( INFO.EQ.0 )  THEN
               IF( JOBN.EQ.0 .OR. JOBN.EQ.2 ) THEN
                  NDIM = N1+N2
               ELSE
                  NDIM = N1
               END IF
C
C              Transpose Q in-situ.
C
               DO I = 2, N
                  CALL DSWAP( I-1, Q(1,I), 1, Q(I,1), LDQ )
               END DO
C
C              Separate spectrum.
C
               CALL TG01NX( 'D', N, M, P, NDIM, A, LDA, E, LDE, B, LDB,
     $                      C, LDC, Q, LDQ, Z, LDZ, IWORK(N+1), INFO )
               IF( INFO.NE.0 ) INFO = 4
            END IF
         END IF
      END IF
C
*      close(1)
C
C Copy output to MATLAB workspace
C
      IF( INFO.EQ.0 ) THEN
         NW = N
         IF( NLHS.GE.1 ) THEN
            PLHS(1) = mxCreateDoubleMatrix( NW, NW, SZERO )
            CALL mxCopyReal8ToPtr( A, mxGetPr( PLHS(1) ), NW*NW )
         END IF
         IF( NLHS.GE.2 ) THEN
            PLHS(2) = mxCreateDoubleMatrix( NW, NW, SZERO )
            CALL mxCopyReal8ToPtr( E, mxGetPr( PLHS(2) ), NW*NW )
         END IF
         IF( NLHS.GE.3 ) THEN
            MW = M
            PLHS(3) = mxCreateDoubleMatrix( NW, MW, SZERO )
            CALL mxCopyReal8ToPtr( B, mxGetPr( PLHS(3) ), NW*MW )
         END IF
         IF( NLHS.GE.4 ) THEN
            PW = P
            PLHS(4) = mxCreateDoubleMatrix( PW, NW, SZERO )
            CALL mxCopyReal8ToPtr( C(1:MAX( 1, P),1:MAX( 1, N )),
     $                             mxGetPr( PLHS(4) ), PW*NW )
         END IF
         IF( NLHS.GE.5 ) THEN
            DUM(1) = N1
            DUM(2) = N2
            DUM(3) = N3
            DUM(4) = N - ND
            PLHS(5) = mxCreateDoubleMatrix( FOURW, ONEW, SZERO )
            CALL mxCopyReal8ToPtr( DUM, mxGetPr( PLHS(5) ), FOURW)
         END IF
         IF( NLHS.GE.6 ) THEN
            IF( JOBN.EQ.0 ) THEN
C              JOB = 'I'
               DO I = 1, NIBLCK
                  DWORK(I) = IWORK(I)
               END DO
               DWORK(NIBLCK+1) = ND
            ELSE
C              JOB = 'F'
               DWORK(1) = ND
               DO I = 2, NIBLCK+1
                  DWORK(I) = IWORK(NIBLCK-I+2)
               END DO
            END IF
            NIBLCKW = NIBLCK
            PLHS(6) = mxCreateDoubleMatrix( NIBLCKW+ONEW, ONEW, SZERO )
            CALL mxCopyReal8ToPtr
     $                 ( DWORK, mxGetPr( PLHS(6) ), NIBLCKW+ONEW )
          END IF
          IF( NLHS.GE.8 ) THEN
             PLHS(8) = mxCreateDoubleMatrix( NW, NW, SZERO )
             IF( SSTYPE .AND. .NOT.METH.EQ.6 ) THEN
                CALL mxCopyReal8ToPtr( Q, mxGetPr( PLHS(8) ), NW*NW )
             ELSE
                CALL mxCopyReal8ToPtr( Z, mxGetPr( PLHS(8) ), NW*NW )
             END IF
          END IF
          IF( NLHS.GE.7 ) THEN
C
C            Transpose Q in-situ.
C
             IF( METH.LE.3 .OR. METH.EQ.5 .OR.
     $           ( SSTYPE .AND. .NOT.METH.EQ.6 ) ) THEN
                DO I = 2, N
                   CALL DSWAP( I-1, Q(1,I), 1, Q(I,1), LDQ )
                END DO
             END IF
             PLHS(7) = mxCreateDoubleMatrix( NW, NW, SZERO )
             CALL mxCopyReal8ToPtr( Q, mxGetPr( PLHS(7) ), NW*NW )
          END IF
C
      END IF
C
C Deallocate local arrays
C !Fortran 90/95
C
      DEALLOCATE ( A, B, C, E, Q, Z, DWORK, IWORK )
C
C Error handling
C
      IF( INFO.NE.0 ) THEN
         IF( METH.LE.2 ) THEN
            IF( INFO.EQ.1 ) THEN
               TEXT = 'The pencil A-sE is not regular'
            END IF
         ELSE IF( METH.EQ.3 ) THEN
            IF( SSTYPE ) THEN
               IF( INFO.EQ.1 ) THEN
                  TEXT = 'QR iteration did not converge'
               END IF
            ELSE
               IF( INFO.EQ.1 ) THEN
                  TEXT = 'The pencil A-sE is not regular'
               ELSE IF( INFO.EQ.2 ) THEN
                  TEXT = 'QZ iteration did not converge'
               END IF
            END IF
         ELSE IF( METH.EQ.4 ) THEN
            IF( SSTYPE ) THEN
               IF( INFO.EQ.1 ) THEN
                  TEXT = 'QR iteration did not converge'
               END IF
            ELSE
               IF( INFO.EQ.1 ) THEN
                  TEXT = 'The pencil A-sE is not regular'
               ELSE IF( INFO.EQ.2 ) THEN
                  TEXT = 'QZ iteration did not converge'
               ELSE IF( INFO.EQ.4 ) THEN
                  TEXT = 'The diagonal blocks have close eigenvalues'
               END IF
            END IF
         ELSE IF( METH.GE.5 ) THEN
            IF( SSTYPE ) THEN
               IF( INFO.EQ.1 ) THEN
                  TEXT = 'QR iteration did not converge'
               ELSE IF( INFO.EQ.2 ) THEN
                  TEXT = 'Failure of eigenvalue reordering'
               END IF
            ELSE
               IF( INFO.EQ.1 ) THEN
                  TEXT = 'The pencil A-sE is not regular'
               ELSE IF( INFO.EQ.2 ) THEN
                  TEXT = 'QZ iteration did not converge'
               ELSE IF( INFO.EQ.3 ) THEN
                  TEXT = 'Failure of generalized eigenvalue reordering'
               ELSE IF( INFO.EQ.4 ) THEN
                  TEXT = 'The diagonal blocks have close eigenvalues'
               END IF
            END IF
         END IF
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of SL_GSEP ***
      END
