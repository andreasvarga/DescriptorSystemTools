C SL_KLF.F - Gateway function for MB04VD via MB04UD.
C
C Matlab call:
C   [A,E,Q,Z,imuk,inuk,imukinf] = sl_klf(A,E,tol,mode,qz)
C
C Purpose
C     To compute orthogonal transformations Q and Z such that the
C     transformed pencil Q'(sE-A)Z is in upper block triangular form.
C     The tolerance TOL is used for rank computations.
C
C     If MODE = 0, then the matrices A and E are transformed into the
C     following generalized Schur form by unitary transformations Q1
C     and Z1 :
C
C                      | sE(eps,inf)-A(eps,inf) |      X     |
C        Q1'(sE-A)Z1 = |------------------------|------------|.   (1)
C                      |            O           | sE(r)-A(r) |
C
C     The pencil sE(eps,inf)-A(eps,inf) is in staircase form, and it
C     contains all Kronecker column indices and infinite elementary
C     divisors of the pencil sE-A. The pencil sE(r)-A(r) contains all
C     Kronecker row indices and elementary divisors of sE-A.
C     Note: X is a pencil.
C
C     If MODE = 1, then the submatrices having full row and column
C     rank in the pencil sE(eps,inf)-A(eps,inf) in (1) are
C     triangularized by applying unitary transformations Q2 and Z2 to
C     Q1'*(sE-A)*Z1.
C
C     If MODE = 2 (default), then the pencil sE(eps,inf)-A(eps,inf) in
C     (1) is separated into sE(eps)-A(eps) and sE(inf)-A(inf) by applying
C     unitary transformations Q3 and Z3 to Q2'*Q1'*(sE-A)*Z1*Z2.
C
C     This gives
C
C                | sE(eps)-A(eps) |        X       |      X     |
C                |----------------|----------------|------------|
C                |        O       | sE(inf)-A(inf) |      X     |
C     Q'(sE-A)Z =|=================================|============| (2)
C                |                                 |            |
C                |                O                | sE(r)-A(r) |
C
C     where Q = Q1*Q2*Q3 and Z = Z1*Z2*Z3.
C
C     QZ specifies which transformation matrices are determined
C     as follows:
C         QZ = 0, Q and Z are not accumulated and
C                 result as empty matrices;
C         QZ = 1, only the matrix Q is accumulated (Z is set to Z=[]);
C         QZ = 2, only the matrix Z is accumulated (Q is set to Q=[]);
C         QZ = 3, both Q and Z are accumulated (default).
C
C     The output vectors IMUK, INUK and IMUKINF specify the structure
C     of computed bloks as follows (see MB04VD)
C
C     IMUK    contains the column dimensions of the submatrices having
C             full column rank in the pencil sE(x)-A(x),
C                where  x = eps,inf  if MODE = 0 or 1,
C                or     x = eps      if MODE = 2.
C
C     INUK    contains the row dimensions of the submatrices having
C             full row rank in the pencil sE(x)-A(x),
C                where  x = eps,inf  if MODE = 0 or 1,
C                or     x = eps      if MODE = 2.
C
C     IMUKINF If MODE = 2, then the elements of this vector contain
C             the dimensions of the square diagonal submatrices in the pencil
C             sE(inf)-A(inf).
C
C Contributor
C   A. Varga, 28-12-1997.
C
C Revisions:
C   A. Varga, March 2002; 13-01-2016.
C   V. Sima, Dec. 2016.
C
C **********************************************************************
C
#include "fintrf.h"
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C .. Parameters ..
      IMPLICIT          NONE
      DOUBLE PRECISION  ONE, TWO, THREE, ZERO
      PARAMETER         ( ONE   = 1.0D0, TWO  = 2.0D0,
     $                    THREE = 3.0D0, ZERO = 0.0D0 )
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
      CHARACTER         JOBQ, JOBZ, MODE
      INTEGER           INFO, LDA, LDE, LDQ, LDZ, LWORK, M,
     $                  MNEI(4), N, NBLCKS, NBLCKI, RANKE
      MWSIZE            KW, MW, NW, NBLCKSW
      DOUBLE PRECISION  TOL
C
C ..Local variables
      CHARACTER*120     TEXT
      INTEGER           I, K, M1, MN, N1, QZ
      MWSIZE            HMIN1W, M1W, N1W
      MWSIZE            ONEW
      LOGICAL           MODES, WITHQ, WITHZ
      DOUBLE PRECISION  XMODE, XQZ
C
      DOUBLE PRECISION, ALLOCATABLE:: A(:,:), E(:,:), Q(:,:), WORK(:),
     $                                Z(:,:)
      INTEGER, ALLOCATABLE:: IMUK(:), INUK(:), ISTAIR(:),
     $                       IMUK0(:), IWORK(:)
C
      EXTERNAL          MB04UD, MB04VD
C     ..Intrinsic Functions..
      INTRINSIC         MAX, MIN
C
C     .. matlab integer*4 flags
C
      INTEGER*4         SZERO
      PARAMETER         ( SZERO = 0 )
C
      LOGICAL  LSAME
      EXTERNAL LSAME
C
      ONEW = 1
C
      IF( NRHS.EQ.0 ) RETURN
C
C Check for proper number of arguments
C
      IF( NRHS.LT.2 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SL_KLF requires at least 2 input arguments' )
      ELSE IF( NLHS.GT.7 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SL_KLF requires at most 7 output argument' )
      END IF
C
C CHECK THE DIMENSIONS OF A, E, TOL, MODE, and QZ.
C They must be N X N, N X N, 1 X 1, 1 X 1, 1 X 1.
C
      MW = mxGetM( PRHS(1) )
      M   = MW
      NW = mxGetN( PRHS(1) )
      N   = NW
      M1W = mxGetM( PRHS(2) )
      M1  = M1W
      N1W = mxGetN( PRHS(2) )
      N1  = N1W
*
* Debuging by writing in a file
*
*      OPEN (1, FILE = 'test.txt', STATUS = 'REPLACE')
*
C
      IF( (M.NE.M1) .OR. (N.NE.N1) ) THEN
         CALL mexErrMsgTxt
     $        ( 'A and E must have the same dimensions' )
      END IF
      IF( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a real matrix' )
      END IF
      IF( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $    mxIsComplex( PRHS(2) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'E must be a real matrix' )
      END IF
      TOL = ZERO
      IF( NRHS.GT.2 ) THEN
         IF( mxGetM( PRHS(3) ).NE.1 .OR. mxGetN( PRHS(3) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'TOL must be a scalar' )
         END IF
         IF( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(3) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'TOL must be a real scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), TOL, ONEW )
      END IF
C
      MODE = 'S'
      IF( NRHS.GT.3 ) THEN
         IF( mxGetM( PRHS(4) ).NE.1 .OR. mxGetN( PRHS(4) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'MODE must be a scalar' )
         END IF
         IF( mxIsNumeric( prhs(4) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(4) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'MODE must be a real scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), XMODE, ONEW )
         IF( XMODE.NE.ZERO .AND. XMODE.NE.1 .AND. XMODE.NE.TWO ) THEN
            CALL mexErrMsgTxt
     $           ( 'MODE has 0, 1 and 2 the only admissible values' )
         END IF
         IF( XMODE.EQ.ZERO ) THEN
            MODE = 'B'
         ELSE IF ( XMODE.EQ.ONE ) THEN
            MODE = 'T'
         END IF
      END IF
C
      IF( NRHS.GT.4 ) THEN
         IF( mxGetM( PRHS(5) ).NE.1 .OR. mxGetN( PRHS(5) ).NE.1 ) THEN
            CALL mexErrMsgTxt( 'QZ must be a scalar' )
         END IF
         IF( mxIsNumeric( PRHS(5) ).EQ.0 .OR.
     $       mxIsComplex( PRHS(5) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'QZ must be a real scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), XQZ, ONEW )
         IF( XQZ.NE.ZERO .AND. XQZ.NE.ONE .AND. XQZ.NE.TWO .AND.
     $       XQZ.NE.THREE ) THEN
            CALL mexErrMsgTxt
     $           ( 'QZ has 0, 1, 2 and 3 the only admissible values' )
         END IF
         QZ = XQZ
      ELSE
         QZ = 3
      END IF
      WITHQ = NLHS.GE.3 .AND. (QZ.EQ.1 .OR. QZ.EQ.3)
      WITHZ = NLHS.GE.4 .AND. QZ.GE.2
      MODES = LSAME( MODE, 'S' )
C
C  SET PROBLEM PARAMETERS
C
      LDA = MAX( 1, M )
      LDE = MAX( 1, M )
      LDQ = MAX( 1, M )
      LDZ = MAX( 1, N )
C
      IF( WITHQ ) THEN
         JOBQ = 'I'
      ELSE
         JOBQ = 'N'
      END IF
      IF( WITHZ ) THEN
         JOBZ = 'I'
      ELSE
         JOBZ = 'N'
      END IF
C
C ALLOCATE WORKING ARRAYS
C
      LWORK = MAX( 1, M, N )
      MN    = MAX( M+1, N+1 )
      ALLOCATE ( A(LDA,N), E(LDE,N), IMUK(MN), INUK(MN),
     $           ISTAIR(MN),  IMUK0(MN), IWORK(N),
     $           WORK(LWORK) )
      ALLOCATE ( Q(LDQ,M), Z(LDZ,N) )
C
C Copy inputs from MATLAB workspace to locally allocated arrays
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), A, MW*NW )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), E, MW*NW )
C
C Do the actual computations
C
      CALL MB04UD( JOBQ, JOBZ, M, N, A, LDA, E, LDE, Q, LDQ, Z, LDZ,
     $             RANKE, ISTAIR, TOL, WORK, INFO )
      IF( WITHQ ) JOBQ = 'U'
      IF( WITHZ ) JOBZ = 'U'
C
      CALL MB04VD( MODE, JOBQ, JOBZ, M, N, RANKE, A, LDA, E, LDE,
     $             Q, LDQ, Z, LDZ, ISTAIR, NBLCKS, NBLCKI, IMUK, INUK,
     $             IMUK0, MNEI, TOL, IWORK, INFO )
*      close(1)
C
C Copy output to MATLAB workspace
C
      IF( INFO.EQ.0 ) THEN
         MW = M
         NW = N
         PLHS(1) = mxCreateDoubleMatrix( MW, NW, SZERO )
         PLHS(2) = mxCreateDoubleMatrix( MW, NW, SZERO )
         CALL mxCopyReal8ToPtr(A, mxGetPr( PLHS(1) ), MW*NW )
         CALL mxCopyReal8ToPtr(E, mxGetPr( PLHS(2) ), MW*NW )
C
         IF( NLHS.GE.3 ) THEN
            IF( WITHQ ) THEN
               K = M
            ELSE
               K = 0
            END IF
            KW = K
            PLHS(3) = mxCreateDoubleMatrix( KW, KW, SZERO )
            CALL mxCopyReal8ToPtr( Q, mxGetPr( PLHS(3) ), KW*KW )
         END IF
C
         IF( NLHS.GE.4 ) THEN
            IF( WITHZ ) THEN
               K = N
            ELSE
               K = 0
            END IF
            KW = K
            PLHS(4) = mxCreateDoubleMatrix( KW, KW, SZERO )
            CALL mxCopyReal8ToPtr( Z, mxGetPr( PLHS(4) ), KW*KW )
         END IF
C
         IF( NLHS.GE.5 ) THEN
            DO 10 I = 1, NBLCKS
               WORK(I) = IMUK(I)
   10       CONTINUE
            NBLCKSW = NBLCKS
            HMIN1W  = MIN( 1, NBLCKS )
            PLHS(5) = mxCreateDoubleMatrix( NBLCKSW, HMIN1W, SZERO )
            CALL mxCopyReal8ToPtr( WORK, mxGetPr( PLHS(5) ), NBLCKSW )
         END IF
C
         IF( NLHS.GE.6 ) THEN
            DO 20 I = 1, NBLCKS
               WORK(I) = INUK(I)
   20       CONTINUE
            NBLCKSW = NBLCKS
            HMIN1W  = MIN( 1, NBLCKS )
            PLHS(6) = mxCreateDoubleMatrix( NBLCKSW, HMIN1W, SZERO )
            CALL mxCopyReal8ToPtr( WORK, mxGetPr( PLHS(6) ), NBLCKSW )
         END IF
C
         IF( NLHS.GE.7 ) THEN
            IF( MODES) THEN
               M1 = NBLCKI
               N1 = MIN(1,NBLCKI)
               DO 30 I = 1, NBLCKI
                  WORK(I) = IMUK0(I)
   30          CONTINUE
            ELSE
               M1 = 0
               N1 = 0
            END IF
            M1W = M1
            N1W = N1
            PLHS(7) = mxCreateDoubleMatrix( M1W, N1W, SZERO )
            CALL mxCopyReal8ToPtr( WORK, mxGetPr( PLHS(7) ), M1W*N1W )
         END IF
C
      END IF
C
C DEALLOCATE LOCAL ARRAYS
C
      DEALLOCATE (A, E, WORK, IMUK, INUK, ISTAIR, IMUK0, IWORK)
      DEALLOCATE (Q, Z)
C
C  ERROR HANDLING
C
      IF( INFO.NE.0 ) THEN
         write(TEXT, '( " INFO =", I4, " ON EXIT FROM MB04VD" )' ) INFO
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of SL_KLF ***
      END
