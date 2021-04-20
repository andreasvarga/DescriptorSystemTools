C SL_GZERO.F - Gateway function for SLICOT descriptor system zeros routine
C              AG08BD.F and AB08NW.F.
C
C MATLAB call:
C   [zf,ni,infdiv,kronr,infe,kronl] = sl_gzero(A,E,B,C,D,tol,scale)
C
C Purpose:
C   To compute the zeros of a descriptor system and the Kronecker
C   structure of the system pencil. Set E = [] for standard state-spase
C   systems. Set D = [] for a system with no direct feedthrough matrix.
C
C Input parameters:
C   A,E,
C   B,C,D   - descriptor system matrices
C   tol     - (optional) tolerance parameter for rank test
C   scale   - (optional) scaling option parameter
C                = 0 - perform scaling
C                = 1 - no scaling (default)
C
C Output parameters:
C   zf      - finite zeros
C   ni      - 2x1 vector
C             ni(1) contains the number of infinite zeros
C             ni(2) contains the normal rank of system pencil
C   infdiv  - multiplicities of infinite zeros
C   kronr   - right Kronecker indices
C   infe    - elementary infinite blocks
C   kronl   - left Kronecker indices
C
C Contributor
C   A. Varga, May 1999.
C
C Revisions:
C   A. Varga, March, 2002; 13-01-2016.
C   V. Sima, Dec. 2016, June 2017.
C
C **********************************************************************
C
#include "fintrf.h"
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
      IMPLICIT NONE
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
      INTEGER*4         mxIsNumeric
C
C .. Scalar parameters used by SLICOT subroutines ..
      CHARACTER         EQUIL
      INTEGER           DINFZ, INFO, L, LDA, LDAF, LDB, LDC, LDD, LDE,
     $                  LDWORK, M, N, NFZ, NINFE, NIZ, NKROL, NKROR,
     $                  NRANK, P
      DOUBLE PRECISION  TOL
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      DOUBLE PRECISION, ALLOCATABLE:: A(:,:), ALPHI(:), ALPHR(:),
     $                                B(:,:), BETA(:), C(:,:), D(:,:),
     $                                DWORK(:), E(:,:)
      INTEGER, ALLOCATABLE::          INFE(:), INFZ(:), IWORK(:),
     $                                KRONL(:), KRONR(:)
C
C .. Local variables and constant dimension arrays ..
      CHARACTER*120     TEXT
      LOGICAL           SSTYPE
      INTEGER           I, INFOFZ, ISCALE, L1, L2, LA, M4, MAXNZ, MP,
     $                  N1, N3, P4
      MWSIZE            HMIN1W, ONEW, TWOW
      DOUBLE PRECISION  DUM(2)
C
C     .. matlab integer*4 flags
C
      INTEGER*4         SONE, SZERO
      PARAMETER         ( SONE = 1, SZERO = 0 )
C
C .. External functions ..
      LOGICAL           CHCKID
      EXTERNAL          CHCKID
C
C .. External Subroutines ..
      EXTERNAL          DLASET, AB08NW, AG08BD
C
C .. Intrinsic functions ..
      INTRINSIC         MAX, MIN
C
      ONEW = 1
      TWOW = 2
*
* Debuging by writing in a file
*
*      OPEN (1, FILE = 'test.txt', STATUS = 'REPLACE')
*
      IF( NRHS.EQ.0 ) RETURN
C
C Check for proper number of arguments
C
      IF( NRHS.LT.5 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SL_GZERO requires at least 5 input arguments' )
      ELSE IF( NLHS.GT.6 ) THEN
         CALL mexErrMsgTxt
     $        ( 'SL_GZERO has at most 6 output argument' )
      END IF
C
C Check dimensions of input parameters and read/set scalar parameters
C
C   A(LxN), E(LxN), B(LxM), C(PxN), D(PxM)
C
      L  = mxGetM( PRHS(1) )
      N  = mxGetN( PRHS(1) )
      L1 = mxGetM( PRHS(2) )
      N1 = mxGetN( PRHS(2) )
      L2 = mxGetM( PRHS(3) )
      M  = mxGetN( PRHS(3) )
      P  = mxGetM( PRHS(4) )
      N3 = mxGetN( PRHS(4) )
      P4 = mxGetM( PRHS(5) )
      M4 = mxGetN( PRHS(5) )
C
      SSTYPE = L1.EQ.0 .AND. N1.EQ.0
      IF( SSTYPE ) THEN
         IF( L.NE.N ) THEN
            CALL mexErrMsgTxt
     $           ( 'A must be square for a standard system' )
         END IF
      ELSE
         IF( N1.NE.N .OR. L1.NE.L ) THEN
            CALL mexErrMsgTxt
     $           ( 'A and E must have the same dimensions' )
         END IF
      END IF
      IF( MIN( M, L2 ).GT.0 .AND. L2.NE.L ) THEN
         CALL mexErrMsgTxt
     $        ( 'B must have the same row dimension as A' )
      END IF
      IF( MIN( P, N3 ).GT.0 .AND. N3.NE.N ) THEN
         CALL mexErrMsgTxt
     $        ( 'C must have the same column dimension as A' )
      END IF
      IF( MIN( P4, M4 ).GT.0 .AND. MIN(L,N).GT.0 .AND.
     $    P4.NE.P .AND. M4.NE.M ) THEN
         CALL mexErrMsgTxt
     $        ( 'D must have compatible dimensions with B and C' )
      END IF
      IF( mxIsNumeric( PRHS(1) ).EQ.0 ) THEN
         CALL mexErrMsgTxt( 'A must be a numeric matrix' )
      END IF
      IF( mxIsNumeric( PRHS(2) ).EQ.0 ) THEN
         CALL mexErrMsgTxt( 'E must be a numeric matrix' )
      END IF
      IF( mxIsNumeric( PRHS(3) ).EQ.0 ) THEN
         CALL mexErrMsgTxt( 'B must be a numeric matrix' )
      END IF
      IF( mxIsNumeric( PRHS(4) ).EQ.0 ) THEN
         CALL mexErrMsgTxt( 'C must be a numeric matrix' )
      END IF
      IF( mxIsNumeric( PRHS(5) ).EQ.0 ) THEN
         CALL mexErrMsgTxt( 'D must be a numeric matrix' )
      END IF
      M = MAX( M, M4 )
      P = MAX( P, P4 )
C
C   tol
C
      TOL = ZERO
      IF( NRHS.GE.6 ) THEN
         IF( mxGetM( PRHS(6) )*mxGetN( PRHS(6) ).GT.1 ) THEN
            CALL mexErrMsgTxt( 'TOL must be scalar' )
         END IF
         IF( mxIsNumeric( PRHS(6) ).EQ.0 ) THEN
            CALL mexErrMsgTxt( 'TOL must be a numeric scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(6) ), TOL, ONEW )
      END IF
C
C     scale
C
      ISCALE = 0
      IF( NRHS.GE.7 ) THEN
         IF( mxGetM( PRHS(7) )*mxGetN( PRHS(7) ).GT.1 ) THEN
            CALL mexErrMsgTxt( 'SCALE must be scalar' )
         END IF
         IF( mxIsNumeric( PRHS(7) ).EQ.0 ) THEN
            CALL mexErrMsgTxt( 'SCALE must be a numeric scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(7) ), DUM, ONEW )
         ISCALE = DUM(1)
         IF( .NOT.( ISCALE.EQ.0 .OR. ISCALE.EQ.1 ) ) THEN
            CALL mexErrMsgTxt
     $           ('SCALE has 0 and 1 the only admissible values')
         END IF
      END IF
C
C  Set options.
C
      IF( ISCALE.EQ.0 ) THEN
         EQUIL = 'N'
      ELSE
         EQUIL = 'S'
      END IF
C
C  Set Problem Parmeters
C
      MP   = MAX( M, P )
      LDA  = MAX( 1, L )
      LDE  = LDA
      LDB  = LDA
      IF( MIN( L2, M ).EQ.0 ) LDB = 1
      LDC  = MAX( 1, P )
      LDD  = LDC
      LDAF = MAX( 1, L+P, N+M )
      LA   = MAX( 1, L, N )
      MAXNZ = 1 + MIN( L+P, N+M )
C
C Determine the lengths of working arrays
C
      IF( SSTYPE ) THEN
         CALL AB08NW( EQUIL, N, M, P, A, LDA, B, LDB, C, LDC, D, LDD,
     $                NFZ, NRANK, NIZ, DINFZ, NKROR, NINFE, NKROL,
     $                INFZ, KRONR, INFE, KRONL, E, LDE, TOL, IWORK,
     $                DUM, -1, INFO )
         CALL DGEEV( 'No vectors', 'No vectors', N, A, LDA, DUM, DUM,
     $               DUM, 1, DUM, 1, DUM(2), -1, INFO )
      ELSE
         CALL AG08BD( EQUIL, L, N, M, P, A, LDA, E, LDE, B, LDB, C, LDC,
     $                D, LDD, NFZ, NRANK, NIZ, DINFZ, NKROR, NINFE,
     $                NKROL, INFZ, KRONR, INFE, KRONL, TOL, IWORK,
     $                DUM, -1, INFO )
         CALL DGEGV( 'No vectors', 'No vectors', N, A, LDA, E, LDE, DUM,
     $               DUM, DUM, DUM, 1, DUM, 1, DUM(2), -1, INFO )
      END IF
      LDWORK = MAX( DUM(1), DUM(2) )
C
C Allocate variable dimension local arrays
C !Fortran 90/95
C
      ALLOCATE ( A(LDA,LA), B(LDB,MAX( 1, M )), C(LDC,LA),
     $           D(LDD,MAX( 1, M )), E(LDE,LA),
     $           ALPHR(LA), ALPHI(LA), BETA(LA),
     $           DWORK(LDWORK), INFE(MAXNZ), INFZ(LA),
     $           IWORK(LA+MP+1), KRONL(L+P+1), KRONR(N+M+1) )
C
C Copy inputs from MATLAB workspace to locally allocated arrays
C
      CALL mxCopyPtrToReal8( mxGetPR( PRHS(1) ), A, L*N )
      IF( .NOT.SSTYPE ) THEN
         CALL mxCopyPtrToReal8( mxGetPR( PRHS(2) ), E, L*N )
         IF( L.EQ.N) SSTYPE = SSTYPE .OR. CHCKID( N, E, LDE )
      END IF
      CALL mxCopyPtrToReal8( mxGetPR( PRHS(3) ), B, L*M )
      CALL mxCopyPtrToReal8( mxGetPR( PRHS(4) ), C, P*N )
      IF( MIN( P4, M4 ).EQ.0 ) THEN
         CALL DLASET( 'F', P, M, ZERO, ZERO, D, LDD )
      ELSE
         CALL mxCopyPtrToReal8( mxGetPR( PRHS(5) ), D, P*M )
      END IF
C
C Do the actual computations
C
      INFOFZ = 0
      IF( SSTYPE) THEN
         CALL AB08NW( EQUIL, N, M, P, A, LDA, B, LDB, C, LDC, D, LDD,
     $                NFZ, NRANK, NIZ, DINFZ, NKROR, NINFE, NKROL,
     $                INFZ, KRONR, INFE, KRONL, E, LDE, TOL, IWORK,
     $                DWORK, LDWORK, INFO )
      ELSE
         CALL AG08BD( EQUIL, L, N, M, P, A, LDA, E, LDE, B, LDB, C, LDC,
     $                D, LDD, NFZ, NRANK, NIZ, DINFZ, NKROR, NINFE,
     $                NKROL, INFZ, KRONR, INFE, KRONL, TOL, IWORK,
     $                DWORK, LDWORK, INFO )
      END IF
      IF( INFO.EQ.0 .AND. NFZ.GT.0 ) THEN
         IF( CHCKID( NFZ, E, LDE ) ) THEN
            CALL DGEEV( 'No vectors', 'No vectors', NFZ, A, LDA,
     $                  ALPHR, ALPHI, DUM, 1, DUM, 1, DWORK, LDWORK,
     $                  INFOFZ )
         ELSE
            CALL DGEGV( 'No vectors', 'No vectors', NFZ, A, LDA, E, LDE,
     $                  ALPHR, ALPHI, BETA, DUM, 1, DUM, 1, DWORK,
     $                  LDWORK, INFOFZ )
            I = 1
   10       IF( I.LE.NFZ ) THEN
               IF( ALPHI(I).EQ.ZERO ) THEN
                  ALPHR(I) = ALPHR(I) / BETA(I)
                  I = I + 1
               ELSE
                  ALPHR(I) = ALPHR(I) / BETA(I)
                  ALPHI(I) = ALPHI(I) / BETA(I)
                  ALPHR(I+1) = ALPHR(I)
                  ALPHI(I+1) = -ALPHI(I)
                  I = I + 2
               END IF
               GO TO 10
            END IF
         END IF
      END IF
C
C Copy output to MATLAB workspace
C
      IF( INFO.EQ.0 .AND. INFOFZ.EQ.0 ) THEN
         IF( NLHS.GE.1 ) THEN
            HMIN1W = MIN( 1, NFZ )
            PLHS(1) = mxCreateDoubleMatrix( NFZ, HMIN1W, SONE )
            CALL mxCopyReal8ToPtr( ALPHR, mxGetPR( PLHS(1) ), NFZ )
            CALL mxCopyReal8ToPtr( ALPHI, mxGetPI( PLHS(1) ), NFZ )
      END IF
         IF( NLHS.GE.2 ) THEN
            DUM(1)  = NIZ
            DUM(2)  = NRANK
            PLHS(2) = mxCreateDoubleMatrix( TWOW, ONEW, SZERO )
            CALL mxCopyReal8ToPtr( DUM, mxGetPR( PLHS(2) ), TWOW )
      END IF
         IF( NLHS.GE.3 ) THEN
            DO 20 I = 1, DINFZ
               DWORK(I) = INFZ(I)
   20       CONTINUE
            HMIN1W = MIN( 1, DINFZ )
            PLHS(3) = mxCreateDoubleMatrix( DINFZ, HMIN1W, SZERO )
            CALL mxCopyReal8ToPtr( DWORK, mxGetPR( PLHS(3) ), DINFZ )
         END IF
         IF( NLHS.GE.4 ) THEN
            DO 30 I = 1, NKROR
               DWORK(I) = KRONR(I)
   30       CONTINUE
            HMIN1W = MIN( 1, NKROR )
            PLHS(4) = mxCreateDoubleMatrix( NKROR, HMIN1W, SZERO )
            CALL mxCopyReal8ToPtr( DWORK, mxGetPR( PLHS(4) ), NKROR )
         END IF
         IF( NLHS.GE.5 ) THEN
            DO 40 I = 1, NINFE
               DWORK(I) = INFE(I)
   40       CONTINUE
            HMIN1W = MIN( 1, NINFE )
            PLHS(5) = mxCreateDoubleMatrix( NINFE, HMIN1W, SZERO )
            CALL mxCopyReal8ToPtr( DWORK, mxGetPR( PLHS(5) ), NINFE )
         END IF
         IF( NLHS.GE.6 ) THEN
            DO 50 I = 1, NKROL
               DWORK(I) = KRONL(I)
   50       CONTINUE
            HMIN1W = MIN( 1, NKROL )
            PLHS(6) = mxCreateDoubleMatrix( NKROL, HMIN1W, SZERO )
            CALL mxCopyReal8ToPtr( DWORK, mxGetPR( PLHS(6) ), NKROL )
         END IF
      END IF
C
C Deallocate local arrays
C !Fortran 90/95
C
      DEALLOCATE (A, B, C, D, E, ALPHR, ALPHI, BETA, DWORK,
     $            INFE, INFZ, KRONL, KRONR, IWORK, STAT = INFO )
*      close(1)
C
C Error handling
C
      IF( INFO.NE.0 .OR. INFOFZ.NE.0 ) THEN
         IF( INFO.NE.0 ) THEN
            write( TEXT,
     $             '( " INFO =", I4, " ON EXIT FROM AG08BD" )' ) INFO
         END IF
         IF( INFOFZ.NE.0 ) THEN
            write( TEXT,
     $             '( " INFO =", I4, " ON EXIT FROM DGEEV/DGEGV" )' )
     $      INFOFZ
         END IF
         CALL mexErrMsgTxt( TEXT )
      END IF
C
      RETURN
C *** Last line of SL_GZERO ***
      END
