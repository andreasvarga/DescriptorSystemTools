C Various additional Mex stuff
C
#include "fintrf.h"

      SUBROUTINE DISP( STR )
C
C     Print a character string STR in MATLAB window.
C
      CHARACTER*( * )   STR
      INTEGER           LCHR
C
      INTEGER*4       mexCallMATLAB
      MWPOINTER       mxCreateString
      MWPOINTER       plhs, prhs
      INTEGER         info
C
C     .. matlab integer*4 flags
C
      INTEGER*4 SONE,SZERO
      PARAMETER(SZERO=0,SONE=1)
C
      LCHR = LEN( STR )
      PRHS = mxCreateString(STR(1:LCHR))
C      call mxSetM(prhs,1)
C      call mxSetN(prhs,LCHR)
      info = mexCallMATLAB( SZERO, plhs, SONE, prhs, 'disp' )
      call mxDestroyArray(prhs)
      RETURN
      END

      INTEGER FUNCTION STR2NUM( STR )
C
C     Numerical value of a character string STR
C
      CHARACTER*( * )  STR
      INTEGER          I, I0, N, TP
C
      TP = 1
      N = 0
      I0 = IACHAR( '0' )
      DO I = LEN( STR ), 1, -1
         N = N + TP * (IACHAR( STR(I:I) ) - I0)
         TP = TP*10
      END DO
C
      STR2NUM = N
C
      RETURN
      END


      SUBROUTINE XERBLA( SRNAME, INFO )
*
*  -- LAPACK auxiliary routine (preliminary version) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
*     ..
*
*  Purpose
*  =======
*
*  XERBLA  is an error handler for the LAPACK routines.
*  It is called by an LAPACK routine if an input parameter has an
*  invalid value.  A message is printed and execution stops.
*
*  Installers may consider modifying the STOP statement in order to
*  call system-specific exception-handling facilities.
*
*  Arguments
*  =========
*
*  SRNAME  (input) CHARACTER*6
*          The name of the routine which called XERBLA.
*
*  INFO    (input) INTEGER
*          The position of the invalid parameter in the parameter list
*          of the calling routine.
*
*
      CHARACTER*72 ERRTXT
      WRITE( ERRTXT, FMT = 9999 ) SRNAME, INFO
*
      call mexErrMsgTxt( ERRTXT )
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ',
     $      'an illegal value' )
      RETURN
*
*     End of XERBLA
*
      END

      LOGICAL FUNCTION CHCKID( N, A, LDA )
C
C     To check if A = I.
C
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         ( ONE = 1.0D0, ZERO = 0.0D0 )
C
      INTEGER           LDA, N
      DOUBLE PRECISION  A(LDA,*)
C
      INTEGER           I, J
      LOGICAL           IDENT
      DOUBLE PRECISION  DIAG
C
      IDENT = .FALSE.
      DO 20 J = 1, N
         DIAG = A(J,J)
         A(J,J) = ZERO
         DO 10 I = 1, N
            IF( A(I,J).NE.ZERO ) THEN
               A(J,J) = DIAG
               GO TO 30
            END IF
   10    CONTINUE
         A(J,J) = DIAG
         IF( DIAG.NE.ONE ) GO TO 30
   20 CONTINUE
      IDENT = .TRUE.
   30 CHCKID = IDENT
C
      RETURN
C *** Last line of CHCKID ***
      END
