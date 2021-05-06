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
