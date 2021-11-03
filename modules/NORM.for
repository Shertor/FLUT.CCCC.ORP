
      SUBROUTINE NORM (Y,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Y(N)
      P=0.D0
      DO 1 I=1,N
      IF (ABS(Y(I)).GT.P) P=ABS(Y(I))
1     CONTINUE
      DO 2 I=1,N
2     Y(I)=Y(I)/P
      RETURN
      END
