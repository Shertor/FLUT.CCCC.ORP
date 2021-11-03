
      SUBROUTINE DIFX1 (D,M,A,B)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION D(M,M)
      PI=3.141592653589D0
      DO 1 NU=1,M
      PSINU=(2.*NU-1.)*PI/2./M
      DO 1 MU=1,M
      PSIMU=(2.*MU-1.)*PI/2./M
      P=0.
      DO 2 K=1,M-1
    2 P=P+K*COS(K*PSIMU)*SIN(K*PSINU)/SIN(PSINU)
    1 D(NU,MU)=4.*P/M/(B-A)
      RETURN
      END
