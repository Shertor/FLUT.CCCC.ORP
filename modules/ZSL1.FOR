
      SUBROUTINE ZSL1(D,K,A,B)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*4 Q
      DIMENSION D(K,K)
      PI=3.141592653589D0
      DO 1 I=1,K
      PSI=(2*I-1)*PI/2./K
      DO 1 J=1,K
      PSJ=(2*J-1)*PI/2./K
      S=1.D0
      DO 2 Q=1,K-1
    2 S=S+COS(Q*PSJ)*((2+Q**2)*COS(Q*PSI)+
     +3.*Q*COS(PSI)*SIN(Q*PSI)/SIN(PSI))
    1 D(I,J)=-2.*(2./(B-A))**2*S/K/SIN(PSJ)**2
      RETURN
      END
