
      SUBROUTINE DIFT (X,M,T,T1,T2,T3,T4)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION T(M),T1(M),T2(M),T3(M),T4(M)
      T(1)=1.D0
      T(2)=X
      T1(1)=0.D0
      T1(2)=1.D0
      T2(1)=0.D0
      T2(2)=0.D0
      T3(1)=0.D0
      T3(2)=0.D0
      T4(1)=0.D0
      T4(2)=0.D0
      DO K=2,M-1
      T(K+1)=2.D0*X*T(K)-T(K-1)
      T1(K+1)=2.D0*T(K)+2.D0*X*T1(K)-T1(K-1)
      T2(K+1)=4.D0*T1(K)+2.D0*X*T2(K)-T2(K-1)
      T3(K+1)=6.D0*T2(K)+2.D0*X*T3(K)-T3(K-1)
      T4(K+1)=8.D0*T3(K)+2.D0*X*T4(K)-T4(K-1)
      ENDDO
      RETURN
      END
