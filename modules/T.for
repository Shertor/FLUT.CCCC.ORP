
      SUBROUTINE T( X,Z,K)                                              T    711
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Z(K)                                                    T    712
      Y=X                                                               T    714
      Z1=1.0d+0                                                         T
      Z2=Y                                                              T    716
      Z(1)=Z1                                                           T    717
      Z(2)=Z2                                                           T    718
      DO 1 I=3,K                                                        T    719
      Z3=2.0D0 *Y*Z2-Z1                                                 T
      IF ( ABS(Z3).LE.1.D-19) Z3=0. D0                                  T
      Z(I)=Z3                                                           T    722
      Z1=Z2                                                             T    723
1     Z2=Z3                                                             T    724
      RETURN                                                            T    725
      END                                                               T    726
