
      DOUBLE PRECISION FUNCTION EIGEN (X,Y,Z,N,A,B)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(N),Z(N)                                               EIGE 664
      X1=(2.*X-B-A)/(B-A)                                               EIGE 665
      PI =3.141592653589 D 0                                            EIGE 666
      EIGEN=0.D0                                                        EIGE 667
      DO 1 I=1,N                                                        EIGE 668
1     EIGEN=EIGEN +Y(I)                                                 EIGE 669
      CALL T(X1,Z,N)                                                    EIGE 670
      NN=N-1                                                            EIGE 671
      DO 2 K=1,NN                                                       EIGE 672
      BE=0.D0                                                           EIGE 673
      DO 3 J=1,N                                                        EIGE 674
      P=(2.*J-1.)*PI/2./N                                               EIGE 675
3     BE=BE+2.*Y(J)*COS(K*P)                                            EIGE 676
2     EIGEN=EIGEN+Z(K+1)*BE                                             EIGE 677
      EIGEN=EIGEN/N
      RETURN                                                            EIGE 678
      END                                                               EIGE 679
