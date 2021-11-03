
C                                                                       86440001
C     ------------------------------------------------------------------86440002
C                                                                       86440003
C     SUBROUTINE HQRM(NM,N,LOW,IGH,H,WR,WI,IERR)                        86440004
      SUBROUTINE HQRM2(NM,N,LOW,IGH,H,WR,WI,IERR,NITER)
C                                                                       86440005
      INTEGER I,J,K,L,M,N,EN,LL,MM,NA,NM,IGH,ITS,LOW,MP2,ENM2,IERR      86440006
      REAL*8 H(NM,N),WR(N),WI(N)                                        86440007
      REAL*8 P,Q,R,S,T,W,X,Y,ZZ,NORM,MACHEP                             86440008
      REAL*8 DSQRT,DABS,DSIGN                                           86440009
      INTEGER MIN0                                                      86440010
      LOGICAL NOTLAS                                                    86440011
C                                                                       86440012
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE HQR,      86440013
C     NUM. MATH. 14, 219-231(1970) BY MARTIN, PETERS, AND WILKINSON.    86440014
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 359-371(1971).   86440015
C                                                                       86440016
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A REAL                   86440017
C     UPPER HESSENBERG MATRIX BY THE QR METHOD.                         86440018
C                                                                       86440019
C     ON INPUT:                                                         86440020
C                                                                       86440021
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         86440022
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          86440023
C          DIMENSION STATEMENT;                                         86440024
C                                                                       86440025
C        N IS THE ORDER OF THE MATRIX;                                  86440026
C                                                                       86440027
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING           86440028
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,          86440029
C          SET LOW=1, IGH=N;                                            86440030
C                                                                       86440031
C        H CONTAINS THE UPPER HESSENBERG MATRIX.  INFORMATION ABOUT     86440032
C          THE TRANSFORMATIONS USED IN THE REDUCTION TO HESSENBERG      86440033
C          FORM BY  ELMHES  OR  ORTHES, IF PERFORMED, IS STORED         86440034
C          IN THE REMAINING TRIANGLE UNDER THE HESSENBERG MATRIX.       86440035
C                                                                       86440036
C     ON OUTPUT:                                                        86440037
C                                                                       86440038
C        H HAS BEEN DESTROYED.  THEREFORE, IT MUST BE SAVED             86440039
C          BEFORE CALLING  HQR  IF SUBSEQUENT CALCULATION AND           86440040
C          BACK TRANSFORMATION OF EIGENVECTORS IS TO BE PERFORMED;      86440041
C                                                                       86440042
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,                86440043
C          RESPECTIVELY, OF THE EIGENVALUES.  THE EIGENVALUES           86440044
C          ARE UNORDERED EXCEPT THAT COMPLEX CONJUGATE PAIRS            86440045
C          OF VALUES APPEAR CONSECUTIVELY WITH THE EIGENVALUE           86440046
C          HAVING THE POSITIVE IMAGINARY PART FIRST.  IF AN             86440047
C          ERROR EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT        86440048
C          FOR INDICES IERR+1,...,N;                                    86440049
C                                                                       86440050
C        IERR IS SET TO                                                 86440051
C          ZERO       FOR NORMAL RETURN,                                86440052
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN               86440053
C                     DETERMINED AFTER 30 ITERATIONS.                   86440054
C                                                                       86440055
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        86440056
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         86440057
C                                                                       86440058
C     ------------------------------------------------------------------86440059
C                                                                       86440060
C     :::::::::: MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING     86440061
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.   86440062
C                MACHEP = 16.0D0**(-13) FOR LONG FORM ARITHMETIC        86440063
C                ON S360 ::::::::::                                     86440064
C     DATA MACHEP/Z3410000000000000/                                    86440065
C     MACHEP = EPS(1.D0)
      MACHEP = 16.0D0**(-13)
C                                                                       86440066
      IERR = 0                                                          86440067
      NORM = 0.0D0                                                      86440068
      K = 1                                                             86440069
C     :::::::::: STORE ROOTS ISOLATED BY BALANC                         86440070
C                AND COMPUTE MATRIX NORM ::::::::::                     86440071
      DO 50 I = 1, N                                                    86440072
C                                                                       86440073
         DO 40 J = K, N                                                 86440074
   40    NORM = NORM + DABS(H(I,J))                                     86440075
C                                                                       86440076
         K = I                                                          86440077
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 50                      86440078
         WR(I) = H(I,I)                                                 86440079
         WI(I) = 0.0D0                                                  86440080
   50 CONTINUE                                                          86440081
C                                                                       86440082
      EN = IGH                                                          86440083
      T = 0.0D0                                                         86440084
C     :::::::::: SEARCH FOR NEXT EIGENVALUES ::::::::::                 86440085
   60 IF (EN .LT. LOW) GO TO 1001                                       86440086
      ITS = 0                                                           86440087
      NA = EN - 1                                                       86440088
      ENM2 = NA - 1                                                     86440089
C     :::::::::: LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT             86440090
C                FOR L=EN STEP -1 UNTIL LOW DO -- ::::::::::            86440091
   70 DO 80 LL = LOW, EN                                                86440092
         L = EN + LOW - LL                                              86440093
         IF (L .EQ. LOW) GO TO 100                                      86440094
         S = DABS(H(L-1,L-1)) + DABS(H(L,L))                            86440095
         IF (S .EQ. 0.0D0) S = NORM                                     86440096
         IF (DABS(H(L,L-1)) .LE. MACHEP * S) GO TO 100                  86440097
   80 CONTINUE                                                          86440098
C     :::::::::: FORM SHIFT ::::::::::                                  86440099
  100 X = H(EN,EN)                                                      86440100
      IF (L .EQ. EN) GO TO 270                                          86440101
      Y = H(NA,NA)                                                      86440102
      W = H(EN,NA) * H(NA,EN)                                           86440103
      IF (L .EQ. NA) GO TO 280                                          86440104
C     IF (ITS .EQ. 30) GO TO 1000
      IF (ITS .EQ. NITER) GO TO 1000                                    86440105
      IF (ITS .NE. 10 .AND. ITS .NE. 20) GO TO 130                      86440106
C     :::::::::: FORM EXCEPTIONAL SHIFT ::::::::::                      86440107
      T = T + X                                                         86440108
C                                                                       86440109
      DO 120 I = LOW, EN                                                86440110
  120 H(I,I) = H(I,I) - X                                               86440111
C                                                                       86440112
      S = DABS(H(EN,NA)) + DABS(H(NA,ENM2))                             86440113
      X = 0.75D0 * S                                                    86440114
      Y = X                                                             86440115
      W = -0.4375D0 * S * S                                             86440116
  130 ITS = ITS + 1                                                     86440117
C     :::::::::: LOOK FOR TWO CONSECUTIVE SMALL                         86440118
C                SUB-DIAGONAL ELEMENTS.                                 86440119
C                FOR M=EN-2 STEP -1 UNTIL L DO -- ::::::::::            86440120
      DO 140 MM = L, ENM2                                               86440121
         M = ENM2 + L - MM                                              86440122
         ZZ = H(M,M)                                                    86440123
         R = X - ZZ                                                     86440124
         S = Y - ZZ                                                     86440125
         P = (R * S - W) / H(M+1,M) + H(M,M+1)                          86440126
         Q = H(M+1,M+1) - ZZ - R - S                                    86440127
         R = H(M+2,M+1)                                                 86440128
         S = DABS(P) + DABS(Q) + DABS(R)                                86440129
         P = P / S                                                      86440130
         Q = Q / S                                                      86440131
         R = R / S                                                      86440132
         IF (M .EQ. L) GO TO 150                                        86440133
         IF (DABS(H(M,M-1)) * (DABS(Q) + DABS(R)) .LE. MACHEP * DABS(P) 86440134
     X    * (DABS(H(M-1,M-1)) + DABS(ZZ) + DABS(H(M+1,M+1)))) GO TO 150 86440135
  140 CONTINUE                                                          86440136
C                                                                       86440137
  150 MP2 = M + 2                                                       86440138
C                                                                       86440139
      DO 160 I = MP2, EN                                                86440140
         H(I,I-2) = 0.0D0                                               86440141
         IF (I .EQ. MP2) GO TO 160                                      86440142
         H(I,I-3) = 0.0D0                                               86440143
  160 CONTINUE                                                          86440144
C     :::::::::: DOUBLE QR STEP INVOLVING ROWS L TO EN AND              86440145
C                COLUMNS M TO EN ::::::::::                             86440146
      DO 260 K = M, NA                                                  86440147
         NOTLAS = K .NE. NA                                             86440148
         IF (K .EQ. M) GO TO 170                                        86440149
         P = H(K,K-1)                                                   86440150
         Q = H(K+1,K-1)                                                 86440151
         R = 0.0D0                                                      86440152
         IF (NOTLAS) R = H(K+2,K-1)                                     86440153
         X = DABS(P) + DABS(Q) + DABS(R)                                86440154
         IF (X .EQ. 0.0D0) GO TO 260                                    86440155
         P = P / X                                                      86440156
         Q = Q / X                                                      86440157
         R = R / X                                                      86440158
  170    S = DSIGN(DSQRT(P*P+Q*Q+R*R),P)                                86440159
         IF (K .EQ. M) GO TO 180                                        86440160
         H(K,K-1) = -S * X                                              86440161
         GO TO 190                                                      86440162
  180    IF (L .NE. M) H(K,K-1) = -H(K,K-1)                             86440163
  190    P = P + S                                                      86440164
         X = P / S                                                      86440165
         Y = Q / S                                                      86440166
         ZZ = R / S                                                     86440167
         Q = Q / P                                                      86440168
         R = R / P                                                      86440169
C     :::::::::: ROW MODIFICATION ::::::::::                            86440170
         DO 210 J = K, EN                                               86440171
            P = H(K,J) + Q * H(K+1,J)                                   86440172
            IF (.NOT. NOTLAS) GO TO 200                                 86440173
            P = P + R * H(K+2,J)                                        86440174
            H(K+2,J) = H(K+2,J) - P * ZZ                                86440175
  200       H(K+1,J) = H(K+1,J) - P * Y                                 86440176
            H(K,J) = H(K,J) - P * X                                     86440177
  210    CONTINUE                                                       86440178
C                                                                       86440179
         J = MIN0(EN,K+3)                                               86440180
C     :::::::::: COLUMN MODIFICATION ::::::::::                         86440181
         DO 230 I = L, J                                                86440182
            P = X * H(I,K) + Y * H(I,K+1)                               86440183
            IF (.NOT. NOTLAS) GO TO 220                                 86440184
            P = P + ZZ * H(I,K+2)                                       86440185
            H(I,K+2) = H(I,K+2) - P * R                                 86440186
  220       H(I,K+1) = H(I,K+1) - P * Q                                 86440187
            H(I,K) = H(I,K) - P                                         86440188
  230    CONTINUE                                                       86440189
C                                                                       86440190
  260 CONTINUE                                                          86440191
C                                                                       86440192
      GO TO 70                                                          86440193
C     :::::::::: ONE ROOT FOUND ::::::::::                              86440194
  270 WR(EN) = X + T                                                    86440195
      WI(EN) = 0.0D0                                                    86440196
      EN = NA                                                           86440197
      GO TO 60                                                          86440198
C     :::::::::: TWO ROOTS FOUND ::::::::::                             86440199
  280 P = (Y - X) / 2.0D0                                               86440200
      Q = P * P + W                                                     86440201
      ZZ = DSQRT(DABS(Q))                                               86440202
      X = X + T                                                         86440203
      IF (Q .LT. 0.0D0) GO TO 320                                       86440204
C     :::::::::: REAL PAIR ::::::::::                                   86440205
      ZZ = P + DSIGN(ZZ,P)                                              86440206
      WR(NA) = X + ZZ                                                   86440207
      WR(EN) = WR(NA)                                                   86440208
      IF (ZZ .NE. 0.0D0) WR(EN) = X - W / ZZ                            86440209
      WI(NA) = 0.0D0                                                    86440210
      WI(EN) = 0.0D0                                                    86440211
      GO TO 330                                                         86440212
C     :::::::::: COMPLEX PAIR ::::::::::                                86440213
  320 WR(NA) = X + P                                                    86440214
      WR(EN) = X + P                                                    86440215
      WI(NA) = ZZ                                                       86440216
      WI(EN) = -ZZ                                                      86440217
  330 EN = ENM2                                                         86440218
      GO TO 60                                                          86440219
C     :::::::::: SET ERROR -- NO CONVERGENCE TO AN                      86440220
C                EIGENVALUE AFTER 30 ITERATIONS ::::::::::              86440221
 1000 IERR = EN                                                         86440222
 1001 RETURN                                                            86440223
C     :::::::::: LAST CARD OF HQR ::::::::::                            86440224
      END                                                               86440225
C                                                                       20210001
