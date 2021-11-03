
      SUBROUTINE C05AZF(X, Y, FX, TOLX, IR, C, IND, IFAIL)
C     MARK 8 RE-ISSUE. NAG COPYRIGHT 1979.
C     .. SCALAR ARGUMENTS ..
      DOUBLE PRECISION FX, TOLX, X, Y
      INTEGER IFAIL, IND, IR
C     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION C(17)
C     ..
C     .. LOCAL SCALARS ..
C$P 1
      DOUBLE PRECISION SRNAME
      DOUBLE PRECISION AB, DIFF1, DIFF2, DIFF, MAX, REL, TOL
      INTEGER I
      LOGICAL T
C     .. FUNCTION REFERENCES ..
      DOUBLE PRECISION DSQRT, X02AAF, X02ABF
C     *** IMPLEMENTATION DEPENDENT DECLARATION ***
C     DOUBLE PRECISION DFLOAT
      INTEGER P01AAF
C     ..
      DATA SRNAME /8H C05AZF /
      I = 0
      IF ((IND.GT.0 .AND. IND.LE.4) .OR. IND.EQ.-1) GO TO 20
C     USER NOT CHECKED IND OR CHANGED IT
      I = 2
      IND = 0
      GO TO 640
   20 IF (TOLX.GT.0.D0 .AND. (IR.EQ.0 .OR. IR.EQ.1 .OR. IR.EQ.2))
     * GOTO 40
      I = 3
      IND = 0
      GO TO 640
   40 REL = 1.D0
      AB = 1.D0
      IF (IR.EQ.1) REL = 0.D0
      IF (IR.EQ.2) AB = 0.D0
      IF (IND.EQ.-1) GO TO 80
      GO TO (60, 100, 180, 480), IND
   60 C(3) = X
      IND = 2
      RETURN
   80 C(3) = X
  100 IF (FX.NE.0.D0) GO TO 140
  120 Y = X
      IND = 0
      I = 0
      GO TO 640
  140 C(4) = FX
      C(15) = DABS(FX)
      C(16) = 0.D0
      X = Y
      Y = C(3)
      C(2) = C(4)
      C(5) = X
      IF (IND.EQ.-1) GO TO 160
      IND = 3
      RETURN
  160 FX = C(1)
      IND = 3
  180 IF (FX.EQ.0.D0) GO TO 120
      IF (DSIGN(1.D0,FX).NE.DSIGN(1.D0,C(2))) GO TO 200
      IND = 0
      I = 1
      GO TO 640
  200 C(6) = FX
      C(13) = DSQRT(X02AAF(0.0D0))
      C(15) = DMAX1(C(15),DABS(FX))
      C(14) = X02ABF(0.0D0)
      C(16) = 0.0D0
  220 C(1) = C(5)
      C(2) = C(6)
      C(17) = 0.D0
  240 IF (DABS(C(2)).GE.DABS(C(4))) GO TO 280
      IF (C(1).EQ.C(5)) GO TO 260
      C(7) = C(5)
      C(8) = C(6)
  260 C(5) = C(3)
      C(6) = C(4)
      X = C(1)
      C(3) = X
      C(4) = C(2)
      C(1) = C(5)
      C(2) = C(6)
  280 TOL = 0.5D0*TOLX*DMAX1(AB,REL*DABS(C(3)))
      DIFF2 = 0.5D0*(C(1)-C(3))
      C(12) = DIFF2
      DIFF2 = DIFF2 + C(3)
      IF (C(12).EQ.0.D0) GO TO 340
      IF (DABS(C(12)).LE.TOL) GO TO 580
      IF (C(17).LT.2.5D0) GO TO 300
      C(11) = C(12)
      GO TO 460
  300 TOL = TOL*DSIGN(1.D0,C(12))
      DIFF1 = (C(3)-C(5))*C(4)
      IF (C(17).GT.1.5D0) GO TO 320
      DIFF = C(6) - C(4)
      GO TO 380
  320 IF (C(7).NE.C(3) .AND. C(7).NE.C(5)) GO TO 360
  340 IND = 0
      I = 5
      GO TO 640
  360 C(9) = (C(8)-C(4))/(C(7)-C(3))
      C(10) = (C(8)-C(6))/(C(7)-C(5))
      DIFF1 = C(10)*DIFF1
      DIFF = C(9)*C(6) - C(10)*C(4)
  380 IF (DIFF1.GE.0.D0) GO TO 400
      DIFF1 = -DIFF1
      DIFF = -DIFF
  400 IF (DABS(DIFF1).GT.C(14) .AND. DIFF1.GT.DIFF*TOL) GO TO 420
      C(11) = TOL
      GO TO 460
  420 IF (DIFF1.GE.C(12)*DIFF) GO TO 440
      C(11) = DIFF1/DIFF
      GO TO 460
  440 C(11) = C(12)
  460 C(7) = C(5)
      C(8) = C(6)
      C(5) = C(3)
      C(6) = C(4)
      C(3) = C(3) + C(11)
      X = C(3)
      Y = C(1)
      IND = 4
      RETURN
  480 IF (FX.EQ.0.D0) GO TO 120
      C(4) = FX
      MAX = DABS(FX)
      IF (C(13)*MAX.LE.C(15)) GO TO 500
      IF (C(16).EQ.1.D0) C(16) = -1.D0
      IF (C(16).EQ.0.D0) C(16) = 1.D0
      GO TO 520
  500 C(16) = 0.D0
  520 IF (C(2).GE.0.D0) GO TO 540
      T = C(4).LE.0.D0
      GO TO 560
  540 T = C(4).GE.0.D0
  560 IF (T) GO TO 220
      I = IFIX(SNGL(C(17)+0.1D0))
      I = I + 1
      IF (C(11).EQ.C(12)) I = 0
      C(17) = DFLOAT(I)
      GO TO 240
  580 IF (C(16).GE.0.D0) GO TO 600
      I = 4
      GO TO 620
  600 Y = C(1)
      I = 0
  620 IND = 0
  640 IFAIL = P01AAF(IFAIL,I,SRNAME)
      RETURN
      END
