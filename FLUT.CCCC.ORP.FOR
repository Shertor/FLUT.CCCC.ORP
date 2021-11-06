      PROGRAM FLUTBP14
C     FLUTTER OF ORTHOTROPIC RECTANGULAR PLATE

      IMPLICIT REAL*8 (A-H, O-Z)
      CHARACTER*26 PATH 
      REAL*8 K
      REAL*8 NUx, NUy
      CHARACTER*3 KEY
      INTEGER :: i
      CHARACTER(len=32) :: arg
      COMMON /DATA/ K,H,BL,RO,DB
      COMMON /DEBUG/ ISDEBUG
      COMMON /D/ Dx,Dxy,Dy
      COMMON //PATH
      DATA PATH /'C:/reports/'/

      ! IF WRITES ARE NEEDED SET ISDEBUG, ELSE IT PRINTS ONLY RESULT SPEED
      ISDEBUG = .FALSE.
      
      IF (ISDEBUG) THEN
          WRITE(*,*)'-----DEBUG PRINT ALL MODE-----'
      ENDIF
      
      IF (ISDEBUG) THEN
          CONTINUE
      ELSE
          IF(COMMAND_ARGUMENT_COUNT().NE.8) THEN
              WRITE(*,*)'BL;HL;Ex;Ey;NUx;NUy;Gxy;RO'
              STOP
          ENDIF
      ENDIF
      
C     ХАРАКТЕРНОЕ ДАВЛЕНИЕ
C     P0 = 1 АТМ (ФИЗИЧЕСКАЯ) = 1.0133E5 Н/М**2 = 1.0333 КГ/СМ**2
C     CHARACTERISTIC PRESSURE
C     P0 = 1 ATMOSPHERE (PHYSICAL) = 1.0133E5 H/M**2 = 1.0333 KG/CM**2
      P0 = 1.0333D0

      
C     КОЭФФИЦИЕНТ ПУАССОНА PNU=0.33
C     POISSONS RATIO
      !PNU = 0.33D0
      NUx = 0.15D0
      IF (ISDEBUG) THEN
          CONTINUE
      ELSE          
          CALL get_command_argument(5, arg)
          READ(arg,*)NUx
      ENDIF
      NUy = 0.12D0
      IF (ISDEBUG) THEN
          CONTINUE
      ELSE          
          CALL get_command_argument(6, arg)
          READ(arg,*)NUy
      ENDIF
      
      
C     ПОКАЗАТЕЛЬ ПОЛИТРОПЫ K = 1.4
C     ЭТО СООТВЕТСТВУЕТ СКОРСТИ ЗВУКА 331.26 М/СЕК
C     POLYTROPY INDICATOR K = 1.4
C     CORRESPONDS TO THE SONIC SPEED 331.26 M/SEC
      K = 1.4D0

      
C     БЕЗРАЗМЕРНАЯ ТОЛЩИНА ПЛАСТИНКИ H
C     DIMENSIONLESS PLATE THICKNESS
      HL = 0.015D0
      IF (ISDEBUG) THEN
          CONTINUE
      ELSE          
          CALL get_command_argument(2, arg)
          READ(arg,*)HL
      ENDIF
      H = HL
      
      IF (ISDEBUG) WRITE(*,*)'H =',H
      IF (ISDEBUG) PAUSE

      
C     МОДУЛЬ ЮНГА В КГ/СМ**2
C     YOUNG'S MODULE IN KG/CM**2
      E = 305914.8D0 !ALUMINIUM
      IF (ISDEBUG) THEN
          CONTINUE
      ELSE          
          CALL get_command_argument(3, arg)
          READ(arg,*)E
      ENDIF
      Ex = E
      Ey = 203943.2D0
      IF (ISDEBUG) THEN
          CONTINUE
      ELSE          
          CALL get_command_argument(4, arg)
          READ(arg,*)Ey
      ENDIF
      

C     МОДУЛЬ СДВИГА Н/М**2
C     SHEAR MODULE H/M**2
      Gxy = 152957.4D0
      IF (ISDEBUG) THEN
          CONTINUE
      ELSE          
          CALL get_command_argument(7, arg)
          READ(arg,*)Gxy
      ENDIF
      

C     ПЛОТНОСТЬ ВОЗДУХА 1.2928 (КГ/М**3)
C     AIR DENSITY 1.2928 KG/M**3
      R = 1.2928D0

      
C     ПЛОТНОСТЬ МАТЕРИАЛА ПЛАСТИНЫ 2.7E3 (КГ/М**3)
C     PLATE MATERIAL DENSITY KG/M**3
      RO = 2.75D3
      IF (ISDEBUG) THEN
          CONTINUE
      ELSE          
          CALL get_command_argument(8, arg)
          READ(arg,*)RO
      ENDIF
      
C     ПОЛОВИНА ШИРИНЫ ПЛАСТИНЫ
C     HALF OF THE PLATE WIDTH (SET 1 FOR SQUARE PLATE)
      BL = 2.4D0
      IF (ISDEBUG) THEN
          CONTINUE
      ELSE          
          CALL get_command_argument(1, arg)
          READ(arg,*)BL
      ENDIF
      
      IF (ISDEBUG) WRITE(*,*)'BL =', BL
      IF (ISDEBUG) WRITE(*,*)' '
      
      
C     SONIC SPEEDS
   10 FORMAT(A)
      C0 = SQRT(K*9.81523D4*P0/R)
      IF (ISDEBUG) WRITE(*,15) C0
   15 FORMAT(' SONIC SPEED C0 = ',E12.5)
      C1 = 1.D2*SQRT(9.81523D4*E/RO)
      C1 = SQRT((1.D6/10.19762D0)*E/RO)
      IF (ISDEBUG) WRITE(*,16) C1
   16 FORMAT(' SONIC SPEED IN THE PLATE C1 = ',E12.5)
      C2 = C1/C0
      IF (ISDEBUG) WRITE(*,17) C2
   17 FORMAT(' DIMENSIOLENS SONIC SPEED IN THE PLATE C2 = ',E12.5)  
      IF (ISDEBUG) PAUSE

C     ОБЕЗРАЗМЕРИВАНИЕ
C     DIMENSIONLESS  
      Ex = Ex/P0
      Ey = Ey/P0
      Gxy = Gxy/P0
      
      
C     ЦИЛИНДРИЧЕСКАЯ ЖЕСТКОСТЬ
C     CYLINDRICAL STIFFNESS
      Dx = Ex*HL**3/12.0/(1.0-NUx*NUy)
      Dy = Ey*HL**3/12.0/(1.0-NUx*NUy)
      Dxy = 2.0*Gxy*HL**3/12.0+Dx*NUy
      D = Dx
      DB = D
      
      
      Dx = Dx/D
      Dy = Dy/D
      Dxy = Dxy/D
      
      RO=K*RO/R


      CALL FLUT34C

      CALL FLUT33D

C     POINTS WRITE TO FILE
      IF (ISDEBUG) WRITE(*,*) 'EXPORT RE OF EIGENFUNCTION? (Y/N)'

      KEY='Y'
      IF (ISDEBUG) READ (*,10)KEY
      IF (KEY.NE.'Y') STOP

      CALL FLUT16C

      CALL FLUT16D
      
      IF (ISDEBUG) WRITE(*,*) 'EXPORT FINISHED'

      STOP
      END


C     ---------------------------------------------------------------------------
      SUBROUTINE FLUT34C

      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*1 KEY,KEY1
      REAL*8 K
      EXTERNAL F
      COMMON /DATA/ K,H,BL,RO,DB
      COMMON /DEBUG/ ISDEBUG
      COMMON /COM1/ IPR,NUMB
      COMMON /COM2/ NLAMDA,TJ,M,N,NT
      COMMON /COM3/ V,IN
      COMMON /FL/ FLG(20)
      COMMON /ABH0H1/ AA,BB,H0,H1
      DIMENSION F1(50),F2(50),V19(2),IN19(2)
      DIMENSION AA(2500),BB(2500),H0(6250000),H1(6250000)

C     THIS IS NUMBER OF BOUNDARY CONDITIONS
C     1 - CCCC    2 - SSSS
      NUMB = 1

      
      IF (ISDEBUG) WRITE(*,*) ' '
      IF (ISDEBUG) WRITE(*,*) 'TJ = ? (IN DEGREES)'
      
      TJ2 = 180.0D0
      TJ1 = 0.0D0
      IF (ISDEBUG) READ (*,*) TJ1
      TJ=TJ1/TJ2
      
 700  IF (ISDEBUG) WRITE(*,*) ' '

   10 FORMAT(A)
      KEY1 = 'A'
      IF(KEY1.EQ.'Q') RETURN
      IF(KEY1.EQ.'A') GO TO 25

      IF (ISDEBUG) WRITE(*,*) ' '
      IF (ISDEBUG) WRITE(*,*) 'M = ?, N = ?'
      IF (ISDEBUG) READ (*,*) M,N
      IF (ISDEBUG) THEN
          IF(M.GT.40.OR.N.GT.40) THEN
              WRITE(*,*) 'M OR N MORE THAN 40. STOP'
              STOP
          ENDIF
      ENDIF

      NT=M*N

  150 IF (ISDEBUG) WRITE(*,*) ' '
      IF (ISDEBUG) WRITE(*,*) ' NLAMDA = ? '
      IF (ISDEBUG) READ (*,*) NLAMDA

      IF (ISDEBUG) WRITE(*,*) ' '
      IF (ISDEBUG) WRITE(*,*) ' A,B = ?'
      IF (ISDEBUG) READ (*,*) A,B

      CALL MATRH2 (M,N,NT,H0,AA,BB,H1,BL,NUMB)

      IF(F(A)*F(B).GT.0.) THEN
          IF (ISDEBUG) WRITE(*,*) 'ON THE ENDS OF THE
     *SET INTERVAL FUNCTION OF STABILITY HAS IDENTICAL SIGNS'
          GO TO 150
      ENDIF
      GO TO 200
      
   25 CONTINUE
C     DEFAULT PARAMETERS
      IPR = 1
      M = 21
      N = 21
      NT = M*N
      NLAMDA = 5
      V0 = 0.
      STEP = 0.01
      V = 8.0D0
      NSTEP = 20
      
      
      IF (ISDEBUG) WRITE(*,26) NLAMDA
   26 FORMAT(' NLAMDA =', I3)
      IF (ISDEBUG) WRITE(*,27) V0
   27 FORMAT(' V0 =', E9.2)
      IF (ISDEBUG) WRITE(*,28) STEP
   28 FORMAT(' STEP =', E9.2)
      IF (ISDEBUG) WRITE(*,29) V
   29 FORMAT(' V =', E9.2)
      IF (ISDEBUG) WRITE(*,30) NSTEP
   30 FORMAT(' NSTEP =', I3)
  310 IF (ISDEBUG) WRITE(*,*) 'GRID 19X19;
     * BEGINNING THE SEARCH OF THE ROOT.'
      XX = V0

      !
      CALL MATRH2 (M,N,NT,H0,AA,BB,H1,BL,NUMB)
      !
      
      FF = F(XX)
      DO I = 1,NLAMDA
          F1(I) = FLG(I)
      ENDDO
      
  400 XX = XX+STEP
      IF(XX.LE.V) THEN
          FF = F(XX)
          DO I = 1,NLAMDA
              F2(I) = FLG(I)
          ENDDO
      ELSE
          GO TO 500
      ENDIF
      
      DO I = 1,NLAMDA
          IF (F1(I)*F2(I).LT.0.) THEN
              IN = I
              GO TO 300
          ENDIF
      ENDDO

      DO I = 1,NLAMDA
          F1(I) = F2(I)
      ENDDO
      
      GO TO 400
  500 IF (ISDEBUG) WRITE(*,*)' END OF THE SEARCH OF THE ROOT'
      IF (ISDEBUG) WRITE(*,*)' IN THE INTERVAL'
      IF (ISDEBUG) WRITE(*,526) V0,V
      
  526 FORMAT(E10.3,1H ,E10.3,1H )
      IF (ISDEBUG) WRITE(*,*)'CHANGES OF SIGNS AT FUNCT
     * OF STABILITY ARE NOT REV'
      WRITE (*,70) 0,0
      STOP
      
  300 IF (ISDEBUG) WRITE(*,15) IN
   15 FORMAT(5H IN =,I3)
      B = XX
      A = XX-STEP
      IF (ISDEBUG) WRITE(*,14) A,B
   14 FORMAT (4H A =,E10.3,4H B =,E10.3)

      NLAMD = NLAMDA
      NLAMDA = IN

  200 EPS1 = 1.E-8
      ETA = 0.0
      IFAIL = 0
      IPR = 2
      
      CALL C05ADF (A,B,EPS1,ETA,F,X,IFAIL)
      
      IF (ISDEBUG) WRITE(*,*)'IFAIL =',IFAIL
      IF (IFAIL.NE.0) THEN
          WRITE(*,70) 0,0
          STOP
      ENDIF
      
      IF (ISDEBUG) WRITE(*,13) M,N
   13 FORMAT(' ON THE GRID M =',I3,'N =',I3,'CRITICAL SPEED IS FOUND:')
      IF (ISDEBUG) WRITE(*,12) X
   12 FORMAT (1X,E15.8)
      IF(KEY1.NE.'Y'.AND.KEY1.NE.'A') THEN
          V = X
          IN = NLAMDA
          IF (ISDEBUG) WRITE(*,*)' REZULT'
          IF (ISDEBUG) WRITE(*,70) V, IN
          RETURN
      ENDIF
      V9 = X
      IN9 = IN
      IF (ISDEBUG) WRITE(*,*) 'GRID 21X21;
     * BEGINNING THE SEARCH OF THE ROOT.'
      M = 21
      N = 21
      NT = M*N
C
      CALL MATRH2 (M,N,NT,H0,AA,BB,H1,BL,NUMB)
C
      FF1 = F(A)
      FF2 = F(B)
      IF (FF1*FF2.GT.0.) THEN
          NLAMDA = NLAMD
          XX = A
          FF = F(XX)
          DO I=1,NLAMD
              F1(I)=FLG(I)
          ENDDO
          
  401     XX=XX+STEP
          IF(XX.LE.V) THEN
              FF = F(XX)
              DO I=1,NLAMD
                  F2(I)=FLG(I)
              ENDDO
          ELSE
              GO TO 501
          ENDIF
          DO I=1,NLAMD
              IF (F1(I)*F2(I).LT.0.) THEN
                  IN=I
                  GO TO 301
              ENDIF
          ENDDO
          DO I=1,NLAMD
              F1(I)=F2(I)
          ENDDO
          GO TO 401
  501     CONTINUE
          NLAMDA = NLAMD
          XX = A
          FF = F(XX)
          DO I=1,NLAMD
              F1(I) = FLG(I)
          ENDDO
  402     XX=XX-STEP

          IF (XX.LE.0.D0) GO TO 502

          IF(XX.GT.VL) THEN
              FF=F(XX)
              DO I=1,NLAMD
                  F2(I)=FLG(I)
              ENDDO
          ELSE
          GO TO 502
          ENDIF
          DO I=1,NLAMD
              IF (F1(I)*F2(I).LT.0.) THEN
                  IN=I
                  GO TO 302
              ENDIF
          ENDDO

          DO I=1,NLAMD
              F1(I)=F2(I)
          ENDDO
          GO TO 402

  302     IF (ISDEBUG) WRITE(*,15) IN
          A1 = XX
          B1 = XX+STEP
          IP = 1
          GO TO 303
  301     IF (ISDEBUG) WRITE(*,15) IN
          B1 = XX
          A1 = XX - STEP
          IP = 2
          IF (ISDEBUG) WRITE(*,14) A1,B1
  303     NLAMDA = IN
          
          CALL C05ADF (A1,B1,EPS1,ETA,F,X,IFAIL)
          
          IF (ISDEBUG) WRITE(*,*)'IFAIL =',IFAIL
          IF (IFAIL.NE.0) THEN
              WRITE(*,70) 0,0
              STOP
          ENDIF
          IF (ISDEBUG) WRITE(*,33) M,N
   33 FORMAT(' ON THE GRID OF M =', I3,'N =', I3,' CRITICAL SPEED:')

          IF (ISDEBUG) WRITE(*,12) X
          IF (ISDEBUG) WRITE(*,15) NLAMDA
          V19(IP) = X
          IN19(IP) = NLAMDA
          IF(IP.EQ.2) THEN
              VL = A - NSTEP*STEP
              IF (ISDEBUG) WRITE (*,34) A
   34 FORMAT(' THE ROOT TO THE LEFT OF THE POINT',E10.3)
              GO TO 501
          ENDIF
          
      ELSE
          CALL C05ADF (A,B,EPS1,ETA,F,X,IFAIL)
          IF (ISDEBUG) WRITE (*,*)' IFAIL =',IFAIL
          IF (IFAIL.NE.0) THEN
              WRITE(*,70) 0,0
              STOP
          ENDIF
          IF (ISDEBUG) WRITE(*,12) X
          IF (ISDEBUG) WRITE(*,15) NLAMDA
          V19(1) = X
          IN19(1) = NLAMDA
          V = X
          IN = NLAMDA
          NLAMDA = NLAMD
          IF (ISDEBUG) WRITE(*,*) 'REZULT'
          WRITE(*,70) V,IN
      ENDIF
      
      
      IF (ISDEBUG) WRITE(*,*) ' '
      IF (ISDEBUG) WRITE(*,50) V9,IN9
      IF (ISDEBUG) WRITE(*,60) V19(1), IN19(1)
      GO TO 800

  502 CONTINUE
      IF (ISDEBUG) WRITE(*,35) A
   35 FORMAT(' SLEVA FROM THE POINT, CRITICAL SPEED IS NOT FOUND',E10.3,
     *'CRITICAL SPEED IS NOT FOUND')
      IF (ISDEBUG) WRITE(*,50) V9,IN9
      IF (ISDEBUG) WRITE(*,60) V19(1),IN19(1)
      IF (ISDEBUG) WRITE(*,60) V19(2),IN19(2)
   50 FORMAT(' V9 =',E13.6,'IN9 =',I3)
   60 FORMAT(' V19 =',E13.6,'IN19 =',I3)
      V = V19(2)
      IN = IN19(2)
      IF(V19(1).GT.0.0.AND.V19(1).LT.V19(2)) THEN
          V = V19(1)
          IN = IN19(1)
      ENDIF
      IF (ISDEBUG) WRITE(*,*) ' '
      IF (ISDEBUG) WRITE(*,*) 'REZULT'
      WRITE(*,70) V,IN
      NLAMDA = NLAMD
   70 FORMAT(E13.6,';',I3)
      GO TO 700
  800 RETURN
      END
      
      
C     ---------------------------------------------------------------------------
      REAL*8 FUNCTION F(XX)
      IMPLICIT REAL*8 (A-H,O-Z)
      character*26 PATH 
      COMMON /COM1/ IPR,NUMB
      COMMON /COM2/ NLAMDA,TJ,M,N,NT
      COMMON /DATA/ K,H,BL,RO,DB
      COMMON /DEBUG/ ISDEBUG
      COMMON /FL/ FLG(20)
      COMMON /ABH0H1/ A,B,H0,H1
      COMMON //PATH

      DIMENSION A(2500),B(2500),DX(2500),DY(2500),H0(6250000),
     1H1(6250000),X(2500),Y(2500),IANA(2500),IJ(2500),SCALE(2500)

      REAL*8 K
      COMPLEX*16 CZ,LAMDA(40)

      CALL MATRH3 (M,N,NT,H0,DX,DY,H1,TJ,BL,K,DB,XX)
      
C     ВЫЧИСЛЕНИЕ СОБСТВЕННЫХ ЗНАЧЕНИЙ МАТРИЦЫ H0
C     CALCULATING THE EIGENVALUES OF A MATRIX H0
      
      NOUT = 4
      OPEN(UNIT=NOUT,FILE=TRIM(ADJUSTL(PATH))//'NOUT')

      CALL BALANC(NT,NT,H0,LOW,IGH,SCALE)
      CALL ELMHES (NT,NT,LOW,IGH,H0,IANA)
      CALL HQRM2(NT,NT,LOW,IGH,H0,X,Y,IERR,100)
      
      IF (IERR.NE.0) THEN
          IF (ISDEBUG) WRITE(*,*)' HQRM2 IERR = ',IERR
          STOP
      ENDIF

  12  FORMAT (4E18.11)
      IL = 0
      
      DO 90 I = 1,NT
   90 IJ(I) = 1
      
C     ВЫЧИСЛЕНИЕ МИНИМАЛЬНОГО ПО МОДУЛЮ СОБСТВЕННОГО ЗНАЧЕНИЯ
C     SEARCH OF THE ABS(MAX) EIGENVALUE
      
  100 IL = IL+1
      P = 1.D25
      
      DO 33 I = 1,NT
          IF(IJ(I).EQ.0) GO TO 33
          Q = SQRT(X(I)**2+Y(I)**2)
          IF(Y(I).LT.0.D0) GO TO 33
          IF(Q.LT.P) THEN
              P = Q
              II = I
          ENDIF
   33 CONTINUE
      
      CZ = CMPLX(X(II),Y(II))
      LAMDA(IL) = CZ
      
C     ПРОВЕРКА УСЛОВИЯ УСТОЙЧИВОСТИ
C     CHECKING THE STABILITY CONDITION
      
      
      FL = REAL(CZ)*K**2-RO*H*DIMAG(CZ)**2
      FLG(IL) = FL
      IF(IL.LT.NLAMDA) THEN
          IJ(II) = 0
          IF(Y(II).NE.0.D0) IJ(II+1) = 0
              GO TO 100
      ENDIF
   52 FORMAT (9E8.1)
      F = FL
      IF(IPR.EQ.2) THEN
          IF (ISDEBUG) WRITE (*,120) F,XX
      ENDIF
  120 FORMAT ('F =',E12.4,'V =',E12.4)
      RETURN
      END



C     ---------------------------------------------------------------------------
      SUBROUTINE FLUT33D
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(160000),B(2500),DX(2500),DY(2500),H0(6250000),
     1H1(6250000),X(2500),Y(2500),INT(2500),IJ(2500),URR(2500),UI(2500),
     2ER(50),YY(50),Z(50),Y2(51),SCALE(2500),ZZ(6250000),RM1(6250000),
     3RV1(2500),RV2(2500),XGRID(50),YGRID(50)
      DATA PI/3.141592653589D0/
      LOGICAL SELECT(1600)
      LOGICAL WERESELECTED
      REAL*8 K
      COMPLEX*16 CZ,LAMDA(40)
      character*26 PATH 
      CHARACTER*1 KEY,KEY1
      EQUIVALENCE (H1(1),ZZ(1)),(A(1),RM1(1))
      COMMON /COM1/ IPR,NUMB
      COMMON /DATA/ K,H,BL,RO,DB
      COMMON /COM2/ NLAMDA,TJ,M,N,NT
      COMMON /COM3/ V,IN
      COMMON /DEBUG/ ISDEBUG
      COMMON //PATH      
      
      WERESELECTED = .FALSE.

      IF (ISDEBUG) WRITE (*,*) 'AUTOMATIC FINDING OF EIGEN VECTOR (Y/Q)'
      
      KEY1 = 'Y'
      IF (ISDEBUG) READ (*,10) KEY1
   10 FORMAT(A)
      IF(KEY1.EQ.'Q') RETURN
      IF(KEY1.EQ.'Y') GO TO 110

  110 CONTINUE
      
      CALL MATRH1 (M,N,NT,H0,A,B,DX,DY,H1,TJ,BL,K,DB,V,NUMB)
      
C     Вычисление собственных значений матрицы H0
C     
      
      MMI = NT
      MI = 2
      NOUT = 4
      OPEN(UNIT=NOUT,FILE=TRIM(ADJUSTL(PATH))//'NOUT')
      IF (ISDEBUG) WRITE (*,75) V,IN
      WRITE (NOUT,75) V,IN
   75 FORMAT(5H V = ,E13.6,5H IN =,I3)
      OPEN (8,FILE=TRIM(ADJUSTL(PATH))//'MATRA')
      
      DO 170 I=1,NT
  170 SELECT(I)=.FALSE.
      
      CALL BALANC(NT,NT,H0,LOW,IGH,SCALE)
      CALL ELMHES(NT,NT,LOW,IGH,H0,INT)
      
      WRITE(8,12) (H0(I),I=1,NT*NT)
      ENDFILE 8

      CALL HQRM2(NT,NT,1,NT,H0,X,Y,IERR,200)

      IF (ISDEBUG) WRITE(*,*) 'IERR = ',IERR
      
  12  FORMAT (4E18.11)
      
      WRITE(NOUT,*) 'REAL PART OF EIGENVALUES'
      WRITE(NOUT,12) (X(I),I=1,NT)
      WRITE(NOUT,*) 'IMAGINARY PART OF EIGENVALUES'
      WRITE(NOUT,12) (Y(I),I=1,NT)
      
      IS=0
      
      DO I=1,NT
          IF (Y(I).GT.0.) IS=IS+1
      ENDDO
      
      IF (ISDEBUG) WRITE(*,120) IS
  120 FORMAT(I3,' COMPLEX COUPLES')
      WRITE(NOUT,120) IS
      IL=0

      DO 90 I=1,NT
   90 IJ(I)=1
      
C     Вычисление минимального по модулю собственного значения
C     
      
  100 IL=IL+1
      WRITE(NOUT,*) 'IL = ',IL
      IF (ISDEBUG) WRITE(*,*) 'IL = ',IL

      P = 1.D25
      
      DO 33 I=1,NT
          IF(IJ(I).EQ.0) GO TO 33
          Q = SQRT(X(I)**2+Y(I)**2)
          IF(Y(I).LT.0.D0) GO TO 33
          IF(Q.LT.P) THEN
              P = Q
              II = I
          ENDIF
   33 CONTINUE
      
      IF (ISDEBUG) WRITE(*,*) 'II = ',II
      CZ = CMPLX(X(II),Y(II))

      LAMDA(IL) = CZ
      WRITE(NOUT,115) IL
  115 FORMAT(I3,' - eigenvalue')
      WRITE(NOUT,*) ' '
      WRITE(NOUT,12) CZ
      IF (ISDEBUG) WRITE(*,115) IL
      IF (ISDEBUG) WRITE(*,*) ' '
      IF (ISDEBUG) WRITE(*,12) CZ
      
C     Проверка условия устойчивости
C     
      
      FL = REAL(CZ)*K**2-RO*H*DIMAG(CZ)**2
      IF (ISDEBUG) WRITE(*,*)' CHECK OF STABILITY CONDITION'
      IF (ISDEBUG) WRITE(*,12) FL
      WRITE(NOUT,*)' CHECK OF STABILITY CONDITION'
      WRITE(NOUT,12) FL

      IF (.FALSE.) THEN
          WRITE(*,*)' CALCULATE THIS EIGENVECTOR (Y/N)'
          READ(*,15) KEY
  15  FORMAT(A)
          IF(KEY.EQ.'Y') SELECT(II) = .TRUE.
          IF(KEY.EQ.'N') GO TO 400
      ELSE
          IF (WERESELECTED) THEN
              GO TO 400
          ELSE
              IF (IL.EQ.IN) THEN
                  SELECT(II) = .TRUE.
                  WERESELECTED = .TRUE.
                  !WRITE(*,*)IL, IN, 'TRUE'
              ELSE
                  GO TO 400
              ENDIF
          ENDIF
      ENDIF
            
      REWIND 8
      READ(8,12) (H0(I),I=1,NT*NT)
      
      CALL INVIT(NT,NT,H0,X,Y,SELECT,MMI,MI,ZZ,IERR,RM1,RV1,RV2)
      
      IF (ISDEBUG) WRITE(*,*) 'IERR = ',IERR
      
      IF (IERR.NE.0) THEN
          IF (ISDEBUG) WRITE(*,*) 'MISTAKE IN THE SUBPROGRAMME OF
     * CALCULATION OF THE EIGENVECTOR'
          IF (ISDEBUG) WRITE(*,*) ' IERR = ', IERR
          STOP
      ENDIF
      
      REWIND 8
      READ(8,12) (H0(I),I=1,NT*NT)
      
      CALL ELMBAK(NT,LOW,IGH,H0,INT,MI,ZZ)
      CALL BALBAK(NT,NT,LOW,IGH,SCALE,MI,H0)
      
      M11 = 2*M
      DO 22 I = 1,NT                                                    ZAP03500
      URR(I) = ZZ(I)
  22  UI(I) = ZZ(NT+I)
      OPEN (9,FILE=TRIM(ADJUSTL(PATH))//'FDATA1')
      NREZ1=18
      OPEN(UNIT=NREZ1,FILE=TRIM(ADJUSTL(PATH))//'REZULT.txt')
      
      DO  I=1,M
         T=(2.D0*I-1.D0)*PI/2.D0/M
         XGRID(I)=COS(T)
      ENDDO
      DO  I=1,N
         T=(2.D0*I-1.D0)*PI/2.D0/N
         YGRID(I)=BL*COS(T) 
      ENDDO    
      
      IJ = 0
            
      DO I = 1,M
          DO J = 1,N
              IJ=IJ+1
              WRITE (NREZ1,150) XGRID(I),YGRID(J),URR(IJ)
          ENDDO
      ENDDO
                  
150   FORMAT(1X,F9.5,F9.5,F9.5)                 
      WRITE (NOUT,*) 'EIGENVECTOR RE'
      WRITE (NOUT,12) (URR(I),I=1,NT)
      WRITE (9,12) (URR(I),I=1,NT)
      WRITE (NOUT,*) 'EIGENVECTOR IM'
      WRITE (NOUT,12) (UI(I),I=1,NT)
      
C     Обязательно M нечетное
C     
      
      DO 40 I=1,N
40    YY(I)=URR(N*(((M-1)/2))+I)
      DO 20 LL=1,21
      X2=-1.+0.1*(LL-1)
20    Y2(LL)=EIGEN(X2,YY,Z,N,-1.D0,+1.D0)
      CALL NORM (Y2,21)
      WRITE(NOUT,*) 'EIGENVECTOR RE ON THE AXIS X IN 21 POINTS'
      WRITE(NOUT,12) Y2
      IF (ISDEBUG) PRINT 12,Y2
      IF(DIMAG(CZ).EQ.0.D0) GO TO 200
      DO 44 I=1,N
44    YY(I)=UI(N*(((M-1)/2))+I)
      DO 21 LL=1,21
      X2=-1.+0.1*(LL-1)
21    Y2(LL)=EIGEN(X2,YY,Z,N,-1.D0,+1.D0)
      CALL NORM (Y2,21)
      WRITE(NOUT,*) 'EIGENVECTOR IM ON THE AXIS X IN 21 POINTS'
      WRITE(NOUT,12) Y2
      IF (ISDEBUG) PRINT 12,Y2
  200 CONTINUE
      IF (ISDEBUG) WRITE (*,*) 'CHECK OF SYMMETRY OF THE EIGENVECTOR'
      WRITE (NOUT,*) 'CHECK OF SYMMETRY OF THE EIGENVECTOR'
      NM=(M-1)/2
      DO 50 I=1,NM
      I1=(I-1)*N
      I2=(M-I)*N
      P=0.D0
      DO 51 J=1,N
      Q=URR(I1+J)-URR(I2+J)
      IF(ABS(Q).GT.P) P=ABS(Q)
   51 CONTINUE
   50 ER(I)=P
      IF (ISDEBUG) WRITE(*,52) (ER(I),I=1,NM)
      WRITE(NOUT,52) (ER(I),I=1,NM)
      IF(DIMAG(CZ).EQ.0.D0) GO TO 300
      DO 60 I=1,NM
      I1=(I-1)*N
      I2=(M-I)*N
      P=0.D0
      DO 61 J=1,N
      Q=UI(I1+J)-UI(I2+J)
      IF(ABS(Q).GT.P) P=ABS(Q)
   61 CONTINUE
   60 ER(I)=P
      IF (ISDEBUG) WRITE(*,52) (ER(I),I=1,NM)
      WRITE(NOUT,52) (ER(I),I=1,NM)
  300 CONTINUE
  400 CONTINUE
      IF(IL.LT.NLAMDA) THEN
      IJ(II)=0
      IF(Y(II).NE.0.D0) IJ(II+1)=0
      GO TO 100
      ENDIF
      WRITE (NOUT,*) 'LAMDA'
      IF (ISDEBUG) WRITE (*,*) 'LAMDA'
      WRITE (NOUT,70) (LAMDA(I),I=1,NLAMDA)
      IF (ISDEBUG) WRITE (*,70) (LAMDA(I),I=1,NLAMDA)
   70 FORMAT(1X,2E18.11)
      RETURN
   52 FORMAT (9E8.1)
      END
      
      
C     ---------------------------------------------------------------------------
      SUBROUTINE FLUT16C
C 23.09.96 Флаттер тетрадь N2,стр. 129
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION U(2500),Y(2500),Z(2500),Y2(101),X2(101)
      REAL*8 K
      character*26 PATH
      COMMON /DATA/ K,H,BL,RO,DB
      COMMON /COM2/ NLAMDA,TJ,M,N,NT
      COMMON //PATH
      
      PI=  3.14159265359D0
      
      OPEN (9,FILE=TRIM(ADJUSTL(PATH))//'FDATA1')
      REWIND 9
      
      NG = 101
      READ(9,12) (U(I),I=1,NT)

  12  FORMAT (4E18.11)
      
      DO 40 I=1,N
40    Y(I) = U(N*(((M-1)/2))+I)
      DO 20 L = 1,NG
          X = -1.+0.02*(L-1)
          X2(L) = X
  20  Y2(L) = EIGEN(X,Y,Z,N,-1.D0,+1.D0)
      CALL NORM1(Y2,NG)
      OPEN(5,FILE=TRIM(ADJUSTL(PATH))//'RDATA1')
      WRITE (5,13) X2
      WRITE (5,13) Y2
  13  FORMAT(1X,4E19.11)

      I1 = (N+1)/2
      DO 41 I=1,M
          Y(I) = U(I1)
  41  I1 = I1+N
      STEP = 2.*BL/(NG-1)
      DO 21 L=1,NG
          X = -BL+STEP*(L-1)
          X2(L) = X
  21  Y2(L) = EIGEN(X,Y,Z,M,-BL,BL)
      CALL NORM1(Y2,NG)
      WRITE (5,13) X2
      WRITE (5,13) Y2
      RETURN
      END                                                               ZAP01380
      

C     ---------------------------------------------------------------------------
      SUBROUTINE NORM1(Y,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Y(1)
      P=0.D0
      DO 1 I=1,N
      IF (ABS(Y(I)).GT.P) IP=I
      IF (ABS(Y(I)).GT.P) P=ABS(Y(I))
  1   CONTINUE
      P=Y(IP)
      DO 2 I=1,N
  2   Y(I)=Y(I)/P
      RETURN
      END
      
      
C     ---------------------------------------------------------------------------
      SUBROUTINE FLUT16D
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION U(2500),Y(2500),Z(2500),X2(101),W2(101,101),
     *Z2(4040)
      REAL*8 K
      character*26 PATH
      COMMON /DATA/ K,H,BL,RO,DB
      COMMON /COM2/ NLAMDA,TJ,M,N,NT
      COMMON /DEBUG/ ISDEBUG
      COMMON //PATH
      NG=101
      PI=3.14159265359D0
      OPEN (9,FILE=TRIM(ADJUSTL(PATH))//'FDATA1')
      REWIND 9
      READ(9,12) (U(I),I=1,NT)
      
  12  FORMAT (4E18.11)

      DO 50 J=1,M
      DO 40 I=1,N
40    Y(I)=U(N*(J-1)+I)
      DO 20 L=1,NG
      X=-1.+0.02*(L-1)
      X2(L)=X

  20  Z2((L-1)*M+J)=EIGEN (X,Y,Z,N,-1.D0,+1.D0)
  50  CONTINUE

      OPEN(7,FILE=TRIM(ADJUSTL(PATH))//'RDATA2')
      WRITE (7,13) X2
  13  FORMAT(1X,4E19.11)

      DO 60 J = 1,NG
      DO 41 I = 1,M
          
  41  Y(I) = Z2((J-1)*M+I)

      STEP = 2.*BL/(NG-1)
      DO 21 L = 1,NG
      X = -BL+STEP*(L-1)
      X2(L) = X

  21  W2(L,J) = EIGEN(X,Y,Z,M,-BL,BL)
  60  CONTINUE

      CALL NORM1(W2,NG*NG)

      WRITE (7,13) X2
      WRITE (7,13) W2

      RETURN
      END
         
                                               
C     ---------------------------------------------------------------------------
      SUBROUTINE MATRH2 (M,N,NT,H0,A,B,H1,BL,NUMB)
C  Вычисление матрицы для свободно опертой и защемленной прямоугольной пластины
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N,N),B(M,M),H0(NT,NT),H1(NT,NT)
      IF (NUMB.NE.1.AND.NUMB.NE.2) THEN
      WRITE (*,*) 'Wrong number of the boundary-value problem. STOP.'
      STOP
      ENDIF
      IF (NUMB.EQ.2) THEN
      CALL BOXBIG (H1,N,M,NT,BL)
      GO TO 110
      ENDIF
C     Построение матрицы H0
      CALL ZSL1 (A,N,-1.D0,+1.D0)
      CALL ZSL1 (B,M,-BL,BL)
      DO 1 I=1,NT
      DO 1 J=1,NT
    1 H0(I,J)=0.D0
      DO 2 I=1,M
      DO 2 J=1,M
      DO 2 NU=1,N
      I1=(I-1)*N+NU
      DO 2 MU=1,N
      J1=(J-1)*N+MU
      IF(I.EQ.J) THEN
      H0(I1,J1)=A(NU,MU)
      IF(NU.EQ.MU) H0(I1,J1)=H0(I1,J1)+B(I,I)
      ENDIF
      IF(I.NE.J) THEN
      IF(NU.EQ.MU) H0(I1,J1)=B(I,J)
      ENDIF
    2 CONTINUE
C Вычисление квадрата матрицы H0
      DO 3 I=1,NT
      DO 3 J=1,NT
      P=0.D0
      DO 4 L=1,NT
    4 P=P+H0(I,L)*H0(L,J)
    3 H1(I,J)=P
C
  110 CONTINUE
      RETURN
      END
      
      
C     ---------------------------------------------------------------------------
      SUBROUTINE MATRH3 (M,N,NT,H0,DX,DY,H1,TJ,BL,K,DB,V)
C Вычисление матрицы для свободно опертой и защемленной прямоугольной пластины
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DX(N,N),DY(M,M),H0(NT,NT),H1(NT,NT)
      REAL*8 K
C
      PI=3.141592653589D0
C
      PSI=PI*TJ
      VX= V*COS(PSI)
      VY= V*SIN(PSI)
C
      CALL DIFX1(DX,N,-1.D0,1.D0)
      CALL DIFX1(DY,M,-BL,BL)
C
      DO 5 I=1,NT
      DO 5 J=1,NT
    5 H0(I,J)=0.D0
      DO 6 I=1,M
      DO 6 J=1,M
      DO 6 NU=1,N
      I1=(I-1)*N+NU
      DO 6 MU=1,N
      J1=(J-1)*N+MU
      IF(I.EQ.J) THEN
      H0(I1,J1)=VX*DX(NU,MU)
      IF(NU.EQ.MU) H0(I1,J1)=H0(I1,J1)+VY*DY(I,I)
      ENDIF
      IF(I.NE.J) THEN
      IF(NU.EQ.MU) H0(I1,J1)=VY*DY(I,J)
      ENDIF
    6 CONTINUE
      DO 7 I=1,NT
      DO 7 J=1,NT
    7 H0(I,J)=DB*H1(I,J)-K*H0(I,J)
      RETURN
      END


C     ---------------------------------------------------------------------------
      SUBROUTINE MATRH1(M,N,NT,H0,A,B,DX,DY,H1,TJ,BL,K,DB,V,NUMB)
C Вычисление матрицы для свободно опертой и защемленной прямоугольной пластины
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N,N),B(M,M),DX(N,N),DY(M,M),H0(NT,NT),H1(NT,NT)
      REAL*8 K
      IF (NUMB.NE.1.AND.NUMB.NE.2) THEN
      WRITE (*,*) 'ЌҐЇа ўЁ«м­л© ­®¬Ґа Єа Ґў®© § ¤ зЁ. STOP.'
      STOP
      ENDIF
      IF (NUMB.EQ.2) THEN
      CALL BOXBIG (H1,N,M,NT,BL)
      GO TO 110
      ENDIF
C     Построение матрицы H0
      CALL ZSL1 (A,N,-1.D0,+1.D0)
      CALL ZSL1 (B,M,-BL,BL)
      DO 1 I=1,NT
      DO 1 J=1,NT
    1 H0(I,J)=0.D0
      DO 2 I=1,M
      DO 2 J=1,M
      DO 2 NU=1,N
      I1=(I-1)*N+NU
      DO 2 MU=1,N
      J1=(J-1)*N+MU
      IF(I.EQ.J) THEN
      H0(I1,J1)=A(NU,MU)
      IF(NU.EQ.MU) H0(I1,J1)=H0(I1,J1)+B(I,I)
      ENDIF
      IF(I.NE.J) THEN
      IF(NU.EQ.MU) H0(I1,J1)=B(I,J)
      ENDIF
    2 CONTINUE
C Вычисление квадрата матрицы H0
      DO 3 I=1,NT
      DO 3 J=1,NT
      P=0.D0
      DO 4 L=1,NT
    4 P=P+H0(I,L)*H0(L,J)
    3 H1(I,J)=P
C
  110 CONTINUE
C
      PI=3.141592653589D0
C
      PSI=PI*TJ
      VX= V*COS(PSI)
      VY= V*SIN(PSI)
C
      CALL DIFX1(DX,N,-1.D0,1.D0)
      CALL DIFX1(DY,M,-BL,BL)
C
      DO 5 I=1,NT
      DO 5 J=1,NT
    5 H0(I,J)=0.D0
      DO 6 I=1,M
      DO 6 J=1,M
      DO 6 NU=1,N
      I1=(I-1)*N+NU
      DO 6 MU=1,N
      J1=(J-1)*N+MU
      IF(I.EQ.J) THEN
      H0(I1,J1)=VX*DX(NU,MU)
      IF(NU.EQ.MU) H0(I1,J1)=H0(I1,J1)+VY*DY(I,I)
      ENDIF
      IF(I.NE.J) THEN
      IF(NU.EQ.MU) H0(I1,J1)=VY*DY(I,J)
      ENDIF
    6 CONTINUE
      DO 7 I=1,NT
      DO 7 J=1,NT
    7 H0(I,J)=DB*H1(I,J)-K*H0(I,J)
      RETURN
      END


C     ---------------------------------------------------------------------------
      SUBROUTINE BOXBIG (H,N,M,NT,BL)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER P
      DIMENSION H(NT,NT),TP(50),TP1(50),TP2(50),TP3(50),TP4(50),
     *TK(50),TK1(50),TK2(50),TK3(50),TK4(50)
      COMMON /D/ Dx,Dxy,Dy
      IF (N.GT.50.OR.M.GT.50) THEN
          !WRITE (*,*) 'M OR N IS MORE THAN 50. STOP.'
          STOP
      ENDIF
      
C     ПОСТРОЕНИЕ МАРИЦЫ H
C     BUILDING MATRIX H
      !WRITE(*,*) 'BL =',BL
      !PAUSE
      PI = 3.141592653589D0
      I1 = 0
      DO 2 MU = 1,M
      Z = COS((2.*MU-1.)*PI/2./M)
      Z2 = Z*Z
      Z3 = Z2*Z
      Z4 = Z3*Z
      CALL DIFT (Z,M,TP,TP1,TP2,TP3,TP4)
      DO 2 NU = 1,N
      X = COS((2.*NU-1.)*PI/2./N)
      X2 = X*X
      X3 = X2*X
      X4 = X3*X
      I1 = I1+1
      CALL DIFT (X,N,TK,TK1,TK2,TK3,TK4)
      I2 = 0
      DO 2 I = 1,M
      TI = (2.*I-1.)*PI/2./M
      STI = SIN(TI)**4
      DO 2 J = 1,N
      TJ = (2.*J-1.)*PI/2./N
      STJ = SIN(TJ)**4
      I2 = I2+1
      SUM = 0.D0
      DO 3 K = 0,N-1
      TKJ = COS(K*TJ)
      IF (K.EQ.0) TKJ = 0.5*TKJ
      D2DX2 = (12.*X2-4.)*TK(K+1)+(8.*X3-8.*X)*TK1(K+1)+
     +(X4-2.*X2+1.)*TK2(K+1)
      D4DX4 = 24.*TK(K+1)+96.*X*TK1(K+1)+(72.*X2-
     -24.)*TK2(K+1)+(16.*X3-16.*X)*TK3(K+1)+
     +(X4-2.*X2+1.)*TK4(K+1)
      DO 3 P=0,M-1
      TPI=COS(P*TI)
      IF (P.EQ.0) TPI=0.5*TPI
      D2DZ2=(12.*Z2-4.)*TP(P+1)+(8.*Z3-8.*Z)*TP1(P+1)+
     +(Z4-2.*Z2+1.)*TP2(P+1)
      D4DZ4=24.*TP(P+1)+96.*Z*TP1(P+1)+(72.*Z2-
     -24.)*TP2(P+1)+(16.*Z3-16.*Z)*TP3(P+1)+
     +(Z4-2.*Z2+1.)*TP4(P+1)
      DELTA=Dx*D4DX4*(Z*Z-1.)**2*TP(P+1)+2.0*Dxy*D2DX2*D2DZ2/BL/BL+
     +(X*X-1.)**2*TK(K+1)*Dy*D4DZ4/BL**4
    3 SUM=SUM+TKJ*TPI*DELTA/STI/STJ
    2 H(I1,I2)=4.*SUM/M/N
      DO 5 I=1,NT
      DO 6 J=I+1,NT
C     IF(J.GE.NT) GO TO 5
      H(J,I)=(H(J,I)+H(I,J))/2.
    6 H(I,J)=H(J,I)
    5 CONTINUE
      RETURN
      END
     