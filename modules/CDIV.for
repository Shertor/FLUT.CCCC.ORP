
      SUBROUTINE CDIV(AR,AI,BR,BI,CR,CI)
      DOUBLE PRECISION AR,AI,BR,BI,CR,CI
C
C     COMPLEX DIVISION, (CR,CI) = (AR,AI)/(BR,BI)
C
      DOUBLE PRECISION S,ARS,AIS,BRS,BIS
      S = DABS(BR) + DABS(BI)
      ARS = AR/S
      AIS = AI/S
      BRS = BR/S
      BIS = BI/S
      S = BRS**2 + BIS**2
      CR = (ARS*BRS + AIS*BIS)/S
      CI = (AIS*BRS - ARS*BIS)/S
      RETURN
      END
