C     Included with permission from W. Bachman.
      SUBROUTINE INFN(ICALL,THETA,DATREC,INDXS,NEWIND)
      INTEGER ICALL,INDXS,NEWIND
      DOUBLE PRECISION THETA
      REAL DATREC
      DIMENSION THETA(*),DATREC(*),INDXS(*)
        IF(ICALL.EQ.3)THEN 
        OPEN(42,FILE='PAR.TAB')
        CALL FILES(42)
        OPEN(43,FILE='COR.TAB')
        CALL FILES(43)
        OPEN(44,FILE='COV.TAB')
        CALL FILES(44)
        CALL WRTAB
        CLOSE(42)
        CALL FILES(42)
        CLOSE(43)
        CALL FILES(43)
        CLOSE(44)
        CALL FILES(44)
      ENDIF
      RETURN
      END
