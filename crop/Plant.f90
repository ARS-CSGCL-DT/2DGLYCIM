!****************************( ACTGRO )*********************************
!*                                                                     *
!* Uses the initial carbon allocation plus any limitations occasioned  *
!* by water stress or nutrient stress to calculate the actual growth   *
!* in size and dry weight of all organs on the plant. Root growth      *
!* occurs in the soil cells where conditions are most favorable and    *
!* is shared with adjacent cells. Initiates branches, pods, and seeds  *
!* according to the plant's ability to support them.                   *
!*                                                                     *
!***********************************************************************
      SUBROUTINE ACTGRO
      INCLUDE 'common.h'
	  INCLUDE 'plant.h'
	  !DEC$ATTRIBUTES DLLEXPORT :: crop, /ShootR/,/shtR_public/, &
		/Weath/, /grid_public/, /nodal_public/, /elem_public/,   &
	    /bound_public/, /time_public/, /module_public/,          &
	     /DataFilenames/   
      REAL dPopt, dP
	  COMMON/ACT1/D17, dPopt, dP

      IF ( INIT ) THEN
         INITR1 = 0
         INITR3 = 0
         INITR4 = 0
         MBRNCH = 0
         SCRTS2 = 0.0
         RSTPCS = 0.0
      END IF

      CBAL1 = SCPOOL + RSTPCS + RCPOOL + USEDCS + XTRACS +		    &
             (PERIOD*(VEGSRC + FRUITC + RSTPC))

! NUMBER OF TOP LEAF ON MAIN STEM.

      IF (ITIME.LT.2) THEN
      NTOPLF = IFIX(VSTAGE + 0.5)
      M11 = INT(VSTMAH)
      IF (NTOPLF.GT.M11) NTOPLF = M11
      END IF

! INITIATE FLOWERS AT R1, PODS AT R3 AND SEEDS AT R4
! FROM R1 TO R3 PLANTS OPEN AN AVERAGE OF 4 FLOWERS PER NODE

      IF (RSTAGE.GE.1.0.AND.INITR1.EQ.0) THEN
         NODES = INT(VSTMAH)
         DO J = 1,NBRNCH
            NODES = NODES + INT(TBMAH(J))
         END DO
! Note 12.01
         NFLWRS = (4*NODES)
         INITR1 = 1
      END IF

! AT R3 PLANTS INITIATE PODS (1 PER 2 NODES AT FIRST)

      IF (RSTAGE.GE.3.0.AND.INITR3.EQ.0) THEN
         NPODS = 1
         SPCSP = 0.0001
         PDWT(1) = 0.001
         PPDWT(1) = 0.001
! PODS = NODES/2.0
! POD(1) = PODS
         INITR3 = 1
      ENDIF

! AT R4 PLANTS INITIATE SEEDS (2 PER POD AT FIRST)

      IF (RSTAGE.GE.4.0.AND.INITR4.EQ.0) THEN
       NSEEDS =  2.5*NPODS
       SPCSS = 0.001
       DO I=1,NSEEDS
       SDWT(I) = 0.001
       PSDWT(I)= 0.001
       ENDDO
!	   SPCSS = 0.001
! NSEEDS = PODS*2.0
! SEED(1) = NSEEDS
      INITR4 = 1
      ENDIF

! GROW EXISTING PODS AND SEEDS
!YA ADDED <6.0 
      IF (RSTAGE.GE.3.0.AND.RSTAGE.LT.6.0) THEN

! FROM R5 ONWARDS SCPOOL IS ACTIVELY DEPLETED TO GROW PODS AND SEEDS ON DEMAND

      IF (RSTAGE.GT.5.0) THEN
        D29 = PDSC + PDPODC
       IF (FRUITC.LT.D29) THEN
         D11 = SCPOOL - ((D29 - FRUITC)*PERIOD)
        IF (D11.LT.0.0) THEN
          FRUITC = FRUITC + (SCPOOL/PERIOD)
          SCPOOL = 0.0
        ELSE
          FRUITC = D29
          SCPOOL = D11
        END IF
      END IF
      END IF

! GROW EXISTING PODS AND SEEDS

      USEDCS = USEDCS + (FRUITC*PERIOD)
      SDWTR  = FRUITC * PERIOD/CONVSE
      IF(NSEEDS.GT.0.0) THEN
      SEEDCD = 0.0
      SEEDWT = 0.0
      D120 = PDSC/NSEEDS*PERIOD/CONVSE
        DO 13 J = 1, NSEEDS
          PSDWT(J) = PSDWT(J) + D120
          D121 = PSDWT(J) - SDWT(J)
          IF(SDWTR.GE.D121)THEN
            SDWT(J) = PSDWT(J)
            SDWTR = SDWTR- D121
          ELSE IF (SDWTR.GT.0.0)THEN
             SDWT(J) = SDWT(J) + SDWTR
             SDWTR = 0.0
             SEEDCD = SEEDCD+(PSDWT(J)-SDWT(J))*CONVSE
          ELSE
             SEEDCD = SEEDCD+(PSDWT(J)-SDWT(J))*CONVSE
          END IF
          SEEDWT = SEEDWT + SDWT(J)
13      CONTINUE
      YLD = SEEDWT * POPAREA * 0.171
      END IF
! Note 12.02
         IF (SDWTR .GT. 0.0) THEN
            PDWTR =SDWTR*CONVSE/CONVPO
            PODCD = 0.0
            PODWT = 0.0
            D122 = PDPODC/NPODS*PERIOD/CONVPO
            DO J = 1, NPODS
               PPDWT(J)=PPDWT(J)+D122
               D123 = PPDWT(J) - PDWT(J)
               IF (PDWTR .GE. D123) THEN
                  PDWT(J) = PPDWT(J)
                  PDWTR = PDWTR - D123
               ELSE IF (PDWTR .GT. 0.0) THEN
                  PDWT(J) = PDWT(J)+PDWTR
                  PDWTR = 0.0
                  PODCD=PODCD+(PPDWT(J)-PDWT(J))*CONVPO
               ELSE
                  PODCD=PODCD+(PPDWT(J)-PDWT(J))*CONVPO
               END IF
               PODWT = PODWT + PDWT(J)
            END DO
         END IF
!YA
       PODWT6=PODWT
       RSTAGEold=RSTAGE
! ENDIF
       ELSEIF(RSTAGE.GT.6.0.AND.RSTAGE.LT.7.0) THEN
         dPopt = 7.0/(POPAREA*0.171)
         dPopt= AMIN1(dPopt,0.5*PODWT6)
         dP=(RSTAGE-RSTAGEold)*SCF*dPopt
         RSTAGEold=RSTAGE
         PODWT=PODWT-dP
         SEEDWT=SEEDWT+dP
         YLD = SEEDWT * POPAREA * 0.171
       ENDIF

! CALCULATE CARBON AVAILABILITY AND LIMITATIONS TO THE CHANGE IN
! SIZE AND WEIGHTS OF VEGETATIVE ORGANS ON THE SHOOT

! IF THERE IS CARBON LEFT AFTER GROWING THE PODS AND SEEDS,
! CALCULATE HOW MUCH IS AVAILABLE TO GROW THE REST OF THE SHOOT.
	
      REMC = 0.0
      SCRATO = 0.0
! Note 12.03
      CASG = (1.0 - PARTRT)*VEGSRC*PERIOD
      IF (SGT .LT. 1.0) CASG = CASG*SGT 
      IF (PDV .LT. 1E-10) THEN
      SCRATO = 1.0
      REMC = CASG
      ELSE
      IF (CASG .GT. 0.0) THEN
      REMC = CASG
      SCRATO = 1.0
      IF ((PCST.GT.0.).OR.(PCST.LT.0.)) THEN
         REMC = 0.0
         USEDCS = USEDCS + CASG

! CALCULATE MAXIMUM AND MINIMUM AMOUNTS OF CARBON NEEDED TO
! GROW SHOOTS AT POTENTIAL RATE LIMITED ONLY BY TEMPERATURE
! AND SHOOT TURGOR.                                               *** Note 12.04
         CMAXSG = PCST*SGT*PERIOD
         CMINSG = QCST*SGT*PERIOD

! DETERMINE IF CARBON AVAILABILITY FURTHER LIMITS SHOOT GROWTH
! AND CALCULATE A SHOOT CARBON SUPPLY/DEMAND RATIO.

         IF (CASG.LT.CMINSG) THEN
! Note 12.05
            SGTR = SGT*CASG/CMINSG
            SCRATO = 0.0
         ELSE
            SGTR = SGT
            IF (CASG .GE. CMAXSG) THEN
               SCRATO = 1.0
               USEDCS = USEDCS - CASG + CMAXSG
               REMC = CASG - CMAXSG
            ELSE
! Note 12.06
            SCRATO = (CASG - CMINSG)/(CMAXSG - CMINSG)
            END IF
		 END IF
! ------------------------------------------------------------------
! CALCULATE ACTUAL CHANGES IN SIZE AND DRY WEIGHT OF ORGANS
! ON THE SHOOT DURING THIS PERIOD.

! CALCULATE CHANGES IN MAINSTEM LEAF AND PETIOLE SIZES AND DRY WEIGHTS.
!YA        IF (VSTAGE.LT.VSTMAH) THEN
           IF (RSTAGE.LE.PARM(16)) THEN
         DO I = 1,3
            D15 = VSTAGE + 1.5 - I
            IF (D15 .LT. 0.0) EXIT
            ITRIF = IFIX(D15)
            ADLA = PDLAM(I)*SGTR*PERIOD
            ADLW = ADLA*TNL/9.0*(1.0 + (0.6*SCRATO))
            ADPL = PDPLM(I)*SGTR*PERIOD
            ADPW = PDPWM(I)*SGTR*PERIOD*(1.0 + (0.6*SCRATO))
            IF (ITRIF .LT. 1) THEN
               ULAREA = ULAREA + ADLA
               ULEAFW = ULEAFW + ADLW
               UPETL = UPETL + ADPL
               UPETW = UPETW + ADPW
               LEAFWT = LEAFWT + ADLW
               PETWT = PETWT + ADPW
               LAREAM = LAREAM + ADLA
               EXIT
            END IF
            MLAREA(ITRIF) = MLAREA(ITRIF) + ADLA
            MLEAFW(ITRIF) = MLEAFW(ITRIF) + ADLW
            MPETL(ITRIF) = MPETL(ITRIF) + ADPL
            MPETW(ITRIF) = MPETW(ITRIF) + ADPW
            LEAFWT = LEAFWT + ADLW
            PETWT = PETWT + ADPW
            LAREAM = LAREAM + ADLA
         END DO
         END IF

! CALCULATE CHANGES IN BRANCH LEAF AND PETIOLE DIMENSIONS AND
!    DRY WEIGHTS.

         IF (NBRNCH.NE.0) THEN
            DO J = 1,NBRNCH
               IF (RSTAGE.LE.PARM(16)) THEN
                  DO II = 1,3
                     I = IFIX(TB(J) + 1.5 - II)
                     IF (I.LE.0) EXIT
		             ADLA = PDLAB(I,J)*SGTR*PERIOD
					 ADLW = ADLA*TNL/12.0*(1.0 + (0.6*SCRATO))
		             ADPW = PDPWB(I,J)*SGTR*PERIOD*(1.0 + (0.6*SCRATO))
					 BLAREA(I,J) = BLAREA(I,J) + ADLA
		             BLEAFW(I,J) = BLEAFW(I,J) + ADLW
					 BPETL(I,J) = BPETL(I,J) + (PDPLB(I,J)*SGTR*PERIOD)
		             BPETW(I,J) = BPETW(I,J) + ADPW
					 LEAFWT = LEAFWT + ADLW
		             PETWT = PETWT + ADPW
				     LAREAB = LAREAB + ADLA
                  END DO
		       END IF
            END DO
         END IF
         LAREAT = LAREAM + LAREAB
         LAI = LAREAT*POPAREA*1.0E-4
! CALCULATE CHANGES IN MAINSTEM AND BRANCH LENGTHS AND DRY WEIGHTS.
         IF((MG.GT.4.AND.RSTAGE.LT.2.0).OR.                            &
            (MG.EQ.4.AND.RSTAGE.LE.2.0)) THEN
            ADMH = PDMH*PERIOD
         ELSE
            ADMH = 0.0
         ENDIF

         IF (RSTAGE.LE.PARM(16)) THEN
            ADMW = PDMW*SGTR*PERIOD*(1.0 + (1.5*SCRATO))
            MSTEMH = MSTEMH + ADMH
            MSTEMW = MSTEMW + ADMW
            STEMWT = STEMWT + ADMW
         END IF

         IF (NBRNCH.GT.0) THEN
            DO J = 1,NBRNCH
               IF (TB(J).LT.(TBMAH(J) - 1.0)) THEN
                  ADBW = PDBRW(J)*SGTR*PERIOD*(1.0 +( 1.5*SCRATO))
                  BRNCHH(J) = BRNCHH(J) + ADMH
                  BRNCHW(J) = BRNCHW(J) + ADBW
                  STEMWT = STEMWT + ADBW
               END IF
            END DO
         END IF

! CALCULATE HEIGHT OF TOP LEAVES ABOVE SOIL.

        IF (NTOPLF.LE.0) NTOPLF = 1
        HTNEW = MSTEMH + MPETL(NTOPLF) + 5.0
        IF (HTNEW.GE.HEIGHT) HEIGHT = HTNEW

! CALCULATE CHANGE IN DRY WEIGHT OF FLOWERS.

        FLWRWT = FLWRWT + PDFW*SGTR*PERIOD
      END IF
      END IF
      END IF

! CALCULATE AN ACCUMULATED SHOOT CARBON SUPPLY/DEMAND RATIO

      SCRTS2 = SCRTS2 + (SCRATO*PERIOD)

! DEPENDING ON THE STAGE OF GROWTH AND AVAILABILITY OF SITES
! AND SUBSTRATES, INITIATE BRANCHES, PODS AND SEEDS
! INITIATE BRANCHES IF CARBON AND SITES ARE AVAILABLE.

      IF (ITIME.EQ.IPERD) THEN
        IF (RSTAGE.LT.3.0) THEN
          IF (RSTAGE.LT.0.0) THEN
! NOTE 12.08
            NBSITE = IFIX(VSTAGE) + 1
            D16 = SCRTS2/24.0

! DATA FOR CALCULATION OF BRNCHF CAME FROM 1992 AND 1993 SOYBEAN
! PLANT DENSITY STUDY AT BARC-WEST, BELTSVILLE, MD.

            BRNCHF = (10.81 * EXP(-0.02854*POPAREA)) * PARM(23)
!      IF (MBRNCH.LT.NBSITE.AND.D16.GT.PARM13) THEN
            IF(MBRNCH.LT.NBSITE.AND.MBRNCH.LT.BRNCHF)THEN
! PARM13 WAS ORIGINALLY 0.5
              MBRNCH = MBRNCH + 1
! NOTE 12.09
              INITBR(MBRNCH) = NTOPLF + 5
            END IF
          END IF
          IF (INITBR(NBRNCH + 1).EQ.NTOPLF) THEN
             NBRNCH =  NBRNCH + 1
          END IF

! INITIATE PODS IF CARBON IS AVAILABLE

        ELSEIF (RSTAGE.LT.4.0) THEN
          FPODS = 0
          IF (SVEGC.GT.0.0) THEN
! NOTE 12.10
            FPODS = INT(NPODS*SVEGC/SPCSP*0.30)
!YA         NPODS = NPODS + FPODS
!YA    Fourty pods minimum
            NPODS = MAX(NPODS+FPODS,40)
          END IF
          SVEGC = 0.0
          SPCSP = 0.0

! INITIATE SEEDS IF CARBON IS AVAILABLE
!
        ELSE IF (RSTAGE.LT.5.0) THEN
      FSEEDS = 0
      IF (SVEGC.GT.0.0) THEN
! NOTE 12.11
        FSEEDS = INT(NSEEDS* SVEGC/SPCSS*0.30)
        NSEEDS = NSEEDS + FSEEDS
        IF(NSEEDS.GT.(NPODS*2.5))NSEEDS=NPODS*2.5
      END IF
        SVEGC = 0.0
        SPCSS = 0.0
      END IF
      SCRTS2 = 0.0
      END IF

! BALANCE CARBON ACCOUNTS

      SCPOOL = SCPOOL + REMC
      SHTWT = LEAFWT + PETWT + STEMWT

! GROW ROOTS IN THE MOST FAVORABLE CELLS UNTIL CARBON AVAILABLE
! IS USED.
!wsun Below is root growth which was in 2d soil

      !PCRTS = 0.0
      !SW40 = 0
      !DO J = 1,NCR
         !L = PCRPAR(J)/1000
         !K = PCRPAR(J) - L*1000
         !IF (PDWR(L,K).LE.0.0) THEN
		    !aux = 0.0
!           !ADWR(L,K) = 0.0
         !ELSE
            !PCRTS = PCRTS + PDWR(L,K)*CONVR
            !IF (PCRS.GE.PCRTS) THEN
			   !aux = PDWR(L,K)
!              !ADWR(L,K) = PDWR(L,K)
             !ELSE IF ((SW40.GT.1).OR.(SW40.LT.1)) THEN
               !SW40 = 1
               !aux = (PDWR(L,K)*CONVR - PCRTS + PCRS)/CONVR
!              !ADWR(L,K) = (PDWR(L,K)*CONVR - PCRTS + PCRS)/CONVR
            !ELSE
			   !aux = 0.0
!              !ADWR(L,K) = 0.0
            !END IF
         !END IF

! CALCULATE ROOT GROWTH RATE.

		 !AWR(L,K) = aux * PERIOD
!        !AWR(L,K) = ADWR(L,K)*PERIOD
         !ROOTWT = ROOTWT + AWR(L,K)/POPSLB
      !END DO
!     DO n = 1,NumNP
!	       ROOTWT = ROOTWT + AWR(n)/POPSLB/24
!     END DO
      !PCRS=AMIN1(PCRTS,PCRS)                                      !YAP

! MINIMUM ROOT LENGTH TO CROSS A SOIL CELL DIAGONALLY

      !IF ( INIT ) RTMINW = RTWL*SQRT(DEPTH*DEPTH + WIDTH*WIDTH)

! REDISTRIBUTE ROOT GROWTH IN THE SOIL PROFILE AND CALCULATE
! THE ROOT WEIGHT,NEW ROOT LENGTH AND CO2 RESPIRED IN EACH CELL.

      !NLRP = NLR
      !NKRP = NKR
      !DO L = 1,NLR
         !DO K = 1,NKR
            !ADRL(L,K) = 0.0
            !IF ((AWR(L,K).GT.0.0).OR.(AWR(L,K).LT.0.0)) THEN
! Note 12.13
               !IF (RTWT(L,K) .GT. RTMINW) THEN
                  !LP1 = L + 1
                  !IF (L.EQ.NL) LP1 = NL
                  !IF (LP1.GT.NLR) NLRP = LP1
                  !KP1 = K + 1
                  !IF (K.EQ.NKH) KP1 = NKH
                  !IF (KP1.GT.NKR) NKRP = KP1
                  !LM1 = L - 1
                  !IF (L.EQ.1) LM1 = 1
                  !KM1 = K - 1
                  !IF (K.EQ.1) KM1 = 1
! Note 12.14
                  !D19 = RGCF(L,K) + RGCF(L,KM1) + 1.5*RGCF(L,KP1) +    &
                        !RGCF(LM1,K) + 5.0*RGCF(LP1,K)
                  !IF (.NOT.((D19.GT.0.).OR.(D19.LT.0.))) GOTO 1112
                  !D20 = AWR(L,K)*RGCF(L,K)/D19
                  !RTWT(L,K) = RTWT(L,K) + D20
! Note 12.15
                  !RESPS(L,K) = RESPS(L,K) + D20*CONRES
				  !ADRL(L,K)  = ADRL(L,K) + D20/RTWL
                  !D21 = AWR(L,K)*RGCF(L,KM1)/D19
                  !RTWT(L,KM1)  = RTWT(L,KM1) + D21
                  !RESPS(L,KM1) = RESPS(L,KM1) + D21*CONRES
                  !ADRL(L,KM1) = ADRL(L,KM1) + D21/RTWL
                  !D22 = AWR(L,K)*1.5*RGCF(L,KP1)/D19
                  !RTWT(L,KP1) = RTWT(L,KP1) + D22
                  !RESPS(L,KP1) = RESPS(L,KP1) + D22*CONRES
                  !ADRL(L,KP1) = ADRL(L,KP1) + D22/RTWL
                  !D23 = AWR(L,K)*RGCF(LM1,K)/D19
                  !RTWT(LM1,K) = RTWT(LM1,K) + D23
                  !RESPS(LM1,K) = RESPS(LM1,K) + D23*CONRES
                  !ADRL(LM1,K) = ADRL(LM1,K) + D23/RTWL
                  !D24 = AWR(L,K)*1.5*RGCF(LP1,K)/D19
                  !RTWT(LP1,K) = RTWT(LP1,K) + D24
                  !RESPS(LP1,K) = RESPS(LP1,K) + D24*CONRES
                  !ADRL(LP1,K) = ADRL(LP1,K) + D24/RTWL
                  !CYCLE
               !END IF
 !1112          !RTWT(L,K) = RTWT(L,K) + AWR(L,K)
               !RESPS(L,K) = RESPS(L,K) + AWR(L,K)*CONRES
               !ADRL(L,K) = ADRL(L,K) + AWR(L,K)/RTWL
            !END IF
         !END DO
      !END DO
      !NLR = NLRP
      !NKR = NKRP
!RDenY(n)=0.0D0		density of young root in soil cell cm/cm2 (new version)
! wsun CALCULATE LENGTH OF YOUNG ROOTS IN CELLS
! ADRL(L,K)	cm	Actual change in root length in soil cell L,K for past period (cm).
      ROOTWT = TotalRootWeight/POPSLB
      !DO K = 1,NKR
		DO n = 1,NumNP
! 0.021 = 1/2/24; **** Note 12.16
           YRL(n) = (YRL(n)*(1.0 - (0.021*PERIOD))) + ADRL(n)/24 !ADRL:  cm/day
           !RUTDEN(L,K) = RUTDEN(L,K) +                                 &
		                 !YRL(L,K)*0.021*PERIOD/(DEPTH*WIDTH)
			!YRL(n) = RDenY(n)*nodeArea(n)
         END DO
      !END DO

! CALCULATE AMOUNT OF CARBON IN THE ROOT STORAGE POOL

      D25 = SGT
      IF (D25.GT.1.0) D25 = 1.0
! Note 12.17
      RCPOOL = RCPOOL + (((VEGSRC*(PARTRT +                            &
	          (1.0 - PARTRT)*(1.0 - D25)))                             &
                  - (PCRS/POPSLB/24))*PERIOD)
 !     USEDCS = USEDCS + (PCRS/POPSLB*PERIOD)
 ! wsun PCRS come from 2d soil g slab-1 day-1
       USEDCS = USEDCS + (PCRS/POPSLB/24*PERIOD)
! FROM R3 TO R4, CARBON EQUIVALENT TO THAT REQUIRED TO GROW
! SEEDS IS PLACED IN AN IRRETRIEVABLE SHOOT CARBON STORAGE POOL
! This has been replaced: RSTPCS = RSTPCS + (RSTPC*PERIOD)
!BA
      IF (INITR4.LT.1) THEN
         RSTPCS = RSTPCS + (RSTPC*PERIOD)
      ELSE
         SCPOOL=SCPOOL+RSTPCS
         RSTPCS=0.0
      ENDIF
!BAEND
! REDISTRIBUTE STORED CARBON BETWEEN ROOT AND SHOOT POOLS
      IF ( ITIME .EQ. IPERD .AND. (SCPOOL+RCPOOL).GT.0.0001) THEN
         RCPMAX = 0.4*ROOTWT
         SCPMAX = 0.4*SHTWT
         XTRAC = RCPOOL - RCPMAX + SCPOOL - SCPMAX
         IF (XTRAC.LE.0.0) THEN
            D26 = RCPOOL + SCPOOL
            D27 = D26/(RCPMAX + SCPMAX)*RCPMAX
            IF (D27 .LT. RCPOOL) THEN
               RCPOOL = D27
               SCPOOL = D26/(RCPMAX + SCPMAX)*SCPMAX
			END IF
         ELSE
            RCPOOL = RCPMAX
            SCPOOL = SCPMAX
            XTRACS = XTRACS + XTRAC
         END IF
      END IF
!YAP
! IF THE CROP HAS REACHED RSTAGE GE. 8.0 THEN TERMINATE THE SIMULATION
! AT FIND R8 DATE FOR SUMMARY FILE.
      IF (RSTAGE.GE.8.0) THEN
         MATURE = 1
         R8DATE = CROP_DATE1
      END IF

! RECORD WHICH FACTOR LIMITED GROWTH THIS PERIOD

      D28 = QCST*SGT + PDPODC + PDSC + RSTPC + PDNODC
      Cstress = TRANSC/D28 !Only for carbon stress output
      IF (SGT.LT.1.0) THEN
         LIMITF(ITIME) = 'W'
      ELSE IF (TRANSC.LT.D28) THEN
         LIMITF(ITIME) = 'C'
      ELSE IF (NRATIO.LT.0.99) THEN
         LIMITF(ITIME) = 'N'
      ELSE
         LIMITF(ITIME) = 'T'
      END IF

END