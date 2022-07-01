! ------------------------------( SHTGRO )------------------------------
! Calculates potential shoot growth assuming that carbohydrate, nitrogen,
! and water are not limiting but allowing for temperature limitation.
! ------------------------------------------------------------------
      SUBROUTINE SHTGRO
      INCLUDE 'common.h'
	  INCLUDE 'plant.h'
	  !DEC$ATTRIBUTES DLLEXPORT :: crop, /ShootR/,/shtR_public/,&
		/Weath/, /grid_public/, /nodal_public/, /elem_public/,  &
	    /bound_public/, /time_public/, /module_public/,         &
	    /error_public/, /DataFilenames/    
      COMMON/SHT1/ALFINT
      INTEGER VSTMAX,TBMAX
      DIMENSION TBMAX(20)

      IF ( INIT ) THEN
         ALFINT = 10000.
         LIMIT = 0
         INITR0 = 0
      END IF
      PDMPC = 0.0
      PDMLC = 0.0
      PDMC = 0.0
      PDBPC = 0.0
      PDBLC = 0.0
      PDBC = 0.0
      PDFC = 0.0
      PDFW = 0.0
      PDSC = 0.0
      PDSW = 0.0
      PDPODC = 0.0

! WHEN R0 IS REACHED, INITIATION OF FLOWERS IN ALL NODES OF
! DETERMINATE PLANTS FIXES THE MAXIMUM NUMBER OF NODES.

      IF (INITR0.EQ.0.AND.RSTAGE.GT.0.0) THEN
! Note 6.06
!VRR NOTE PARM(16) WAS 0.38 FOR INDETERMINATE AND 0.15 FOR DETERMINATE
!YA      VSTMAX = IFIX(8.0*VSTAGE**PARM(16))
         VSTMAX=PARM(3)
         VSTMAH = VSTMAX + 0.5

         DO J = 1,NBRNCH
            TBMAX(J) = IFIX(TB(J) + VSTMAX - VSTAGE + 0.5)
            TBMAH(J) = TBMAX(J) + 0.5
         END DO
         NBP1 = NBRNCH + 1
         DO J = NBP1,20
            TBMAX(J) = 5
            TBMAH(J) = 5.5
         END DO
! END IF
         INITR0 = 1
      END IF

! REDUCE RATE OF VEGETATIVE DEVELOPMENT AFTER FLOWERING.

      SLOW = 1.0

! THICKNESS OF NEW LEAVES

      IF (ITIME .LT. 2) THEN                                    ! *** Note 8.01
         TNL = 5.76E-4*WATWK/(1.0 + (0.021*WATWK)) + (7.E-6*(CO2-300))
         IF (VSTAGE.LT.5.0) THEN
            D11 = 7.5E-3 + (2.5E-3*VSTAGE)
            IF (TNL.GT.D11) TNL = D11
         END IF

! NUMBER OF DAYS SINCE POD GROWTH STARTED
! NUMBER OF DAYS SINCE SEED FILL STARTED

         IF (RSTAGE.GE.4.0) FILDAY = FILDAY + 1
      END IF

      IF (RSTAGE.GE.4.0.AND.FILDAY.EQ.0) FILDAY = 1
! ------------------------------------------------------------------
! CALCULATE POTENTIAL RATES OF CHANGE IN MAINSTEM LEAF AND
! PETIOLE DIMENSIONS AND DRY WEIGHTS.  CALCULATE MINIMUM
! RATE OF CARBON SUPPLY NEEDED TO SYNTHESIZE THESE.
! ------------------------------------------------------------------
! LIMIT THE VEGETATIVE GROWTH OF DETERMINATE PLANTS.

      IF (LIMIT.LT.1) THEN

! CALCULATE A CARBON/DRY WEIGHT CONVERSION FACTOR FOR EACH
! ORGAN DEPENDING ON ITS PROTEIN, (I.E., N), CONTENT.             *** Note 8.02
         D12 = 6.25*NRATIO/100
         D13 = 5.0*D12
         CONVL = (D13*0.65) + ((1.0-D13)*0.48)
         D14 = 2.0*D12
         CONVPE = (D14*0.65) + ((1.0-D14)*0.48)
         CONVST = CONVPE

! CALCULATE TEMPERATURE LIMITED RATES OF LEAF AREA AND PETIOLE EXTENSION.
! Note 8.03
         PDLAT = 10.0*(TAIR(ITIME) - 9.0)*PDV/3.0
         IF (PDLAT.LT.0.) PDLAT = 0.0
         PDPLT = 1.5*(TAIR(ITIME) - 9.0)*PDV/3.0
         IF (PDPLT.LT.0.) PDPLT = 0.0
! POTENTIAL RATE OF CHANGE IN AREA OF EACH GROWING MAINSTEM LEAF

         PDMLA = 0.
         PDMPW = 0.
         DO I = 1,3
            D15 = VSTAGE + 1.5 - I
            IF (D15 .LT. 0.) EXIT
            ITRIF = IFIX(D15)
!           IF (ITRIF.GT.1) PDLAM(I) = ((ITRIF*19.5) - 3.0)*PDV/3.0
            IF (ITRIF.GT.1) PDLAM(I) = PDLAT
            IF (ITRIF.LT.2) PDLAM(I) = 20.0*PDV/3.0

!  TEMPERATURE LIMITED RATE OF MAINSTEM LEAF AREA GROWTH.

            IF (PDLAT.LT.PDLAM(I)) PDLAM(I) = PDLAT
            PDLAM(I) = PDLAM(I)*SLOW * PARM(25)

!  RATE OF CHANGE IN MAINSTEM LEAF AREA.

            PDMLA = PDMLA + PDLAM(I)

!  POTENTIAL RATE OF CHANGE IN LENGTH OF EACH GROWING MAINSTEM PETIOLE.

            PDPLM(I) = ((ITRIF*4.0) + 1.0)*PDV/3.0
            IF (ITRIF .LE. 0) PDPLM(I) = 0.

!  TEMPERATURE LIMITED RATE OF MAINSTEM PETIOLE EXTENSION.

            IF (PDPLT .LT. PDPLM(I)) PDPLM(I) = PDPLT
            PDPLM(I) = PDPLM(I)*SLOW

!  RATE OF CHANGE IN MAINSTEM PETIOLE DRY WEIGHT AND MINIMUM
!  RATE OF CARBON SUPPLY NEEDED FOR ITS SYNTHESIS.

            PDPWM(I) = 0.0
            IF (ITRIF.NE.0) THEN
! Note 8.04
!VRR PARM(17) WAS 0.0024
               PDPWM(I) = PDPLM(I)*(PARM(17) +                         &
			              8.7E-5*(PDPLM(I)+2.0*MPETL(ITRIF)))
            END IF
            PDMPW = PDMPW + PDPWM(I)
         END DO
         PDMPC = PDMPW*CONVPE

! RATE OF CHANGE IN MAINSTEM LEAF VOLUME AND DRY WEIGHT, AND
! MINIMUM RATE OF CARBON SUPPLY NEEDED FOR ITS SYNTHESIS.

         PDMLV = PDMLA*TNL
         PDMLW = PDMLV/16.0
! Parameter 16.0 was 9.0 in original c-e. data
         PDMLC = PDMLW*CONVL
! ------------------------------------------------------------------
! CALCULATE POTENTIAL RATES OF CHANGE IN BRANCH LEAF AND PETIOLE
! DIMENSIONS AND DRY WEIGHT. CALCULATE MINIMUM RATE OF CARBON
! SUPPLY NEEDED TO SYNTHESIZE THESE.
! ------------------------------------------------------------------
! POTENTIAL RATE OF CHANGE IN AREA OF EACH GROWING BRANCH LEAF

         PDBLA = 0.0
         PDBPW = 0.0
         IF (NBRNCH .GE. 1) THEN
            DO J = 1,NBRNCH
               IF (TB(J).LT.TBMAH(J)) THEN
                  DO II = 1,3
                     I = IFIX(TB(J) + 1.5 - II)
                     IF (I .LE. 0) CYCLE
                     PDLAB(I,J) = ((I*43.3) - 13.0)*PDV/3.0
! TEMPERATURE LIMITED RATE OF BRANCH LEAF AREA GROWTH.

                     IF ((PDLAT).LT.PDLAB(I,J)) PDLAB(I,J)= PDLAT
                     PDLAB(I,J) = PDLAB(I,J)*SLOW * PARM(25)

! RATE OF CHANGE IN BRANCH LEAF AREA.

                     PDBLA = PDBLA + PDLAB(I,J)

! POTENTIAL RATE OF CHANGE IN LENGTH OF EACH GROWING BRANCH PETIOLE.

                     PDPLB(I,J) = I*12.*PDV/3.0

! TEMPERATURE LIMITED RATE OF BRANCH PETIOLE EXTENSION.
                     IF (PDPLT.LT.PDPLB(I,J)) PDPLB(I,J) = PDPLT
                     PDPLB(I,J) = PDPLB(I,J)*SLOW

! RATE OF CHANGE IN BRANCH PETIOLE DRY WEIGHT AND MINIMUM RATE
! OF CARBON SUPPLY NEEDED FOR ITS SYNTHESIS.
!VRR PARM(17) WAS 0.0024
                     PDPWB(I,J) = PDPLB(I,J)*(PARM(17) +               &
                          8.7E-5*(PDPLB(I,J) + 2.0*BPETL(I,J)))
                     PDBPW = PDBPW + PDPWB(I,J)
                  END DO
               END IF
            END DO
            IF (TB(NBRNCH).GT.TBMAH(NBRNCH)) LIMIT = 1
         END IF
         PDBPC = PDBPW*CONVPE

! RATE OF CHANGE IN BRANCH LEAF VOLUME AND DRY WEIGHT, AND
! MINIMUM RATE OF CARBON SUPPLY NEEDED FOR ITS SYNTHESIS.

         PDBLV = PDBLA*TNL
         PDBLW = PDBLV/21.0
!   Parameter 21.0 was 12.0 in original c-e. data
         PDBLC = PDBLW*CONVL
! ------------------------------------------------------------------
! CALCULATE POTENTIAL RATE OF CHANGE IN STEM AND BRANCH LENGTH and dry
! WEIGHT AND MINIMUM RATE OF CARBON SUPPLY NEEDED FOR THEIR SYNTHESIS.
! Note 8.05
! THE DATA FOR CALCULATION OF HTFAIR CAME FROM 1992 AND 1993 FIELD 
! STUDIES AT BARC-WEST, BELTSVILLE, MD. V. R. REDDY. HTFAIR IS ADJUSTMENT
! FOR HEIGHT DUE TO INCREASE IN PLANT DENSITY
! OLD CODE BEGIN
!    IF(ITIME.EQ.8)ALFINT = (RI/POPROW)/(LAREAT/(100/POPROW))
!    HTFATR = 1.0
!    IF(ALFINT.LE.4000)HTFATR=3.0
!    PDMH=(PARM(21)+(PARM(22)*(VSTAGE+PDV/2.0)**1.37))*PDV*SLOW*HTFATR
! OLD CODE END
         HTFATR = 1.87 - EXP(-0.0140 * POPAREA)
!YA      PDMH=(PARM(21)+(PARM(22)*(VSTAGE+PDV/2.0)**1.37))*PDV*SLOW*HTFATR
         IF (VSTAGE.GT.0.0) THEN
            F1=PARM(21)*PARM(22)*(VSTAGE)**(PARM(22)-1.)
         ELSE
            F1=0.0
         END IF
         IF (VSTAGE+PDV.GT.0.0) THEN
            F2=PARM(21)*PARM(22)*(VSTAGE+PDV)**(PARM(22)-1.)
         ELSE
            F2=0.0
         END IF
         PDMH=0.5*(F1+F2)*PDV
!YAEND
! PARM11 was 1.5 in original c-e. data
! PARM12 WAS 0.23
!YA      IF (VSTAGE.LT.(VSTMAH - 1.0)) THEN
         IF (RSTAGE .LE. PARM(16)) THEN
            PDMW = PDMH * (0.011 + 8.7E-5*PARM(24)*(PDMH+2.0*MSTEMH))
! PARM14 WAS 0.007
! Parameter 0.007 was 0.011 in original c-e. data
            PDMC = PDMW*CONVST
         END IF
         PDBW = 0.
         DO J = 1, NBRNCH
            IF (TB(J).LT.(TBMAH(J) - 1.0)) THEN
               PDBRW(J) = PDMH * (0.011 +                              &
			              8.7E-5*PARM(24)*(PDMH + 2.0*BRNCHH(J)) )
! Parameter 0.007 was 0.011 in original c-e. data
               PDBW = PDBW + PDBRW(J)
            END IF
         END DO
         PDBC = PDBW*CONVST

! MAXIMUM AND MINIMUM RATES OF CARBON USE BY VEGETATIVE SHOOT
! EXPANDING AT POTENTIAL RATE.                                    *** Note 8.06
         PCST = ((PDMPC + PDMLC + PDBPC + PDBLC)*1.6)                  &
              + ((PDMC + PDBC)*2.5)
         QCST = PDMPC + PDMLC + PDBPC + PDBLC + PDMC + PDBC
	  END IF
	  
! PCST Maximum rate of carbon use by shoot expanding at potential rate (g plant-1 hr-1).
! QCST Minimum rate of carbon use by shoot expanding at potential rate (g plant-1 hr-1).
! ------------------------------------------------------------------
! CALCULATE POTENTIAL RATE OF CHANGE IN DRY WEIGHT OF FLOWERS
! PODS AND SEEDS AND THE CARBON SUPPLY NEEDED FOR THEIR SYNTHESIS
! CALCULATE A CARBON/DRY WEIGHT CONVERSION FACTOR FOR FLOWERS AND PODS.

      D16 = 4.0*NRATIO*6.25/100
      CONVPO = (D16*0.58) + ((1.0 - D16)*0.46)

! IF RSTAGE IS BETWEEN 1 AND 3,
! CALCULATE POTENTIAL RATE OF CHANGE IN FLOWER DRY WEIGHT AND
! RATE OF CARBON SUPPLY NEEDED FOR ITS SYNTHESIS.

      IF (RSTAGE.GE.1.0) THEN
         IF (RSTAGE.LT.3.0) THEN
            PDFW = NFLWRS*0.001*PDR
! 0.001 = 0.002/2
            PDFC = PDFW*CONVPO

! IF RSTAGE IS BETWEEN 3 AND 6,
! CALCULATE POTENTIAL RATE OF CHANGE IN POD DRY WEIGHT AND
! RATE OF CARBON SUPPLY NEEDED FOR ITS SYNTHESIS.

         ELSE
            IF (RSTAGE .LE. 6.0) THEN
               PDPODW = (NPODS*0.06)*PDR * PARM(19)
               IF (RSTAGE.GE.4.0)                                      &
			       PDPODW = PDPODW*(6.0 - RSTAGE)*0.5*PARM(19)
               PDPODC = PDPODW*CONVPO
			END IF
! CALCULATE POTENTIAL RATE OF CHANGE IN SEED DRY WEIGHT AND
! RATE OF CARBON SUPPLY NEEDED FOR ITS SYNTHESIS.  SEEDS DO NOT
! GROW BEFORE R4 BUT POTENTIAL GROWTH RATE IS DETERMINED FROM
! R3 TO DECIDE HOW MANY PODS THE PLANT CAN CARRY.                 *** Note 8.07
            IF (RSTAGE.LT.4.0) THEN
               OIL = 0.22
               NSEEDS = NPODS*2.5
            ELSE
               OIL = (FILDAY + 4.)/100.
               IF (OIL.GT.0.22) OIL = 0.22
            END IF
            PROTIN = 0.4063
! 0.4063 = 6.5*6.25/100
            CONVSE = OIL*0.95 + PROTIN*0.58 + (1.0 - OIL - PROTIN)*0.46
! The following line is modified from the original code.
            IF (NSEEDS.GT.0.0) THEN
               FILL = 2.683 + (16.94*SDWMAX)                    ! *** Note 8.08
               IF (FILL.LT.(SDWMAX*5.82E-2)) FILL = SDWMAX*5.82E-2
! Note 8.09
               PDSW =NSEEDS*FILL/24./1000.                             &
                     *(1.0 - (0.025*ABS(TAIR(ITIME) - 24.0)))*PARM(20)
               IF (PDSW.LE.0.) PDSW = 0.
               PDSC = PDSW*CONVSE

! BETWEEN R3 AND R4 CALCULATE THE CARBON THAT WOULD BE NEEDED BY
! SEEDS IF THEY WERE GROWING.

               IF (RSTAGE.LT.4.0) THEN                          ! *** Note 8.10
                  EQSC = PDSC + (PDSW*0.13)
                  RSTPC = EQSC
                  PDSC = 0.0
                  PDSW = 0.0
                  NSEEDS = 0.0
               ELSE
                  RSTPC = 0.0
               END IF
            END IF
         END IF
      END IF
      RETURN
      END