!*****************************( PARTIT )********************************
!** This module calculates an initial partitioningof carbon to various**
!** organs based upon the priorities that change with the stage of    **
!** growth and assuming that water and nutrients are plentiful.       **
!***********************************************************************
      SUBROUTINE PARTIT
      INCLUDE 'common.h'
	  INCLUDE 'plant.h'
	  !DEC$ATTRIBUTES DLLEXPORT :: crop, /ShootR/,/shtR_public/,&
		/Weath/, /grid_public/, /nodal_public/, /elem_public/,  &
	    /bound_public/, /time_public/, /module_public/,         &
	    /error_public/, /DataFilenames/    
      IF ( INIT ) SCRTS1 = 0.0
      NODC = 0.0
      FRUITC = 0.0
      VEGSRC = 0.0
      IF (RSTAGE.LT.3.0) THEN

! PRIOR TO R1 THE PRIORITY FOR CARBON USE IS: NODULES > VEGETATIVE PARTS

         IF (RSTAGE .LT. 1.0) THEN
            NODC = PDNODC
            IF (TRANSC.LT.NODC) THEN
               NODC = TRANSC
            ELSE
               VEGSRC = TRANSC-NODC
            END IF

! FROM R1 TO R3 THE PRIORITY FOR CARBON USE IS: NODULE GROWTH IN ANTICIPATION
! OF FRUITING > NORMAL NODULE ACTIVITY > VEGETATIVE PARTS

         ELSE
            PRENOD = TRANSC/6.0
            D11 = TRANSC - PRENOD
            NODC = PDNODC
            IF (D11.LT.NODC) THEN
               NODC = D11
            ELSE
               VEGSRC = D11 - NODC
            END IF
         END IF

! CALCULATE A SHOOT CARBON SUPPLY/DEMAND RATIO FOR THE DAY
! AND SET A FLAG (SINKC) TO ADJUST PHOTOSYNTHETIC RATE

         IF (VEGSRC.LT.QCST) THEN
            SCRATO = 0.0
         ELSE
            SCRATO = 1.0
            IF (VEGSRC.LT.PCST) THEN
               D12 = PCST - QCST
               IF (D12.GT.0.0) SCRATO = (VEGSRC - QCST)/D12
            END IF
         END IF
         SCRTS1 = SCRTS1 + (SCRATO*PERIOD)
         IF (ITIME.EQ.IPERD) THEN
            D13 = SCRTS1/24.0
            IF (D13.GT.0.9) THEN
               SINKC = 1.0 !wsun Change SINK TO SINKC
            ELSE
               SINKC = -1.0
            END IF
            SCRTS1 = 0.0
         END IF

! FROM R3 TO R4 THE PRIORITY FOR CARBON USE IS:
! FRUITS = STORAGE OF CARBON EQUIVALENT TO THAT REQUIRED FOR
! POTENTIAL SEEDS = NODULES > VEGETATIVE PARTS

      ELSE IF (RSTAGE .LT. 4.0) THEN
         SPCSP = SPCSP+((PDPODC+(PDPODW*0.08)+EQSC)*PERIOD)
         PRENOD = 0.0
!VRR
! IF(SGT.LT.1.0) THEN
!   VEGSRC = TRANSC
!   FRUITC = 0.0
!   NODC = 0.0
!   RSTPC = 0.0
! ELSE
         FRUITC = PDPODC + (PODCD/PERIOD)
         NODC = PDNODC
         D14 = FRUITC + RSTPC + NODC
         IF (TRANSC.LT.D14) THEN
            IF (ITIME.EQ.IPERD) SINKC = -1.0
            FRUITC = TRANSC*FRUITC/D14
            NODC = TRANSC*NODC/D14
            RSTPC = TRANSC*RSTPC/D14
         ELSE
            IF (ITIME.EQ.IPERD) SINKC = 1.0
            VEGSRC = TRANSC - FRUITC - RSTPC - NODC
            SVEGC = SVEGC + (VEGSRC * PERIOD)
         END IF
! FROM R4 ONWARDS THE PRIORITY FOR CARBON USE IS:
! FRUITS = NODULES > VEGETATIVE ORGANS.
! PHOTOSYNTHETIC RATE IS ADJUSTED (WITHIN LIMITS) TO JUST MEET
! THE CARBON DEMAND OF THE FRUITS PLUS NODULES

      ELSE
         SPCSS = SPCSS + ((PDSC+(PDSW*0.13)*PERIOD))
         FRUITC = PDPODC + (PODCD/PERIOD) + PDSC +(SEEDCD/PERIOD)
         NODC = PDNODC
         IF (TRANSC.LT.(FRUITC + NODC)) THEN
            IF (ITIME.EQ.IPERD) SINKC = -1.0
            D15 = FRUITC + NODC
            FRUITC = TRANSC*FRUITC/D15
            NODC = TRANSC*NODC/D15
         ELSE
            IF (ITIME.EQ.IPERD) SINKC = 1.0
            VEGSRC = TRANSC - FRUITC - NODC
            SVEGC = SVEGC + (VEGSRC*PERIOD)
         END IF
      END IF
      
      END