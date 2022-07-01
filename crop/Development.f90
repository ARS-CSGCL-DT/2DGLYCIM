!****************************** (PHEN) *********************************
! This subroutine calculates the vegetative and reproductive stages
! of growth as defined by FEHR and CAVINESS, 1977. 
! There are no other subroutines in this module.
      SUBROUTINE PHEN
      INCLUDE 'common.h'
	  INCLUDE 'plant.h'   
	  !DEC$ATTRIBUTES DLLEXPORT :: CROP, /ShootR/, /shtR_public/, &
	    /Weath/, /grid_public/, /nodal_public/, /elem_public/,	  &
	    /bound_public/, /time_public/, /module_public/,			  &
	    /error_public/, /DataFilenames/ 
      COMMON /ACCU/ DEFICIT, DDAE, DAE, U0, DR2, U5, DDU5, DU6, DR8,&
                    IRFLAG1, IRFLAG4, UMAX
!      INTEGER CurYear, julday
      IF ( INIT ) THEN
         DEFICIT= 0.
         DDAE   = 0.
         DAE = 0.
         DU6 = 0.
         IRFLAG1 = 0
         IRFLAG4 = 0
! COMPUTING AUXILIARY VARIABLES
         U0=4E+04 * Parm(2) -255.0
!		 WSUN change BEGIN_DATE TO Emerge		 
!		 write(myStr, '("01/01/",i4)') CurYear(BEGIN_DATE)
!		 write(myStr, '("01/01/",i4)') CurYear(Emerge)
		 
!        WSUN change JDFRST to emergeDay
!		 tmpday = (JDFRST- julday(myStr))+1
!		 tmpday = (emergeDay - julday(myStr))+1		 
		 DR2=FLOAT(tmpday)*Parm(9)+Parm(10)-tmpday
         ! BY GR DR2=FLOAT(JDFRST)*Parm(9)+Parm(10)-JDFRST
         UMAX=U0+Parm(3)/Parm(2)
         VMAX=Parm(3)
         DDU5=Parm(12)	 
         !IF (JDFRST .GT. 177) DDU5=10.
		 IF (tmpday .GT. 177) DDU5=10.
	  END IF
!YA
!To find evaporation 
!wsun 2d soil already calculated the evaporation
! ------------------------------------------------------------------
! CALCULATE POTENTIAL WATER EVAPORATION RATE FROM THE SOIL SURFACE
! FOR DAYLIGHT PERIODS: CALCULATE NET UPWARD LONG-WAVE RADIATION

      !IF (WATTSM(ITIME).GT.0.0) THEN
! Note 5.04
         !RNLU = (1.6E-3*WATRAT*(100.0 - TAIR(ITIME)))*697.6

! CALCULATE THE ALBEDO OF THE EXPOSED SOIL CELLS                  *** Note 5.05
! THE TOTAL RADIATION FALLING ON THE SOIL IS ASSUMED TO BE
! CONCENTRATED ON THE EXPOSED SOIL CELLS AND AN EQUIVALENT
! TOTAL RADIATION IS CALCULATED.                                  *** Note 5.07
         !D11 = WATTSM(ITIME)

! CALCULATE NET RADIATION ON THE EXPOSED SOIL CELLS.              *** Note 5.06
!YA
         !RNS = D11 - RNLU
         !IF (RNS.LE.0.0) RNS = 0.0

! CALCULATE POTENTIAL EVAPORATION RATE FROM EXPOSED SOIL CELLS    *** Note 5.08
         !D12 = 1.0
! Note 5.09
         !ESO = (                                                       &
		 !          DEL(ITIME)/GAMMA*RNS*3600.0 /                       &
		 !          (2500.8-2.3668*TAIR(ITIME))                         &
         !        + VPD(ITIME)*109.375*(1.0+0.149*WIND(ITIME)*D12)      &
		 !	   ) / (DEL(ITIME)/GAMMA + 1.0)

! FOR NIGHT PERIODS:
! CALCULATE POTENTIAL EVAPORATION RATE FROM EXPOSED SOIL CELLS

      !ELSE
      !D12 = 1.0
      !ESO = (VPD(ITIME)*109.375*(1.0 + (0.149*WIND(ITIME)*D12)))       &
      !     /((DEL(ITIME)/GAMMA) + 1.0)
      !END IF
! POTENTIAL EVAPORATION RATE ADJUSTED BASED ON NOAA MAPS.

      ESO = ESO * PARM(1)
!YA
      SUMESO=SUMESO+ESO/1000.           
      SUMRAIN=SUMRAIN+RAIN*PERIOD/24.  !RAIN mm/day
!wsun in 2D soil, Precipitation and irrigation together in RAIN
!      IF ( .NOT. INIT ) THEN
!         DO JIRRI = 1, NIRRI
!            IF (IRRDAY(JIRRI).EQ.JDAY)                                 &
!              SUMRAIN=SUMRAIN+WAPR(JIRRI)*10.*PERIOD/24.
!         END DO
!     END IF
!YAEND
! CALCULATE PROGRESS BETWEEN VSTAGES. DEGREE-DAYS, DAYS AFTER EMERGENCE, DEFICIT
      TAIRL=TAIR(ITIME)
      DDAE=DDAE+TAIRL*PERIOD/24.
      DAE=DAE+PERIOD/24.
      DEFICIT=SUMESO-SUMRAIN
! CALCULATE PROGRESS BETWEEN VEGETATIVE STAGES
      VOLD=VSTAGE
      IF (DDAE.LE.U0) THEN
         VSTAGE=0.
      ELSE IF (DDAE.GE.UMAX) THEN
         VSTAGE=VMAX
      ELSE
         DV=Parm(2)*TAIRL*PERIOD/24.
!wsun change clay(1) to fracclay(1) as it transferred FROM 2DSOIL to crop
!		 Clay(1) = 0.2
!		 IF (CLAY(1).GT.0.4 .AND. VSTAGE.LT.2.0) DV=DV*Parm(4)
		 IF (FracClay(1).GT.0.4 .AND. VSTAGE.LT.2.0) DV=DV*Parm(4)
  
         VSTAGE=VSTAGE+DV
      END IF
      DV=VSTAGE-VOLD
      IF (VSTAGE.LT.VMAX) PDV=DV

      DO J = 1, NBRNCH
         TB(J)=TB(J)+DV                                         ! *** Note 6.03
      ENDDO

! CALCULATE PROGRESS IN REPRODUCTIVE STAGES

      ROLD=RSTAGE
      IF (RSTAGE .LT. 0.) THEN
!wsun calculate DLNGMAX in 2d soil. 
! PROGRESS TOWARDS R0 (FLORAL INDUCTION)

         IF (DayOfYear .LT. 173) THEN
! BEFORE SOLSTICE
            DR=(Parm(5)+Parm(6)*(DAYLNG-DLNGMAX))*PERIOD/24.
         ELSE
! AFTER SOLSTICE
            DR=(Parm(5)+Parm(7)*(DAYLNG-DLNGMAX))*PERIOD/24.
         END IF
         RSTAGE=RSTAGE+DR
         IF (RSTAGE.GT.0.0) THEN
! PARTIAL PROGRESS TOWARDS R1
            D11=PERIOD*(1.0+ROLD/DR)
            RSTAGE=D11*Parm(8)
         END IF
      ELSE IF (RSTAGE .LT. 2.0) THEN

! PROGRESS TOWARDS R2

         DR=Parm(8)*PERIOD/24.
         RSTAGE=AMIN1(RSTAGE+DR,2.0)
         IF (DAE .GT. DR2) THEN
!  PARTIAL PROGRESS TOWARDS R3
            D12=(DAE-DR2)/PERIOD
            RSTAGE=2.+(DAE-DR2)*Parm(11)*TAIRL
         END IF
!  AT FIND DATE OF R1 FOR SUMMARY
         IF (RSTAGE.GE.1.0 .AND. IRFLAG1.EQ.0) THEN
            R1DATE = CROP_DATE1
            IRFLAG1 = 1
			R1DAY = DayOfYear        
         END IF
      ELSE IF (RSTAGE .EQ. 2.0) THEN

! PLATEAU R2

         IF (DAE.GT.DR2) THEN
!  PARTIAL PROGRESS TOWARDS R3
            D12=(DAE-DR2)/PERIOD
            RSTAGE=2.+(DAE-DR2)*Parm(11)*TAIRL
         END IF
      ELSE IF (RSTAGE .LT. 5.0) THEN

! PROGRESS TOWARDS R5

         DR=Parm(11)*TAIRL*PERIOD/24.
         RSTAGE=RSTAGE+DR
         IF (RSTAGE .GT. 5.0) THEN
            U5=DDAE-((RSTAGE-5.0)/DR)*TAIRL*PERIOD/24.
            RSTAGE=5.0
            IF (DDAE .GT. U5+DDU5) THEN
! PARTIAL PROGRESS TOWARDS R6
               RSTAGE=5+(DDAE-U5-DDU5)*Parm(11)
            END IF
         END IF
! AT FIND R4 DATE
         IF (RSTAGE.GE.4.0.AND.IRFLAG4.EQ.0) THEN
            R4DATE = CROP_DATE1
            IRFLAG4 = 1
			R4DAY = DayOfYear  
         END IF
      ELSE IF (RSTAGE .EQ. 5.0) THEN
         R5DATE = CROP_DATE1         !For regional simulation
		 R5DAY = DayOfYear  
! PLATEAU R5

         IF (DDAE.GT.U5+DDU5) THEN
!  PARTIAL PROGRESS TO R6
            RSTAGE=5+(DDAE-U5-DDU5)*Parm(11)
         END IF
      ELSE IF(RSTAGE.LT.6.0) THEN

! PROGRESS TOWARDS R6

         DR=Parm(11)*TAIRL*PERIOD/24.
         RSTAGE=RSTAGE+DR
         IF (RSTAGE.GT.6.0) THEN
! PARTIAL PROGRESS ALONG PLATEAU R6
            U6=DDAE-((RSTAGE-6.0)/DR)*TAIRL*PERIOD/24.
            DDU6= (DDAE - U6) /                                        &
                  MIN ( Parm(13), MAX(0.001,Parm(13)-Parm(14)*DEFICIT))
            DU6=DU6+DDU6
            RSTAGE=6.0
            IF (DU6.GT.1.0) THEN
! PARTIAL PROGRESS TOWARDS R7
               RSTAGE=6.0+Parm(15)*TAIRL*(DU6-1.0)/DDU6
            END IF
         END IF
      ELSE IF (RSTAGE .EQ. 6.0) THEN

! PLATEAU R6

         DDU6=TAIRL*PERIOD/24./                                        &
              MIN (Parm(13), MAX(0.001,Parm(13)-Parm(14)*DEFICIT))
         DU6=DU6+DDU6
         IF (DU6 .GT. 1.0) THEN
! PARTIAL PROGRESS TOWARDS R7
            RSTAGE=6.0+Parm(15)*TAIRL*(DU6-1.0)/DDU6
         END IF
      ELSE IF (RSTAGE.LT.7.0) THEN

!  PROGRESS TOWARDS R7

         DR=Parm(15)*TAIRL*PERIOD/24.
         RSTAGE=RSTAGE+DR
         J8=NTOPLF-ILOW
         IF (J8 .LE. 0) THEN
            DR8=DR
         ELSE
            DR8=1.0/(NTOPLF-ILOW)/IPERD
         END IF
      ELSE

!  PROGRESS TOWARDS R8

         DROP=1
         RSTAGE=RSTAGE+DR8
         IF (RSTAGE .GE. 8.0) then
				MATURE=1
				R8DAY = DayOfYear  !ONLY FOR REGIONAL SIMULATION
		 end if
! THIS IS IT
      END IF
      PDR=RSTAGE-ROLD
      RETURN
      END