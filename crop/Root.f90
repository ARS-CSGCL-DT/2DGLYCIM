 SUBROUTINE ROOT
      INCLUDE 'common.h'
	  INCLUDE 'plant.h'
	  !DEC$ATTRIBUTES DLLEXPORT :: crop, /ShootR/,/shtR_public/,&
		/Weath/, /grid_public/, /nodal_public/, /elem_public/,  &
	    /bound_public/, /time_public/, /module_public/,         &
	    /error_public/, /DataFilenames/  
	

	  
! CALCULATE A CARBON TO DRY WEIGHT CONVERSION FACTOR FOR THE ROOTS.

      D15 = NRATIO*0.1563
! 0.1563 = 2.5*6.25/100
! Note 10.20
      CONVR  = D15*0.650 + (1.0 - D15)*0.480
      CONRES = D15*0.416 + (1.0 - D15)*0.124


	  
! CALCULATE RATE AT WHICH CARBON IS SUPPLIED TO THE GROWING
! ROOTS IN A HALF SOIL SLAB WHEN:
! (1) NO CARBON IS TRANSLOCATED FROM THE SHOOT
! (2) SHOOTS DO NOT GROW AND ROOTS GET ALL TRANSLOCATED CARBON

      PARTRT = 1.0 - (VSTAGE*0.9/5.0)
      IF (PARTRT.LT.0.1) PARTRT = 0.1
      RCPAV = AMAX1(RCPOOL*(1.0 - (0.75**PERIOD))/PERIOD,0.)
!Rate at which carbon becomes available from RCPOOL for root growth (g plant-1 hr-1).
! Note 10.18
       PCRL = (RCPAV + PARTRT*VEGSRC)*POPSLB*24
	   
!PCRL Minimum carbon needed to grow roots g C day-1
! Note 10.19
!YAP
! AMULT=PARTRT
! IF (RSTAGE.GT.1.0) AMULT=PARTRT
! PCRQ = (RCPAV + AMULT*VEGSRC)*POPSLB 
         PCRQ = (RCPAV + 0.95*VEGSRC)*POPSLB*24  !keep this one right now 12/11/2020
!		 PCRQ = (FIXC)*POPSLB*24   
!PCRQ Maximum carbon provided to grow roots g C day-1
end