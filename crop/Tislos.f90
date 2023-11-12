!***************************( TISLOS )**********************************
!** This module removes from the plant a tissue that is lost for any  **
!** reason. If the lower LEAVES ARE PARASITIC OR PLANT IS SHORT OF    **
!** NITROGEN OR A BRANCH IS GROWING IN THE LEAF AXIL, ABSCISE LEAVES  **
!** AND PETIOLES STARTING FROM BASE OF PLANT.                         **
!***********************************************************************
      SUBROUTINE TISLOS
      INCLUDE 'common.h'
      INCLUDE 'plant.h'
	  !DEC$ATTRIBUTES DLLEXPORT :: crop, /ShootR/,/shtR_public/, &
		/Weath/, /grid_public/, /nodal_public/, /elem_public/,   &
	    /bound_public/, /time_public/, /module_public/,          &
	     /DataFilenames/  
      IF ( INIT ) THEN
         NPA = 3
         LONDAY = 0
         NRATS = 0.0
         PDWT(1) = 0.0
		 ABSDW = 0.0  !WSUN put ABSDW to calculate the total lost weight
      END IF
       iday=int(time)
!YAP Number of days after which seeds that have stopped growing will be aborted
      MSEEDS=0
      DWBAL1 = LEAFWT + PETWT + ROOTWT + PODWT + SEEDWT
	  
! CALCULATE ACCUMULATED NITROGEN SUPPLY/DEMAND RATIO

      NRATS = NRATS + (NRATIO*PERIOD)

! ONCE A DAY CHECK FOR ORGANS THAT MUST BE ABSCISED               *** Note 13.01
      IF (ITIME .EQ. IPERD) THEN

! IF THE TEMPERATURE IS TOO LOW THE WHOLE PLANT DIES.

         IF (TMIN .LT. 0.0) THEN
            WRITE (10,"(' PLANTS ARE FROZEN TO DEATH')")
            MATURE = 1
         ELSE

! IS A BRANCH GROWING IN THE LEAF AXIL ?

            IF (VSTAGE.GE.5.) THEN
! Note 13.02
               DO J = 1,NBRNCH
                  DAYBR(J) = DAYBR(J) + 1
                  DAYLF = 21.0
                  IF (ILOW.LT.3.0) DAYLF = ILOW*7.0
                  IF (DAYBR(J).GT.DAYLF.AND.J.GE.ILOW) THEN
                     ABSCIS = 1
!                     WRITE (10, 112) IDAY, ILOW
                     GO TO 1111              !YA
                  END IF
               END DO

! IS THE PLANT SHORT OF NITROGEN ?

               AVNRAT = NRATS/24.0
               NRATS = 0.0
               IF (RSTAGE.LE.1.0) THEN
! NOTE 13.03
                  IF (AVNRAT.LT.0.99 .AND. (VSTAGE - 5.0).GT.ILOW)     &
                  LONDAY = LONDAY + 1
                  IF (LONDAY.GE.5) THEN
                     LONDAY = 0
!                     WRITE (10, 113) IDAY, ILOW
                     GOTO 1111               !YAP
                  END IF
               END IF

! ARE THE LOWER LEAVES PARASITIC ?                                *** NOTE 13.04
!YA
! IF (DROP.EQ.0.OR.RSTAGE.LT.5.0) GOTO 1112
               IF (RSTAGE.LE.6.0) GOTO 1112
               IF (DROP.EQ.0) GOTO 1112
!YAEND
!               WRITE (10, 114) IDAY,ILOW

! ABSCISE THE UNIFOLIOLATES

 1111          IF (ILOW.LT.1) THEN
! NOTE 13.05
                  ABSLPN = 0.02*ULEAFW + 0.0004 + UPETW*0.008
                  LEAFWT = LEAFWT - ULEAFW - 0.02
                  PETWT  = PETWT  - UPETW
                  ABSDW  = ABSDW  + ULEAFW + 0.02 + UPETW
                  LAREAM = LAREAM - ULAREA - 4.00
                  ULEAFW = 0.0
                  UPETW  = 0.0
                  ULAREA = 0.0
                  UPETL  = 0
                  ILOW = ILOW + 1
                  IF (DROP .EQ. 1 .OR. ABSCIS .EQ. 1) GOTO 1113
               END IF

!              IF (RSTAGE.LE.6.0) GOTO 1112          !YA
! ABSCISE TRIFOLIOLATES ON THE MAINSTEM

               IF (ILOW.LE.NTOPLF) THEN
                  I = ILOW
                  ABSLPN = ABSLPN + (MLEAFW(I)*2.0 + MPETW(I)*0.8)/100
                  LEAFWT = LEAFWT - MLEAFW(I)
                  PETWT = PETWT - MPETW(I)
                  ABSDW = ABSDW + MLEAFW(I) + MPETW(I)
                  LAREAM = LAREAM - MLAREA(I)
                  MLEAFW(I) = 0.0
                  MPETW(I) = 0.0
                  MLAREA(I) = 0.0
                  MPETL(I) = 0.0
                  ILOW = ILOW + 1
               END IF

! ABSCISE TRIFOLIOLATES ON THE BRANCHES

			   IF (NBRNCH .GT. 0) THEN
               DO J = 1, NBRNCH
! Note 13.06
                  IF (INITBR(J) .GE. ILOW ) GOTO 1113
                  IF (ILOWB(J).LE.(IFIX(TBMAH(J) - 0.5))) THEN
                     JJ = ILOWB(J)
                     ABSLPN = ABSLPN +									&
							  0.02*BLEAFW(JJ,J) + 0.008*BPETW(JJ,J)
                     LEAFWT = LEAFWT - BLEAFW(JJ,J)
                     PETWT  = PETWT - BPETW(JJ,J)
                     ABSDW  = ABSDW + BLEAFW(JJ,J) + BPETW(JJ,J)
                     LAREAB = LAREAB - BLAREA(JJ,J)
                     BLEAFW(JJ,J) = 0.0
                     BPETW(JJ,J) = 0.0
                     BLAREA(JJ,J) = 0.0
                     BPETL(JJ,J) = 0.0
                     ILOWB(J) = ILOWB(J) + 1
                  END IF
               END DO
               IF (LAREAT .LE. 1.0E-3 .AND. RSTAGE .GT. 7.9) GOTO 1114
!              IF (LAREAT .LE. 1.0E-3) GOTO 1114
			   END IF

! RECALCULATE TOTAL LEAF AREA AND RESET THE FLAGS CONTROLLING ABSCISSION.
!
 1113          LAREAT = LAREAM + LAREAB
               LAI = LAREAT * POPAREA * 1.0E-4

             
               
!DT Changed this to 1.0e-3 because it sometimes became 1.0E-7
               IF (LAREAT .LE. 1.0E-3 .AND. RSTAGE .GT. 7.9) GOTO 1114
!              IF (LAREAT .LE. 1.0E-3) GOTO 1114
               DROP = 0
               ABSCIS = 0
               GOTO 1112

! WHEN ALL LEAVES HAVE ABSCISED, GROWTH STOPS.

 1114          MATURE = 1
               R8DATE = CROP_DATE1
               LAREAM = 0.
               LAREAB = 0.
               LAREAT = 0.
               LAI    = 0.
               LEAFWT = 0.
               PETWT  = 0.
            END IF

 1112       CONTINUE
         END IF
      END IF
      RETURN
!112  FORMAT (1X,'DAY',I3,' LEAF',I3,' BRANCH PUSH')
!113  FORMAT (1X,'DAY',I3,' LEAF',I3,' MILD NITROGEN STRESS')
!114  FORMAT (1X,'DAY',I3,' LEAF',I3,' PARASITIC') 
      END
