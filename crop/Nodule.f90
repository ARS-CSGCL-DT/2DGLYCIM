! ----------------------------( NODULE )----------------------------
! THIS ROUTINE CALCULATES THE NITROGEN FIXING ACTIVITY, GROWTH AND
! DEATH OF NODULES IN EACH SOIL CELL.  IT STARTS NITROGEN FIXATION
! BY NODULES 3 DAYS AFTER THEY ARE INITIATED AND IT KILLS THOSE THAT
! HAVE NOT FIXED NITROGEN FOR 10 DAYS.
! ------------------------------------------------------------------
      SUBROUTINE NODULE
      use common_block
      INCLUDE 'common.h'
	  
	 
      REAL      NFCF1, NFCF2, NODSIT, NCPM
      DIMENSION VALUE(200)

! CALCULATE INITIAL NODULE DISTRIBUTION.

!      IF ( INIT ) THEN
!         NLN = NLR
!         NKN = NKR
!         DO L = 1,NLN
!            DO K = 1,NKN
!               NFIXWT(L,K) = RTWT(L,K)/100.
!               NFIXNU(L,K) = NFIXWT(L,K)/CONVN/1.58E-4
!               NODNUM(L,K) = NFIXNU(L,K)
!            END DO
!         END DO
!      END IF
! CALCULATE INITIAL NODULE DISTRIBUTION.
! WSUN use 2dsoil cell based	  
		IF ( INIT ) THEN
            DO n=1,NumNP
               NFIXWT(n) = RTWT(n)/100.
               NFIXNU(n) = NFIXWT(n)/CONVN/1.58E-4
               NODNUM(n) = NFIXNU(n)
            END DO
		END IF

! AT THE START OF EACH DAY, AGE RECENTLY INITIATED NODULES BY ONE DAY.

!      IF ( ITIME .EQ. 1) THEN
!         DO L = 1,NLN
!            DO K = 1,NKN
!               IF ((NODNUM(L,K).GT.0.0).OR.(NODNUM(L,K).LT.0.0)) THEN
!                  NFIXNU(L,K) = NFIXNU(L,K) + NODNEW(3,L,K)
!                  NFIXWT(L,K) = NFIXWT(L,K) + NNEWWT(3,L,K)
!                  NODNEW(3,L,K) = NODNEW(2,L,K)
!                  NNEWWT(3,L,K) = NNEWWT(2,L,K)
!                  NODNEW(2,L,K) = NODNEW(1,L,K)
!                  NNEWWT(2,L,K) = NNEWWT(1,L,K)
!                  NODNEW(1,L,K) = 0.0
!                  NNEWWT(1,L,K) = 0.0
!               END IF
!            END DO
!         END DO
!      END IF
		
! WSUN use 2dsoil cell based
! AT THE START OF EACH DAY, AGE RECENTLY INITIATED NODULES BY ONE DAY.
      IF ( ITIME .EQ. 1) THEN
			DO n=1,NumNP
               IF ((NODNUM(n).GT.0.0).OR.(NODNUM(n).LT.0.0)) THEN
                  NFIXNU(n) = NFIXNU(n) + NODNEW(3,n)
                  NFIXWT(n) = NFIXWT(n) + NNEWWT(3,n)
                  NODNEW(3,n) = NODNEW(2,n)
                  NNEWWT(3,n) = NNEWWT(2,n)
                  NODNEW(2,n) = NODNEW(1,n)
                  NNEWWT(2,n) = NNEWWT(1,n)
                  NODNEW(1,n) = 0.0
                  NNEWWT(1,n) = 0.0
               END IF
            END DO
      END IF
		
! CALCULATE CARBON AVAILABLE TO NODULES IN A HALF SOIL SLAB
! DURING THE PRESENT PERIOD.                                     *** Note 11.06
      RNODC = (NODC*POPSLB*PERIOD) + NCPOOL
      NCPOOL = 0.0

! SINCE THERE IS NO CARBON FOR NODULES, ALL EXISTING NODULES ARE
! MARKED FOR DEATH IN 10 DAYS.

!      IF (RNODC.LE.0.0) THEN
!         IF (RSTAGE.GE.4.0) THEN
!            NLNP = 1
!            NKNP = 1
!		    DO L = 1,NLN
!			   DO K = 1,NKN
!                  IF ((NODNUM(L,K).GT.0.0).OR.(NODNUM(L,K).LT.0.0)) THEN
!	                 NODIE(L,K) = NODIE(L,K) + (1.0/IPERD)
!			         IF (NODIE(L,K).GE.10) THEN
!			            NODWT = NODWT - (NFIXWT(L,K)/POPSLB)
!			            NFIXWT(L,K) = 0.0
!			            NFIXNU(L,K) = 0.0
!			            NODNUM(L,K) = 0.
!				     ELSE
!			            IF (L.GT.NLNP) NLNP = L
!				        IF (K.GT.NKNP) NKNP = K
!			         END IF
!				  END IF
!			   END DO
!		    END DO
!            NLN = NLNP
!            NKN = NKNP
!         END IF
!      END IF
!      IF(RNODC.LE.0.0.AND.PRENOD.LE.0.0) RETURN

! WSUN use 2dsoil cell based
! SINCE THERE IS NO CARBON FOR NODULES, ALL EXISTING NODULES ARE
! MARKED FOR DEATH IN 10 DAYS.
	  IF (RNODC.LE.0.0) THEN
         IF (RSTAGE.GE.4.0) THEN
             m = 1
		       DO n=1,NumNP
                  IF ((NODNUM(n).GT.0.0).OR.(NODNUM(n).LT.0.0)) THEN
	                 NODIE(n) = NODIE(n) + (1.0/IPERD)
			         IF (NODIE(n).GE.10) THEN
			            NODWT = NODWT - (NFIXWT(n)/POPSLB)
			            NFIXWT(n) = 0.0
			            NFIXNU(n) = 0.0
			            NODNUM(n) = 0.0
				     ELSE
			            IF (n.GT.m) m = n
			         END IF
				  END IF
			   END DO
            n = m
         END IF
      END IF
      IF(RNODC.LE.0.0.AND.PRENOD.LE.0.0) RETURN
	  
	  
	  
! THE REDUCTION IN NODULE NITROGEN-FIXING ACTIVITY CAUSED BY
! TEMPERATURE AND OXYGEN CONCENTRATION IS CALCULATED FOR EACH
! SOIL CELL AND THE VALUES ARE SORTED INTO AN ARRAY TO
! INDICATE WHICH CELLS HAVE THE MOST FAVORABLE CONDITIONS

! CALCULATE PROPORTIONAL REDUCTION OF NODULE NITROGEN FIXATION
! BY SOIL TEMPERATURE.

!      NCR = 0
!      DO 80 L = 1,NLR
! Note 11.08
!        DO 90 K = 1,NKR
!			IF (TS(L,K).LE.9.0) NFCF1 = 0.
!	        IF ((TS(L,K).GT.9.0).AND.(TS(L,K).LT.26.0))                &
!               NFCF1 = (TS(L,K) - 9.0)*0.0588
!	       IF ((TS(L,K).GE.26.0).AND.(TS(L,K).LE.34.0)) NFCF1 = 1.0
!      IF ((TS(L,K).GT.34.0).AND.(TS(L,K).LT.44.0))					   &
!               NFCF1 = (44.0 - TS(L,K))*0.1
!	              IF (TS(L,K).GE.44.0) NFCF1 = 0.


! CALCULATE PROPORTIONAL REDUCTION OF NODULE NITROGEN FIXATION
! BY SOIL OXYGEN CONTENT.                                *** Note 11.09
!            NFCF2 = (8.0*CXT(L,K)) - (15.0*CXT(L,K)*CXT(L,K))

! CALCULATE PROPORTIONAL REDUCTION OF NODULE NITROGEN FIXATION
! FROM ALL PHYSICAL CAUSES.

!            NFCF(L,K) = MIN(NFCF1, NFCF2)

! MAKE A SINGLE DIMENSIONED ARRAY FROM THE TWO-DIMENSIONED

!            NCR = NCR + 1
!            NODWAR(NCR) = L*1000 + K
!            VALUE(NCR) = NFCF(L,K)
!90         CONTINUE
!80	CONTINUE

! WSUN use 2dsoil cell based
! THE REDUCTION IN NODULE NITROGEN-FIXING ACTIVITY CAUSED BY
! TEMPERATURE AND OXYGEN CONCENTRATION IS CALCULATED FOR EACH
! SOIL CELL AND THE VALUES ARE SORTED INTO AN ARRAY TO
! INDICATE WHICH CELLS HAVE THE MOST FAVORABLE CONDITIONS

! CALCULATE PROPORTIONAL REDUCTION OF NODULE NITROGEN FIXATION
! BY SOIL TEMPERATURE.
! WSun TS(n) Temperature of soil in cell. This need tranfer from 2dsoil
! WSun Change TS to Tmpr
	  NCR = 0
      DO 80 n=1,NumNP            ! Note 11.08
			IF (tmpr(n).LE.9.0) NFCF1 = 0.
	        IF ((tmpr(n).GT.9.0).AND.(tmpr(n).LT.26.0))                &
               NFCF1 = (tmpr(n) - 9.0)*0.0588
	       IF ((tmpr(n).GE.26.0).AND.(tmpr(n).LE.34.0)) NFCF1 = 1.0
      IF ((tmpr(n).GT.34.0).AND.(tmpr(n).LT.44.0))                     &
               NFCF1 = (44.0 - tmpr(n))*0.1
	              IF (tmpr(n).GE.44.0) NFCF1 = 0.


! CALCULATE PROPORTIONAL REDUCTION OF NODULE NITROGEN FIXATION
! BY SOIL OXYGEN CONTENT.                                *** Note 11.09
! wsun make CXT(n) to be constant TRANSFER IT FROM 2DSOIL
			CXT(n) = 0.21	
            NFCF2 = (8.0*CXT(n)) - (15.0*CXT(n)*CXT(n))

! CALCULATE PROPORTIONAL REDUCTION OF NODULE NITROGEN FIXATION
! FROM ALL PHYSICAL CAUSES.

            NFCF(n) = MIN(NFCF1, NFCF2)

! MAKE A SINGLE DIMENSIONED ARRAY FROM THE TWO-DIMENSIONED
!			 L = NODWAR(J)/1000
!			 K = NODWAR(J) - L*1000
		     NCR = NCR + 1
!			 NODWAR(NCR) = L*1000 + K
			 NODWAR(NCR) = iFavRoot(NCR)
!             VALUE(NCR) = NFCF(NCR)
!		     iFavRoot(NumNPD) which came from 2dsoil

			 
80	CONTINUE
! SORT SOIL CELLS IN ORDER HIGHEST TO LOWEST TO SHOW WHERE
! CONDITIONS ARE MOST FAVORABLE FOR ROOT GROWTH.
!        CALL SORT(NCR, iFavRoot, NFCF)
		 CALL SORT(NCR, NODWAR, NFCF)
! 	     CALL SORT(NumNP, iFavRoot, NFCF)
! ------------------------------------------------------------------
! NITROGEN FIXING ACTIVITY,GROWTH AND INITIATION OF NODULES ARE
! CALCULATED FOR EACH SOIL CELL
! CALCULATE MAXIMUM POTENTIAL CHANGE IN DRY WEIGHT OF EACH
! NODULE DURING THIS PERIOD.                                      *** Note 11.10
      PDNWM = 5.27E-5/24*PERIOD

! CALCULATE POTENTIAL CHANGE IN DRY WEIGHT OF EXISTING NODULES
! AS LIMITED BY TEMPERATURE AND OXYGEN CONCENTRATION IN EACH
! SOIL CELL.  START WITH SOIL CELL WHERE CONDITIONS ARE MOST
! FAVORABLE FOR NITROGEN FIXATION.

      PDNWS = 0.0
      CNFIXS = 0.0
      !NLNP = 1
      !NKNP = 1
	  m=1
! Note 11.11
      DO J = 1,NumNP
         !L = NODWAR(J)/1000
         !K = NODWAR(J) - L*1000
		   n = NODWAR(J)
		 ! n= iFavRoot(J)
         IF ((NODNUM(n).GT.0.0).OR.(NODNUM(n).LT.0.0)) THEN
            PDNW(n) = NODNUM(n)*NFCF(n)*PDNWM
            PDNWS = PDNWS + PDNW(n)

! CALCULATE CARBON CONSUMPTION DURING NITROGEN FIXATION IN EACH
! CELL UNTIL ALL CARBON IS USED.

            IF ((NFIXWT(n).GT.0.0).OR.(NFIXWT(n).LT.0.0)) THEN
               IF (RNODC.GT.0.0) THEN
! Note 11.12
                  CNFIX = NFIXWT(n)*NFCF(n)*2.147E-3*PERIOD*2.0
                  NODIE(n) = 0.
                  IF (RNODC.LE.CNFIX) CNFIX = RNODC
                  RNODC = RNODC - CNFIX
                  CNFIXS = CNFIXS + CNFIX
                  RESPS(n) = RESPS(n) + (CNFIX*3.67)
                  CYCLE
               END IF

! NODULES IN CELLS NOT RECEIVING CARBON FOR NITROGEN FIXATION
! ARE MARKED FOR DEATH IN 10 DAYS.

               IF (RSTAGE.LT.4.0) CYCLE
! Note 11.13
               NODIE(n) = NODIE(n) + (1.0/IPERD)
               IF (NODIE(n).LT.10) GO TO 1112
               NODWT = NODWT - (NFIXWT(n)/POPSLB)
               NFIXWT(n) = 0.0
               NFIXNU(n) = 0.0
               NODNUM(n) = 0.
               CYCLE
            END IF

 1112        IF (n.GT.m) m = n
!            IF (K.GT.NKNP) NKNP = K
         END IF
	END DO
!      NLN = NLNP
!      NKN = NKNP
	  n=m
      USEDCS = USEDCS + (CNFIXS/POPSLB)

! CALCULATE AMOUNT OF NITROGEN FIXED PER PLANT DURING THIS PERIOD.

      NODN = CNFIXS/2.0/POPSLB

! IF ANY CARBON IS AVAILABLE FILL THE NODULES CARBON STORES.

      IF (RNODC.GT.0.0) THEN
! Note 11.14
         NCPM = NODWT*0.03*POPSLB
         IF (RNODC.GT.NCPM) THEN
            NCPOOL = NCPM
            RNODC = RNODC - NCPM
         ELSE
            NCPOOL = RNODC
            RNODC = 0.0
         END IF
      END IF

! CALCULATE CARBON NEEDED TO GROW EXISTING NODULES AT POTENTIAL
! RATE IN A HALF SOIL SLAB.                                       *** Note 11.15
      PDNCS = PDNWS*CONVN

! ADD IN THE CARBON PARTITIONED TO NODULE GROWTH

      RNODC = RNODC + (PRENOD*POPSLB*PERIOD)

! CALCULATE HOW MUCH CARBON IS AVAILABLE TO INITIATE NEW NODULES
! IN A HALF SOIL SLAB. IF THERE IS NONE, JUMP TO APPROPRIATE
! PART OF ROUTINE.                                                *** Note 11.16
      CNINT = RNODC - PDNCS

! SINCE THERE IS CARBON AVAILABLE FOR INITIATING NEW NODULES,
! GROW ALL EXISTING NODULES AT POTENTIAL RATE AND INITIATE NEW
! ONES UNTIL ALL CARBON IS USED. START WITH SOIL CELL WHERE
! CONDITIONS ARE MOST FAVORABLE FOR NITROGEN FIXATION.
! RESPS CO2 respired in soil cell L,K in this period (g).
      IF (CNINT.GT.0.0) THEN
         DO n=1,NumNP
            !L = NODWAR(J)/1000
            !K = NODWAR(J) - L*1000
			! n1=iFavRoot(n)
			  n1 = NODWAR(n)
            IF ((NFCF(n1).GT.0.).OR.(NFCF(n1).LT.0.)) THEN
               D11 = NFCF(n1)*PDNWM
               IF ((NODNUM(n1).GT.0).OR.(NODNUM(n1).LT.0)) THEN
                  D12 = D11*NFIXNU(n1)
                  NFIXWT(n1) = NFIXWT(n1) + D12
                  NODWT = NODWT + (D12/POPSLB)
! Note 11.17
                  RESPS(n1) = RESPS(n1) + (D12*0.212)
                  DO I = 1,10
                     II = 11 - I
                     D13 = D11*NODNEW(II,n1)
                     NNEWWT(II,n1) = NNEWWT(II,n1) + D13
                     NODWT = NODWT + (D13/POPSLB)
                     RESPS(n1) = RESPS(n1) + (D13*0.212)
                  END DO
               END IF
               IF (CNINT.LE.0.0) CYCLE
! Note 11.18
! YRL Length of young root in soil cell  (cm).
! CDT D14 is only given a value in the else loop.  thus it will be undefined for the RESPS calc if the 
! else loop is not entered. In the default loop, D12 is used like D14. Thus it seems logical that D12 can
!  be used for respiration calcs. I moved the resp calculation into the if statements using appropriate
! D12 or D14
               NODSIT = YRL(n1)
               PNONEW = CNINT/CONVN/PDNWM/NFCF(n1)
               IF (PNONEW.LT.NODSIT) THEN
                  NODNEW(1,n1) = NODNEW(1,n1) + PNONEW
                  NODNUM(n1) = NODNUM(n1) + PNONEW
                  D12 = CNINT/CONVN
                  NNEWWT(1,n1) = NNEWWT(1,n1) + D12
                  NODWT = NODWT + (D12/POPSLB)
                  CNINT = 0.0
                  RESPS(n1) = RESPS(n1) + (D12*0.212)
			   ELSE
                  NODNEW(1,n1) = NODNEW(1,n1) + NODSIT
                  NODNUM(n1) = NODNUM(n1) + NODSIT
                  D14 = NODSIT*D11
                  NNEWWT(1,n1) = NNEWWT(1,n1) + D14
                  CNINT = CNINT - D14*CONVN
                  NODWT = NODWT + (D14/POPSLB)
                  RESPS(n1) = RESPS(n1) + (D14*0.212)
               END IF

!			   IF (L.GT.NLN) NLN = L
            END IF
         END DO
         SCPOOL = SCPOOL + (CNINT/POPSLB)
         USEDCS = USEDCS + ((RNODC - CNINT)/POPSLB)
      ELSE
         USEDCS = USEDCS + (RNODC/POPSLB)

! SINCE THERE IS NO CARBON AVAILABLE FOR INITIATING NEW NODULES,
! GROW EXISTING NODULES UNTIL ALL CARBON USED. START WITH THE
! SOIL CELL WHERE CONDITIONS ARE MOST FAVORABLE FOR NITROGEN FIXATION.

         DO n=1,NumNP
            !L = NODWAR(J)/1000
            !K = NODWAR(J) - L*1000
			! n1=iFavRoot(n)
			  n1 = NODWAR(n)
            IF ((NODNUM(n1).GT.0.).OR.(NODNUM(n1).LT.0.)) THEN
               IF (RNODC.LE.0.0) RETURN
               PDNFW = NFIXNU(n1)*NFCF(n1)*PDNWM
               D15 = PDNFW*CONVN
               IF (RNODC.LE.D15) THEN
                  D16 = RNODC/CONVN
                  NFIXWT(n1) = NFIXWT(n1) + D16
                  NODWT = NODWT + (D16/POPSLB)
                  RESPS(n1) = RESPS(n1) + (D16*0.212)
                  RNODC = 0.0
                  RETURN
               END IF
               NFIXWT(n1) = NFIXWT(n1) + PDNFW
               NODWT = NODWT + (PDNFW/POPSLB)
               RESPS(n1) = RESPS(n1) + (PDNFW*0.212)
               RNODC = RNODC - D15
               DO I = 1,5
                  II = 6 - I
                  PNNWWT = NODNEW(II,n1)*NFCF(n1)*PDNWM
                  D20 = PNNWWT*CONVN
                  IF (RNODC.LE.D20) THEN
                     D21 = RNODC/CONVN
                     NNEWWT(II,n1) = NNEWWT(II,n1) + D21
                     NODWT = NODWT + (D21/POPSLB)
                     RESPS(n1) = RESPS(n1) + (D21*0.212)
                     RNODC = 0.0
                     RETURN
                  END IF
                  NNEWWT(II,n1) = NNEWWT(II,n1) + PNNWWT
                  NODWT = NODWT + (PNNWWT/POPSLB)
                  RESPS(n1) = RESPS(n1) + (PNNWWT*0.212)
                  RNODC = RNODC - D20
               END DO
            END IF
         END DO
      END IF
      RETURN
      END