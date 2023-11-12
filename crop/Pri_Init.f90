 SUBROUTINE PRI_INIT
	  Include 'Common.h'
	  Include 'plant.h'
	  !DEC$ATTRIBUTES DLLEXPORT :: crop, /ShootR/,/shtR_public/, &
		/Weath/, /grid_public/, /nodal_public/, /elem_public/,   &
	    /bound_public/, /time_public/, /module_public/,          &
	     /DataFilenames/
!      Pi     = 3.1415926535897932D0
! 	  TwoPi  = 6.2831853071795864D0
!	  HlfPi  = 1.5707963267948966D0
!	  DEGRAD = Pi/180.0D0              ! 0.0174533...
!	  cm2in  = 0.393700787D0		   ! convert ratio for cm to in

      DO I=1,3
         PDLAM(I) = 0.0
         PDPLM(I) = 0.0
         PDPWM(I) = 0.0
	  END DO

      !DO J=1,5
         !ABPOD(J) = 0.0  !wsun Not used in GLYCIM
         !ABPWT(J) = 0.0  !wsun Not used in GLYCIM
         !ABSEED(J)= 0.0  !wsun Not used in GLYCIM
!wsun Below is original GLYCIM soil part	  
         !AIRDR(J) = 0.0
         !BD(J)    = 0.0
         !BETA(J)  = 0.0
         !CLAY(J)  = 0.0
         !DIFF0(J) = 0.0
         !ETA(J)   = 0.0
         !FCI(J)   = 0.0
         !PSISAT(J)= 0.0
         !SAND(J)  = 0.0
         !SCOND(J) = 0.0
         !SDEPTH(J)= 0.0
         !SILT(J)  = 0.0
         !SSIND(J) = 0.0
         !THETAS(J)= 0.0
         !THETA0(J)= 0.0
      !END DO

      !DO I=1,7
         !CLIMAT(I)=0.0
      !END DO

      !DO I=1,10
         !SURFLU(I)=0.0                 !YAP
         !Cum_FL_Inf(i)=0.0             !DT
         !EWU(I)=0.0
         !TNAIRD(I)=0.0
      !END DO

      DO I=1,20
         BRNCHH(I)=0.0
         BRNCHW(I)=0.0
         DAYBR(I)=0.0
		!FC(I)=0.0
         ILOWB(I)=1
         INITBR(I)=0
		!NBTPLF(I)=0 !wsun Not used in GLYCIM
		!NHORIZ(I)=0
         PDBRW(I)=0.0
         TB(I)=-1
         TBMAH(I)=30.0
      END DO

      DO I = 1, 24
! Iintial the LIMITF to blank (space)
         LIMITF(I)= ' '
! wsun calulated in 2dsoil		
		!EO(I)    = 0.0		
        !PAR(I)   = 0.0     
        !PARINT(I)= 0.0     
        !PSIL(I)  = 0.0     
        !RADINT(I)= 0.0
        !RHS(I)   = 0.0
        !TAIR(I)  = 0.0
        !WATTSM(I)= 0.0
        !WIND(I)  = 0.0
      END DO

      DO I = 1, 35
         MLAREA(I) = 0.0
         MLEAFW(I) = 0.0
         MPETL(I)  = 0.0
         MPETW(I)  = 0.0
      END DO

      !DO I = 1, 60
         !POD(I)  = 0.0
      !END DO
      DO I = 1, 4000
         NODWAR(I)=0.0
         !PCRPAR(I)=0.0
      END DO

      DO I = 1,500
         PDWT(I)=0.0
         PPDWT(I)=0.0
      END DO
      DO I = 1,1000
         PSDWT(I)=0.0
         SDWT(I)=0.0
      END DO

      DO n=1,4000
         !DO J=1, 20
            !!FCAP(J,K) = 0.0            !YAP
            !!ADRL(J,K) = 0.0
!           ADWR(J,K) = 0.0
            !!AWR(J,K)  = 0.0
            !!AWUP(J,K) = 0.0
            !!BDC(J,K)  = 0.0
            !!COND(J,K) = 0.0
            CXT(n)  = 0.0
            NFCF(n) = 0.0
            NFIXNU(n) = 0.0
            NFIXWT(n) = 0.0
            NODIE(n)  = 0.0
            NODNUM(n) = 0.0
            PDNW(n)   = 0.0
           !PDWR(n)   = 0.0
		   !PERTWT(n) = 0.0
		   !PPDRL(n)  = 0.0
		   !PSIRD(J,K)  = -0.1
		   !PSIS(J,K)   = -0.3
           RESPS(n)  = 0.0
		   !RGCF(J,K)   = 0.0
           !RTWT(J,K)   = 0.0
           !RUTDEN(J,K) = 0.0
		   !RVR(J,K) = 0.0
           !TOPORS(J,K) = 0.0
           !TPRD(J,K) = 2.0
            TS(n)     = 0.0
            !UPTH2O(J,K) = 0.0
            !VCLAY(J,K)  = 0.0
            !VH2OC(J,K)  = 0.0
            !VMUCK(J,K)  = 0.0
            !VNC(J,K)    = 0.0
            !VNH4C(J,K)  = 0.0
            !VNO3C(J,K)  = 0.0
            !VSIND(J,K)  = 0.0
            YRL(n)    = 0.0
         END DO
      !END DO

      DO J=1,20
         DO I=1,20
            BLAREA(I,J)=0.0
            BLEAFW(I,J)=0.0
            BPETL(I,J)=0.0
            BPETW(I,J)=0.0
            PDLAB(I,J)=0.0
            PDPLB(I,J)=0.0
            PDPWB(I,J)=0.0
         END DO
      END DO

      DO n=1,4000
            DO I=1,10
               NNEWWT(I,n)=0.0
               NODNEW(I,n)=0.0
            END DO
      END DO

      ABSCIS = 0
      ABSLPN = 0.0
      ABSPON = 0.0
     !AVP = 0.0
     !AWUPS = 0.0
     !AWUPSS = 0.0
     !CEC = .55    !wsun READ IN 2DSOIL INITIAL FILE
     !CLOUD = 0.0  !wsun 2DSOIL 
      CNO3 = 0.0
     !CONRES = 0.0 !wsun CALCULATE A CARBON TO DRY WEIGHT CONVERSION FACTOR FOR THE ROOTS.
      CONVN = 0.53 !wsun Constant value
      CONVPO = 0.0
      CONVSE = 0.0
	 !Cum_Flood=0.0                    !DT
      DEADRN = 0.0
     !DIFIMG = 0.0 !wsun not used in GLYCIM crop
     !DRAIN = 0.0
      DROP = 0
      DRPDAY = 0
     !EORS = 0.0
     !EORSCS = 0.0
     !EWUS = 0.0
      FILDAY = 0
      FIXCS = 0.0
      FLWRWT = 0.0
      FSEEDS = 0.0
     !GAMMA = 0.0
      HEIGHT = 5.5
     !HVALUES = 2
     !IDAY = 0
      ILOW = 0
      ILOWP = -1
      IPERD = 24
     !KEXPOS = 0
      LAREAB = 0.0
      LAREAT = 7.33
     !LCAI = 0.0
      LEAFWT = 0.03
     !LKH = 0
     !LKRP1 = 0
      LLUE = 0.012                 !wsun gasexchange initial
      LTC = 0.17				   !wsun gasexchange initial
      MATURE = 0
      MINPCD = 999.0
      MINSCD = 999.0
      MSTEMH = 0.5
      MSTEMW = 0.0
      NBRNCH = 0
      !NCD = 0
      NCPOOL = 0.0
      NFLWRS = 0
	 !NKH = 0
     !NKN = 1
     !NKR = 2
     !NLN = 1
     !NLR = 6  !wsun not used 
      NODES = 0
      NODN = 0.0
      NODWT = 0.0
      NRATIO = 1.0
      NTOPLF = 0
     !OSMFAC  = 0.0   !wsun defined in 2dsoil
     !OSMREG = -0.25
     ! PCRS = 0.0      !wsun need to transferred to 2dsoil
      PCST = 0.0
      PDFW = 0.0
      PDMH = 0.0
      PDMW = 0.0
      PDNODC = 0.0
      PDPODC = 0.0 
      PDPODW = 0.0
      PDR = 0.0
      PDSC = 0.0
      PDSW = 0.0
      PDV = 0.0
     !PERIOD = 0.0
      PETWT = 0.001
      PHOTO = 0.0
     !PILD = -9.0
     !PILOSM = -9.0
     ! PLANTN = 0.0065
      PNLLFS = 0.0
      PODCD = 0.0
      NPODS = 0.0
      PODWT = 0.0
      PRENOD = 0.0
	 !PSILD = -2.0
	 !PSISM = 0.0
      SDWT(1) = 0.0
      QCST = 0.0
     !RAIN = 0.0
      RCPOOL = 0.0
     ! ROOTWT = 0.035
      ROUGH = 0.0
     !RRRM = 4000.0   ! wsun was used in 2d soil
     !RRRY = 750.0
      RSTAGE = -1.0
      RSTPC = 0.0
      !RTEMP = 0.0
      !RTMINW = 0.0
      !RTWL = 9.0E-5 !wsun already defined in 2dsoil
      !RVRL = 17.50  !wsun not used
      !SCF = 0.01
      SCPOOL = 0.0124     !wsun gasexchange initial
      SCRAT0 = 0.0
      SEEDCD = 0.0
      NSEEDS = 0.0
      SEEDWT = 0.0
      !SGT = 0.0
      !SGTLI = 0.0
      SGTR = 0.0
      SINKC = 0.0 !wsun Change SINK TO SINKC
      SPCSP = 0.0
      SPCSS = 0.0
	 !SPDPC = 0.0
	 !SPDSC = 0.0
      STEMWT = 0.002
      SUPNO3 = 0.0
      SUMESO=0.0
	  SUMRAIN=0.0
      SVEGC = 0.0
      TAUN = 0.0
	 !TCADS = 0.0
      TNL = 0.0
	 !TNYT = 0.0
	 !TPL = 0.0
	 !TPLD = 5.0
	 !TRAIN = 0.0
      TRANSC = 0.0
      ULAREA = 3.33
      ULEAFW = 0.01
      UPETL = 0.17
      UPETW = 0.001
	 !UPH2OS = 0.0
      USEDCS = 0.0
      VEGSRC = 0.0
	 !VLOST = 0.0
      VSTAGE = 0.0
      VSTMAH = 30.0
      WATWK = 0.0
      XTRACS = 0.0
      YLD=0.0                          !YAP
     ! DO J = 1,5
        ! DO I = 1,5
            !TDSX(I,J) = 0.0
         !END DO
      !END DO

RETURN
END
