
!***********************************************************************
! This subroutine outputs to a graphics file.
! dt added leafwt and Stem weight 12/6/2006
!-----------------------------------------------------------------------
 SUBROUTINE Initialize
	 INCLUDE 'common.h'
	 Include 'plant.h'
	 Character InString1*120, location*120, location1*12,       &
	 nitrogencrp*156, nitrogencrp1*156, location2*15,		    &
	 plantstresscrp*156, plantstresscrp1*156
	 INTEGER index, declen, len, last, trim
	 integer*4 i,ls1,ls2,ps1,ps2
    
	 !DEC$ATTRIBUTES DLLEXPORT :: crop, /ShootR/, /shtR_public/,&
	    /Weath/, /grid_public/, /nodal_public/, /elem_public/,	&
	    /bound_public/, /time_public/, /module_public/,         &
	    /error_public/, /DataFilenames/  
	
   Write(*,*)'************************ GLYCIM ***********************'  
   Write(*,*)'*                    Version 08.4.30                  *'  
   Write(*,*)'*                                                     *'  
   Write(*,*)'*  A DYNAMIC SIMULATOR FOR SOYBEAN CROPS ORIGINALLY   *'  
   Write(*,*)'*  CREATED BY BASIL ACOCK AND SUBSEQUENTLY IMPROVED   *'  
   Write(*,*)'*  BY YAKOV PACHEPSKY, V.R. REDDY, DENNIS TIMLIN, AND *'  
   Write(*,*)'*  ANTHONY TRENT.                                     *'  
   Write(*,*)'*  USDA-ARS,REMOTE SENSING AND MODELING LABORATORY    *'  
   Write(*,*)'*  BELTSVILLE, MD  20705.   TEL:(301)504-5872         *'   
   Write(*,*)'*******************************************************'
   
   !Read variety file 

   Open(40,file = VarietyFile, status='old', ERR=13)

      im=220
      il=0
      Read(40,*, ERR=13)
      im=im+1
      il=il+1
      Read(40,*, ERR=13)
      im=im+1
      il=il+1
      Read(40,*, ERR=13)
      im=im+1
      il=il+1
      Read(40,*, ERR=13) MG, SEEDLB, FILL,						  &
	  PARM(2), PARM(3), PARM(4), PARM(5),                         &
	  PARM(6), PARM(7), PARM(8), PARM(9), PARM(10),               &
	  PARM(11), PARM(12), PARM(13), PARM(14), PARM(15),           &
	  PARM(16), PARM(17), PARM(18), PARM(19), PARM(20),           &
	  PARM(21), PARM(22), PARM(23), PARM(24), PARM(25)
      close(40)
	  
	  Open(85,file=PlantGraphics)  
	  Open(86,file=LeafGraphics)
	  
	
      declen = len(PlantGraphics)
	  last = index(PlantGraphics, '\', back=.true.)
	  location=PlantGraphics(1:last)
	  location=trim(location)
	  location1='nitrogen.crp'
	  location2='plantstress.crp'
	  nitrogencrp1=location//location1
	  plantstresscrp1=location//location2
	  
      ls1 = len_trim(nitrogencrp1)
      ls2 = 0
      do i = 1,ls1
        if(nitrogencrp1(i:i).ne.' ') then
           ls2 = ls2 + 1
           nitrogencrp(ls2:ls2) = nitrogencrp1(i:i)
        endif
	  enddo
	  
	  ps1 = len_trim(plantstresscrp1)
      ps2 = 0
      do i = 1,ps1
        if(plantstresscrp1(i:i).ne.' ') then
           ps2 = ps2 + 1
           plantstresscrp(ps2:ps2) = plantstresscrp1(i:i)
        endif
      enddo
	   
	  Open(87,file=nitrogencrp)
	  Open(88,file=plantstresscrp)
	  
	  Write(85,5) "date", "jday", "time", "RSTAGE", "VSTAGE",         &
	  "PFD", "SolRad", "Tair", "Tcan", "Pgross", "Pnet", "gs", "PSIL",&
	  "LAI", "LAREAT", "totalDM", "rootDM", "stemDM", "leafDM",	      &
	  "seedDM", "podDM", "DeadDM", "Tr_pot", "Tr_act", "wstress", 	  &
	  "Nstress", "Limit"
	 
	  !Write(86,6) "date", "DayOfYear", "Time", "CWAD", "GWAD",		  &
	  !"PN", "Ags", "Canopy T"
	  
	  Write(87,7) "date", "jday","time","total_N","leaf_N","stem_N",  &
		          "pod_N", "seed_N", "root_N", "dead_N", "plant_N_C", &
	              "N_uptake", "N_demand", "Nstress"
      Write(88,8) "date", "jday","time", "wstress", "Nstress",        &
	              "Cstress", "NEffect_veg", "wstress2"
	 
5	  FORMAT (A10, 26(",",A10))
6	  FORMAT (6(A10))
7	  FORMAT (A10, 13(",",A10))
8	  FORMAT (A10, 7(",",A10))
	  
	  RETURN
13	  Write(*,*) 'Error in VARIETY FILE'
	  GOTO 11
11	  CONTINUE	         

   END