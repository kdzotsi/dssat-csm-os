!======================================================================================
! MZ_GM_LeafArea
! Calculate leaf area expansion based on Thijs Tollenaar' model:
!    (1) Leaf area distribution by leaf number (Dwyer & Stewart, 1987)
!    (2) Rate and duration of expansion of individual leaves (Stewart & Dwyer, 1994)
!    (2) Effect of temperature on LAI from 4 to 12 leaf stage (Tollenaar, 1989)
!    (3) Effect of plant density on LAI (Grignon 1994 dataset)
!----------------------------------------------------------------------
! Revision History
! 08/15/2012 TT/SK/KAD Translated from MS Basic (originally written by Tollenaar et al.)
!----------------------------------------------------------------------
! Called from: Main program: MZ_AG_AGMAIZE
! Calls      : None
!----------------------------------------------------------------------
! TO DO: 
! - Improve nitrogen stress effect on leaf senescence
!----------------------------------------------------------------------
!======================================================================================
subroutine MZ_AG_Leafarea(control, eop, ep, tlu, tnleaf, lfnc, tmpa, rla,       & !Input 
           dvs, olg, gstdyrdoySim, swfac, turfac, agefac,                       & !Input 
           pla, prlsen, xlai, xhlai, latf, lapot, greenla, greenlaSilk,         & !Output
           lamaxpot, lftips, lfcols, lfn, maxlarf, lflon, tluladur)               !Output    

use ModuleDefs
use ModuleData
use MZ_AG_ModuleDefs
implicit  none
save      

!----------------------------------------------------------------------
! Define Variables
!----------------------------------------------------------------------
integer, parameter :: d=5, lenla=50, m2tocm2=1.E4                          
real, intent(in)  :: tlu, tnleaf, tmpa, rla, dvs, olg
type (ControlType):: control                          
type(FileioType)  :: datafileio      
type(SpeciesType) :: dataspecies                             
real, intent(out) :: pla, xlai, xhlai, latf, lftips, lfcols, maxlarf, tluladur
real, dimension(lenla) :: lapot, greenla, greenlaSilk, lamaxpot, lflon, lflontlu

! Other variables
character :: files*12, pathsr*80
integer :: dynamic, yrdoy, idelt, iidelt, i, ii, icol, itip, lfn, gstdyrdoySim(20)
real :: pltpop, lftop, lalfx, ampli, asymp, agefac, swfac, turfac, tlud, lfnc, lnl, wha
real :: tempCoef(3,3), decimal, coltlu, eop, ep, rlastress, rlawatstres, rlanitstres, prlsen
real, dimension(lenla) :: rlepot, rle, rlspot, rls, lamax, tlulabeg, tlulaend, tlustress, collartlu, lfcollars

! Functions
real :: getcoltlu, tempeffect

! Transfer values from constructed data types into local variables
dynamic = control % dynamic
yrdoy   = control % yrdoy
     
!----------------------------------------------------------------------
! Dynamic = runinit or dynamic = seasinit
!----------------------------------------------------------------------
if(dynamic==runinit .OR. dynamic==seasinit) then
 
  !Initialize some variables
  latf        = 1.0        
  pla         = 0.0
  prlsen      = 0.0
  xlai        = 0.0
  xhlai       = 0.0
  greenla     = 0.0
  greenlaSilk = 0.0
  lapot       = 0.0
  lamaxpot    = 0.0
  lamax       = 0.0
  tlulabeg    = 0.0
  tlulaend    = 0.0
  lftips      = 0.0
  lfcols      = 0.0
  icol        = 0
  itip        = 0
  coltlu      = 0.0
  maxlarf     = 1.0
  lflon       = 0.0
  lflontlu    = 0.0
  tluladur    = 0.0
  lfn         = 0
  lfcollars   = 0.0
  tlustress   = 0.0
  rlastress   = 0.0
  rlawatstres = 0.0
  rlawatstres = 0.0
  rle         = 0.0
  rlepot      = 0.0
  rls         = 0.0
  rlspot      = 0.0  
  
  !Read all sections of fileio and transfer variables
  call readfileio(control, 'ALLSEC', datafileio)
  pltpop = datafileio % pltpop 
  lftop  = datafileio % lftop
  lalfx  = datafileio % lalfx * m2tocm2
  ampli  = datafileio % ampli
  asymp  = datafileio % asymp
  files  = datafileio % files
  pathsr = datafileio % pathsr
  
  !Read temperature coefficients parameters from species file and transfer variables
  call readspecies(files, pathsr, '*LEAF ', dataspecies)
  tempCoef = dataspecies % tempCoef
   
!----------------------------------------------------------------------
! Dynamic = rate 
!----------------------------------------------------------------------
else if(dynamic == rate) then


!-----------------------------------------------------------------------
! Dynamic = integr
!-----------------------------------------------------------------------
else if(dynamic == integr) then
      
  !- Temperature effect on leaf expansion
  latf = tempeffect(tlu, tmpa, tempCoef)
  if(latf <= 0.0) latf = 0.0    !To prevent negative latf when tmpa is not yet computed
        
  !- Preliminary calculations
  lfn = ceiling(tnleaf)         !This is equivalent to the Basic code fix(tnleaf+1) since tnleaf > 0
  
  !- Leaf senescence parameters (see Lisazo et al. 2003)
  lnl = lfnc
  wha = tnleaf/3.0              !Equation 12, Width of the relationship between leaf longevity and nodal position at half amplitude

  !- Maximum leaf area reduction factor; based on Lizaso's Grignon dataset 1994; 2/01/2013
  !  For low plant lai (< 2) reduction factor depends on plant density (plants sense each other)
  !  If lai >= 2, reduction factor depends upon lai itself because of shading.
  if(xhlai < 2.0) then
     maxlarf = exp(-0.0085*(max(pltpop,4.0) - 4.0))     !Attn: this equation is again used later in the season when senescence brings LAI down to < 2
  else  
     maxlarf = min(1.1589 - 0.0873*xhlai, exp(-0.0085*(max(pltpop,4.0)-4.0)))
  end if  
  
  !- Effect of water and nitrogen stress on leaf senescence (rlastress)
  if(swfac < 1.0 .OR. agefac < 1.0) then
     rlawatstres = rla*(-4.0*swfac  + 4.0)
     !rlanitstres = rla*(-1.0*agefac + 1.0)
     rlanitstres = 0.0
  else
     rlawatstres = 0.0
     rlanitstres = 0.0
  end if
  rlastress = max(rlawatstres, rlanitstres)
  
  !- Leaf expansion loop
  tlud = tlu - rla               !tlu is cumulative rla so tlud is yesterday's tlu
  do idelt = 1, d                
     tlud = tlud + rla/d         !Increment tlud by (1/d)th in this loop
    
     EXPANDLEAF: do i = 1, lfn                
        if(tlud < i) exit EXPANDLEAF                 !Not enough tlu for tip of leaf i to appear
        if(tlulabeg(2) /= 0.0) tlulabeg(1) = tlulabeg(2)
        if(tlulabeg(i) == 0.0) tlulabeg(i) = tlud    !Tip of leaf i has already appeared, save corresponding tlu
      
        !Compute collar tlu for each expanding leaf depending on their position relative to the largest leaf (lfnc)
  	    !collar tlu beyond lfnc increases linearly to a maximum of tnleaf+3
  	    if(real(i) < tnleaf) then
  	       coltlu = getcoltlu(real(i), lfnc)
  	    else 
  	       coltlu = coltlu + tnleaf - real(i-1)
        end if
        collartlu(i) = coltlu
  	  
        EXPANDING: if(tlud <= coltlu) then         !Leaf is still expanding
           !Compute potential area of leaf i (Stewart and Dwyer 1994, equation 9)
           lamaxpot(i) = lalfx * exp(-0.0359*(i-lfnc)**2. + 0.00117*(i-lfnc)**3.)   !lfnc = 12.4, to be checked against data
              
           !Rate of leaf area expansion (unit area/day/leaf). Area of an expanding leaf between the date
           !of leaf-tip and leaf-collar appearance = lamaxpot(i) * (tlud - i) / (coltlu – i); adapt this to have a logistic type growth?
           rlepot(i) =  maxlarf * latf * lamaxpot(i) * (rla/real(d)) / (coltlu-real(i))    
           rle(i) = rlepot(i)*min(turfac, agefac)
  	       
 	       !Integrate for each leaf (move to integrate); lapot(i) is cumulative leaf area for leaf i (m2/leaf)
  	       if(real(i) <= tnleaf) then 
              lapot(i) = lapot(i) + rlepot(i)
              greenla(i) = greenla(i) + rle(i)                         
           else 
           !Since lfn is ceiling(tnleaf), add growth for the decimal part of tnleaf
              lapot(i) = lapot(i) + rlepot(i)*decimal(tnleaf)        
              greenla(i) = greenla(i) + rle(i)*decimal(tnleaf)
           end if
                      
        else EXPANDING                            !Expansion is completed for current leaf
            if(tlulaend(i) == 0.0) then
              tlulaend(i) = tlud  
              lamax(i) = greenla(i)               !Maximum leaf area after accounting for stresses
              icol = i                            !Node number of collar that has appeared (or leaf that has completed expansion)
              lfcols = real(icol)
              tluladur = tlulaend(i) - tlulabeg(i)              
           end if
           
        end if EXPANDING                  !Leaf is still epanding
  
     end do EXPANDLEAF    !Leaf expansion loop
     
     !- Leaf senescence
     !The rate of leaf senescence can also increase due to a higher leaf area at the beginning
     !of senescence
     SENESCELEAF: do ii = 1, lfn
        if(tlulabeg(ii) == 0.0) exit                                   !Leaf(ii) has not started expanding
        lflon(ii) = asymp + ampli*exp(-(ii-lnl)**2.0/(2.0*wha**2.0))   !Longevity in GDD (Eq. 10, Lizaso et al. 2003 FCR)
        lflontlu(ii) = (0.0264*lflon(ii) - 0.2467) - tlustress(ii)     !Longevity in TLU and adjusted for stress
        lflon(ii) = (lflontlu(ii) + 0.2467)/0.0264                     !Longevity in GDD adjusted for stress
        if(tlud <= (tlulabeg(ii)+lflontlu(ii))) then                   !Senescence has not started yet   
           rlspot(ii) = 0.0
           rls(ii) = 0.0                                     
           tlustress(ii) = tlustress(ii) + rlastress/d * (tlud-tlulabeg(ii))/lflontlu(ii)
           if(rlastress == 0.) tlustress(ii) = 0.8 * tlustress(ii)
        else if(greenla(ii) == 1.E-4) then     !Leaf has fully senesced
           cycle SENESCELEAF 
        else
           !Rate of leaf senescence (negative rate)
           rlspot(ii) = (-1) * lamax(ii) * (rla/d) / (tlulaend(ii)-tlulabeg(ii))
           rls(ii) = rlspot(ii) 
        end if
             
        !Integrate; greenla getting negative, senescence rate might be too high. 
        !We need to address this specifically instead of truncating at 0.0001
        greenla(ii) = max(1.E-4, greenla(ii) + rls(ii))           
     end do SENESCELEAF
    
  end do                  !tlud increment loop

  !- Check upper bound: some problems related to photoperiod sensitivity:
  !lamaxpot dependency on a varying tnleaf
  !lapot(1:lfn) = (/ (min(lapot(i), lamaxpot(i)), i=1,lfn) /)
  !greenla(1:lfn) = (/ (min(greenla(i), lamaxpot(i)), i=1,lfn) /)
 
  !- Plant leaf area and LAI 
  if(tlu < 2.0) then 
     pla = 0.0
     prlsen = 0.0
  else   
     pla = sum(greenla)/m2tocm2           !in m2/plant
     prlsen = -1.0 * sum(rls)/m2tocm2     !Leaf senescence (m2/plt, positive value)
  end if
  xhlai = pla*pltpop               !in m2/m2
  xlai = xhlai
  if(yrdoy == gstdyrdoySim(7)) greenlaSilk = greenla
  
  !- Set leaf tips and leaf ligules to be at most equal to total leaf number
  lftips = min(tlu,    tnleaf)
  lfcols = min(lfcols, tnleaf)
  
  !- Save variables for use in SPAM
  call put('PLANT','LFN',      lfn)
  call put('PLANT','GREENLA',  greenla)
  call put('PLANT','LFLON',    lflon)
  call put('PLANT','LAMAXPOT', lamaxpot)
  call put('PLANT','LAPOT',    lapot)
	   
!-----------------------------------------------------------------------
! End of dynamic if structure
!-----------------------------------------------------------------------
end if  

!----------------------------------------------------------------------
! End of subroutine
!----------------------------------------------------------------------
return      
end subroutine MZ_AG_LeafArea
 
!==============================================================================================================================
! Functions used in this leaf area subroutine
!-----------------------------------------------------------------------      
! Function to compute temperature effect on leaf expansion
! Depends only on average daily temperature and replaces the previous
! data based interpolation procedure. 9/14/2012
!-----------------------------------------------------------------------    
function tempeffect(leafNum, tAvg, tCoef) result(latf)
   real, intent(in) :: leafNum, tavg, tCoef(3,3)
   real :: latf
   real :: latf1, latf2
   if(leafNum < 4.0) then
       latf = tCoef(1,1)*tAvg**2 +tCoef(1,2)*tAvg +tCoef(1,3)
   else if(leafNum >= 4.0 .and. leafNum < 8.0) then
       latf1 = tCoef(1,1)*tAvg**2 +tCoef(1,2)*tAvg +tCoef(1,3)
       latf2 = tCoef(2,1)*tAvg**2 +tCoef(2,2)*tavg +tCoef(2,3)
       latf = latf1 + (latf2-latf1)/4*(leafNum-4.0)
   else if(leafNum >= 8.0 .and. leafNum < 12.0) then
       latf1 = tCoef(2,1)*tAvg**2 +tCoef(2,2)*tAvg +tCoef(2,3)
       latf2 = tCoef(3,1)*tAvg**2 +tCoef(3,2)*tAvg +tCoef(3,3)
       latf = latf1 + (latf2-latf1)/4*(leafNum-8.0)
   else
       latf = tCoef(3,1)*tAvg**2 +tCoef(3,2)*tAvg +tCoef(3,3)
   end if
end function tempeffect                  

      
!==============================================================================================================================
! Variable definitions                                                                                     Unit 
!------------------------------------------------------------------------------------------------------------------------------
! agefac        !Nitrogen stress factor affecting expansion (0-1)
! lalfx         Area of the largest leaf on the plant under non-stressed conditions                        cm2
! ampli         Amplitude of the leaf longevity vs. leaf nodal position relationship                       degree-day
! asymp         Asymptote of the leaf longevity vs. leaf nodal position relationship                       degree-day
! coltlu        Thermal leaf units from leaf tip appearance to full leaf expansion                         tlu
! control       Constructed type for control variables
! datafileio    Constructed type for variables read from DSSAT's input/output file (DSSAT45.INP)
! d             Number of steps for calculation of leaf area expansion and senescence
! dvs           Development stage (0 is planting, 1 is silking, and 2 is maturity) 
! dynamic       Modular control (runinit=1 for run initialization, seasinit=2 for seasonal initialization,
!               rate=3 for rate calculations, integr=4 for integration of state variables, output=5 for
!               writing daily outputs, seasend=6 for closing output files                                                                            
! getcoltlu     [Function] For computing thermal leaf unit at which leaf is fully expanded 
! greenla(i)    Cumulative area of leaf i during expansion and senescence                                  cm2/leaf
! greenlaSilk(i)
! i             Integer loop counter
! ii            Counter of leaf number during leaf senescence
! idelt         Counter of D steps during leaf expansion
! iidelt        Counter of D steps during leaf senescence
! lamax(i)      Maximum area actually reached by leaf i at full expansion                                  m2/leaf
! lamaxpot(i)   Maximum area of leaf i under non-stressed conditions                                       m2/leaf
! lapot(i)      Cumulative area of leaf i during expansion (tip to collar) under non-stressed conditions   m2/leaf
! latf          Effect of temperature on leaf growth
! lenla         Arbitrary definition of the length of leaf area vector
! lfcols        Leaf node position of leaf that has completed expansion
! lfn           Total leaf number rounded up
! lfnc          Node position of largest leaf 
! lflon(i)      Longevity of leaf i in GDD                                                                 degree-day
! lflon(i)      Longevity of leaf i in TLU                                                                 TLU
! lnl           Node position of leaf with the largest longevity
! maxlarf       Reduction factor of maximum leaf area due to plant density
! m2tocm2       Constant for converting areas from m2 to cm2
! olg           Onset linear dry matter accumulation of grain                                              dvs
! pla           Plant leaf area                                                                            m2/plant
! pltpop        Plant density                                                                              plants/m2
! prlsen
! rla           Mean daily rate of leaf appearance                                                         tlu/day
! rlastress
! rle(i)        Actual rate of area expansion of leaf i                                                    m2/leaf/day
! rlepot(i)     Rate of area expansion of leaf i under non-stressed conditions                             m2/leaf/day
! rls(i)        Actual rate of area senescence of leaf i                                                   m2/leaf/day
! rlspot(i)     Rate of area senescence of leaf i under non-stressed conditions                            m2/leaf/day
! tempCoef      Coefficient for estimating temperature effect on leaf expansion
! tempeffect    [Function] For calculating temperature effect on leaf expansion
! tlu           Cumulative thermal leaf unit                                                               tlu
! tlud          Accumulated thermal leaf unit until yesterday                                              tlu
! tlulabeg(i)   Thermal leaf unit from planting to the appearance of leaf i's tip                          tlu
! tluladur      Duration of leaf expansion                                                                 tlu
! tlulaend(i)   Thermal leaf unit from planting to the appearance of leaf i's collar                       tlu
! tlustress
! tmpa          Mean air temperature for current day                                                       degree C
! tnleaf        Total number of initiated leaves on the plant                                              leaves
! turfac        Water stress factor affecting leaf expansion
! wha           Width of the leaf nodal position vs. leaf longevity relationship at half amplitude         leaves   
! xhlai         Whole-plant leaf area index                                                                m2[leaf]/m2[ground]
!==============================================================================================================================

