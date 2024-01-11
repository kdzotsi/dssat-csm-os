!======================================================================================
! MZ_AG_AGMAIZE
! Interface between DSSAT and AgMaize model
! This subroutine is called from plant.for
!-----------------------------------------------------------------------
! Revision History
! 04/17/2012 CHP/JZW Interface routine for Global-Maize model
! 08/23/2012 KAD Renamed and re-structured to add new subroutines
! 02/20/2013 KAD Created a new model "case" in plant.for
!-----------------------------------------------------------------------
! Called from: Plant subroutine: Plant
! Calls      : MZ_AG_Phenol, MZ_AG_Leafarea, MZ_AG_Photsynt, MZ_AG_Respiration,
!              MZ_AG_PartCoefs, MZ_AG_Growth, MZ_AG_RootGro, MZ_AG_Opgrow, MZ_AG_Opharv
! TODO:
! - Need to add yrend to drive harvest
! - harvres (harvest residue)
! - If harvesting at specific growth stages, send gstdyrdoySim to plant.for
! - Check N stress factors on different processes
! - Eventually use composite variables (like "growth") to transfer variables
! - Refactor GET_Real_Array_NL and PUT_Real_Array_NL in ModuleDefs.for. Separate AgMaize variables 
!   so that SPAM variables are not affected.
!======================================================================================

subroutine MZ_AG_AGMAIZE(control, iswitch, soilprop, weather, st,       & !Input
           sw, yrplt, yrend, snow, no3, nh4, eop, ep, trwup, harvfrac,  & !Input
           mdate, xlai, xhlai, rlv, nstres, senesce, unh4, uno3)          !Output

!----------------------------------------------------------------------
use ModuleDefs
use MZ_AG_ModuleDefs
implicit none
save
!----------------------------------------------------------------------
integer, intent(in) :: yrplt, yrend
real,    intent(in) :: snow, eop, trwup, harvfrac(2)
real, dimension(nl), intent(in) :: st, sw, no3, nh4
integer, intent(out) :: mdate

! Constructed variables based on definitions in ModuleDefs.for      
type (ControlType), intent(in) :: control
type (SoilType),    intent(in) :: SoilProp
type (SwitchType),  intent(in) :: iswitch
type (WeatherType), intent(in) :: weather
type(SpeciesType) dataspecies
type(FileioType)  datafileio
type(ResidueType) senesce
type(GrowthType)  growth
integer :: dynamic, yrdoy
     
! Phenology variables
integer :: growthStage, gstdyrdoySim(20), gstddoySim(20), gstdyearSim(20)
real :: gti, tlu, tnleafrm, tnleaf, lftips, anthes, silk, olg, dvs, dvsrate, tmpa, rla, gddcer, dtt, gddae
real :: gddten, gddapsim, growthtlu(10), tnleafphot, tnleaftemp, stLayer1, rlaSoil, rlaAir, lfnc, xstage 

! Leaf area variables
integer :: lfn
real :: lapot(50), greenla(50), greenlaSilk(50), lamaxpot(50), lflon(50)
real :: pla, xlai, xhlai, latf, lfcols, coltlu, maxlarf, tluladur, prlsen

! Photosynthesis and respiration variables
real :: pgross, maint, weolg, cvf

! Growth and partitioning
real :: pstres, fr, fsts, flv, fvear, kn, grt, grtplt, gsh, glv, gst, egr, ggr, wrt, wlv
real :: wlvtop, wlvbot, wsts, wstr, wst, we, grain, kw, tadrw, cgr, cgrf, cgrw, swfac, turfac
real :: tadrwAnth, tadrwSilk, seedno, skerwt, stover, wtmain, ep

! Root growth
real :: cumdep, rtdep, rlv(nl), rwuep1
character :: iswwat*1, iswnit*1, files*12, pathsr*80, pathex*80

! Nitrogen
real :: unh4(nl), uno3(nl), rootn, stovn, ranc, tanc, trnu, nfac, agefac, ndef3, nstres
real :: wtnvg, cannaa, wtnsd, wtncan, wtnup, pcngrn, pcnl

! Output
real :: bwah, sdwtah

! For checking radiation
real :: srad, parday, parsunday, parshday, radtotday
   
! Transfer values from constructed data types into local variables
dynamic = control % dynamic
yrdoy   = control % yrdoy
iswwat  = iswitch % iswwat
iswnit  = iswitch % iswnit

!----------------------------------------------------------------------
! Dynamic = runinit or dynamic = seasinit
!----------------------------------------------------------------------
if(dynamic==runinit .or. dynamic==seasinit) then   
  
  swfac  = 1.0
  turfac = 1.0
  nfac   = 1.0
  agefac = 1.0
  ndef3  = 1.0
  nstres = 1.0  
  
 !Read all sections of fileio and transfer variables
 call readfileio(control, 'ALLSEC', datafileio)
 files  = datafileio % files
 pathsr = datafileio % pathsr
 pathex = datafileio % pathex
  
 !Read root parameters from species file and transfer variables
 call readspecies(files, pathsr, '*ROOT ', dataspecies)
 rwuep1 = dataspecies % rwuep1
   
  call MZ_AG_Phenol(control, iswitch, SoilProp, weather, yrplt, sw, st, snow, cgr, nfac,         & !Input
       gti, tlu, tnleafrm, tnleaf, anthes, silk, olg, dvs, dvsrate, growthStage, tmpa, rla,      & !Output
       gddcer, dtt, gddae, gddten, gddapsim, mdate, gstdyrdoySim, gstddoySim, gstdyearSim,       & !Output 
       growthtlu, tnleafphot, tnleaftemp, stLayer1, rlaSoil, rlaAir, lfnc, coltlu, xstage)         !Output
       
  call MZ_AG_Leafarea(control, eop, ep, tlu, tnleaf, lfnc, tmpa, rla,   & !Input 
       dvs, olg, gstdyrdoySim, swfac, turfac, agefac,                   & !Input 
       pla, prlsen, xlai, xhlai, latf, lapot, greenla, greenlaSilk,     & !Output
       lamaxpot, lftips, lfcols, lfn, maxlarf, lflon, tluladur)           !Output 
       
  call MZ_AG_Photsynt(control, soilprop, weather, sw,          & !Input
       gddae, dvs, lfn, greenla, lflon, lamaxpot, lapot,       & !Input
       pgross, srad, parday, parsunday, parshday, radtotday)     !Output  
        
  call MZ_AG_Respiration(control, weather,                     & !Input
       dvs, wlvtop, wlvbot, wst, weolg, wrt, cgrw,             & !Input
       maint)                                                    !Output

!  call MZ_AG_RespirIXIM(control, weather,                      & !Input
!       pgross, wtmain, glv, grt, gst, egr, ggr,                & !Input
!       maint, cvf)                                               !Output
       
  call MZ_AG_PartCoefs(control, weather, yrplt,                & !Input
       silk, tlu, dvs, olg, pstres,                            & !Input
       fr, fsts, flv, fvear)                                     !Output 
       
  call MZ_AG_Growth(control, weather, iswitch, yrplt, eop, ep, trwup, dvs, tlu, olg, & !Input
       gstdyrdoySim, xhlai, prlsen, pgross, maint, fr, fsts, flv, fvear, swfac, nstres,      & !Input
       ndef3, pstres, kn, grt, grtplt, gsh, glv, gst, egr, ggr, wrt, wlv, wlvtop,    & !Output
       wlvbot, wsts, wstr, wst, we, weolg, grain, kw, tadrw, cgr, cgrf, cgrw, &        !Output
       tadrwAnth, tadrwSilk, seedno, skerwt, stover, wtmain, growth)                   !Output 
   
  call MZ_AG_RootGro(control, SoilProp, iswitch, no3, nh4, sw, & !Input
       growthStage, gstdyrdoySim, gddcer, dtt, grtplt, swfac,  & !Input                                   
       cumdep, rtdep, rlv)                                       !Output  
   
  if(iswnit /= 'N') then
     call MZ_AG_NPlant(control, iswitch, weather, soilprop, sw, yrplt, mdate,     & !Input
          xstage, growth, rlv, no3, nh4, gstdyrdoySim, growthStage, turfac, xhlai,       & !Input
          senesce, rootn, stovn, ranc, tanc, trnu, nfac, agefac, ndef3,           & !Output
          nstres, uno3, unh4, wtnvg, cannaa, wtnsd, wtncan, wtnup, pcngrn, pcnl)    !Output       
   else
     nfac   = 1.0
     agefac = 1.0
     ndef3  = 1.0
     nstres = 1.0               
  end if      
       
  call MZ_AG_Opgrow(control, iswitch, soilprop, yrplt, mdate, growthStage, gti, tlu, tnleafrm, tnleaf,  & !Input
       anthes, silk, olg, dvs, gddcer, rla, pla, xhlai, latf, greenla, gddten, gddapsim,                & !Input
       lftips, lfcols, coltlu, tnleafphot, tnleaftemp, maxlarf, lflon, tluladur, stLayer1, rlaSoil,     & !Input
       rlaAir, pgross, maint, fr, fsts, flv, fvear, kn, grt, gsh, glv, gst, egr, ggr,                   & !Input
       wrt, wlv, wlvtop, wlvbot, wsts, wstr, wst, we, grain, kw, tadrw, cgr, cgrf,cgrw, eop,            & !Input
       trwup, swfac, nstres, rtdep, rlv, srad, parday, parsunday, parshday, radtotday)    
      
  call MZ_AG_Opharv(control, iswitch, gstddoySim, gstdyearSim, growthtlu)
  
  call MZ_AG_Opharv2(control, iswitch, yrplt, harvfrac, swfac, turfac, senesce,   & !Input
       agefac, nstres, wtnvg, cannaa, wtnsd, wtncan, wtnup, pcngrn, growthStage,  & !Input 
       mdate, gstdyrdoySim, xhlai, tlu, tnleaf, greenla, greenlaSilk,             & !Input
       tadrwAnth, tadrwSilk, kn, seedno, we, skerwt, stover, tadrw, grain,        & !Input
       bwah, sdwtah)                                                                !Output

!----------------------------------------------------------------------
! Dynamic = rate
!----------------------------------------------------------------------
else if(dynamic == rate) then
  
  swfac = 1.0
  turfac = 1.0
  WATERON: if(iswwat /= 'N') then
     if(eop > 0.) then
       if(trwup*10. < eop) then
       !if(ep < eop) then
       ! swfac = ep/eop
        swfac = trwup*10./eop
        if(trwup*10./eop < rwuep1) turfac = trwup*10./(eop*rwuep1)
       !if(swfac < rwuep1) turfac = swfac/rwuep1
       end if
     end if    
  end if WATERON    
  swfac  = real(int(swfac*1000))/1000
  turfac = real(int(turfac*1000))/1000   
  
  call MZ_AG_Phenol(control, iswitch, SoilProp, weather, yrplt, sw, st, snow, cgr, nfac,         & !Input
       gti, tlu, tnleafrm, tnleaf, anthes, silk, olg, dvs, dvsrate, growthStage, tmpa, rla,      & !Output
       gddcer, dtt, gddae, gddten, gddapsim, mdate, gstdyrdoySim, gstddoySim, gstdyearSim,       & !Output 
       growthtlu, tnleafphot, tnleaftemp, stLayer1, rlaSoil, rlaAir, lfnc, coltlu, xstage)         !Output

if(growthStage >=1 .AND. growthStage <=10) then       
  call MZ_AG_Leafarea(control, eop, ep, tlu, tnleaf, lfnc, tmpa, rla,   & !Input 
       dvs, olg, gstdyrdoySim, swfac, turfac, agefac,                   & !Input 
       pla, prlsen, xlai, xhlai, latf, lapot, greenla, greenlaSilk,     & !Output
       lamaxpot, lftips, lfcols, lfn, maxlarf, lflon, tluladur)           !Output 
       
  call MZ_AG_Photsynt(control, soilprop, weather, sw,          & !Input
       gddae, dvs, lfn, greenla, lflon, lamaxpot, lapot,       & !Input
       pgross, srad, parday, parsunday, parshday, radtotday)     !Output  
        
  call MZ_AG_Respiration(control, weather,                     & !Input
       dvs, wlvtop, wlvbot, wst, weolg, wrt, cgrw,             & !Input
       maint)                                                    !Output

!  call MZ_AG_RespirIXIM(control, weather,                      & !Input
!       pgross, wtmain, glv, grt, gst, egr, ggr,                & !Input
!       maint, cvf)                                               !Output
       
  call MZ_AG_PartCoefs(control, weather, yrplt,                & !Input
       silk, tlu, dvs, olg, pstres,                            & !Input
       fr, fsts, flv, fvear)                                     !Output 
       
  call MZ_AG_Growth(control, weather, iswitch, yrplt, eop, ep, trwup, dvs, tlu, olg, & !Input
       gstdyrdoySim, xhlai, prlsen, pgross, maint, fr, fsts, flv, fvear, swfac, nstres,      & !Input
       ndef3, pstres, kn, grt, grtplt, gsh, glv, gst, egr, ggr, wrt, wlv, wlvtop,    & !Output
       wlvbot, wsts, wstr, wst, we, weolg, grain, kw, tadrw, cgr, cgrf, cgrw, &        !Output
       tadrwAnth, tadrwSilk, seedno, skerwt, stover, wtmain, growth)                   !Output 
  
  if(iswwat /= 'N') then     
  call MZ_AG_RootGro(control, SoilProp, iswitch, no3, nh4, sw,  & !Input
       growthStage, gstdyrdoySim, gddcer, dtt, grtplt, swfac,   & !Input                                   
       cumdep, rtdep, rlv)                                       !Output  
  end if   
  
  if(iswnit /= 'N') then
     call MZ_AG_NPlant(control, iswitch, weather, soilprop, sw, yrplt, mdate,     & !Input
          xstage, growth, rlv, no3, nh4, gstdyrdoySim, growthStage, turfac, xhlai,       & !Input
          senesce, rootn, stovn, ranc, tanc, trnu, nfac, agefac, ndef3,           & !Output
          nstres, uno3, unh4, wtnvg, cannaa, wtnsd, wtncan, wtnup, pcngrn, pcnl)    !Output       
  else
     nfac   = 1.0
     agefac = 1.0
     ndef3  = 1.0
     nstres = 1.0               
  end if     
     
end if

!----------------------------------------------------------------------
! Dynamic = integr
!----------------------------------------------------------------------
else if(dynamic == integr) then
  
  call MZ_AG_Phenol(control, iswitch, SoilProp, weather, yrplt, sw, st, snow, cgr, nfac,         & !Input
       gti, tlu, tnleafrm, tnleaf, anthes, silk, olg, dvs, dvsrate, growthStage, tmpa, rla,      & !Output
       gddcer, dtt, gddae, gddten, gddapsim, mdate, gstdyrdoySim, gstddoySim, gstdyearSim,       & !Output 
       growthtlu, tnleafphot, tnleaftemp, stLayer1, rlaSoil, rlaAir, lfnc, coltlu, xstage)         !Output

if(growthStage >=1 .AND. growthStage <=10) then       
  call MZ_AG_Leafarea(control, eop, ep, tlu, tnleaf, lfnc, tmpa, rla,   & !Input 
       dvs, olg, gstdyrdoySim, swfac, turfac, agefac,                   & !Input 
       pla, prlsen, xlai, xhlai, latf, lapot, greenla, greenlaSilk,     & !Output
       lamaxpot, lftips, lfcols, lfn, maxlarf, lflon, tluladur)           !Output 
       
  call MZ_AG_Photsynt(control, soilprop, weather, sw,            & !Input
       gddae, dvs, lfn, greenla, lflon, lamaxpot, lapot,         & !Input
       pgross, srad, parday, parsunday, parshday, radtotday)       !Output  
        
  call MZ_AG_Respiration(control, weather,                       & !Input
       dvs, wlvtop, wlvbot, wst, weolg, wrt, cgrw,               & !Input
       maint)                                                      !Output

!  call MZ_AG_RespirIXIM(control, weather,                      & !Input
!       pgross, wtmain, glv, grt, gst, egr, ggr,                & !Input
!       maint, cvf)                                               !Output
       
  call MZ_AG_PartCoefs(control, weather, yrplt,                  & !Input
       silk, tlu, dvs, olg, pstres,                              & !Input
       fr, fsts, flv, fvear)                                       !Output 
       
  call MZ_AG_Growth(control, weather, iswitch, yrplt, eop, ep, trwup, dvs, tlu, olg, & !Input
       gstdyrdoySim, xhlai, prlsen, pgross, maint, fr, fsts, flv, fvear, swfac, nstres,      & !Input
       ndef3, pstres, kn, grt, grtplt, gsh, glv, gst, egr, ggr, wrt, wlv, wlvtop,    & !Output
       wlvbot, wsts, wstr, wst, we, weolg, grain, kw, tadrw, cgr, cgrf, cgrw,        & !Output
       tadrwAnth, tadrwSilk, seedno, skerwt, stover, wtmain, growth)                   !Output 

  if(iswwat /= 'N') then       
  call MZ_AG_RootGro(control, SoilProp, iswitch, no3, nh4, sw, & !Input
       growthStage, gstdyrdoySim, gddcer, dtt, grtplt, swfac,   & !Input                                   
       cumdep, rtdep, rlv)                                       !Output  
  end if   
  
  if(iswnit /= 'N') then
     call MZ_AG_NPlant(control, iswitch, weather, soilprop, sw, yrplt, mdate,     & !Input
          xstage, growth, rlv, no3, nh4, gstdyrdoySim, growthStage, turfac, xhlai,       & !Input
          senesce, rootn, stovn, ranc, tanc, trnu, nfac, agefac, ndef3,           & !Output
          nstres, uno3, unh4, wtnvg, cannaa, wtnsd, wtncan, wtnup, pcngrn, pcnl)    !Output       
  else
     nfac   = 1.0
     agefac = 1.0
     ndef3  = 1.0
     nstres = 1.0               
  end if     
     
end if

!----------------------------------------------------------------------
! Dynamic = output
!----------------------------------------------------------------------
else if(dynamic == output) then

  if(iswnit /= 'N') then
     call MZ_AG_NPlant(control, iswitch, weather, soilprop, sw, yrplt, mdate,     & !Input
          xstage, growth, rlv, no3, nh4, gstdyrdoySim, growthStage, turfac, xhlai,       & !Input
          senesce, rootn, stovn, ranc, tanc, trnu, nfac, agefac, ndef3,           & !Output
          nstres, uno3, unh4, wtnvg, cannaa, wtnsd, wtncan, wtnup, pcngrn, pcnl)    !Output       
  else
     nfac   = 1.0
     agefac = 1.0
     ndef3  = 1.0
     nstres = 1.0               
  end if     
  
  call MZ_AG_Opgrow(control, iswitch, soilprop, yrplt, mdate, growthStage, gti, tlu, tnleafrm, tnleaf,  & !Input
       anthes, silk, olg, dvs, gddcer, rla, pla, xhlai, latf, greenla, gddten, gddapsim,                & !Input
       lftips, lfcols, coltlu, tnleafphot, tnleaftemp, maxlarf, lflon, tluladur, stLayer1, rlaSoil,     & !Input
       rlaAir, pgross, maint, fr, fsts, flv, fvear, kn, grt, gsh, glv, gst, egr, ggr,                   & !Input
       wrt, wlv, wlvtop, wlvbot, wsts, wstr, wst, we, grain, kw, tadrw, cgr, cgrf, cgrw, eop,           & !Input
       trwup, swfac, nstres, rtdep, rlv, srad, parday, parsunday, parshday, radtotday)    
       
  call MZ_AG_Opharv(control, iswitch, gstddoySim, gstdyearSim, growthtlu)
  
  call MZ_AG_Opharv2(control, iswitch, yrplt, harvfrac, swfac, turfac, senesce,   & !Input
       agefac, nstres, wtnvg, cannaa, wtnsd, wtncan, wtnup, pcngrn, growthStage,  & !Input 
       mdate, gstdyrdoySim, xhlai, tlu, tnleaf, greenla, greenlaSilk,             & !Input
       tadrwAnth, tadrwSilk, kn, seedno, we, skerwt, stover, tadrw, grain,        & !Input
       bwah, sdwtah)                                                                !Output
                                                                             
  
!----------------------------------------------------------------------
! Dynamic = seasend
!----------------------------------------------------------------------
else if(dynamic == seasend) then

  if(iswnit /= 'N') then
     call MZ_AG_NPlant(control, iswitch, weather, soilprop, sw, yrplt, mdate,     & !Input
          xstage, growth, rlv, no3, nh4, gstdyrdoySim, growthStage, turfac, xhlai,      & !Input
          senesce, rootn, stovn, ranc, tanc, trnu, nfac, agefac, ndef3,           & !Output
          nstres, uno3, unh4, wtnvg, cannaa, wtnsd, wtncan, wtnup, pcngrn, pcnl)    !Output       
  else
     nfac   = 1.0
     agefac = 1.0
     ndef3  = 1.0
     nstres = 1.0               
  end if       
 
  call MZ_AG_Opgrow(control, iswitch, soilprop, yrplt, mdate, growthStage, gti, tlu, tnleafrm, tnleaf,  & !Input
       anthes, silk, olg, dvs, gddcer, rla, pla, xhlai, latf, greenla, gddten, gddapsim,                & !Input
       lftips, lfcols, coltlu, tnleafphot, tnleaftemp, maxlarf, lflon, tluladur, stLayer1, rlaSoil,     & !Input
       rlaAir, pgross, maint, fr, fsts, flv, fvear, kn, grt, gsh, glv, gst, egr, ggr,                   & !Input
       wrt, wlv, wlvtop, wlvbot, wsts, wstr, wst, we, grain, kw, tadrw, cgr, cgrf,cgrw, eop,            & !Input
       trwup, swfac, nstres, rtdep, rlv, srad, parday, parsunday, parshday, radtotday)    
       
  call MZ_AG_Opharv(control, iswitch, gstddoySim, gstdyearSim, growthtlu)
  
  call MZ_AG_Opharv2(control, iswitch, yrplt, harvfrac, swfac, turfac, senesce,   & !Input
       agefac, nstres, wtnvg, cannaa, wtnsd, wtncan, wtnup, pcngrn, growthStage,  & !Input 
       mdate, gstdyrdoySim, xhlai, tlu, tnleaf, greenla, greenlaSilk,             & !Input
       tadrwAnth, tadrwSilk, kn, seedno, we, skerwt, stover, tadrw, grain,        & !Input
       bwah, sdwtah)                                                                !Output
       
  ! Copy daily output file to filex's directory
  !call copy_file('*.OUT', trim(pathex))
  !call copy_file(outaltmz, trim(pathex)//outaltmz)

  
!-----------------------------------------------------------------------
! End of dynamic if structure
!-----------------------------------------------------------------------
end if
      
!----------------------------------------------------------------------
! End of subroutine
!----------------------------------------------------------------------
return
end subroutine MZ_AG_AGMAIZE
     

!==============================================================================================================================
! Subroutine and module definitions
!------------------------------------------------------------------------------------------------------------------------------
! ModuleDefs          Contains definitions of global constants and variables and contructed variables
! ModuleData          Contains definitions of constructed variables and interface subroutines for storing and retrieving data
! MZ_AG_ModuleData    Contains definitions of AgMaize constants and parameters that are genotype-files candidates
! MZ_AG_ModuleDefs    Contains definitions of AgMaize constructed variables and subroutines for reading input files
! MZ_AG_Phenol        Phenology subroutine
! MZ_AG_Leafarea      Leaf area subroutine
! MZ_AG_Photsynt      Subroutine for radiation absorption and photosynthesis
! MZ_AG_Respiration   Subroutine for calculating respiration
! MZ_AG_PartCoefs     Subroutine for calculating coefficients for partitioning assimilates
! MZ_AG_Growth        Subroutine for partitioning net assimilates to different plant components
! MZ_AG_RootGro       Root growth subroutine. Calculates the root length density
! MZ_AG_Opgrow        Subroutine for writing daily outputs to file
! MZ_AG_Opharv        Subroutine for writing end-of-season output to file
!==============================================================================================================================
     
!==============================================================================================================================
! Variable definitions                                                                                        Unit 
!------------------------------------------------------------------------------------------------------------------------------
! anthes           Leaf stage at anthesis                                                                     tlu   
! cgr              Crop growth rate (aboveground dry matter)                                                  kg[dm]/ha/day
! cgrw             Mean crop growth rate (aboveground dry matter) during the previous week                    kg[dm]/ha/day
! coltlu           Thermal leaf units from leaf tip appearance to full leaf expansion                         tlu
! control          Constructed type for control variables
! cumdep           Cumulative depth, to bottom of soil layer                                                  cm
! dtt              Daily thermal time (for current day)                                                       degree-day
! dvs              Development stage (0 is planting, 1 is silking, and 2 is maturity) 
! dvsrate          Rate of development stage
! dynamic          Modular control (runinit=1 for run initialization, seasinit=2 for seasonal initialization,
!                  rate=3 for rate calculations, integr=4 for integration of state variables, output=5 for
!                  writing daily outputs, seasend=6 for closing output files                                 
! egr	           Ear growth rate	                                                                          kg[dm]/ha/day
! eop              Potential plant transpiration                                                              mm/day
! flv              Proportion of assimilates partitioned to leaves                                            fraction
! fvear	           Effect of stage of development on partitioning to non-grain ear part                       fraction
! fr               Proportion of assimilates partitioned to roots                                             fraction
! fsts             Proportion of assimilates partitioned to structural stems                                  fraction
! gddae            Cumulative daily thermal time after emergence (CERES method)                               degree-day
! gddapsim         Cumulative daily thermal time after planting using APSIM's approach                        degree-day
! gddcer           Cumulative daily thermal time after planting (CERES method)                                degree-day                                                            
! gddten           Cumulative daily thermal time after planting using GDD[10,30]                              degree-day
! ggr	           Grain growth rate	                                                                      kg[dm]/ha/day
! glv	           Growth rate of leaves' dry matter	                                                      kg[dm]/ha/day
! grain	           Grain dry matter	                                                                          kg[dm]/ha
! greenla(i)       Cumulative area of leaf i during expansion and senescence                                  m2/leaf
! growthStage      Integer value of growth stage
! growthtlu(i)     Leaf stage at which growth stage i occurred                                                tlu
! grt	           Growth rate of root dry matter	                                                          kg[dm]/ha/day
! grtplt            Root growth rate                                                                           g/plant/day
! gsh	           Growth rate of vegetative components of shoot	                                          kg[dm]/ha/day
! gst	           Growth of stem dry matter	                                                              kg[dm]/ha/day
! gstddoySim(i)    Component day of year of gstdyrdoySim                                                      ddd
! gstdyearSim(i)   Component year of gstdyrdoySim                                                             yyyy
! gstdyrdoySim(i)  Simulated year and day of year for growth stage i                                          yyyyddd
! gti              General Thermal Index heat unit accumulation                                               gti (Q.what is the real unit of gti?)
! iswitch          Constructed type for control switches
! kn               Kernel number                                                                              kernel/plant
! kw	           Kernel weight 	                                                                          mg/kernel
! lamaxpot(i)      Maximum area of leaf i under non-stressed conditions                                       m2/leaf
! lapot(i)         Cumulative area of leaf i during expansion (tip to collar) under non-stressed conditions   m2/leaf
! latf             Effect of temperature on leaf growth
! lfcols           Leaf node position of leaf that has completed expansion
! lftips
! lflon(i)         Longevity of leaf i                                                                        degree-day
! lfn              Total leaf number rounded up
! maint            Maintenance respiration	                                                                  kg[glucose]/ha/day 
! maxlarf          Reduction factor of maximum leaf area due to plant density
! mdate            Maturity date                                                                              yyyyddd
! mulch            Constructed type for mulch variables
! nh4(l)           Ammonium in soil layer l                                                                   ppm
! no3(l)           Nitrate in soil layer l                                                                    ppm
! olg              Onset linear dry matter accumulation of grain                                              dvs
! pgross           Daily canopy gross assimilation (glucose equivalent in kg/ha/day)                          kg[glucose]/ha/day
! pla              Plant leaf area                                                                            m2/plant
! pstres           Fraction of stem reserves                                                                  fraction
! rla              Rate of Leaf tip Appearence                                                                leaves/day
! rlaAir           Rate of leaf tip appearance using air temperatures                                         leaves/day
! rlaSoil          Rate of leaf tip appearance uing soil temperature of layer 1                               leaves/day
! rlv(l)           Root length density of layer l                                                             cm[root]/cm3[soil]
! rtdep            Rooting depth, initially set at emergence                                                  cm
! silk             Leaf stage at silking                                                                      tlu 
! snow             Snow depth                                                                                 mm
! SoilProp         Constructed type for variables related to soil properties
! st(l)            Soil temperature in soil layer l                                                           degree C
! stLayer1         Soil temperature of layer 1                                                                degree C
! sw(l)            Volumetric soil water content in layer l                                                   cm3[water]/cm3[soil]
! swfac            Soil water stress effect on growth (0-1), 1 is no stress, 0 is full 
! tadrw	           Total above ground dry weight	                                                          kg[dm]/ha
! tlu              Thermal leaf unit (cumulative rate of leaf tip appearance)                                 tlu
! tluladur         Duration of leaf expansion                                                                 tlu
! tmpa             Daily mean air temperature                                                                 degree C
! tnleaf           Total number of initiated leaves                                                           leaves
! tnleafphot       Change in leaf number due to photoperiod                                                   leaves 
! tnleaftemp       Change in leaf number due to temperature                                                   leaves    
! trwup            Total potential daily root water uptake                                                    cm/day
! we	           Weight of ears	                                                                          kg[dm]/ha
! weather          Constructed type for weather variables
! weolg            Weight of ears at onset of linear grain filling                                            kg[dm]/ha
! wlv	           Total weight of leaves	                                                                  kg[dm]/ha
! wlvbot           Weight of leaves below top layer (wlv - wlvtop)                                            kg[dm]/ha               
! wlvtop           Weight of leaves in top layer (lai <= 2)                                                   kg[dm]/ha
! wrt              Weight of roots                                                                            kg[dm]/ha
! wst	           Weight of stems	                                                                          kg[dm]/ha
! wstr	           Weight of (temporarily-stored) stem reserves                                               kg[sucrose]/ha
! wsts	           Weight of structural stem 	                                                              kg[dm]/ha
! wtmain
! xhlai            Whole-plant leaf area index                                                                m2[leaf]/m2[ground]
! yrplt            Planting date                                                                              yyyyddd
!==============================================================================================================================
