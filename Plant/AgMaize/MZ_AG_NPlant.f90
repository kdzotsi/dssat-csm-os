!======================================================================
!  MZ_AG_NPlant, Subroutine
!  Determines Nitrogen transformations in the plant
!----------------------------------------------------------------------
!  Revision history
!  07/14/2014 Subroutine created from N model extracted from CERES-Maize
!----------------------------------------------------------------------
!  Called : MZ_AG_AGMAIZE  
!  Calls  : MZ_AG_NUptake, MZ_AG_OpNPlant
! 
!  TO DO:
!  - Enable N concentrations to start at optimum
!  - Calculate leaf senescence and add to litter
!  - Check units when cumulative N senesced is calculated
!  - Senesced N not substracted from shoot N
!----------------------------------------------------------------------

subroutine MZ_AG_NPlant(control, iswitch, weather, soilprop, sw, yrplt, mdate,      & !Input
           xstage, growth, rlv, no3, nh4, gstdyrdoySim, growthStage, turfac, xhlai, & !Input
           senesce, rootn, stovn, ranc, tanc, trnu, nfac, agefac, ndef3,            & !Output
           nstres, uno3, unh4, wtnvg, cannaa, wtnsd, wtncan, wtnup, pcngrn, pcnl)    !Output       

use ModuleDefs
use MZ_AG_ModuleDefs
implicit  none
save

!----------------------------------------------------------------------
! DEFINITION OF VARIABLES
!----------------------------------------------------------------------
real,parameter :: frleafn=0.75
real,dimension(nl) :: rlv, no3, nh4, uno3, unh4, sw, dlayr
character :: iswnit*1, files*12, pathsr*80
integer :: l, nlayr, dynamic, yrdoy, yrplt, mdate, gstdyrdoySim(20), growthStage
real :: pltpop, tmnc2,  ctcnp1, ctcnp2, tance, rance, pligrt, pliglf, tmin, tmax, tempm, xstage
real :: tcnp, tanc, tmnc, rcnp, ranc, rmnc, gnp, frstemn
real :: pcnrt, pcnst, pcnl, pcnveg, pcngrn, wtnrt, wtnst, wtnlf, wtnvg, wtnsd, wtncan, wtnup, cannaa
real :: drootn, rootn, dstovn, stovn, dstemn, stemn, dleafn, leafn, nsink, grainn, ptf
real :: agefac, ndef3, nstres, nfac, sfac, tfac, turfac, CumLeafSenes, CumLeafSenesY, CumLfNSenes
real :: npool, npool1, npool2, xnf, rnlab, tnlab, trnu, dgsh, rndem, tndem, ndem, nsdr
real :: pgrort, gshplt, grogrn, rtwt, stmwt, lfwt, stovwt, grnwt, rnout
real :: trlv, rnloss, trnlos, rootnloss, shutnloss, stemnloss, leafnloss, glvsenplt, leafndead, slw, sln, xhlai

! Constructed variables based on definitions in modules
type(ControlType) control
type(SwitchType)  iswitch
type(WeatherType) weather
type(SoilType)    soilprop
type(ResidueType) senesce
type(GrowthType)  growth
type(FileioType)  datafileio
type(SpeciesType) dataspecies

! Transfer values from constructed data types into local variables
dynamic   = control  % dynamic
yrdoy     = control  % yrdoy
iswnit    = iswitch  % iswnit
tmax      = weather  % tmax
tmin      = weather  % tmin
dlayr     = soilprop % dlayr
nlayr     = soilprop % nlayr  
rtwt      = growth   % wrtplt
lfwt      = growth   % wlvplt
stmwt     = growth   % wstplt
stovwt    = growth   % stoverplt
grnwt     = growth   % grainplt
pgrort    = growth   % grtplt
gshplt    = growth   % gshplt
glvsenplt = growth   % glvsenplt
grogrn    = growth   % ggrplt
ptf       = growth   % ptf
slw       = growth   % slw


!----------------------------------------------------------------------
! DYNAMIC = RUNINIT OR SEASINIT
!----------------------------------------------------------------------
if(dynamic==runinit .OR. dynamic==seasinit) then

   !Read fileio (DSSATxx.INP) and transfer variables
   call readfileio(control, 'ALLSEC', datafileio)
   pltpop = datafileio % pltpop 
   files  = datafileio % files
   pathsr = datafileio % pathsr
   
   !Read species file and transfer variables
   call readspecies(files, pathsr, 'ALLSEC', dataspecies)
   tmnc2  = dataspecies % tmnc
   rcnp   = dataspecies % rcnp
   ctcnp1 = dataspecies % ctcnp1
   ctcnp2 = dataspecies % ctcnp2
   tance  = dataspecies % tance
   rance  = dataspecies % rance
   pligrt = dataspecies % pligrt
   pliglf = dataspecies % pliglf
   
   !Calculate stem N partitioning coefficient
   frstemn = 1.0 - frleafn
   
   !Critical concentrations (g[N]/g[DM])
   tcnp = tance       !Shoots
   rcnp = rance       !Roots
   
   !Minimum concentrations (g[N]/g[DM])
   tmnc = 0.0         !Shoots
   rmnc = 0.0         !Roots
   
   !Actual concentrations (g[N]/g[DM])
   tanc = tance       !Shoots
   ranc = rance       !Roots
   
   !Actual concentrations after partitioning to different tissues (%, i.e. g[N]/100g[DM])
   pcnrt  = 0.0       !Roots
   pcnveg = 0.0       !Vegetative tissues
   pcnl   = 0.0       !Leaves
   pcnst  = 0.0       !Stems
   pcngrn = 0.0       !Grains
   
   !Demand (g N/plant)
   dgsh   = 0.0       !N demand for new shoot growth
   tndem  = 0.0       !Shoot total demand
   rndem  = 0.0       !Root total demand
   ndem   = 0.0       !Total demand
   nsink  = 0.0       !Grain N demand (only for new growth)
   
   !N uptake/mass in plant parts (g N/plt)
   stemn  = 0.0
   leafn  = 0.0
   grainn = 0.0
   stovn  = 0.0
   rootn  = 0.0
   
   !Mass of N in different tissues (g N/m2)
   wtnrt  = 0.0       !Roots
   wtnsd  = 0.0       !Grains
   wtnlf  = 0.0       !Leaves
   wtnst  = 0.0       !Stems
   wtnvg  = 0.0       !Stover (leaves + stems)
   cannaa = 0.0       !Stover at anthesis
   wtncan = 0.0       !Canopy (stover + grains)
   wtnup  = 0.0       !Total
   
   !Stress factors
   nstres = 1.0
   agefac = 1.0
   ndef3  = 1.0
   nfac   = 1.0 
   
   !Senescence variables
   senesce % ResWt  = 0.0
   senesce % ResLig = 0.0
   senesce % ResE   = 0.0
   CumLeafSenes     = 0.0
   CumLeafSenesY    = 0.0
   CumLfNSenes      = 0.0
   trlv             = 0.0
   rnloss           = 0.0
   trnlos           = 0.0
   rootnloss        = 0.0
   shutnloss        = 0.0
   stemnloss        = 0.0
   leafnloss        = 0.0
   leafndead        = 0.0
   sln              = 0.0   
   
   !- Initialize N uptake variables
   call MZ_AG_NUptake(control, soilprop, sw, growth, rlv, no3, nh4,       & !Input 
                      rndem, tndem, ndem, frstemn, frleafn,               & !Input
                      drootn, dstemn, dleafn, dstovn, trnu, uno3, unh4)     !Output
   
   call MZ_AG_OpNPlant(control, iswitch, yrplt, mdate, pltpop, nstres,   & !Input
                       agefac, rcnp, tcnp, rmnc, tmnc, pcnrt, pcnveg,    & !Input
                       pcnst, pcnl, pcngrn, wtnrt, wtnvg, wtnst, wtnlf,  & !Input
                       wtnsd, wtncan, rndem, tndem, ndem, sln)             !Input


!----------------------------------------------------------------------
! DYNAMIC = RATE
!----------------------------------------------------------------------
else if(dynamic==rate) then


!----------------------------------------------------------------------
! DYNAMIC = INTEGR
!----------------------------------------------------------------------
else if(dynamic==integr) then

   if(iswnit == 'N') return
   if(growthStage > 0 .AND. growthStage < 2) return    !Before emergence
   
   !- Initialize N variables at emergence
   if(yrdoy==gstdyrdoySim(2)) then
      !Critical concentrations
      rcnp = rance
      tcnp = tance
      
      !Actual concentrations
      ranc = rance
      tanc = tance
      
      !Mass of N in shoots and roots (g N/plant)
      rootn = ranc * rtwt
      stovn = tanc * stovwt
      stemn = frstemn * stovn
      leafn = frleafn * stovn
   end if
   
   !- N at anthesis and silking. Check if yrdoy==gstdyrdoySim(7) is YRDOY.EQ.STGDOY(4) in CERES-Maize
   if(yrdoy==gstdyrdoySim(6)) cannaa = stovn * pltpop
   
   !- N mobilization during grain filling period
   GRAINFILL: if(growthStage>=7 .AND. growthStage<=9) then
      !Grain N allowed to vary between .01 and .018. 
      !High temperature, low soil water, and high N increase grain N  
      tempm = 0.5 * (tmin + tmax)
      tfac = 0.690 + 0.0125*tempm
      sfac = 1.125 - 0.1250*turfac                      
      gnp = (0.004 + 0.0130*nfac)*amax1(sfac,tfac)      
      nsink = grogrn*gnp                                 
      
      !Calculate root and shoot mobilization pools. Grain N demand (nsink) is
      !preferentially satisfied using the shoot mobilization pool first.                                                                   
      if(nsink /= 0.0) then                           
         rmnc = 0.75*rcnp                              
         if(ranc < rmnc) ranc = rmnc                                 
         tanc = stovn/stovwt             !Calculate shoot concentration based on today's stover dry weight                            
         if(tanc < tmnc) tanc = tmnc                        
         npool1 = stovwt * (tanc-tmnc)                     
         npool2 = rtwt   * (ranc-rmnc)                     
         xnf = 0.15  + 0.25*nfac                      
         tnlab = xnf * npool1                         
         rnlab = xnf * npool2                         
         npool = tnlab + rnlab                          

         nsdr = npool/nsink                            
         if(nsdr < 1.0) nsink = nsink*nsdr                          
         if(nsink > tnlab) then                      
            stovn = stovn - tnlab                       
            rnout = nsink - tnlab                       
            rootn = rootn - rnout                       
         else 
            stovn = stovn - nsink                       
         end if                                                                         
      end if                                               
   end if GRAINFILL
   
   !- Calculate critical and minimum N concentrations
   tcnp = exp(ctcnp1 - ctcnp2 * xstage) / 100.0 
   if(xstage < 4.0) then
      tmnc = (1.25 - 0.200 * xstage)/100.0
   else
      tmnc = tmnc2
   end if
   
   !- Re-calculate N mass in stems and leaves due to N mobilitization out of shoots
   if(stovwt > 0.0) stemn = frstemn * stovn
   if(stovwt > 0.0) leafn = frleafn * stovn
   
   !- Calculate N concentrations based on increased dry matter (that includes today's addition)
   !  but before N uptake (dilution) in g[N]/g[DM]
   if(rtwt   > 0.0) ranc = rootn / rtwt
   if(stovwt > 0.0) tanc = stovn / stovwt
   
   !- Calculate shoot N demand for new growth (dgsh in g[N]/plt)
   if(gshplt == 0.0) gshplt = 1.0 ; dgsh = gshplt * tcnp
   if(xstage <= 1.2) dgsh = 0.0
   
   !- Calculate N demand in g[N]/plt
   rndem = rtwt   * (rcnp-ranc) + pgrort*rcnp
   tndem = stovwt * (tcnp-tanc) + dgsh
   ndem  = tndem + rndem
   
   !- Calculate daily N uptake
   call MZ_AG_NUptake(control, soilprop, sw, growth, rlv, no3, nh4,       & !Input 
                      rndem, tndem, ndem, frstemn, frleafn,               & !Input
                      drootn, dstemn, dleafn, dstovn, trnu, uno3, unh4)     !Output
   
   !- Calculate root senescence losses @ 0.5%/day
   trlv = sum( (/ (rlv(l)*dlayr(l),  l=1,nlayr) /) )
   senesce % ResWt  = 0.0
   senesce % ResLig = 0.0
   senesce % ResE   = 0.0  
   rnloss = 0.0
   trnlos = 0.0
   if(dstovn > 0.0 .AND. drootn > 0.0) then
      do l = 1, nlayr
         if(rlv(l) /= 0.0) then
            rnloss = ranc * rtwt * 0.005 * pltpop * rlv(l)*dlayr(l) / trlv
         
            !Calculate N in senesced roots (kg[N]/ha)
            senesce % ResE(l,1) = rnloss * 10.0
            !kg[N]/ha           =  g[N]/m2  * 10.
      
            !Back calculate senesced root mass from N senesced
            if(ranc > 0.0) then
               senesce % ResWt(l) = senesce % ResE(l,1) / ranc   
            else                               !kg[root DM]/ha
               senesce % ResWt(l) = senesce % ResE(l,1) * 10.0 / 0.40   
               !kg[root DM]/ha    =            kg[N]/ha * kg[C]/kg[N] / kg[C]/kg[dry matter]
            end if
      
            !Compute lignin, cellulose and carbohydrate portions
            senesce % ResLig(l) = senesce % ResWt(l) * pligrt
            
            !Total N loss by root exudation (g[N]/m2)
            trnlos = trnlos + rnloss
         end if
      end do
   end if
   
   ! Partition root N losses (g[N]/plt/day)
   rootnloss = 0.0
   shutnloss = 0.0
   if(pltpop > 0.0) then
      rootnloss = (1.0-ptf)*trnlos/(pltpop)
      shutnloss = ptf*trnlos/(pltpop)
   end if   
   stemnloss = frstemn * shutnloss
   leafnloss = frleafn * shutnloss
   
   !- Calculate N loss from leaf senescence
   leafndead = tmnc * glvsenplt
   
   !- Integrate: calculate N mass in plant parts (g[N]/plt)
   !  Adjust dstovn and drootn to compensate for N lost to FON
   grainn = grainn + nsink             !Grain N demand (nsink) was computed prior to root and shoot uptake
   rootn = max(0.0, rootn + drootn - rootnloss)
   stemn = max(0.0, stemn + dstemn - stemnloss - leafndead)
   leafn = max(0.0, leafn + dleafn - leafnloss)
   stovn = stemn + leafn

   !- Compute plant nitrogen mass variables
   wtnrt  = rootn  * pltpop           !Nitrogen in roots (g N/m2)
   wtnst  = stemn  * pltpop           !Nitrogen in stems (g N/m2)
   wtnlf  = leafn  * pltpop           !Nitrogen in leaves (g N/m2)
   wtnsd  = grainn * pltpop           !Nitrogen in grains (g N/m2)
   wtnvg  = wtnlf  + wtnst            !Nitrogen in vegetative tissue (stems+leaves, g[N]/m2)
   wtncan = wtnvg  + wtnsd            !Nitrogen in canopy (stems+leaves+grains, g N/m2)
   wtnup  = wtncan + wtnrt            !Total (roots+stems+leaves+grains, g N/m2)
   
   !- Compute specific leaf nitrogen
   sln = 0.0 ; if(xhlai > 0.0) sln = wtnlf / xhlai
   
   !- Update N concentrations in plant parts and roots after uptake (accumulation, g[N]/g[dry matter])
   if(stovwt > 0.0) tanc = stovn / stovwt
   if(rtwt > 0.0 .AND. rtwt > 0.1*rtwt) ranc = rootn / (rtwt-0.01*rtwt)
   
   !- Plant nitrogen concentration variables (%, i.e. g[N]/100g[dry matter])
   pcnrt  = 0.0 ; if(rtwt > 0.0)                            pcnrt  = wtnrt/(rtwt  *pltpop)*100.0
   pcnl   = 0.0 ; if(lfwt  > 0.0 .AND. pltpop > 0.0)        pcnl   = wtnlf/(lfwt  *pltpop)*100.0
   pcnst  = 0.0 ; if(stmwt > 0.0 .AND. pltpop > 0.0)        pcnst  = wtnst/(stmwt *pltpop)*100.0
   pcngrn = 0.0 ; if(grnwt > 0.0 .AND. pltpop > 0.0)        pcngrn = wtnsd/(grnwt *pltpop)*100.0
   pcnveg = 0.0 ; if((lfwt+stmwt) > 0.0 .AND. pltpop > 0.0) pcnveg = wtnvg/(stovwt*pltpop)*100.0
   
   !- Compute N lost in leaf senescence
   !  N in senesced leaves - calculate based on today's N content
   if(stovwt > 0.0) CumLfNSenes = CumLfNSenes + (CumLeafSenes - CumLeafSenesY) * stovn / stovwt
   
   !- Compute N stress factors --------------------------------------------------------------------------------
   ! If the actual nitrogen concentration in shoots (tanc) decreases below a critical level (tcnp), then 
   ! compute N stress based on proportion of nitrogen below critical level (tcnp-tanc) and total nitrogen 
   ! between the critical level and minimum level (tcnp-tmnc). 
   !Compute nitrogen stress factor affecting grain N concentration
   nfac = 1.0 - (tcnp-tanc)/(tcnp-tmnc)
   nfac = amin1(amax1(0.001, nfac), 1.0)
   
   !Compute nitrogen stress factor for reducing expansion
   agefac = amin1(2*nfac, 1.0)
   
   !Compute nitrogen stress factor affecting growth 
   if(nfac > 0.5) then
      nstres = 0.4*nfac + 0.6
   else   
      nstres = 1.2*nfac + 0.2
   end if   
   
   !Make grain per plant calculation less sensitive to N stress
   if(nfac < 0.8) then
      ndef3 = amin1(0.2 + nfac, 1.0)
   else
      ndef3 = 1.0
   end if   
   ! ----------------------------------------------------------------------------------------------------------


!----------------------------------------------------------------------
! DYNAMIC = OUTPUT
!----------------------------------------------------------------------
else if(dynamic==output) then

   call MZ_AG_OpNPlant(control, iswitch, yrplt, mdate, pltpop, nstres,   & !Input
                       agefac, rcnp, tcnp, rmnc, tmnc, pcnrt, pcnveg,    & !Input
                       pcnst, pcnl, pcngrn, wtnrt, wtnvg, wtnst, wtnlf,  & !Input
                       wtnsd, wtncan, rndem, tndem, ndem, sln)             !Input
   

!----------------------------------------------------------------------
! DYNAMIC = SEASEND
!----------------------------------------------------------------------
else if(dynamic==seasend) then

   !- Senesced leaves do not fall to the ground and so are added to surface litter only at harvest.
   senesce % ResWt(0)  = CumLeafSenes
   senesce % ResLig(0) = CumLeafSenes * 0.4 * pliglf
   senesce % ResE(0,1) = CumLfNSenes                   !N in senesced leaves
   
   call MZ_AG_OpNPlant(control, iswitch, yrplt, mdate, pltpop, nstres,   & !Input
                       agefac, rcnp, tcnp, rmnc, tmnc, pcnrt, pcnveg,    & !Input
                       pcnst, pcnl, pcngrn, wtnrt, wtnvg, wtnst, wtnlf,  & !Input
                       wtnsd, wtncan, rndem, tndem, ndem, sln)             !Input

!-----------------------------------------------------------------------
! END OF DYNAMIC IF STRUCTURE
!-----------------------------------------------------------------------  
end if    

!-----------------------------------------------------------------------
! END OF SUBROUTINE
!-----------------------------------------------------------------------    
return
end subroutine MZ_AG_NPlant


!--------------------------------------------------------------------------------------------------
! Define Variables
!--------------------------------------------------------------------------------------------------
! agefac      !Nitrogen stress factor affecting expansion (0-1)
! cannaa      !Stover N at anthesis, g N/m2
! ctcnp1      !Maximum value for critical tissue N concentration (in developing seed embryo), g N/100g[DM]
! ctcnp2      !Coefficent for change in N concentration with growth stage, g N/100g[DM]/xstage unit
! CumLeafSenes!Cumulative leaf senescence, kg[leaf DM]/ha
! dng         !N demand of potential new growth of tops (g N/plant)
! gnp         !Nitrogen concentration in new grain growth, gN/g dry matter
! grainn      !Grain nitrogen content, g N/plant
! grnwt       !Grain weight, g/plant
! grogrn      !Daily growth of the grain - g/plt/day
! lfwt        !Leaf weight, g/plant
! ndef3       !Nitrogen stress factor affecting grains per plant (0-1)
! nfac        !Nitrogen stress factor based on actual and critical nitrogen content in vegetative tisue
! npool       !Total plant N available for translocation to grain (g/plant)
! npool1      !Tops N available for translocation to grain (g/plant)
! npool2      !Root N available for translocation to grain (g/plant)
! nsdr        !Plant N supply/demand ratio used to modify grain N content
! nsink       !Demand for N associated with grain filling (g[N]/plant/day)
! nstres      !Nitrogen stress factor affecting growth (0-1)
! pdwi        !Potential increment in new shoot growth, g/plant
! pcngrn      !Percent nitrogen in grain, %
! pcnl        !Percent nitrogen in leaves, %
! pcnrt       !Percent nitrogen in roots, %
! pcnst       !Percent nitrogen in stems, %
! pcnveg      !Percent nitrogen in vegetative tissue (leaf and stem), %
! pgrort      !Potential increment in new root growth, g/plant
! pliglf      !Lignin fraction in leaves
! pligrt      !Lignin fraction in roots
! pltpop      !Plant population, plants/m2
! ptf         !Ratio of above ground biomass to total biomass
! ranc        !Root actual N concentration, g N/g root
! rance       !Root nitrogen concentration at emergence, g N/g root dry weight
! rcnp        !Root critical nitrogen concentration, g N/g root dry weight
! rmnc        !Root minimum nitrogen concentration (g N/g root dry weight)
! rndem       !Plant root demand for nitrogen (g/plant)
! rnlab       !
! rnloss      !Loss of N from the plant via root exudation in one layer (g N/m2)
! rnout       !
! rootn       !Root nitrogen content, g N/plant
! rtwt        !Root weight, g/plant
! sfac        !Drought stress factor for grain nitrogen concentration
! stmwt       !Stem weight, g/plant
! stovn       !Nitrogen content in stover, g N/plant
! stovwt      !Stover weight (Stem + leaf), g/plant
! tanc        !Nitrogen content in above ground biomass, decimal
! tance       !Nitrogen content in above ground biomass at emergence, g N/g dry weight
! tcnp        !Critical nitrogen concentration in tops, g N/g dry wt.
! tempm       !Mean daily temperature, C
! tfac        !Temperature stress factor for grain nitrogen concentration
! tmnc        !Plant top minimum N concentration g N/g dry matter
! tmnc2       !Plant top minimum N concentration g N/g dry matter
! tndem       !Plant tops demand for nitrogen (g N/plant)
! tnlab       !
! trnlos      !Total plant N lost by root exudation (g N/m2)
! trnu        !Total potential root nitrogen uptake, kg N/ha
! unh4(20)    !Plant uptake of ammonium from layer (kg N/ha/day)
! uno3(20)    !Plant uptake of nitrate from a layer (kg N/ha/day)
! wtncan      !Weight of nitrogen in above ground biomass (stem, leaf, grain), g N/m2
! wtnlf       !Weight of nitrogen in leaf tissue, g N/m2
! wtnsd       !Weight of nitrogen in seed, g[N] / m2[ground]
! wtnst       !Weight of nitrogen in stem tissue, g N/m2
! wtnup       !Total N uptake, g/m2
! wtnvg       !Weight of nitrogen in vegetative tissue, g N/m2
! xnf         !Modified nitrogen factor based on critical N concentration in vegetative biomass
! xstage      !Non-integer growth stage indicator
!--------------------------------------------------------------------------------------------------


