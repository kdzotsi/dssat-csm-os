!======================================================================
!  MZ_AG_NUptake, Subroutine
!  Determines N uptake
!----------------------------------------------------------------------
!  Revision history
!                 Written
!  02/08/1993 PWW Header revision and minor changes    
!  06/01/1994 WTB Modified 
!  06/28/1994 JR/BDB Changed water content dependent factor 
!  03/29/2001 WBD Converted to modular format      
!  12/01/2001 WDB Major revision for 2002 release   
!  03/12/2003 CHP Changed senescence variable to composite (SENESCE) as defined in ModuleDefs.for
!  03/02/2005 CHP/JIL Add check for negative RNH4U & RNO3U.
!  11/07/2005 CHP Replaced FAC with SOILPROP variable KG2PPM
!  01/31/2006 CHP/JIL Fixed computation of RNLOSS to match root loss calculation in MZ_GROSUB (from 5% to .5%)
!  03/30/2006 CHP Added composition of senesced matter for SOM modules
!  07/10/2014 KAD/TT/SK Adapted for AgMaize
!----------------------------------------------------------------------
!  Called :  
!  Calls  : None
!  TODO:
!  - Move N in senesced material to NPlant
!  - Part of the N loss due to root senescence substracted from shoot N?
!  - Check for negative shoot demand at emergence (not just total demand)
!----------------------------------------------------------------------

subroutine MZ_AG_NUptake(control, soilprop, sw, growth, rlv, no3, nh4,       & !Input 
                         rndem, tndem, ndem, frstemn, frleafn,               & !Input
                         drootn, dstemn, dleafn, dstovn, trnu, uno3, unh4)     !Output


use ModuleDefs
use MZ_AG_ModuleDefs
implicit  none
save

integer :: dynamic, l, l1, nlayr, yrdoy
real :: rndem, tndem, ndem, andem, xrndem, xtndem, xndem, drootn, dstovn, dstemn, dleafn, frstemn, frleafn
real :: factor, fnh4, fno3, nuf, pltpop, rfac, smdfr, trnu, xmin, stmwt, lfwt, stovwt
real, dimension(nl) :: dlayr, kg2ppm, ll, sat, sw, nh4, no3, rlv, rnh4u, rno3u, shf, snh4, sno3, unh4, uno3

! Constructed variables based on definitions in modules
type(ControlType) control
type(SoilType)    soilprop
type(FileioType)  datafileio
type(GrowthType)  growth

! Transfer values from constructed data types into local variables
dynamic = control % dynamic
yrdoy   = control % yrdoy
dlayr   = soilprop % dlayr  
nlayr   = soilprop % nlayr  
kg2ppm  = soilprop % kg2ppm    
ll      = soilprop % ll     
sat     = soilprop % sat    
shf     = soilprop % wr
stmwt   = growth   % wstplt
lfwt    = growth   % wlvplt
stovwt  = growth   % stoverplt

!----------------------------------------------------------------------
! DYNAMIC = RUNINIT OR SEASINIT
!----------------------------------------------------------------------
if(dynamic==runinit .OR. dynamic==seasinit) then

   !Read all sections of fileio and transfer variables
   call readfileio(control, 'ALLSEC', datafileio)
   pltpop = datafileio % pltpop 
   
   !Shoot and root demand 
   andem  = 0.0
   xndem  = 0.0
   xrndem = 0.0
   xtndem = 0.0
   
   !Supply variables
   sno3   = 0.0
   snh4   = 0.0
   
   !Uptake variables (kg N/ha)
   xmin   = 0.0
   unh4   = 0.0
   uno3   = 0.0
   rno3u  = 0.0
   rnh4u  = 0.0
   
   !Mass of N in shoots and roots (g N/plant)
   drootn = 0.0
   dstovn = 0.0
   dstemn = 0.0
   dleafn = 0.0
   
   !Other variables
   nuf    = 0.0
   factor = 0.0
   trnu   = 0.0
   smdfr  = 0.0
   l1     = 0
 
!----------------------------------------------------------------------
! DYNAMIC = RATE
!----------------------------------------------------------------------
else if(dynamic==rate) then


!----------------------------------------------------------------------
! DYNAMIC = INTEGR
!----------------------------------------------------------------------
else if(dynamic==integr) then

   !- Calculate potential N supply in soil layers with roots (trnu)
   !  Potential ammonia and nitrate uptake from layer l in kg[N]/ha (rnh4u and rno3u)
   l1    = 0
   rno3u = 0.0
   rnh4u = 0.0
   trnu  = 0.0
   do l = 1, nlayr
      if(rlv(l) /= 0.0) then
         l1 = l
         !The following code was written by JTR to simplify the code for the
         !generic model and to make the functions more like the water uptake
         !functions.  Done on June 28, 1994.
         !New water content dependent factor for uptake and new uptake
         !methods that work like the soil water uptake  -- JTR 6/94
         smdfr    = 1.5 - 6.0*((sw(l)-ll(l))/(sat(l)-ll(l))-0.5)**2
         smdfr    = amax1(smdfr, 0.0)
         smdfr    = amin1(smdfr, 1.0)
         rfac     = 1.0 - exp(-8.0*rlv(l))
         fnh4     = shf(l) * 0.075
         fno3     = shf(l) * 0.075
         rnh4u(l) = rfac*smdfr*fnh4*(nh4(l)-0.5)*dlayr(l)
         rno3u(l) = rfac*smdfr*fno3*no3(l)*dlayr(l)
         rnh4u(l) = max(rnh4u(l), 0.0)
         rno3u(l) = max(rno3u(l), 0.0) 
         trnu     = trnu + rno3u(l) + rnh4u(l)
      end if
   end do
   
   !- Calculate factor (nuf) to reduce N uptake to the level of demand
   andem = ndem * pltpop * 10.0       !Total N demand in kg[N]/ha
   if(andem <= 0.0) then
      trnu = 0.0
      nuf  = 0.0
   else
      if(trnu == 0.0) return
      if(andem < trnu) then
         nuf = andem/trnu
      else
         nuf = 1.0
      end if   
      trnu  = 0.0
   end if
   
   !- Calculate N uptake in soil layers with roots based on demand (kg/ha)
   !  Uptake cannot exceed supply minus xmin, i.e. 0.0 <= unh4(l) <= sno3 - xmin
   sno3(1:nlayr) = (/ (no3(l)/kg2ppm(l), l=1,nlayr) /)     !Convert N supply to kg[N]/ha
   snh4(1:nlayr) = (/ (nh4(l)/kg2ppm(l), l=1,nlayr) /)
   unh4 = 0.0
   uno3 = 0.0
   do l = 1, l1
      uno3(l) = rno3u(l)*nuf
      unh4(l) = rnh4u(l)*nuf
      xmin    = 0.25/kg2ppm(l)
      uno3(l) = min(uno3(l), sno3(l) - xmin)
      uno3(l) = max(uno3(l), 0.0)  
      xmin    = 0.5/kg2ppm(l)
      unh4(l) = min(unh4(l), snh4(l) - xmin)
      unh4(l) = max(unh4(l), 0.0) 
      trnu    = trnu + uno3(l) + unh4(l)     !kg[N]/ha
   end do
   
   !- Total potential N uptake in g[N]/plt
   if(pltpop > 0.0) then
     trnu = trnu/(pltpop*10.0)               
   else
     trnu = 0.0
   end if
   
   !- Reduce nitrogen demand to the level of total potential uptake
   if(ndem > trnu) then
      xndem  = trnu
      factor = xndem / ndem
      xrndem = rndem * factor
      xtndem = tndem * factor
   else
      xndem  = ndem
      xrndem = rndem   
      xtndem = tndem
   end if
   
   !- N uptake
   drootn = 0.0
   dstovn = 0.0
   dstemn = 0.0
   dleafn = 0.0
   if(xndem <= 0.0 .OR. trnu <= 0.0) then
      dstemn = 0.0
      dleafn = 0.0
      dstovn = 0.0
      drootn = 0.0
   else
      !Adjust dstovn and drootn to compensate for N lost to FON
      if(xndem > 0.0 .AND. pltpop > 0.0 .AND. xtndem > 0.0) then
         !JIL 06/25/2007
         dstovn = max(0.0, (xtndem/xndem)*trnu)
         drootn = max(0.0, (xrndem/xndem)*trnu)
         
         !Partition stover N uptake to stems and leaves 
         dstemn = frstemn * dstovn
         dleafn = frleafn * dstovn
      end if
   end if


!-----------------------------------------------------------------------
! END OF DYNAMIC IF STRUCTURE
!-----------------------------------------------------------------------  
end if    

!-----------------------------------------------------------------------
! END OF SUBROUTINE
!-----------------------------------------------------------------------    
return
end subroutine MZ_AG_NUptake


!--------------------------------------------------------------------------------------------------
! Define Variables
!--------------------------------------------------------------------------------------------------
! ANDEM       !Crop N demand (kg N/ha)
! ANO3        !Total extractable nitrate N in soil profile (kg N/ha)
! ANH4        !Total extractable ammonium N in soil profile (kg N/ha)
! DLAYR(L)    !Soil thickness in layer L (cm)
! DLEAFN      !Daily change in plant stover nitrogen content (g N/plant)
! DROOTN      !Daily change in plant root nitrogen content (g N/plant)
! DSTEMN      !Daily change in plant stover nitrogen content (g N/plant)
! DSTOVN      !Daily change in plant stover nitrogen content (g N/plant)
! FACTOR      !Ratio of root N uptake to plant N demand
! FNH4        !Unitless soil ammonium supply index
! FNO3        !Unitless soil nitrate supply index
! FON(20)     !Fresh organic nitrogen in soil layer L due to root senescence, kg N/ha
! KG2PPM(20)  !Factor that convergs mg elemental N/kg soil to kg N/ha for soil layer L
! L           !Index counter
! L1          !Lowest soil layer with roots
! LL(20)      !Lower limit of plant extractable water for soil layer L, cm3/cm3
! NDEM        !Plant nitrogen demand, g/plant
! NH4(20)     !Ammonium in soil layer L, mg elemental N/kg soil
! NLAYR       !Number of soil layer
! NO3(20)     !Nitrate in soil layer L (mg elemental N/kg soil)
! NUF         !Plant N supply to demand used to modify N uptake
! PLTPOP      !Plant population, plants/m2
! RFAC        !Interim variable describing the effects of root length density on potential N uptake from a layer
! RLV(20)     !Root length density for soil layer L, cm root/cm3 soil
! RNDEM       !Plant root demand for nitrogen (g/plant)
! RNH4U(20)   !Potential ammonia uptake from layer L, kg N/ha
! RNO3U(20)   !Potential nitrate uptake from layer L, kg N/ha
! ROOTN       !Root nitrogen content, g N/plant
! SHF(20)     !Relative root distribution in soil layer L (0-1)
! SMDFR       !Soil moisture deficit factor affecting N uptake
! SNH4(20)    !Ammonium nitrogen in soil layer L, kg N/ha
! SNO3(20)    !Nitrate content in soil layer L, kg N/ha
! STOVN       !Nitrogen content in stover, g N/plant
! SAT(20)     !Saturated water holding capacity for soil layer L, cm3/cm3
! SW(20)      !Soil water content in layer L, cm3 water/cm3 soil
! TANC        !Nitrogen content in above ground biomass, decimal
! TCNP        !Critical nitrogen concentration in tops, g N/g dry weight
! TNDEM       !Plant tops demand for nitrogen (g N/plant)
! TRLV        !Total root length density, cm root/cm2 soil
! TRNU        !Total potential root nitrogen uptake, g N/plt
! UNH4(20)    !Plant uptake of ammonium from layer (kg N/ha/day)
! UNO3(20)    !Plant uptake of nitrate from a layer (kg N/ha/day)
! XMIN        !
! XNDEM       !Re-scaled total demand limited by potential N uptake (g N/plant)
! XRNDEM      !Re-computed root N demand using re-scaled total demand (g N/plant)
! XTNDEM      !Re-computed canopy N demand using re-scaled total demand (g N/plant)




