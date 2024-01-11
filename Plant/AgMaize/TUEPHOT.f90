!=======================================================================
! Note: Should be able to use today's phenology and leaf area?
!=======================================================================
subroutine TUEPHOT(control, soilprop, weather, sw, & !Input
                   eop)                                !Output
!-----------------------------------------------------------------------
   use ModuleDefs
   use ModuleData
   implicit none
   save
!-----------------------------------------------------------------------
   integer, parameter :: lenla = 50
   real :: sw(nl), cgr, cgrf, tmax, tmin, esatTmax, esatTmin, vpdday, wue, rhmin, vpdmax, gddae
   real :: dvs, srad, parday, parsunday, parshday, radtotday, wlvtop, wlvbot, wst, weolg, wrt, cgrw
   real :: pgross, asg, maint, wrequired, eop
   real, dimension(lenla) :: greenla, lflon, lamaxpot, lapot
   integer :: dynamic, lfn, growthStage

! Constructed types
   type(ControlType) :: control
   type(SoilType)   :: soilprop
   type(WeatherType) :: weather
   dynamic = control%dynamic
   tmax = weather%tmax
   tmin = weather%tmin
   rhmin = weather%rhum

   if (dynamic == runinit .or. dynamic == seasinit) then
      wrequired = 0.0
      eop = 0.0
      gddae = 0.0
      dvs = 0.0
      cgrf = 0.0
      lfn = 0
      growthStage = 0
      greenla = 0.0
      lflon = 0.0
      lamaxpot = 0.0
      lapot = 0.0
      wlvtop = 0.0
      wlvbot = 0.0
      wst = 0.0
      weolg = 0.0
      wrt = 0.0
      cgrw = 0.0

      !Call AgMaize's photosynthesis and respiration
      call MZ_TUE_Photsynt(control, soilprop, weather, sw, & !Input
                           gddae, dvs, lfn, greenla, lflon, lamaxpot, lapot, & !Input
                           pgross, srad, parday, parsunday, parshday, radtotday)     !Output

      call MZ_TUE_Respiration(control, weather, & !Input
                              dvs, wlvtop, wlvbot, wst, weolg, wrt, cgrw, & !Input
                              maint)                                                    !Output

   else if (dynamic == rate) then

      call get('PLANT', 'GSTD', growthStage)

      ! Wait for crop model to start (i.e. germination) before retrieving variables from data storage.
      ! This avoids getting values from the previous run into a subsequent run.
      if (growthStage >= 1 .AND. growthStage < 10) then

         !Obtain values of variables from data storage
         call get('PLANT', 'GDDAE', gddae)
         call get('PLANT', 'DVS', dvs)
         call get('PLANT', 'CGRF', cgrf)
         call get('PLANT', 'LFN', lfn)
         call get('PLANT', 'GREENLA', greenla)
         call get('PLANT', 'LFLON', lflon)
         call get('PLANT', 'LAMAXPOT', lamaxpot)
         call get('PLANT', 'LAPOT', lapot)
         call get('PLANT', 'WLVTOP', wlvtop)
         call get('PLANT', 'WLVBOT', wlvbot)
         call get('PLANT', 'WST', wst)
         call get('PLANT', 'WEOLG', weolg)
         call get('PLANT', 'WRT', wrt)
         call get('PLANT', 'CGRW', cgrw)
         !-----------------------------------------------------------------------

         !Call AgMaize's photosynthesis
         call MZ_TUE_Photsynt(control, soilprop, weather, sw, & !Input
                              gddae, dvs, lfn, greenla, lflon, lamaxpot, lapot, & !Input
                              pgross, srad, parday, parsunday, parshday, radtotday)     !Output

         call MZ_TUE_Respiration(control, weather, & !Input
                                 dvs, wlvtop, wlvbot, wst, weolg, wrt, cgrw, & !Input
                                 maint)                                                    !Output

         !Expected crop growth rate given today's asg
         asg = pgross - maint
         cgr = cgrf*asg

         !Calculate daytime vapor pressure deficit (vpdday). If relative humidity is available (rhmin) use it
         !to calculate maximum vpd then derive vpdday (Stockle's email, Aug 2013); the rhum input in DSSAT is
         !to be interpreted as the rhum recorded at tmax. If rh is not available (i.e. -99, missing value), use
         !saturation vapor pressure at tmax and tmin to calculate vpdday (Basso and Ritchie, 2014).
         !Equations for esat and vpdmax from http://agsys.cra-cin.it/tools/evapotranspiration/help/Maximum_vapour_pressure_deficit.html
         esatTmax = 0.6108*exp(17.27*tmax/(tmax + 237.3))
         if (rhmin > 0.0) then    !Use Stockle's method
            vpdmax = esatTmax*(1 - rhmin/100.)
            vpdday = 0.67*vpdmax
         else                    !rhmin missing, use Basso and Ritchie (2014). Nature Clim. Chg 4.
            esatTmin = 0.6108*exp(17.27*tmin/(tmin + 237.3))
            vpdday = 0.75*(esatTmax - esatTmin)
         end if

         !Calculate water use efficiency (wue), water requirement (wrequired) and actual transpiration (ep)
         !Conversion: 1 kg[water]/m2[ground] = 0.001 m3[water]/m2[ground] = 1 mm water height
         wue = 7.44*vpdday**(-0.42)    !g[dry matter]/kg[water], Kremer et al. 2008 (book chapter)
         wrequired = cgr*0.1/wue         !Convert dry matter from kg/ha to g/m2. Get wrequired in kg[water]/m2/day or mm/day
         eop = wrequired
      end if
   end if

   return
end subroutine TUEPHOT

!==============================================================================================================================
! Variable definitions                                                                                  Unit
!------------------------------------------------------------------------------------------------------------------------------
! esatTmax  Saturation vapor pressure at maximum temperature                                            kPa
! rhmin     Minimum relative humidity                                                                   %
! cgr       Crop growth rate (aboveground dry matter)                                                   kg[dm]/ha/day
! tmax      Maximum daily temperature                                                                   degree C
! tmin      Minimum daily temperature                                                                   degree C
! vpdmax    Maximum vapor pressure deficit                                                              kPa
! vpdday    Daytime vapor pressure deficit                                                              kPa
! weather   Constructed type for weather variables
! wrequired Water requirement                                                                           mm/day
! wue       Water use efficiency                                                                        g[dm]/kg[water]
!==============================================================================================================================

!=======================================================================
!  MZ_TUE_Photsynt
!  Computes canopy daily gross photosynthesis (g biomass/m2 d) by
!  integrating each leaf contribution over hourly time steps in the day.
!  Adapted from ETPHOT by N.B. Pickering
!-----------------------------------------------------------------------
!  REVISION HISTORY
!  02/07/2003 JIL Written
!  09/12/2007 JIL Modified and adapted for IXIM model
!  03/08/2013 TT/SK/KAD Adapted for AgMaize, respiration part removed
!-----------------------------------------------------------------------
!  Questions for Jon:
!  -What is cans (= 4.5) in the equation for canopy height? Equation here different from the one in paper
!  -Canopy width cannot be larger than row spacing?
!
!  Questions for CHP
!  -azir is azimuth relative to North while azzon is relative to South?
! TODO:
! -Error message for rowsp or pltpop = 0.0
!=======================================================================
subroutine MZ_TUE_Photsynt(control, soilprop, weather, sw, &
                           gddae, dvs, &   !From Phenology
                           lfn, greenla, lflon, lamaxpot, lapot, &   !From Leaf area
                           pgross, srad, parday, parsunday, parshday, radtotday) !Outputs

   use ModuleDefs
   use MZ_AG_ModuleDefs
   implicit none
   save

!----------------------------------------------------------------------
! Define variables
!----------------------------------------------------------------------
   integer, parameter :: lenla = 50
   real, parameter :: cans = 4.5
   real, parameter :: scvp = 0.2
   real, parameter :: tincr = 24.0/ts
   real, parameter :: insp = 4.514
   real, parameter :: as = 288.615
   real, parameter :: cvx = -75.774
   real, parameter :: pang = 64.0
   real, parameter :: cm2tom2 = 1.0E-4
   integer, intent(in) :: lfn
   real, intent(in)  :: sw(nl), gddae, dvs
   real, intent(in), dimension(lenla)  :: lflon, lamaxpot
   real, intent(out) :: pgross

! Constructed variables based on definitions in ModuleDefs.for
   type(ControlType), intent(in) :: control
   type(WeatherType), intent(in) :: weather
   type(SoilType), intent(in) :: soilprop
   type(FileioType)  :: datafileio
   type(SpeciesType) :: dataspecies

! Variables from constructed types
   integer :: dynamic, yrdoy
   real, dimension(ts) :: amtrh, tairhr, azzon, beta, frdifp, parhr
   real :: salb, dul(nl), snup, sndn, azir, pltpop, rowspcm, rowspmeter, sgfun, gasfn, tmin
   real :: asmax, canh, xc    !Photosynthesis parameters from species file

! KAD 09/24/2013 - For checking radiation subroutine
   real, dimension(ts) :: parsunHour, parshHour, radtotHour
   real :: srad, parday, parsunday, parshday, radtotday

! Other variables
   character(len=12) :: files
   character(len=80) :: pathsr
   integer :: i, h, leafnum(lenla)
   logical :: light, daytim
   real :: lfln, lflp, canw, canwh, hs, turfac, mult, amplt, froll, palb, palbd
   real :: pg, pghr, betn, canht, hlai, pgday, radtot
   real, dimension(lenla) :: plaisl, plaish, parsh, parsun, lapot, greenla
   real, dimension(lenla) :: arefhr, adifhr, addrhr, addfhr

! Transfer values from constructed data types into local variables
   dynamic = control%dynamic
   yrdoy = control%yrdoy
   azzon = weather%azzon
   beta = weather%beta
   amtrh = weather%amtrh
   tairhr = weather%tairhr
   sndn = weather%sndn
   snup = weather%snup
   frdifp = weather%frdifp
   parhr = weather%parhr
   tmin = weather%tmin
   salb = soilprop%salb
   dul = soilprop%dul

!----------------------------------------------------------------------
! Dynamic = runinit or dynamic = seasinit
!---------------------------------------------------------------------
   if (dynamic == runinit .or. dynamic == seasinit) then
!open(unit=9000, file="PHOTAGMAIZE.OUT")
!Initialize variables
      light = .true.
      pghr = 0.0
      betn = 0.0
      canht = 0.0
      hlai = 0.0
      lfln = 0.0
      lflp = 0.0
      canw = 0.0
      canwh = 0.0
      daytim = .false.
      hs = 0.0
      turfac = 1.0
      mult = 0.0
      amplt = 0.0
      froll = 0.0
      palb = 0.0
      palbd = 0.0
      pg = 0.0
      pgross = 0.0
      srad = 0.
      parday = 0.
      parsunday = 0.
      parshday = 0.
      radtotday = 0.

!Read all sections of fileio and transfer variables
      call readfileio(control, 'ALLSEC', datafileio)
      files = datafileio%files
      pathsr = datafileio%pathsr
      sgfun = datafileio%sgfun
      gasfn = datafileio%gasfn
      pltpop = datafileio%pltpop
      rowspcm = datafileio%rowspc
      rowspmeter = rowspcm/100.0
      azir = datafileio%azir

!Read photosynthesis parameters from species file and transfer variables
      call readspecies(files, pathsr, '*PHOTO', dataspecies)
      asmax = dataspecies%asmax
      xc = dataspecies%xc
      canh = dataspecies%canh

!CHP 2/23/2009 added array index for tairhr -- I know this is wrong, but
!what does it need to be?  Call subroutine in a loop from 1 to 24?
      call MZ_TUE_Iphotsynt(dynamic, tmin, asmax, gddae, greenla, plaisl, plaish, & !Input
                            lapot, lflon, lfn, light, parsh, parsun, tairhr(1), lamaxpot, sgfun, gasfn, & !Input
                            pghr)                                                                            !Output

!----------------------------------------------------------------------
! Dynamic = rate
!----------------------------------------------------------------------
   else if (dynamic == rate) then

! Calculate between plant spacing (m)
      if (rowspmeter > 0.0 .and. pltpop > 0.0) then
         betn = 1.0/(rowspmeter*pltpop)
      else
         betn = 0.0
      end if

!Calculate canopy growth and update LAI using updated green leaf area
      if (dvs <= 1.0) then
         canht = 1.85*canh/(1 + EXP(-cans*(dvs - 0.95)))
      end if

      hlai = 0.0
      do i = 1, lfn
         if (greenla(i) < 0.0) greenla(i) = 0.0
         hlai = hlai + cm2tom2*greenla(i)*pltpop
         if (lapot(i) >= 0.8*lamaxpot(i)) then
            lfln = (((insp*greenla(i) + as) - ((insp*greenla(i) + as)**2.0 - (4.0*insp*greenla(i)*as*cvx))**0.5)/(2.0*cvx))
            lflp = lfln*cos(pang*rad)
            canw = lflp*2.0/100.0
            canwh = min(max(canwh, canw), rowspmeter)
         end if
      end do
      if (canht > 0.0 .and. canht < 0.01) canht = 0.01
      if (canht > 0.01 .and. canwh < 0.01) canwh = 0.01

! Calculate plant albedo to PAR as a function of surface SW
      palb = 0.6*salb
      if (sw(1) < dul(1)) then
         palbd = palb*1.25
         palb = palbd - (palbd - palb)/dul(1)*sw(1)
      end if

!***Begin hourly loop
      pgday = 0.0
      light = .true.

!open(unit=9000, file="RADABS.OUT")
!write(9000,'(9(1X,A6))') 'YRDOY', 'HOUR', 'LEAFNO', 'LAISL', 'LAISH', 'ADIF', 'ADDR', 'ADDF', 'AREF'
      do h = 1, ts
         !Calculate effect of leaf rolling
         !mult  = 5.0 - 4.0*turfac
         !amplt = (25.0+(100.0-25.0)*exp(-3.5*(1.0-amtrh(h))))*mult
         !froll = amin1(turfac+(real(h)-14.0)**2.0/amplt,1.0)
         froll = 1.0

         !Calculate real and solar time
         hs = real(h)*tincr
         if (hs > snup .and. hs < sndn) then
            daytim = .true.
         else
            daytim = .false.
         end if

         !Calculate hourly radiation absorption by canopy/soil
         call MZ_TUE_Radabs(azir, azzon(h), beta(h), betn, canht, canwh, daytim, frdifp(h), &  !Input
                            froll, greenla, h, lfn, palb, parhr(h), pltpop, rowspmeter, scvp, hlai, xc, &  !Input
                            parsh, parsun, plaish, plaisl, radtot, arefhr, adifhr, addrhr, addfhr, leafnum)     !Output

!   !write(9000, '(i7,1x,2(1x,i2),4(1x,f5.0))') (yrdoy, h, leafnum(i), adifhr(i), addrhr(i), addfhr(i), arefhr(i), i=1,25)
!   (write(9000, '(i7,1x,2(1x,i2),4(1x,f5.0))') yrdoy, h, leafnum(i), adifhr(i), addrhr(i), addfhr(i), arefhr(i), i=22,1,-1)
!   do i = 22,1,-1
!      write(9000, '(i7,2(1x,i6),2(1x,f6.3),4(1x,f6.1))') yrdoy, h, leafnum(i), plaisl(i), plaish(i), adifhr(i), addrhr(i), addfhr(i), arefhr(i)
!   end do

         !KAD - Save absorbed radiation values for each hour
         parsunHour(h) = sum(parsun*plaisl)*3600.
         parshHour(h) = sum(parsh*plaish)*3600.
         radtotHour(h) = radtot*3600.

         !Calculate instantaneous gross assimilation
         call MZ_TUE_Iphotsynt(dynamic, tmin, asmax, gddae, greenla, plaisl, plaish, & !Input
                               lapot, lflon, lfn, light, parsh, parsun, tairhr(h), lamaxpot, sgfun, gasfn, & !Input
                               pghr)                                                                           !Output

         !Integrate instantaneous canopy photoynthesis (�mol[CO2]/m2/s) to get daily values (g[CO2]/m2/day)
         !1 mol[CO2] = 44 g[CO2] so 1umol[CO2] = 44.10-6g[CO2] and 1umol[CO2]/s = 44.10-6 x 3600 g[CO2]/hr
         pgday = pgday + tincr*pghr*44.0*0.0036
      end do
!close(9000)
!***End hourly loop

! Daily gross photosynthesis (CH2O): g[CH2O]/m2/day
!This is actually glucose produced: 6 moles[CO2] (44 g/mole) would produce 1 mole[glucose] (180 g/mole)
!also equivalent to 6 moles[CH2O] (30 g/mole)
      pg = pgday*30.0/44.0
      pgross = pg*10.0           !Convert pg to kg[CH2O]/ha/day same as kg[glucose]/ha/day

! write(9000,'(f6.2,1x,4(1x,f6.1),1x,i6,1x,f6.2,1x,f6.0,2(1x,f6.1))')  &
!     hlai,gddae,sum(greenla),lapot(10),lflon(10),lfn,dvs,lamaxpot(10),pgross

! KAD 09/24/2013. Incoming and absorbed radiation
      srad = weather%srad
      parday = sum(parhr*3600.)*1E-6/4.6
      parsunday = sum(parsunHour)*1E-6/4.6
      parshday = sum(parshHour)*1E-6/4.6
      radtotday = sum(radtotHour)*1E-6/4.6

!-----------------------------------------------------------------------
! End of dynamic if structure
!-----------------------------------------------------------------------
   end if

!----------------------------------------------------------------------
! End of subroutine
!----------------------------------------------------------------------
   return
end subroutine MZ_TUE_Photsynt

!==============================================================================================================================
! Variable definitions                                                                                  Unit
!------------------------------------------------------------------------------------------------------------------------------
! amtrh        Hourly atmospheric transmission coefficient or ratio of solar:extraterrestrial radiation -
! as           Asymptote of the leaf length vs. leaf area relationship                                  cm
! asmax        Maximum instantaneous assimilation at 30 degC                                            �mol[CO2]/m2/s
! azir         Row azimuth relative to North                                                            degrees
! azzon(h)     Hourly solar azimuth (+/- from South)                                                    degrees
! beta(h)      Hourly solar elevation (+/- from horizontal)                                             degrees
! betn         Spacing between plants along a row                                                       m/plant
! canh         Potential canopy height                                                                  m
! canht        Canopy height                                                                            m
! canw         Canopy width (general calculation)                                                       m
! canwh        Canopy width with boundaries defined (e.g. cannot exceed row spacing)                    m
! cm2tom2      Constant for converting leaf area from cm2 to m2                                         m2/cm2
! cvx          Curvature of the leaf length vs. leaf area relationship                                  -
! datafileio   Constructed variable containing all variables read from the DSSAT45.INP                  -
! dataspecies  Constructed variable containing all variables read from the species file                 -
! daytim       Logical variable that describes daytime or nighttime conditions                          -
! dvs          Development stage (0 is planting, 1 is silking, and 2 is maturity)                       Unitless
! dynamic      Main control variable to tell each module which section of code to run                   -
! files        Species file name                                                                        -
! froll        Leaf rolling factor associated with soil water stress affecting cell expansion           -
! frdifp(h)    Hourly fraction diffuse photon flux density after correcting                             -
!              for circumsolar radiation (Spitters, 1986)
! gddae        Cumulative growing degree day after emergence                                            degree-days
! greenla(l)   Green leaf area for leaf l                                                               cm2
! h            Hourly iteration counter (1-24)                                                          -
! hlai         "Healthy" leaf area index (excluding senesced parts)                                     m2/m2
! hs           Hourly counter                                                                           hour
! i            Iteration counter of leaf number                                                              -
! insp         Initial slope of the leaf length vs. leaf area relationship                              cm/cm2
! integr       Program control variable to execute code to integrate daily rate variables (value=4)
! lamaxpot(l)  Maximum value of potential leaf area for leaf l (when fully expanded)                    cm2
! lapot(l)     Potential leaf area of leaf l from tip appearance to full expansion                      cm2
! lenla        Dimension (number of elements) of leaf area vector variables                             -
! lfln         Leaf length                                                                              cm
! lflon(l)     Longevity of leaf l                                                                      degree-days
! lfn          Total leaf number rounded up                                                             -
! lflp         Projection of leaves in the horizontal plane                                             cm
! light        Logical variable that differentiates daytime from nighttime
! palb         Plant albedo accounting for soil water in the first soil layer                           -
! palbd        Intermediary plant albedo calculation                                                    -
! pang         Average angle of leaves with the horizontal                                              degrees
! parhr(h)     Hourly photosynthetically active radiation (PAR)                                         �mol[quanta]/m2/s
! parsh(l)     Photosynthetically active radiation absorbed by leaf l in the shaded zone                �mol[quanta]/m2/s
! parsun(l)    Photosynthetically active radiation absorbed by leaf l in the sunlit zone                �mol[quanta]/m2/s
! pathsr       Directory containing (or path to ) the species file                                      -
! pg           Daily canopy gross assimilation (glucose equivalent)                                     g[glucose]/m2/day
! pgday        Daily canopy gross assimilation (CO2 uptake)                                             g[CO2]/m2/day
! pghr         Canopy instantaneous gross assimilation                                                  �mol[CO2]/m2/s
! pgross       Daily canopy gross assimilation (glucose equivalent in kg/ha/day)                        kg[glucose]/ha/day
! pi           Mathematical constant pi                                                                 -
! plaish(l)    Shaded leaf area index for leaf l                                                        m2/m2
! plaisl(l)    Shaded leaf area index for leaf l                                                        m2/m2
! pltpop       Plant density                                                                            plants/m2
! rad          Constant for converting angles from degrees to radians                                   radians/degree
! radtot       Sum of PAR components (energy balance check: should equal parhr)                         �mol[quanta]/m2/s
! rate         Program control variable to execute code to compute daily rate variables (value=3)
! rowspmeter   Row spacing in meter                                                                     m
! salb         Soil albedo                                                                              -
! scvp         Scattering coefficient used to calculate diffuse reflectance
! sndn         Time of sunset                                                                           hour
! snup         Time of sunrise                                                                          hour
! tairhr(h)    Hourly air temperature                                                                   degrees
! tincr        Time increment                                                                           hour
! ts           Number of hourly time steps per day                                                      -
! xc           Parameter X for calculating black layer extinction coefficient                           -
!              according to Campbell (1986)
!==============================================================================================================================

!=======================================================================
! Respiration subroutine for AgMaize by Thijs Tollenaar
! Calculates maintenance respiration
!-----------------------------------------------------------------------
! REVISION HISTORY
! Originally written by Thijs Tollenaar
! 06/17/2013 Translated into Fortran - KAD
!=======================================================================
subroutine MZ_TUE_Respiration(control, weather, & !DSSAT inputs
                              dvs, & !From Phenology
                              wlvtop, wlvbot, wst, weolg, wrt, cgrw, & !From Growth
                              maint)                                    !Outputs

!-----------------------------------------------------------------------
   use ModuleDefs
   implicit none
   save
!-----------------------------------------------------------------------
   real, parameter :: qten = 2.0
   integer :: dynamic
   real :: tmax, tmin, tmpa, dvs, wlvtop, wlvbot, wst, weolg, wrt, cgrw, teff, maint

! Constructed variables based on definitions in ModuleDefs.for
   type(ControlType), intent(in) :: control
   type(WeatherType), intent(in) :: weather

! Transfer values from constructed data types into local variables
   dynamic = control%dynamic
   tmax = weather%tmax
   tmin = weather%tmin

!-----------------------------------------------------------------------
! Dynamic = runinit
!-----------------------------------------------------------------------
   if (dynamic == runinit .or. dynamic == seasinit) then
      maint = 0.0
      wlvtop = 0.0
      wlvbot = 0.0
      wst = 0.0
      weolg = 0.0
      wrt = 0.0
      cgrw = 0.0

!-----------------------------------------------------------------------
! Dynamic = rate
!-----------------------------------------------------------------------
   else if (dynamic == rate) then

! Effect of temperature on maintenance respiration
      tmpa = (tmax + tmin)/2.0       !Daily average air temperature
      teff = qten**(0.1*tmpa - 2.5)  !No temperature effect at 25 degrees C (i.e. teff = 1.0)

! Maintenance respiration
      if (cgrw <= 0.0) then
         maint = 0.0
      else
         maint = teff*(0.03*wlvtop + 0.015*(wlvbot + wst) + 0.01*weolg + 0.01*wrt)
      end if

!-----------------------------------------------------------------------
! End of dynamic if structure
!-----------------------------------------------------------------------
   end if

!-----------------------------------------------------------------------
! End of subroutine
!-----------------------------------------------------------------------
   return
end subroutine MZ_TUE_Respiration

!=======================================================================
!  MZ_TUE_Radabs, Subroutine, J.I. Lizaso
!  Calculates hourly canopy absorption of PAR (J/m2 s) by shaded and
!  sunlit leaf area
!  Adapted from RADABS by N.B. Pickering
!-----------------------------------------------------------------------
!  REVISION HISTORY
!  02/07/2003 JIL Written
!  09/12/2007 JIL Adapted for IXIM model
!  03/13/2013 Adapted for AgMaize
!-----------------------------------------------------------------------
!  TO DO:
!  -Fix reflected radiation
!=======================================================================
subroutine MZ_TUE_Radabs(azir, azzon, beta, betn, canht, canwh, daytim, & !Input
                         frdifp, froll, gla, h, lfn, palb, parh, pltpop, rowspc, scvp, hlai, xc, & !Input
                         parsh, parsun, plaish, plaisl, radtot, & !Output
                         arefhr, adifhr, addrhr, addfhr, leafnum)                                    !Output

   implicit none
   save

   logical :: daytim
   integer :: h, i, lfn, leafnum(50)
   real, dimension(50) :: gla, parsh, parsun, plaish, plaisl, arefhr, adifhr, addrhr, addfhr
   real :: azir, azzon, beta, betn, canht, canwh, fracsh, frdifp, froll, frshv, hlai, kdifbl
   real :: kdirbl, kdrblv, laish, laishv, laisl, laislv, palb, parh, parsl, parss, pcabsp, pcabsr
   real :: pcintp, pcintr, pcrefp, pcrefr, pltpop, rowspc, scvp, xc, radtot

! Initialize.
   if (h == 1) then        !Daily initialization
      frshv = 0.0
      laishv = 0.0
      laislv = 0.0
   end if
   pcabsp = 0.0         !Hourly initialization
   pcintp = 0.0
   pcrefp = 0.0
   parsl = 0.0
   parss = 0.0
   pcabsr = 0.0
   pcintr = 0.0
   pcrefr = 0.0
   parsun = 0.0
   parsh = 0.0
   radtot = 0.0

   if (hlai > 0.000) then
      !Calculate fraction shaded and LAI's for vertical sun position
      if (h == 1) then
         call MZ_TUE_Shadow(azir, azzon, 90.0, betn, canht, canwh, rowspc, & !Input
                            frshv)                                              !Output
         call MZ_TUE_Lfextn(90.0, froll, frshv, h, hlai, xc, & !Input
                            kdifbl, kdrblv, laishv, laislv)                     !Output
      end if

      if (daytim) then
         !Calculate fraction shaded, leaf extinction and LAI's
         call MZ_TUE_Shadow(azir, azzon, beta, betn, canht, canwh, rowspc, & !Input
                            fracsh)                                             !Output
         call MZ_TUE_Lfextn(beta, froll, fracsh, h, hlai, xc, & !Input
                            kdifbl, kdirbl, laish, laisl)                       !Output

         !Calculate PAR absorbed by canopy during day
         call MZ_TUE_Canabs(palb, beta, betn, canht, canwh, fracsh, froll, & !Input
                            gla, lfn, pltpop, frdifp, kdifbl, kdirbl, parh, rowspc, scvp, & !Input
                            pcabsp, pcintp, pcrefp, plaish, plaisl, parsh, parsun, radtot, & !Output
                            arefhr, adifhr, addrhr, addfhr, leafnum)
      else
         !Night time with canopy
         fracsh = frshv
         kdirbl = kdrblv
         laish = laishv
         laisl = laislv
      end if

   else
      !Bare soil (day or night)
      pcabsp = (1.0 - palb)*100.0
      pcrefp = palb*100.0
      parss = (1.0 - palb)*parh
   end if

   return
end subroutine MZ_TUE_Radabs

!==============================================================================================================================
! Variable definitions                                                                                  Unit
!------------------------------------------------------------------------------------------------------------------------------
! azir         Row azimuth relative to North                                                            degrees
! azzon        Hourly solar azimuth (+/- from South)                                                    degrees
! beta         Hourly solar elevation (+/- from horizontal)                                             degrees
! betn         Spacing between plants along a row                                                       m/plant
! canht        Canopy height                                                                            m
! canwh        Canopy width                                                                             m
! daytim       Logical variable for test for daylight hours between sunrise and sunset                  -
! fracsh       Fraction of land shaded by the canopy
! frdifp(h)    Hourly fraction diffuse photon flux density after correcting                             -
!              for circumsolar radiation (Spitters, 1986)
! froll        Leaf rolling factor associated with soil water stress affecting cell expansion           -
! frshv        Fraction of land area shaded with vertical sun position
! gla(l)       Green leaf area for leaf l                                                               cm2
! h            Hourly iteration counter (1-24)                                                          -
! hlai         "Healthy" leaf area index (excluding senesced parts)                                     m2/m2
! kdifbl       Extinction coefficient of black leaves to diffuse radiation
! kdrblv       Extinction coefficient of black leaves to direct radiation for vertical sun position
! kdirbl       Extinction coefficient of black leaves to direct radiation
! laish        Shaded leaf area index                                                                   m2[leaf]/m2[ground]
! laishv       Shaded leaf area index for vertical sun position                                         m2[leaf]/m2[ground]
! laisl        Sunlit leaf area index                                                                   m2[leaf]/m2[ground]
! laislv       Sunlit leaf area index for vertical sun position                                         m2[leaf]/m2[ground]
! lfn          Total leaf number rounded up                                                             -
! palb         Plant albedo accounting for soil water in the first soil layer                           -
! parh         Hourly photosynthetically active radiation (PAR)                                         �mol[quanta]/m2/s
! parsh(l)     Photosynthetically active radiation absorbed by leaf l in the shaded zone                �mol[quanta]/m2/s
! parsun(l)    Photosynthetically active radiation absorbed by leaf l in the sunlit zone                �mol[quanta]/m2/s
! pcabsp       Hourly percent absorbed photosynthetic radiation                                         %
! pcintp       Hourly percent intercepted photosynthetic radiation                                      %
! pcrefp       Hourly percent reflected photosynthetic radiation                                        %
! plaish(l)    Shaded leaf area index for leaf l                                                        m2/m2
! plaisl(l)    Sunlit leaf area index for leaf l                                                        m2/m2
! pltpop       Plant density                                                                            plants/m2
! radtot       Sum of PAR components (energy balance check: should equal parh)                          �mol[quanta]/m2/s
! rowspc       Row spacing                                                                              m
! scvp         Scattering coefficient used to calculate diffuse reflectance
! xc           Parameter X for calculating black layer extinction coefficient                           -
!              according to Campbell (1986)
!==============================================================================================================================

!=======================================================================
!  SHADOW, Subroutine, N.B. Pickering, J.W. Jones
!  Calculates fraction shaded for sun and row geometry using ellipses.
!-----------------------------------------------------------------------
!  REVISION HISTORY
!  03/14/1991 NBP Written
!  11/15/1991 NBP Modified
!  11/23/1993 NBP Included more error checks and limits
!  09/12/2007 JIL Modified and adapted for IXIM model
!=======================================================================
subroutine MZ_TUE_Shadow(azir, azzon, beta, betn, canht, canwh, rowspc, &  !Input
                         fracsh)                                              !Output

   implicit none
   save

   real, parameter :: pi = 3.14159, rad = pi/180.0, zero = 1.0E-6
   real :: a, b, c1, c2, c3, c4, azimd, azir, azzon, beta, betn, canht, canwh, eta
   real :: fracsh, gamma, rbeta, rowspc, shade, shlen, shperp, stouch

! Initialization
   a = 0.0
   b = 0.0
   c1 = 0.0
   c2 = 0.0
   c3 = 0.0
   c4 = 0.0
   azimd = 0.0
   eta = 0.0
   fracsh = 0.0
   gamma = 0.0
   rbeta = 0.0
   shade = 0.0
   shlen = 0.0
   shperp = 0.0
   stouch = 0.0

! Set fraction shaded to 0.0 for zero width or height
   if (canwh <= zero .or. canht <= zero) then
      fracsh = 0.0

! Set fraction shaded to 1.0 for full cover
   else if (canwh >= rowspc) then
      fracsh = 1.0

! Calculate fraction shaded
   else
      !Adjust BETA if sun exactly overhead or at horizon.  Calculate acute positive angle between
      !sun and row azimuths. Initialize other constants
      rbeta = min(max(beta*rad, 1.0e-6), pi/2.0 - 1.0e-6)
      azimd = abs(azzon - azir)*rad
      if (azimd > pi/2.0) azimd = pi - azimd
      a = (canwh/canht)**2
      gamma = atan(a*tan(rbeta))
      c1 = a*(tan(rbeta))**2
      c2 = (a*tan(rbeta))**2

      ! Calculate shadow length assuming ellipsoidal shape for the plant
      shlen = canht*cos(rbeta - gamma)/sin(rbeta)*sqrt((1.0 + c2)/(1.0 + c1))
      b = (shlen/canwh)**2
      c3 = b*(tan(azimd))**2
      c4 = (b*tan(azimd))**2
      stouch = shlen/(cos(azimd)*sqrt(1.+c3))

      ! CALCULATE FRACTION SHADED
      ! Sun parallel to row.  Shadow limited to BETN
      if (azimd <= zero) then
         shlen = min(shlen, betn)
         shade = 0.25*pi*shlen*canwh

         !Sun not parallel to row
      else
         !Calculate perpendicular shadow length
         azimd = max(azimd, 1.0e-6)
         eta = atan(1.0/(b*tan(azimd)))
         shperp = canwh*sin(azimd + eta)*sqrt((1.0 + c4)/(1.0 + c3))

         !Hedgerow (plant shadows overlap)
         if (stouch >= betn) then
            !Shadow length is perpendicular and limited to ROWSPC
            shlen = min(shperp, rowspc)
            shade = shlen*betn

            !Individual plants
         else
            !Limit shadow length to within one ROWSPC
            if (shperp > rowspc) shlen = shlen*rowspc/shperp
            shade = 0.25*pi*shlen*canwh
         end if
      end if

      fracsh = min(shade/(rowspc*betn), 1.0)

   end if

   fracsh = min(max(fracsh, 1.0e-6), 1.0)

end subroutine MZ_TUE_Shadow

!==============================================================================================================================
! Variable definitions                                                                                  Unit
!------------------------------------------------------------------------------------------------------------------------------
! azir         Row azimuth relative to North                                                            degrees
! azzon        Hourly solar azimuth (+/- from South)                                                    degrees
! beta         Hourly solar elevation (+/- from horizontal)                                             degrees
! betn         Spacing between plants along a row                                                       m/plant
! canht        Canopy height                                                                            m
! canwh        Canopy width                                                                             m
! fracsh       Fraction of land shaded by the canopy
! rowspc       Row spacing                                                                              m
! shade        Shadow area per plant                                                                    m2
! shlen        Shadow length per plant                                                                  m
! shperp       Length of the shadow perpendicular to the row                                            m
!==============================================================================================================================

!=======================================================================
!  LFEXTN, Subroutine, N.B. Pickering, K.J. Boote
!  Computes leaf extinction coefficients based on leaf angle distribution
!  (Goudriaan, 1988) and sunlit and shaded leaf area indices.
!-----------------------------------------------------------------------
!  REVISION HISTORY
!  ??/??/???? KJB Written
!  05/14/1991 NBP Removed COMMON and reorganized.
!  08/31/2004 JIL Adapted for Maize photosynthesis
!  09/12/2007 JIL Modified and adapted for IXIM model
!=======================================================================

subroutine MZ_TUE_Lfextn(beta, froll, fracsh, h, hlai, xc, & !Input
                         kdifbl, kdirbl, laish, laisl)                        !Output

   implicit none
   save

   integer :: i, h
   real, parameter :: pi = 3.14159, rad = pi/180.0, expbound = -40.0
   real :: beta, fracsh, kdirbl, kdifbl, laish, froll, laisl, hlai
   real :: bt, kdir, taudiff, taudir, xc, exponent

! Calculate black leaf extinction coefficients for diffuse and direct radiation
! Campbell (1986); Campbell (1990); Campbell and Norman (1998)
! Initialize
   if (h == 1) then     !Daily initialization
      bt = 0.0
      kdifbl = 0.0
      kdir = 0.0
      taudiff = 0.0
      taudir = 0.0
   end if

   i = 1           !Hourly initialization
   kdirbl = 0.0
   laish = 0.0
   laisl = 0.0

! Calculate KDIFBL once a day assuming Uniform OverCast sky (UOC)
   if (h == 1) then
      taudiff = 0.0
      do i = 1, 89, 2
         bt = real(i)
         kdir = ((xc**2.+1./(tan(bt*rad))**2.)**0.5)/(xc + 1.744*(xc + 1.182)**(-0.733))

         !Underflows  chp 8/5/2009
         exponent = -kdir*hlai*froll
         if (exponent > expbound) then
            taudir = exp(-kdir*hlai*froll)
         else
            taudir = 0.0
         end if

         taudiff = taudiff + (taudir*sin(bt*rad)*cos(bt*rad))
      end do
      taudiff = 2.0*taudiff*2.0*rad
      kdifbl = -log(taudiff)/(hlai*froll)
   end if

! Calculate KDIRBL hourly as a function of solar elevation
   kdirbl = ((xc**2.+1./(tan(beta*rad))**2.)**0.5)/(xc + 1.744*(xc + 1.182)**(-0.733))

! Calculate sunlit and shaded leaf area indices
   exponent = -kdirbl*hlai*froll/fracsh
   if (exponent > expbound) then
      laisl = (fracsh/kdirbl)*(1.0 - exp(-kdirbl*hlai*froll/fracsh))
   else
      laisl = fracsh/kdirbl
   end if

   if (hlai*froll > 0.02) then
      laisl = max(laisl, 0.02)
   else
      laisl = hlai*froll
   end if

   laish = hlai*froll - laisl

end subroutine MZ_TUE_Lfextn

!==============================================================================================================================
! Variable definitions                                                                                  Unit
!------------------------------------------------------------------------------------------------------------------------------
! beta         Hourly solar elevation (+/- from horizontal)                                             degrees
! fracsh       Fraction of land shaded by the canopy
! froll        Leaf rolling factor associated with soil water stress affecting cell expansion           -
! h            Hourly iteration counter (1-24)                                                          -
! hlai         "Healthy" leaf area index (excluding senesced parts)                                     m2/m2
! kdifbl       Extinction coefficient of black leaves to diffuse radiation
! kdirbl       Extinction coefficient of black leaves to direct radiation
! laish        Shaded leaf area index                                                                   m2[leaf]/m2[ground]
! laisl        Sunlit leaf area index                                                                   m2[leaf]/m2[ground]
!==============================================================================================================================

!=======================================================================
!  CANABS, Subroutine, K.J. Boote, N.B. Pickering
!  Computes radiation absorbed by soil, sunlit and shaded leaves.
!-----------------------------------------------------------------------
!  REVISION HISTORY
!  ??/??/?? KJB Written
!  05/14/91 NBP Removed COMMON and reorganized.
!  08/31/04 JIL Adapted to photosynthesis in maize
!  09/12/2007 JIL Modified and adapted for IXIM model
! Q. Why not add addrsl to cumasun instead of (1.-scvr)*raddir*kdirbl*ylaisl
!=======================================================================

subroutine MZ_TUE_Canabs(albedo, beta, betn, canht, canwh, fracsh, froll, & !Input
                         gla, lfn, pltpop, frdif, kdifbl, kdirbl, parh, rowspc, scvr, & !Input
                         pctabs, pctint, pctref, plaish, plaisl, radsh, radsun, radtot, & !Output
                         arefhr, adifhr, addrhr, addfhr, leafnum)

   implicit none
   save

   real, parameter :: pi = 3.14159, rad = pi/180.0
   integer :: l, lfn, leafnum(50)
   real, dimension(50) :: gla, plaish, plaisl, radsh, radsun, arefhr, adifhr, addrhr, addfhr
   real :: addf, addfsh, addfsl, addr, addrsl, adif, adifsh, adifsl, adir, adirsl, albedo
   real :: aref, arefsh, arefsl, atot, beta, betn, canht, canwh, cumash, cumasun, delwp
   real :: delwr, difp, difpr, difr, fracsh, frdif, froll, incsoi, intcan, kdifbl, kdirbl
   real :: parh, pathp, pathr, pctabs, pctint, pctref, pltpop, raddif, raddir, radss, radtot
   real :: rdifsl, refdf, refdif, refdir, refdr, refsoi, reftot, rowspc, scvr, sinb
   real :: sqv, ucumash, ucumasun, uylaish, uylaisl, yhlai, ylaish, ylaisl
   real :: yhlaican, ylaislcan, ylaishcan, adircan, adifcan, intcanall, incsoiall, refsoiall

! Initialization.
   addf = 0.0
   addfsh = 0.0
   addfsl = 0.0
   addr = 0.0
   addrsl = 0.0
   adif = 0.0
   adifsh = 0.0
   adifsl = 0.0
   adir = 0.0
   adirsl = 0.0
   aref = 0.0
   arefsh = 0.0
   arefsl = 0.0
   atot = 0.0
   cumash = 0.0
   cumasun = 0.0
   delwp = 0.0
   delwr = 0.0
   difp = 0.0
   difpr = 0.0
   difr = 0.0
   incsoi = 0.0
   intcan = 0.0
   l = 1
   pathp = 0.0
   pathr = 0.0
   pctabs = 0.0
   pctint = 0.0
   pctref = 0.0
   raddif = 0.0
   raddir = 0.0
   radss = 0.0
   radtot = 0.0
   rdifsl = 0.0
   refdf = 0.0
   refdif = 0.0
   refdir = 0.0
   refdr = 0.0
   refsoi = 0.0
   reftot = 0.0
   ucumash = 0.0
   ucumasun = 0.0
   uylaish = 0.0
   uylaisl = 0.0
   yhlai = 0.0
   ylaish = 0.0
   ylaisl = 0.0
   plaish = 0.0
   plaisl = 0.0
   radsh = 0.0
   radsun = 0.0
   yhlaican = 0.0
   adircan = 0.0
   adifcan = 0.0
   intcanall = 0.0
   incsoiall = 0.0
   refsoiall = 0.0
   sqv = sqrt(1.0 - scvr)
   sinb = sin(beta*rad)
   arefhr = 0.0
   adifhr = 0.0
   addrhr = 0.0
   addfhr = 0.0
   leafnum = 0.0

!-Compute reflection coefficients for direct and diffuse light
   refdf = (1.0 - sqv)/(1.0 + sqv)
   refdr = 2.0*kdirbl/(1.0 + kdirbl)*refdf

!-Compute difpr, the fraction of sky "seen" by plants.
!Diffuse skylight is absorbed over an effective area equal to the
!canopy height plus width for an isolated row.  For interfering rows,
!Eqns. 2.130 and 2.128 from Goudriaan (1977) are used.  Concept
!extended for both between plants (P) and rows (R)
   if (canwh < betn) then
      pathp = betn - canwh
      delwp = pathp + canht - sqrt(pathp**2 + canht**2)
      difp = min((canwh + delwp)/betn, 1.0)
   else
      difp = 1.0
   end if

   if (canwh < rowspc) then
      pathr = rowspc - canwh
      delwr = pathr + canht - sqrt(pathr**2 + canht**2)
      difr = min((canwh + delwr)/rowspc, 1.0)
   else
      difr = 1.0
   end if
   difpr = min(max(difp*difr, 1.0e-6), 1.0)

!-Partition total radiation into direct and diffuse components
   raddif = frdif*parh
   raddir = parh - raddif

!-Compute total soil-reflected radiation before entering the loop. The original codes
!recompute soil reflected light in each iteration and then calculate the plant absorption
!from it. In reality, for a given parh and canopy, soil reflected light is constant but it is
!the plant absorption from it that changes depending on the absorbing leaf layer.
!Correction added KAD 6/20/2013
   yhlaican = sum(gla(1:lfn))*froll*pltpop*0.0001
   ylaislcan = (fracsh/kdirbl)*(1.0 - exp(-kdirbl*yhlaican/fracsh))
   ylaishcan = yhlaican - ylaislcan
   adircan = fracsh*(1.0 - refdr)*raddir*(1.0 - exp(-kdirbl*sqv*yhlaican/fracsh))
   adifcan = difpr*(1.0 - refdf)*raddif*(1.0 - exp(-kdifbl*sqv*yhlaican/difpr))
   refdir = fracsh*refdr*raddir
   refdif = difpr*refdf*raddif
   intcanall = refdir + refdif + adircan + adifcan
   incsoiall = parh - intcanall
   refsoiall = albedo*incsoiall

!-Direct radiation absorbed in shaded zone by considering the direct
!(ADDR) and diffuse/scattered (ADDF) components of the direct beam.
! ** JIL Beginning per-leaf loop
   do l = lfn, 1, -1
      if (gla(l) > 0.0) then
         !-Cumulative leaf area index (m2/m2) from top to bottom with
         !sunlit and shaded components separate
         yhlai = yhlai + gla(l)*froll*pltpop*0.0001
         ylaisl = (fracsh/kdirbl)*(1.0 - exp(-kdirbl*yhlai/fracsh))
         ylaish = yhlai - ylaisl

         !-Absorbed from direct (Equations 17-23 in Lizaso et al. 2005)
         !First, compute total absorbed from direct PAR (adir) and its two components
         !(not scattered addr, and scattered addf)
         adir = fracsh*(1.0 - refdr)*raddir*(1.0 - exp(-kdirbl*sqv*yhlai/fracsh))
         addr = fracsh*(1.0 - scvr)*raddir*(1.0 - exp(-kdirbl*yhlai/fracsh))
         addf = adir - addr
         !From the previous three components, compute quantities absorbed by sunlit and shaded LAI.
         !Note that shaded LAI only absorbs the scattered component not absorbed by the sunlit LAI.
         adirsl = fracsh*(1.0 - refdr)*raddir*(1.0 - exp(-kdirbl*sqv*ylaisl/fracsh))
         addrsl = fracsh*(1.0 - scvr)*raddir*(1.0 - exp(-kdirbl*ylaisl/fracsh))
         addfsl = adirsl - addrsl
         addfsh = addf - addfsl

         !-Absorbed from incident diffuse
         adif = difpr*(1.0 - refdf)*raddif*(1.0 - exp(-kdifbl*sqv*yhlai/difpr))
         adifsl = difpr*(1.0 - refdf)*raddif*(1.0 - exp(-kdifbl*sqv*ylaisl/difpr))
         adifsh = adif - adifsl

         !-Light reflected from the soil assumed to be isotropic and diffuse.
         !Absorption handled in the same manner as diffuse skylight
         !refdir = fracsh * refdr * raddir
         !refdif = difpr * refdf * raddif
         intcan = refdir + refdif + adir + adif       !Light intercepted by the canopy
         !incsoi = parh - intcan                      !Light intercepted by the soil
         !refsoi = albedo * incsoi                    !Light reflected by the soil

         !-Absorbed light by the canopy from the soil-reflected light
         !aref = difpr * (1.0-refdf) * refsoi * (1.0-exp(-kdifbl*sqv*yhlai/difpr))
         !arefsh = difpr * (1.0-refdf) * refsoi * (1.0-exp(-kdifbl*sqv*yhlaish/difpr))
         !Check and correct
         aref = difpr*(1.0 - refdf)*refsoiall*(1.0 - exp(-kdifbl*sqv*(yhlaican - yhlai)/difpr))
         arefsh = difpr*(1.0 - refdf)*refsoiall*(1.0 - exp(-kdifbl*sqv*(ylaishcan - ylaish)/difpr))
         arefsl = aref - arefsh

         !-Total absorbed and total reflected light
         atot = adir + adif + aref
         reftot = refdir + refdif + refsoiall - aref

         !-Cumulative light absorbed down the canopy by sunlit and shaded LAI
         !radss = incsoi * (1.0-albedo)
         radss = incsoiall*(1.0 - albedo)
         rdifsl = addfsl + adifsl + arefsl
         cumasun = rdifsl + (1.-scvr)*raddir*kdirbl*ylaisl
         cumash = addfsh + adifsh + arefsh

         !-Determine per-leaf sunlit and shaded LAI and sunlit and shaded PAR
         plaisl(l) = ylaisl - uylaisl
         plaish(l) = ylaish - uylaish

         if (plaisl(l) < 1e-4) then
            radsun(l) = 0.0
         else
            radsun(l) = (cumasun - ucumasun)/plaisl(l)
         end if

         if (plaish(l) == 0.0) then
            radsh(l) = 0.0
         else
            radsh(l) = (cumash - ucumash)/plaish(l)
         end if

         !-Save current LAI and cumulative radiation values for next iteration
         uylaish = ylaish
         uylaisl = ylaisl
         ucumasun = cumasun
         ucumash = cumash

         !-KAD 09/24/2013 - Save different absorbed radiation components for checking
         arefhr(l) = aref
         adifhr(l) = adif
         addrhr(l) = addr
         addfhr(l) = addf
         leafnum(l) = l
      end if
   end do
! JIL End of per-leaf loop

! Set radiation array and calculate ratios of components
   if (parh > 0.0) then
      pctint = 100.0*intcan/parh
      pctabs = 100.0*atot/parh
      pctref = 100.0*reftot/parh
   else
      pctint = 0.0
      pctabs = 0.0
      pctref = 0.0
   end if

! Energy balance check (RADTOT=PARH)
   radtot = atot + reftot + radss

end subroutine MZ_TUE_Canabs

!==============================================================================================================================
! Variable definitions                                                                                  Unit
!------------------------------------------------------------------------------------------------------------------------------
! addf         Absorbed direct PAR that is scattered and becomes diffuse                                �mol[quanta]/m2/s
! addfsh       Absorbed direct PAR in the shaded LAI that is scattered and becomes diffuse              �mol[quanta]/m2/s
! addfsl       Absorbed direct PAR in the sunlit LAI that is scattered and becomes diffuse              �mol[quanta]/m2/s
! addr         Absorbed direct PAR that is not scattered                                                �mol[quanta]/m2/s
! addrsl       Absorbed direct PAR in the sunlit LAI that is not scattered                              �mol[quanta]/m2/s
! adif         Total absorbed diffuse PAR                                                               �mol[quanta]/m2/s
! adifsh       Diffuse PAR captured by shaded LAI                                                       �mol[quanta]/m2/s
! adifsl       Diffuse PAR captured by sunlit LAI                                                       �mol[quanta]/m2/s
! adir         Total absorbed direct PAR                                                                �mol[quanta]/m2/s
! adirsl       Direct PAR captured by sunlit LAI                                                        �mol[quanta]/m2/s
! albedo       Plant albedo accounting for soil water in the first soil layer                           -
! aref         Total canopy absorption from the soil-reflected PAR                                      �mol[quanta]/m2/s
! arefsh       Shaded LAI PAR absorption from the soil-reflected PAR                                    �mol[quanta]/m2/s
! arefsl       Sunlit LAI PAR absorption from the soil-reflected PAR                                    �mol[quanta]/m2/s
! atot         Total PAR absorbed by the canopy                                                         �mol[quanta]/m2/s
! beta         Hourly solar elevation (+/- from horizontal)                                             degrees
! betn         Spacing between plants along a row                                                       m
! canht        Canopy height                                                                            m
! canwh        Canopy width                                                                             m
! cumash       Cumulative PAR absorbed by shaded LAI (down to current leaf)                             �mol[quanta]/m2/s
! cumasun      Cumulative PAR absorbed by sunlit LAI (down to current leaf)                             �mol[quanta]/m2/s
! difpr        Fraction of sky "seen" by plants
! fracsh       Fraction of land shaded by the canopy (during daytime)
! frdif        Fraction diffuse photon flux density after correcting for circumsolar
!              radiation (Spitters, 1986)
! froll        Leaf rolling factor associated with soil water stress affecting cell expansion           -
! gla(l)       Green leaf area for leaf l                                                               cm2
! incsoi       Total PAR intercepted by the soil                                                        �mol[quanta]/m2/s
! intcan       Total PAR intercepted by the canopy                                                      �mol[quanta]/m2/s
! kdifbl       Extinction coefficient of black leaves to diffuse radiation
! kdirbl       Extinction coefficient of black leaves to direct radiation
! l            Integer loop counter
! lfn          Total leaf number rounded up                                                             -
! parh         Hourly photosynthetically active radiation (PAR)                                         �mol[quanta]/m2/s
! pctabs       Percent of total PAR absorbed by the canopy                                              �mol[quanta]/m2/s
! pctint       Percent of total PAR intercepted by the canopy                                           �mol[quanta]/m2/s
! pctref       Percent of total PAR reflected by the canopy and the soil                                �mol[quanta]/m2/s
! plaish(l)    Shaded leaf area index for leaf l                                                        m2[leaf]/m2[ground]
! plaisl(l)    Sunlit leaf area index for leaf l                                                        m2[leaf]/m2[ground]
! pltpop       Plant density                                                                            plants/m2
! raddif       Hourly diffuse PAR                                                                       �mol[quanta]/m2/s
! raddir       Hourly direct PAR                                                                        �mol[quanta]/m2/s
! radsh(l)     Photosynthetically active radiation absorbed by leaf l in the shaded zone                �mol[quanta]/m2/s
! radss        Radiation on soil                                                                        �mol[quanta]/m2/s
! radsun(l)    Photosynthetically active radiation absorbed by leaf l in the sunlit zone                �mol[quanta]/m2/s
! radtot       Sum of PAR components (energy balance check: should equal parh)                          �mol[quanta]/m2/s
! rdifsl       Total diffuse PAR absorbed by sunlit LAI                                                 �mol[quanta]/m2/s
! refdif       Canopy-reflected diffuse PAR                                                             �mol[quanta]/m2/s
! refdir       Canopy-reflected direct PAR                                                              �mol[quanta]/m2/s
! refdr        Reflection coefficient for direct radiation
! refdf        Reflection coefficient for diffuse radiation
! refsoi       Diffuse PAR reflected from the soil                                                      �mol[quanta]/m2/s
! reftot       Total PAR reflected by the canopy and the soil                                           �mol[quanta]/m2/s
! rowspc       Row spacing                                                                              m
! scvr         Scattering coefficient used to calculate diffuse reflectance
! ucumash      Cumulative PAR absorbed by shaded LAI (down to previous leaf)                            �mol[quanta]/m2/s
! ucumasun     Cumulative PAR absorbed by sunlit LAI (down to previous leaf)                            �mol[quanta]/m2/s
! uylaish      Cumulative shaded LAI from top to bottom of canopy (down to previous leaf)               m2[leaf]/m2[ground]
! uylaisl      Cumulative sunlit LAI from top to bottom of canopy (down to previous leaf)               m2[leaf]/m2[ground]
! yhlai        Cumulative leaf area index from top to bottom of canopy                                  m2[leaf]/m2[ground]
! ylaish       Cumulative shaded LAI from top to bottom of canopy (down to current leaf)                m2[leaf]/m2[ground]
! ylaisl       Cumulative sunlit LAI from top to bottom of canopy (down to current leaf)                m2[leaf]/m2[ground]
!==============================================================================================================================

!=======================================================================
!  Iphotsynt, Subroutine, J.I. Lizaso
!  Calculates instantaneous plant photosynthesis (�mol CO2/m2 s) of
!  shaded and sunlit leaf area on a per-leaf basis and integrates for
!  the whole canopy
!  Adapted from CANOPG by K.J. Boote, J.W. Jones, G. Hoogenboom
!-----------------------------------------------------------------------
!  REVISION HISTORY
!  02/07/2003 JIL Written
!  09/12/2007 JIL Adapted for IXIM model
!  03/13/2013 KAD/TT/SK Adapted for AgMaize
!  08/28/2013 KAD Added definitions of variables
!-----------------------------------------------------------------------
! Questions for Jon:
! -What is pgsum? Looks like it is not being used?
! -Parameters for equations to estimate light response curve parameters in MZ_GM_PSParam? What is axx?
! -Data for the equations mentioned above
! -Temperature effect for curvature parameter in MZ_GM_PSParam tair**2 instead of tair*2
!=======================================================================
subroutine MZ_TUE_Iphotsynt(dynamic, & !Control
                            tmin, asmax, gddae, gla, plaisl, plaish, lap, lfl, lfn, & !Input
                            light, parsh, parsun, tairh, yx, sgfun, gasfn, & !Input
                            pghr)                                                       !Output

   use ModuleDefs
   implicit none
   save

   logical :: light
   integer :: dynamic, i, lfn
   real :: asmax, gddae, pghr, pgsum, pgsun, tairh, sgfun, gasfn, tmin
   real, dimension(50) :: assat, cvxty, gla, intslp, lap, lfl, parsh, parsun, pgsh, pgsl, plaish, plaisl, yx

!----------------------------------------------------------------------
! Dynamic = runinit
!----------------------------------------------------------------------
   if (dynamic == runinit .or. dynamic == seasinit) then

! Initialize.
      i = 1
      pghr = 0.0
      pgsum = 0.0
      pgsun = 0.0
      pgsh = 0.0
      pgsl = 0.0

      call MZ_TUE_PSParam(dynamic, & !Control
                          tmin, asmax, gddae, gla, lap, lfl, lfn, light, tairh, yx, sgfun, gasfn, & !Input
                          assat, cvxty, intslp)                                                       !Output

!-----------------------------------------------------------------------
! Dynamic = rate
!-----------------------------------------------------------------------
   else if (dynamic == rate) then

! JIL Calculate light response curve parameters per leaf as a function of leaf age and hourly air temperature
      call MZ_TUE_PSParam(dynamic, tmin, asmax, gddae, gla, lap, lfl, lfn, light, tairh, yx, sgfun, gasfn, & !Input
                          assat, cvxty, intslp)                                                               !Output

! Calculate per-leaf assimilation
      pghr = 0.0
      do i = 1, lfn
         if (gla(i) > 0.05*yx(i)) then
            pgsum = 0.0

            !Compute photosynthesis for sunlit leaves
            call MZ_TUE_PSLeaf(parsun(i), assat(i), intslp(i), cvxty(i), & !Input
                               pgsun)                                         !Output
            pgsl(i) = pgsun

            !Compute photosynthesis for shaded leaves
            call MZ_TUE_PSLeaf(parsh(i), assat(i), intslp(i), cvxty(i), &  !Input
                               pgsh(i))                                       !Output

            !Compute instantaneous canopy gross photosynthesis (�mol[CO2]/m2/s)
            pghr = pghr + pgsl(i)*plaisl(i) + pgsh(i)*plaish(i)
         end if
      end do

   end if       !Endif for DYNAMIC LOOP

   return

end subroutine MZ_TUE_Iphotsynt

!==============================================================================================================================
! Variable definitions                                                                                  Unit
!------------------------------------------------------------------------------------------------------------------------------
! assat(l)     Light-saturated assimilation rate (corresponds to the asymptote of the light             umol[CO2]/m2[leaf]/s
!              response curve) for leaf l
! cvxty(l)     Ratio of the diffusion resistance to the total resistance to CO2 (curvature              -
!              parameter for the light response curve) for leaf l
! dynamic      Main control variable to tell each module which section of code to run                   -
! asmax        Maximum instantaneous assimilation at 30 degC                                            �mol[CO2]/m2[leaf]/s
! gddae        Cumulative growing degree day after emergence                                            degree-days
! gla(l)       Green leaf area for leaf l                                                               cm2
! i            Iteration counter of leaf number                                                         -
! intslp(l)    Quantum efficiency of CO2 assimilation (initial slope of the light response              �mol[CO2]/�mol[quanta]
!              curve) for leaf l
! plaish(l)    Shaded leaf area index for leaf l                                                        m2/m2
! plaisl(l)    Shaded leaf area index for leaf l                                                        m2/m2
! lap(l)       Potential leaf area of leaf l from tip appearance to full expansion                      cm2
! lfl(l)       Longevity of leaf l                                                                      degree-days
! lfn          Total leaf number rounded up                                                             -
! light        Logical variable that differentiates daytime from nighttime
! parsh(l)     Photosynthetically active radiation absorbed by leaf l in the shaded zone                J/m2/s
! parsun(l)    Photosynthetically active radiation absorbed by leaf l in the sunlit zone                J/m2/s
! pghr         Canopy instantaneous gross assimilation                                                  �mol[CO2]/m2[ground]/s
! pgsh(l)      Shaded photosynthesis for leaf l                                                         �mol[CO2]/m2[leaf]/s
! pgsl(l)      Sunlit photosynthesis for leaf l                                                         �mol[CO2]/m2[leaf]/s
! pgsun(l)     Sunlit photosynthesis for leaf l (intermediary variable)                                 �mol[CO2]/m2[leaf]/s
! tairh(h)     Hourly air temperature                                                                   degrees
! yx(l)        Maximum value of potential leaf area for leaf l (when fully expanded)                    cm2
!==============================================================================================================================

!=======================================================================
!  MZ_AG_PSParam, Subroutine, J.I. Lizaso
!  Calculates parameters of light response curve per leaf as affected by
!  leaf age and hourly air temperature
!  Adapted from PGLFEQ by K.J. Boote, N.B. Pickering
!-----------------------------------------------------------------------
!  REVISION HISTORY
!  02/07/2003 JIL Written
!  09/12/2007 JIL Adapted for IXIM model
!=======================================================================
subroutine MZ_TUE_PSParam(dynamic, tmin, asmax, gddae, gla, lap, lfl, lfn, light, tairh, yx, sgfun, gasfn, & !Input
                          assat, cvxty, intslp)                                                                !Output

   use ModuleDefs
   implicit none
   save

   integer :: dynamic, i, lfn
   logical :: light
   real :: asf, asfn, gasfn, tmin, ak, asmax, axx, axz, ayz, ck, cx, cxf, cxz, exponent, gddae, isf, tairh, xx, sgfun
   real, dimension(50) :: apk, asat, assat, cvty, cvxty, gla, insl, intslp, lap, lfl, yx

!----------------------------------------------------------------------
! Dynamic = runinit
!----------------------------------------------------------------------
   if (dynamic == runinit .or. dynamic == seasinit) then

! Initialize
      asf = 0.0
      asfn = 1.0
      ak = 0.0
      axx = 0.0
      axz = 0.0
      ayz = 0.0
      ck = 0.0
      cx = 0.0
      cxf = 0.0
      cxz = 0.0
      i = 1
      isf = 0.0
      xx = 0.0
      apk = 0.0
      asat = 0.0
      assat = 0.0
      cvty = 0.0
      cvxty = 0.0
      insl = 0.0
      intslp = 0.0

! JIL Parameters are calculated once a day and updated hourly with changing hour temperature

!-----------------------------------------------------------------------
! Dynamic = integr
!-----------------------------------------------------------------------
   else if (dynamic == rate) then

! JIL Calculating light response curve parameters for each leaf
      if (light) then
         light = .false.
         do i = 1, lfn
            asat(i) = 0.0
            cvty(i) = 0.0
            insl(i) = 0.0
            if (gla(i) > 0.0) then
               xx = lap(i)/yx(i)       !KAD: Age of leaf while it is expanding
               if (xx < 0.99) then      !Expanding leaf
                  ayz = 0.66
                  axx = 0.34
                  ak = 10.0
                  axz = 0.5
                  cx = 0.95
                  ck = 3.55
                  cxz = 0.186
                  insl(i) = 0.06         !KAD: No effect of leaf age on quantum efficiency during expansion
               else                     !Leaf maturing after complete expansion
                  if (apk(i) < 1.0) then
                     apk(i) = gddae
                  end if
                  ! CHP 11/27/2007 CHECK FOR LFL(I) = 0.
                  if (lfl(i) > 1.E-10) then
                     xx = (gddae - apk(i))/lfl(i)   !KAD: Age of leaf while it is senescing
                  else
                     xx = 1.0
                  end if

                  ayz = 0.18
                  axx = 0.85
                  ak = -7.0
                  axz = 0.47
                  axz = sgfun     !TT 08/06/2014
                  cx = 0.95
                  ck = -16.7
                  cxz = 0.88
                  !KAD: Effect of leaf age on quantum efficiency during leaf senescence
                  exponent = 8.0*(xx - 0.75)
                  if (exponent < -40.0) then
                     insl(i) = 0.06
                  else if (exponent > 40.0) then
                     insl(i) = 0.04
                  else
                     insl(i) = 0.04 + (0.02/(1.0 + exp(exponent)))
                  end if
               end if
               !Leaf age determined. Continue to apply its effect on asat and cvty

               asat(i) = asmax*(ayz + (axx/(1.0 + exp(-ak*(xx - axz)))))
               cvty(i) = cx/(1.0 + exp(-ck*(xx - cxz)))
               if (apk(i) > 1.0 .and. xx < 0.25) then
                  cx = 0.9
                  ck = 0.2
                  cvty(i) = cx + (ck*xx)
               end if
            end if    !gla(i) > 0.0 loop
         end do
      end if      !if(light) loop

! JIL Effect of T on light curve parameters curve parameters are defined at 30 C and scaled using hourly air temp
      asf = 0.0886 - 0.00867*tairh + 0.002840*tairh**2.0 - 0.00005070*tairh**3.0
      isf = 0.6783 + 0.02910*tairh - 0.000756*tairh**2.0 + 0.00000513*tairh**3.0
      cxf = 1.0108 - 0.00050*tairh - 0.000010*tairh**2.0 + 0.00000050*tairh**3.0

! Relative effect of night temperature < 10C on assat (Dwyer and Tollenaar, 1989 - CJPS 69,81; Ying et al. 2000 - FCR 68,87; Ying et al., 2002 - Crop Sci. 42,1164)
      if (tmin < 10.) then
         asfn = 1./exp(gasfn*(10.-tmin))
      else
         asfn = 1.
      end if

      do i = 1, lfn
         assat(i) = asat(i)*asf*asfn*1.0
         intslp(i) = insl(i)*isf*1.0
         cvxty(i) = cvty(i)*cxf*1.0
      end do

   end if       !Endif for DYNAMIC LOOP

   return

end subroutine MZ_TUE_PSParam

!==============================================================================================================================
! Variable definitions                                                                                  Unit
!------------------------------------------------------------------------------------------------------------------------------
! ak           Curvature parameter for the relationship between light-saturated assimilation rate
!              and leaf age
! apk(l)       Cumulative thermal time after emergence when expansion of leaf l was completed           degree-days
! asf          Relative effect of temperature on the light-saturated assimilation rate                  -
! asfn         Relative effect of night temperature on the light-saturated assimilation rate
! asat(l)      Light-saturated assimilation rate for leaf l (no temperature effect applied)             umol[CO2]/m2[leaf]/s
! asmax        Maximum instantaneous assimilation at 30 degC                                            �mol[CO2]/m2[leaf]/s
! assat(l)     Light-saturated assimilation rate (corresponds to the asymptote of the light             �mol[CO2]/m2[leaf]/s
!              response curve) for leaf l
! axx          "Upper end" of the sigmoid relatiohsip between light-saturated assimilation rate
!              and leaf age
! axz          Relative leaf age at which 50% of axx is reached
! ayz          Vertical offset of the relationship between light-saturated assimilation rate
!              and leaf age
! ck           Parameter controlling the shape of the sigmoid relationship between the light response
!              curve curvature parameter (cvxty) and leaf age
! cx           "Upper end" of the sigmoid relationship between the light response curve curvature
!              parameter (cvxty) and leaf age
! cxf          Relative effect of temperature on the curvature parameter of the light response curve    -
! cxz          Relative leaf age when 50% of cx is reached
! cvty(l)      Curvature parameter of light response curve for leaf l (no temperature effect applied)   -
! cvxty(l)     Ratio of the diffusion resistance to the total resistance to CO2 (curvature              -
!              parameter for the light response curve) for leaf l
! dynamic      Main control variable to tell each module which section of code to run                   -
! exponent     Intermediary variable
! gasfn        Genotype coefficient for effect of low (<10C) night temperatutres on leaf photosynthesis
!              (varies from 0.03 to 0.07 for newer to older hybrids)
! gddae        Cumulative growing degree days after emergence                                           degree-days
! gla(l)       Green leaf area for leaf l                                                               cm2
! i            Iteration counter of leaf number                                                         -
! insl(l)      Quantum efficiency of CO2 assimilation for leaf l (no temperature effect applied)        �mol[CO2]/�mol[quanta]
! intslp(l)    Quantum efficiency of CO2 assimilation (initial slope of the light response              �mol[CO2]/�mol[quanta]
!              curve) for leaf l
! isf          Relative effect of temperature on the quantum efficiency of CO2 assimilation             -
! lap(l)       Potential leaf area of leaf l from tip appearance to full expansion                      cm2
! lfl(l)       Longevity of leaf l                                                                      degree-days
! lfn          Total leaf number rounded up                                                             -
! light        Logical variable that differentiates daytime from nighttime
! tairh(h)     Hourly air temperature                                                                   degrees
! xx           Relative leaf age during expansion and senescence
! yx(l)        Maximum value of potential leaf area for leaf l (at full expansion)                      cm2
!==============================================================================================================================

!=======================================================================
!  MZ_AG_PSLeaf, Subroutine, J.I. Lizaso
!  Calculates gross photosynthesis (�mol CO2/m2 s) per unit leaf area as
!  a function of instantaneous PAR (�mol/m2 s)
!  Adapted from PGLEAF by K.J.Boote, J.W.Jones, G.Hoogenboom
!-----------------------------------------------------------------------
!  REVISION HISTORY
!  02/07/03 JIL Written
!  09/12/2007 JIL Adapted for IXIM model
!  08/28/2013 KAD Added variable definitions
!=======================================================================
subroutine MZ_TUE_PSLeaf(parlf, ast, isp, cvt, &  !Input
                         pglf)                                   !Output

   implicit none
   save
   real :: a, b, c, ast, cvt, isp, parlf, pglf

! Initialize
   a = 0.0
   b = 0.0
   c = 0.0
   pglf = 0.0

! JIL Using a non-rectangular hyperbolae (Thornley and Johnson, 1990)
! PGLF is average leaf gross assimilation (umol/m2 s)
! KAD: for details, see equations 39 and 40 in Lizaso et al. 2005 Agron. J. 97
   a = cvt
   b = (isp*parlf) + ast
   c = isp*parlf*ast

! JIL 09/17/09
   if (ast > 0.0 .and. a > 0.0) then
      pglf = (b - sqrt(b**2.0 - 4.0*a*c))/(2.0*a)
   else
      pglf = 0.0
   end if

   return

end subroutine MZ_TUE_PSLeaf

!==============================================================================================================================
! Variable definitions                                                                                  Unit
!------------------------------------------------------------------------------------------------------------------------------
! a, b and c   Parameters of a generalized non-rectangular hyperbolae
! ast          Light-saturated assimilation rate (corresponds to the asymptote of the light             �mol[CO2]/m2[leaf]/s
!              response curve)
! cvt          Ratio of the diffusion resistance to the total resistance to CO2 (curvature              -
!              parameter for the light response curve)
! isp          Quantum efficiency of CO2 assimilation (initial slope of the light response curve)       �mol[CO2]/�mol[quanta]
! parlf        Photosynthetically active radiation absorbed by leaf                                     �mol[quanta]/m2/s
! pglf         Instantaneous CO2 assimilation                                                           �mol[CO2]/m2[leaf]/s
!==============================================================================================================================

