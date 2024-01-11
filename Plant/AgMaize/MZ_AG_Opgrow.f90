!======================================================================================
! MZ_AG_Opgrow
! Produces AgMaize's daily outputs
! Write the following files:
!   PlantGro.OUT  Daily growth output
!   LaindGro.OUT  Daily leaf-level output
!----------------------------------------------------------------------
! Revision History
! 04/17/2012 CHP/JZW Written (based on CSM's template)
! 08/23/2012 KAD General revisions, added leaf area outputs
!----------------------------------------------------------------------
! Called from: MZ_AG_AGMAIZE
! Calls      : None
!======================================================================================
subroutine MZ_AG_Opgrow(control, iswitch, soilprop, yrplt, mdate, growthStage, gti, tlu, & !Input
                        tnleafrm, tnleaf, anthes, silk, olg, dvs, gddcer, rla, pla, xhlai, latf, greenla, gddten, & !Input
                        gddapsim, lftips, lfcols, coltlu, tnleafphot, tnleaftemp, maxlarf, lflon, tluladur, & !Input
                        stLayer1, rlaSoil, rlaAir, pgross, maint, fr, fsts, flv, fvear, kn, grt, gsh, glv, & !Input
                        gst, egr, ggr, wrt, wlv, wlvtop, wlvbot, wsts, wstr, wst, we, grain, kw, tadrw, cgr, & !Input
                        cgrf, cgrw, eop, trwup, swfac, nstres, rtdep, rlv, srad, parday, parsunday, parshday, &
                        radtotday)

   use ModuleDefs
   use MZ_AG_ModuleDefs
   implicit none
   save

!----------------------------------------------------------------------
! Define variables
!----------------------------------------------------------------------
! Constructed variables
   type(ControlType), intent(in):: control
   type(SwitchType), intent(in):: iswitch
   type(SoilType), intent(in):: soilProp
   type(FileioType) :: datafileio

! Other variables
   character(len=15), parameter :: outaltmz = 'PlantGro.OUT', outlai = "LaindGro.OUT"
   integer, parameter :: lenla = 50
   character :: pathex*80
   logical :: fexist, fexistlai
   integer :: i, doy, dap, growthStage, noutmz, noutlai, errnum, errnumlai, year, yrplt, mdate, imax
   real :: gti, gddcer, gddten, gddapsim, tlu, tnleaf, lftips, lfcols, tnleafphot, tnleaftemp, anthes, silk, pgross, maint
   real :: greenla(lenla), lflon(lenla), olg, dvs, rla, pla, xhlai, maxlarf, latf, coltlu, tluladur
   real :: stLayer1, rlaSoil, rlaAir, fr, fsts, flv, fvear, kn, grt, gsh, glv, gst, egr, ggr, wrt, wlv
   real :: wlvtop, wlvbot, wsts, wstr, wst, we, grain, kw, tadrw, cgr, cgrf, cgrw
   real :: eop, trwup, swfac, rtdep, rlv(nl), dlayr(nl), cwad, gwad, lwad, swad, tnleafrm, nstres
   real :: srad, parday, parsunday, parshday, radtotday, lawd, pltpop
   character(len=4) :: lfhead(lenla)
!real, dimension(50) :: lamax, gddlabeg, gddlaend

! Functions/ parameters
   integer timdif
   character :: zero_blanks*30

! Constructed type variables
   character(len=1) :: idetg
   integer :: das, dynamic, run, yrdoy, frop

! Transfer values from constructed data types into local variables
   idetg = iswitch%idetg
   if (idetg == 'N') return
   das = control%das
   dynamic = control%dynamic
   run = control%run
   yrdoy = control%yrdoy
   frop = control%frop
   dlayr = soilprop%dlayr

! Create local variables that are compatible with DSSAT output names
   cwad = tadrw
   gwad = grain
   lwad = wlv
   swad = wst

!----------------------------------------------------------------------
! Dynamic = runinit
!----------------------------------------------------------------------
   if (dynamic == runinit) then

! Get unit number for the output file
      call getlun('outaltmz', noutmz)
      call getlun('outlai', noutlai)

      ! Read fileio case 'FILENP' and transfer variables
      call readfileio(control, 'ALLSEC', datafileio)
      pathex = datafileio%pathex
      pltpop = datafileio%pltpop

!----------------------------------------------------------------------
! Dynamic = seasinit
!----------------------------------------------------------------------
   else if (dynamic == seasinit) then
      lawd = 0.0
      imax = nint(tnleafrm)  !This is the maximum number of leaves for which output will be printed
      !lfhead(1:lenla) = paste('LA', lenla, len('LA'))
      !lfhead = repeat('LA',lenla)
      do i = 1, lenla
         write (lfhead(i), '(I2)') i
         lfhead(i) = 'LA'//lfhead(i)
         lfhead(i) = zero_blanks(lfhead(i))
      end do

      !Indices of leaf numbers to display in output file !Or [(i, i=1,ceiling(tnleafrm))]
      !lfind = [1:size(lfind)];  lfind = lfind(1:(nint(tnleafrm)+4))

      ! Open PlantGro.OUT as new or existing file
      inquire (file=outaltmz, exist=fexist)
      if (fexist) then
         open (unit=noutmz, file=outaltmz, status='old', iostat=errnum, position='append')
      else
         open (unit=noutmz, file=outaltmz, status='new', iostat=errnum)
         write (noutmz, '("*AgMAIZE DAILY OUTPUT FILE")')
      end if

      ! Open LaindGro.OUT as new or existing file
      inquire (file=outlai, exist=fexistlai)
      if (fexistlai) then
         open (unit=noutlai, file=outlai, status='old', iostat=errnumlai, position='append')
      else
         open (unit=noutlai, file=outlai, status='new', iostat=errnumlai)
         write (noutlai, '("*AgMAIZE DAILY LEAF-LEVEL OUTPUT FILE")')
      end if

      ! Generate file heading
      call header(seasinit, noutmz, run)
      call header(seasinit, noutlai, run)

      ! Write variable headers
      write (noutmz, '(A5, 1X,A3, 2(1X,A5), 1X,A4, 2(1X,A6), 12(1X,A6), 13(1X,A6), 30(1X,A6), 7(1X,A6))') &
         '@YEAR', 'DOY', 'DAS', 'DAP', 'GSTD', 'GTI', 'GDDCER', &
         'GDDTEN', 'GDDAPS', 'PHAD', 'MAINRE', 'TLU', 'NTIP', 'NLIG', 'TNLEAF', 'TNLFP', 'TNLFT', 'ANTHES', 'SILK', &
         'OLG', 'DVS', 'RLA', 'PLA', 'LAID', 'LADF', 'LATF', 'COLTLU', 'LGLF10', 'TLUDUR', 'RLASOL', 'RLAAIR', 'STLYR1', &
         'FR', 'FSTS', 'FLV', 'FVEAR', 'GWGD', 'GRT', 'GSH', 'GLV', 'GST', 'EGR', 'GGR', &
         'G#AD', 'RWAD', 'LWAD', 'LWADT', 'LWADB', 'SSAD', 'RSAD', 'SWAD', 'EWAD', 'GWAD', 'CWAD', 'CGRAD', 'CGRW', 'CGRAF', &
         'EOPD', 'TWUPD', 'WSPD', 'RDPD', 'RTLD', 'SRAD', 'PARDAY', 'PARDSL', 'PARDSH', 'RADTOT', 'NSTRES', 'LAWD'
      write (noutlai, '(A5, 1X,A3, 2(1X,A5), 50(1X,A6))') '@YEAR', 'DOY', 'DAS', 'DAP', lfhead(1:imax)

!----------------------------------------------------------------------
! Dynamic = output
!----------------------------------------------------------------------
   else if (dynamic == output) then
      dap = max(0, timdif(yrplt, yrdoy))
      if (dap > das) dap = 0
      if (lwad > 0) lawd = xhlai*1.E4/(lwad/10.)
      ! Write output based on user specified frequency:
      ! Use advance='no' to position file within the current record
      if ((mod(das, frop) == 0) &     !Daily output every FROP days,
          .OR. (yrdoy == yrplt) &     !on planting date, and
          .OR. (yrdoy == mdate)) then    !at harvest maturity

         call yr_doy(yrdoy, year, doy)
         write (noutmz, '(1X,I4, 1X,I3, 2(1X,I5), 1X,I4, 2(1X,F6.1), 2(1X,F6.1), &
    &                    2(1X,F6.1), 1X,F6.2, 2(1X,F6.2), 5(1X,F6.2)             &
    &                    8(1X,F6.2), 1X,F6.0, 1X,F6.2, 3(1X,F6.2),               &
    &                    5(1X,F6.2), 13(1X,F6.0),                                &
    &                    6(1X,F6.0), 5(1X,F6.2), 1X,F6.1,                        &
    &                    7(1X,F6.2))') &
         !<imax>(1X,F6.1),
         !5(1X,F6.2), 5(1X,F6.2))')
         !To plot in GBuild, use correctyear(year) instead of year, and change formatting to character
                    year, doy, das, dap, growthStage, gti, gddcer, gddten, gddapsim, &
                    pgross/10., maint/10., tlu, lftips, lfcols, tnleaf, tnleafphot, tnleaftemp, anthes, silk, &
                    olg, dvs, rla, pla, xhlai, maxlarf, latf, coltlu, lflon(10), tluladur, rlaSoil, rlaAir, stLayer1, &
                    fr, fsts, flv, fvear, kw, grt, gsh, glv, gst, egr, ggr, kn*pltpop, wrt, lwad, wlvtop, wlvbot, wsts, wstr, &
                    swad, we, gwad, cwad, cgr, cgrw, cgrf, eop, trwup*10., (1.-swfac), rtdep/100., sum(rlv*dlayr), &
                    srad, parday, parsunday, parshday, radtotday, (1.-nstres), lawd
         !lamax(1:imax),
         !gddlabeg(1:5), gddlaend(1:5)
         write (noutlai, '(1X,I4, 1X,I3, 2(1X,I5), 50(1X,F6.1))') year, doy, das, dap, greenla(1:imax)
      end if

!----------------------------------------------------------------------
! Dynamic = seasend
!----------------------------------------------------------------------
   else if (dynamic == seasend) then
      ! Close daily output file
      close (noutmz)
      close (noutlai)

!-----------------------------------------------------------------------
! End of dynamic if structure
!-----------------------------------------------------------------------
   end if
   return

!----------------------------------------------------------------------
! End of subroutine
!----------------------------------------------------------------------
end subroutine MZ_AG_Opgrow

!==============================================================================================================================
! Variable definitions                                                                                  Unit
!------------------------------------------------------------------------------------------------------------------------------
! anthes       Leaf stage at anthesis                                                                   tlu
! cgrw         Mean crop growth rate (aboveground dry matter) during the previous week                  kg[dm]/ha/day
! coltlu       Thermal leaf units from leaf tip appearance to full leaf expansion                       tlu
! control      Constructed type for control variables
! cwad               Total above ground dry weight                                                                kg[dm]/ha
! dap          Day after planting
! das          Days after start of simulation
! datafileio   Constructed type for variables read from DSSAT's input/output file (DSSAT45.INP)
! doy          Day of year
! dvs          Development stage (0 is planting, 1 is silking, and 2 is maturity)
! dynamic      Modular control (runinit=1 for run initialization, seasinit=2 for seasonal initialization,
!              rate=3 for rate calculations, integr=4 for integration of state variables, output=5 for
!              writing daily outputs, seasend=6 for closing output files
! egr               Ear growth rate                                                                                kg[dm]/ha/day
! errnum       Integer variable that informs on the data transfer status for the plant growth output file
! errnumlai    Integer variable that informs on the data transfer status for the detailed leaf area output file
! eop          Potential plant transpiration                                                            mm/day
! fexist       Logical variable for testing whether the plant growth output file exist
! fexistlai    Logical variable for testing whether the detailed leaf area output file exist
! flv          Proportion of assimilates partitioned to leaves                                          fraction
! fvear               Effect of stage of development on partitioning to non-grain ear part                     fraction
! fr           Proportion of assimilates partitioned to roots                                           fraction
! frop         Frequency of daily output                                                                days
! fsts         Proportion of assimilates partitioned to structural stems                                fraction
! gddapsim     Cumulative daily thermal time after planting using APSIM's approach                      degree-day
! gddcer       Cumulative daily thermal time after planting (CERES method)                              degree-day
! gddten       Cumulative daily thermal time after planting using GDD[10,30]                            degree-day
! ggr               Grain growth rate                                                                            kg[dm]/ha/day
! glv               Growth rate of leaves' dry matter                                                            kg[dm]/ha/day
! grain               Grain dry matter                                                                                kg[dm]/ha
! greenla(i)   Cumulative area of leaf i during expansion and senescence                                m2/leaf
! growthStage  Integer value of growth stage
! grt               Growth rate of root dry matter                                                                kg[dm]/ha/day
! gsh               Growth rate of vegetative components of shoot                                                kg[dm]/ha/day
! gst               Growth of stem dry matter                                                                    kg[dm]/ha/day
! gti          General Thermal Index heat unit accumulation                                             gti (Q.what is the real unit of gti?)
! gwad               Grain dry matter                                                                                kg[dm]/ha
! imax         Maximum number of leaves for which output will be printed
! idetg        Growth output switch
! iswitch      Constructed type for control switches
! kn           Kernel number                                                                            kernel/plant
! kw               Kernel weight                                                                                 mg/kernel
! lapot(i)     Cumulative area of leaf i during expansion (tip to collar) under non-stressed conditions m2/leaf
! latf         Effect of temperature on leaf growth
! lfcols       Leaf node position of leaf that has completed expansion
! lflon(i)     Longevity of leaf i                                                                      degree-day
! lwad               Total weight of leaves                                                                        kg[dm]/ha
! maint        Maintenance respiration                                                                        kg[glucose]/ha/day
! maxlarf      Reduction factor of maximum leaf area due to plant density
! mdate        Maturity date                                                                            yyyyddd
! noutmz       Integer unit number for plant growth output file
! noutlai      Integer unit number for detailed leaf area output file
! olg          Onset linear dry matter accumulation of grain                                            dvs
! outaltmz     Daily plant growth output file name
! outlai       Leaf-level daily leaf area output file name
! pgross       Daily canopy gross assimilation (glucose equivalent in kg/ha/day)                        kg[glucose]/ha/day
! pla          Plant leaf area                                                                          m2/plant
! rla          Rate of Leaf tip Appearence                                                              leaves/day
! rlaAir       Rate of leaf tip appearance using air temperatures                                       leaves/day
! rlaSoil      Rate of leaf tip appearance uing soil temperature of layer 1                             leaves/day
! rlv(l)       Root length density of layer l                                                           cm[root]/cm3[soil]
! rtdep        Rooting depth, initially set at emergence                                                cm
! silk         Leaf stage at silking                                                                    tlu
! stLayer1     Soil temperature of layer 1                                                              degree C
! swad               Weight of stems                                                                                kg[dm]/ha
! swfac        Soil water stress effect on growth (0-1), 1 is no stress, 0 is full
! tadrw               Total above ground dry weight                                                                kg[dm]/ha
! timdif       [Function] For computing difference between two dates                                    returns number of days
! tlu          Thermal leaf unit (cumulative rate of leaf tip appearance)                               tlu
! tluladur     Duration of leaf expansion of youngest leaf                                              tlu
! tnleaf       Total number of initiated leaves                                                         leaves
! tnleafphot   Change in leaf number due to photoperiod                                                   leaves
! tnleaftemp   Change in leaf number due to temperature                                                   leaves
! trwup        Total potential daily root water uptake                                                  cm/day
! we               Weight of ears                                                                                kg[dm]/ha
! wlv               Total weight of leaves                                                                        kg[dm]/ha
! wlvbot       Weight of leaves below top layer (wlv - wlvtop)                                          kg[dm]/ha
! wlvtop       Weight of leaves in top layer (lai <= 2)                                                 kg[dm]/ha
! wrt          Weight of roots                                                                          kg[dm]/ha
! wst               Weight of stems                                                                                kg[dm]/ha
! wstr               Weight of (temporarily-stored) stem reserves                                             kg[sucrose]/ha
! wsts               Weight of structural stem                                                                     kg[dm]/ha
! xhlai        Whole-plant leaf area index                                                              m2[leaf]/m2[ground]
! year         Four-digit year                                                                          yyyy
! yrdoy        Year day-of-year                                                                         yyyyddd
! yrplt        Planting date                                                                            yyyyddd
!==============================================================================================================================
