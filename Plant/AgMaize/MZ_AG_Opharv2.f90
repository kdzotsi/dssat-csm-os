!=======================================================================
!  MZ_OPHARV, Subroutine C.H. Porter
!  Generates data for use by OPSUM and OVERVIEW for MAIZE.
!-----------------------------------------------------------------------
!  REVISION HISTORY
!  01/25/2002 CHP Written, based on OPHARV
!  03/03/2002 GH  Modified logic for calling of fileA
!  08/12/2003 CHP Added I/O error checking and changed call to READA
!  02/04/2005 CHP Added PODWT, LAI, HI to Summary.out output
!  08/11/2005 CHP Use BWAH as byproduct variable for Summary.OUT
!                 (byproduct harvested), BWAM as byproduct variable
!                 for Overview.OUT (byproduct produced).  Both variables
!                 now include senesced stover.
!  10/24/2005 CHP Added environmental & stress factors to Overview.OUT
!  07/13/2006 CHP Added P model
!  02/09/2007 GH  Add path for FileA
!  08/28/2009 CHP added EDAT, EDAP
!  10/01/2013 KAD Adapted for AgMaize
!
!  TODO:
!  - Is Overview.OUT's interpretation of SNAM 'Stem N at maturity' correct?
!    I think SNAM (observed value) is being compared to stover N (simulated)
!    instead, change simulated value to stem N?
!=======================================================================

subroutine MZ_AG_Opharv2(control, iswitch, & !Control inputs
                         yrplt, harvfrac, swfac, turfac, & !Inputs from MZ_AG_AGMAIZE
                         senesce, agefac, nstres, wtnvg, cannaa, wtnsd, wtncan, wtnup, pcngrn, & !Inputs from NPlant
                         growthStage, mdate, gstdyrdoySim, & !Inputs from Phenology
                         xhlai, tlu, tnleaf, greenla, greenlaSilk, & !Inputs from Leaf Area
                         tadrwAnth, tadrwSilk, kn, seedno, we, skerwt, stover, tadrw, grain, & !Inputs from Growth
                         bwah, sdwtah)                                                                 !Output

!-----------------------------------------------------------------------
   use ModuleDefs
   use MZ_AG_ModuleDefs
   implicit none
   save
!-----------------------------------------------------------------------
   real, parameter :: cm2tom2 = 1E-4
   integer ::  growthStage, run, timdif, lfmax, i, j, igst, iobsdat, ilist(5), glist(5), dsimdat(5), dobsdat(5)
   integer :: dynamic, isens, mdate, trtnum, yrdoy, yrsim, yrplt, trt_rot, gstdyrdoySim(20)
   real :: agefac, wtnvg, bwah, bwam, cannaa, tadrwAnth, tadrwSilk, wtnsd, kn, hi, StovSenes, maxlai, greenla(50)
   real :: greenlaSilk(50), nstres, we, psdwt, Pstres1, Pstres2, sdrate, sdwt, sdwtah, harvfrac(2), grain, tlu, tnleaf, nleaf
   real :: sdwtam, pltpop, laitop, seedno, skerwt, stover, swfac, tadrw, turfac, wtncan, wtnup, pcngrn, xhlai

! Arrays which contain data for printing in SUMMARY.OUT file
   integer, parameter :: sumnum = 18
   character*4, dimension(sumnum) :: label
   real, dimension(sumnum) :: value

! Arrays which contain Simulated and Measured data for printing in OVERVIEW.OUT and EVALUATE.OUT files (OPVIEW subroutine)
   integer acount
   character(len=6), dimension(EvaluateNum) :: olab, olap, X     !OLAP in dap
   character(len=8), dimension(EvaluateNum) :: Simulated, Measured
   character(len=50), dimension(EvaluateNum):: descrip

   type(ControlType) control
   type(SwitchType) iswitch
   type(ResidueType) senesce
   type(FileioType) :: datafileio

! Variables added for environmental and stress factors output
   Type(PlStresType) PlantStres

   character(len=1) :: rnmode, iplti, ideto, idets
   character(len=2) :: crop
   character(len=6), dimension(5) :: desclist
   character(len=10), dimension(20):: stname
   character(len=12) :: filea
   character(len=80) :: pathex

   dynamic = control%dynamic
   crop = control%crop
   yrsim = control%yrsim
   rnmode = control%rnmode
   run = control%run
   yrdoy = control%yrdoy
   iplti = iswitch%iplti
   ideto = iswitch%ideto
   idets = iswitch%idets
   StovSenes = senesce%ResWt(0)

!-----------------------------------------------------------------------
   acount = 26  !Number of FILEA headings.

! Headings in FILEA for Measured data
   data olab/ &
      'EDAT  ', & !1  Emergence date
      'ADAT  ', & !2  DANTH
      'LDAT  ', & !3  Silking date
      'KDAT  ', & !4  Milk line date
      'PD1T  ', & !5  IFPD
      'PDFT  ', & !6  IFSD
      'MDAT  ', & !7  DMAT Maturity date
      'HWAM  ', & !8  XGWT
      'PWAM  ', & !9  XPDW
      'H#AM  ', & !10 XNOGR
      'HWUM  ', & !11 XGWU
      'H#UM  ', & !12 XNOGU
      'CWAM  ', & !13 XCWT
      'BWAM  ', & !14 XSWT
      'LAIX  ', & !15 XLAM
      'LAIT  ', & !16 LAI of top 7 leaves
      'HIAM  ', & !17 XHIN
      'THAM  ', & !18 XTHR
      'GNAM  ', & !19 XNGR
      'CNAM  ', & !20 XNTP
      'SNAM  ', & !21 XNST
      'GN%M  ', & !22 XNPS
      'CWAA  ', & !23 XCWAA
      'CWAS  ', & !24 CWAS
      'CNAA  ', & !25 XCNAA
      'L#SM  ', & !26 XLFNO
      14*'      '/

!-----------------------------------------------------------------------
   data stname/ & !Stage
      'Germinate ', & !1
      'Emergence ', & !2
      'End Juv Ph', & !3
      'Tassel Ini', & !4
      'Topmost Lf', & !5
      'Anthesis  ', & !6
      'Silking   ', & !7
      'Ons Lin GF', & !8
      '50% Mk Lin', & !9
      'Black Layr', & !10
      '          ', & !11
      '          ', & !12
      'Planting  ', & !13
      'Start Sim ', & !14
      'End Sim   ', & !15
      'Harvest   ', & !16
      '          ', & !17
      '          ', & !18
      '          ', & !19
      'Harvest   '/     !20

   nleaf = min(tlu, tnleaf)

!***********************************************************************
!***********************************************************************
!     RUN INITIALIZATION
!***********************************************************************
   if (dynamic == runinit) then
!-----------------------------------------------------------------------

! Variables not currently calculated
      pstres1 = 1.0
      pstres2 = 1.0

!Read fileio case 'FILENP' and transfer variables
      call readfileio(control, 'FILENP', datafileio)
      filea = datafileio%filea
      pathex = datafileio%pathex

!Read fileio case 'INPUTS' and transfer variables
      call readfileio(control, 'INPUTS', datafileio)
      trtnum = datafileio%trtnum
      pltpop = datafileio%pltpop

! Assign descriptions to Measured and Simulated data from DATA.CDE.
      call getdesc(acount, olab, descrip)
      olap = olab

      call Opview(control, tadrw, acount, descrip, ideto, nleaf, &
                  measured, plantstres, simulated, gstdyrdoySim, &
                  stname, wtncan, xhlai, nint(grain), yrplt, growthStage)

!***********************************************************************
!***********************************************************************
!     SEASONAL INITIALIZATION
!***********************************************************************
   elseif (dynamic == seasinit) then
!-----------------------------------------------------------------------
      Simulated = ' '
      Measured = ' '

! Establish # and names of stages for environmental & stress summary
      Plantstres%active = .FALSE.
      Plantstres%nstages = 5

      PlantStres%StageName(0) = 'Planting to Harvest    '
      PlantStres%StageName(1) = 'Emegence-Tassel Init   '
      PlantStres%StageName(2) = 'Tassel Init-Topmost Lf '
      PlantStres%StageName(3) = 'Topmost leaf-Silking   '
      PlantStres%StageName(4) = 'Silking-Beg Lin Gr Fill'
      PlantStres%StageName(5) = 'Beg Lin Gr Fill-Blk Lyr'

      call Opview(control, tadrw, acount, descrip, ideto, nleaf, &
                  measured, plantstres, simulated, gstdyrdoySim, &
                  stname, wtncan, xhlai, nint(grain), yrplt, growthStage)

      maxlai = 0.0

!***********************************************************************
!***********************************************************************
!     DAILY OUTPUT
!***********************************************************************
   else if (dynamic == output) then
!-----------------------------------------------------------------------
      maxlai = amax1(maxlai, xhlai)      !Maximum lai season

      PlantStres%W_grow = turfac
      PlantStres%W_phot = swfac

      PlantStres%N_grow = nstres    !Growth
      PlantStres%N_phot = agefac    !Leaf area expansion

      PlantStres%P_grow = pstres2
      PlantStres%P_phot = pstres1
      PlantStres%active = .FALSE.

!if(growthStage > 0 .AND. growthStage < 6) then
!  Plantstres % active(growthStage) = .TRUE.
!endif
!KAD 10/08/2013- Do not want to add more phases to PlantStres so, here I am manually combining
!AgMaize's phenological phases to match the maximum number of 5 stages shown in OVERVIEW.OUT
      if (growthStage == 2 .OR. growthStage == 3) Plantstres%active(1) = .TRUE.    !Emergence to Tassel Initiation
      if (growthStage == 4) Plantstres%active(2) = .TRUE.    !Tassel Initiation to Appearance of topmost leaf
      if (growthStage == 5 .OR. growthStage == 6) Plantstres%active(3) = .TRUE.    !Appearance of topmost leaf to Silking
      if (growthStage == 7) Plantstres%active(4) = .TRUE.    !Silking to Onset of Linear Grain Filling
      if (growthStage == 8 .OR. growthStage == 9) Plantstres%active(5) = .TRUE.    !Onset of Linear Grain Filling to Black Layer

      if (yrdoy >= yrplt) then
         if (mdate < 0 .OR. (mdate > 0 .AND. yrdoy < mdate)) then
            PlantStres%active(0) = .TRUE.
         end if
      end if

! Send data to Overview.out data on days where stages occur
      call Opview(control, tadrw, acount, descrip, ideto, nleaf, &
                  measured, plantstres, simulated, gstdyrdoySim, &
                  stname, wtncan, xhlai, nint(grain), yrplt, growthStage)

!***********************************************************************
!***********************************************************************
!     Seasonal Output
!***********************************************************************
   else if (dynamic == seasend) then
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Calculate variables for output
! Update nitrogen and residue applications after routines have been
! modified to handle automatic management
!-----------------------------------------------------------------------
      if (seedno > 0.0) then
         psdwt = grain/(seedno*10.)
      else
         psdwt = 0.0
      end if
      if (tadrw > 0.0 .AND. grain >= 0.0) then
         hi = grain/tadrw
      else
         hi = 0.0
      end if

!LAI for the top 7 leaves
      lfmax = ceiling(tnleaf)                                           !Total number of actual leaves
      laitop = sum(greenlaSilk((lfmax - 6):lfmax))*cm2tom2*pltpop   !in m2/m2

!-----------------------------------------------------------------------
! Actual yield harvested (default is 100 %)
!-----------------------------------------------------------------------
      sdwt = grain/10.0
      sdwtam = sdwt
      sdwtah = sdwt*harvfrac(1)

!-----------------------------------------------------------------------
! Actual byproduct harvested (default is 0 %). Byproduct not harvested is incorporated
! 08/11/2005 Senesced leaf and stem stay on plant and are available for by-product harvest.
!-----------------------------------------------------------------------
      bwam = stover + stovsenes
      bwah = (stover + stovsenes)*harvfrac(2)

!-----------------------------------------------------------------------
      if ((ideto == 'Y' .OR. index('IAEBCGDT', rnmode) > 0) .OR. &
          (index('AY', idets) > 0 .AND. crop /= 'FA')) then
         if (index('FQ', rnmode) > 0) then
            trt_rot = control%rotnum
         else
            trt_rot = trtnum
         end if
         call reada(filea, pathex, olab, trt_rot, yrsim, x)

         !-----------------------------------------------------------------------
         ! Convert from YRDOY format to DAP.  Change descriptions to match.
         !Observed dates of phenological events
         ilist = [1, 2, 3, 4, 7]   !List of indices in olab for variables to show in Overview.OUT
         desclist = ['EDAP  ', 'ADAP  ', 'LDAP  ', 'KDAP  ', 'MDAP  ']  !Corresponding codes in DATA.CDE
         do j = 1, size(ilist)
            i = ilist(j)
            call reada_dates(x(i), yrsim, iobsdat)
            if (iobsdat > 0 .AND. iplti == 'R' .AND. isens == 0) then
               dobsdat(j) = timdif(yrplt, iobsdat)
            else
               dobsdat(j) = -99
            end if
            olap(i) = desclist(j)
            call getdesc(1, olap(i), descrip(i))
         end do

         !Simulated dates of phenological events
         !2=Emergence; 6=Anthesis; 7=Silking; 9=Milk line; 10=Black layer
         glist = [2, 6, 7, 9, 10]      !List of growth stages to show in Overview.OUT
         do j = 1, size(glist)
            igst = glist(j)
            dsimdat(j) = timdif(yrplt, gstdyrdoySim(igst))
            if (dsimdat(j) <= 0 .OR. yrplt <= 0) dsimdat(j) = -99
         end do

         write (Simulated(1), '(i8)') dsimdat(1); write (measured(1), '(i8)') dobsdat(1) !EDAT
         write (Simulated(2), '(i8)') dsimdat(2); write (measured(2), '(i8)') dobsdat(2) !ADAT
         write (Simulated(3), '(i8)') dsimdat(3); write (measured(3), '(i8)') dobsdat(3) !LDAT
         write (Simulated(4), '(i8)') dsimdat(4); write (measured(4), '(i8)') dobsdat(4) !KDAT
         write (Simulated(5), '(i8)') - 99; write (measured(5), '(i8)') - 99       !PD1T
         write (Simulated(6), '(i8)') - 99; write (measured(6), '(i8)') - 99       !PDFT
         write (Simulated(7), '(i8)') dsimdat(5); write (measured(7), '(i8)') dobsdat(5) !MDAT
         write (Simulated(8), '(i8)') nint(grain); write (measured(8), '(a8)') x(8)      !HWAM
         write (simulated(9), '(i8)') - 99; write (measured(9), '(i8)') - 99       !PWAM
         write (simulated(10), '(i8)') nint(seedno); write (measured(10), '(a8)') x(10)    !H#AM
         write (simulated(11), '(f8.4)') skerwt; write (measured(11), '(a8)') x(11)    !HWUM
         write (simulated(12), '(f8.1)') kn; write (measured(12), '(a8)') x(12)    !H#UM
         write (simulated(13), '(i8)') nint(tadrw); write (measured(13), '(a8)') x(13)     !CWAM
         write (Simulated(14), '(i8)') nint(bwam); write (measured(14), '(a8)') x(14)     !BWAM
         write (Simulated(15), '(f8.2)') maxlai; write (measured(15), '(a8)') x(15)     !LAIX
         write (Simulated(16), '(f8.2)') laitop; write (measured(16), '(i8)') - 99       !LAIT
         write (Simulated(17), '(f8.3)') hi; write (measured(17), '(a8)') x(17)     !HIAM
         write (Simulated(18), '(i8)') - 99; write (measured(18), '(i8)') - 99       !THAM
         write (Simulated(19), '(i8)') nint(wtnsd*10.); write (measured(19), '(a8)') x(19)     !GNAM
         write (Simulated(20), '(i8)') nint(wtncan*10.); write (measured(20), '(a8)') x(20)     !CNAM
         write (Simulated(21), '(i8)') nint(wtnvg*10.); write (measured(21), '(a8)') x(21)     !SNAM
         write (Simulated(22), '(f8.1)') pcngrn; write (measured(22), '(a8)') x(22)     !GN%M
         write (Simulated(23), '(i8)') nint(tadrwAnth); write (measured(23), '(a8)') x(23)     !CWAA
         write (Simulated(24), '(i8)') nint(tadrwSilk); write (measured(24), '(i8)') - 99       !CWAS
         write (Simulated(25), '(i8)') nint(cannaa*10); write (measured(25), '(a8)') x(25)     !CNAA
         write (Simulated(26), '(f8.2)') tnleaf; write (measured(26), '(a8)') x(26)     !L#SM
      end if

      call Opview(control, tadrw, acount, descrip, ideto, nleaf, &
                  measured, plantstres, simulated, gstdyrdoySim, &
                  stname, wtncan, xhlai, nint(grain), yrplt, growthStage)

!-------------------------------------------------------------------
! Send information to OPSUM to generate SUMMARY.OUT file
!-------------------------------------------------------------------
      psdwt = skerwt
      sdrate = -99.0

! Store Summary.out labels and values in arrays to send to
! OPSUM routines for printing.  Integers are temporarily
! saved as real numbers for placement in real array.
      label(1) = 'ADAT'; value(1) = float(gstdyrdoySim(6))
      label(2) = 'MDAT'; value(2) = float(mdate)
      label(3) = 'DWAP'; value(3) = sdrate
      label(4) = 'CWAM'; value(4) = tadrw
      label(5) = 'HWAM'; value(5) = sdwtam*10.
      label(6) = 'HWAH'; value(6) = sdwtah*10.
! BWAH multiplied by 10.0 in OPSUM - divide by 10. here to preserve units.
      label(7) = 'BWAH'; value(7) = bwah/10.
      label(8) = 'HWUM'; value(8) = psdwt       !*1000.
      label(9) = 'H#AM'; value(9) = seedno
      label(10) = 'H#UM'; value(10) = kn
      label(11) = 'NFXM'; value(11) = 0.0         !wtnfx*10.
      label(12) = 'NUCM'; value(12) = wtnup*10.
      label(13) = 'CNAM'; value(13) = wtncan*10.
      label(14) = 'GNAM'; value(14) = wtnsd*10.
      label(15) = 'PWAM'; value(15) = we
      label(16) = 'LAIX'; value(16) = maxlai
      label(17) = 'HIAM'; value(17) = hi
      label(18) = 'EDAT'; value(18) = float(gstdyrdoySim(2))

! Send labels and values to OPSUM
      call sumvals(sumnum, label, value)

! Send Measured and Simulated datat to OPSUM
      call evaluatedat(acount, measured, simulated, descrip, olap)

!***********************************************************************
!***********************************************************************
!     END OF DYNAMIC IF CONSTRUCT
!***********************************************************************
   end if
!***********************************************************************
   return
end subroutine MZ_AG_Opharv2
!=======================================================================

!==============================================================================================================================
! Variable definitions                                                                                  Unit
!------------------------------------------------------------------------------------------------------------------------------
! control       Constructed type for control variables
! agefac        Nitrogen stress factor affecting cell expansion
! cannaa        Stover N at anthesis, g N/M2
! tadrwAnth     Canopy weight at anthesis, kg/ha
! tadrwSilk     Canopy weight at silking, kg/ha
! kn            Grain number per plant, grains/plant
! seedno        Grain numbers, grains/m2
! harvfrac      Two-element array containing fractions of (1) yield harvested and (2) by-product harvested (fraction)
! ideto         Switch for printing overview.out file
! idets         Code to generate several output files
! isdate        Year and day of year for end of leaf growth
! growthStage   Growth stage
! mdate         Maturity data, year and day of year
! nstres        Nitrogen stress factor affecting growth (0-1)
! we            Pod (ear) weight, kg/ha
! seedno        Seed number, #/m2
! skerwt        Weight per kernel, g/kernel
! gstdyrdoySim(20)    Array storing the dates that different growth stages occurred
! stover        Stover weight (leaf+stem), kg/ha
! swfac         Soil water stress effect on growth (0-1), 1 is no stress, 0 is full stress
! tadrw         Total above ground biomass, kg/ha
! turfac        Soil water stress effecting cell expansion
! wtncan        Weight of nitrogen in above ground biomass (stem, leaf, grain), kg N/ha
! wtnup         Total nitrogen uptake, g/m2
! grain         Yield in kg/ha at 0% moisture content
! xhlai         Leaf area index, m2/m2
! yremrg        Year and day of year of emergence
! yrplt         Year and day of year of planting
!==============================================================================================================================
