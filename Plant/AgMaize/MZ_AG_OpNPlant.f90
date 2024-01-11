!=======================================================================
!  MZ_AG_OpNPlant, Subroutine
!-----------------------------------------------------------------------
!  Generates output file for daily plant nitrogen simulation
!-----------------------------------------------------------------------
!  REVISION       HISTORY
!  07/17/2014 KAD Written
!-----------------------------------------------------------------------
!  Called by: MZ_AG_NPlant
!  Calls:     None
!=======================================================================
subroutine MZ_AG_OpNPlant(control, iswitch, yrplt, mdate, pltpop, nstres, agefac, rcnp, tcnp, rmnc, tmnc, & !Input
                          pcnrt, pcnveg, pcnst, pcnl, pcngrn, wtnrt, wtnvg, wtnst, wtnlf, wtnsd, wtncan, rndem, tndem, & !Input
                          ndem, sln)                                                                                         !Input

!-----------------------------------------------------------------------
   use ModuleDefs
   use ModuleData
   implicit none
   save
!-----------------------------------------------------------------------

   character :: iswnit*1, crop*2, outn*12
   integer :: count, dap, das, doy, dynamic, errnum, frop, noutdn, run, timdif, year, yrdoy, mdate, yrplt
   real :: pltpop, ns1_av, ns2_av, nstres, agefac, rcnp, tcnp, rmnc, tmnc, pcnrt, pcnveg, pcnst, pcnl, pcngrn
   real :: wtnrt, wtnvg, wtnst, wtnlf, wtnsd, wtncan, rndem, tndem, ndem, sln
   logical fexist, first

   type(ControlType) control
   type(SwitchType) iswitch

! Transfer values from constructed data types into local variables
   dynamic = control%dynamic
   crop = control%crop
   das = control%das
   frop = control%frop
   run = control%run
   yrdoy = control%yrdoy
   iswnit = iswitch%iswnit

! No output for fallow crop
   if (crop == 'FA' .OR. iswnit == 'N') return

!----------------------------------------------------------------------
! DYNAMIC = RUNINIT
!----------------------------------------------------------------------
   if (dynamic == runinit) then
      ! Get unit number for the output file
      outn = 'PlantN.OUT'
      call getlun(outn, noutdn)

!----------------------------------------------------------------------
! DYNAMIC = SEASINIT
!----------------------------------------------------------------------
   else if (dynamic == seasinit) then
      ! Initialize daily growth output file
      inquire (file=outn, exist=fexist)
      if (fexist) then
         open (unit=noutdn, file=outn, status='old', iostat=errnum, position='append')
         first = .false.
      else
         open (unit=noutdn, file=outn, status='new', iostat=errnum)
         write (noutdn, '("*Plant Nitrogen Daily Output")')
         first = .true.
      end if

      ! Write file headers
      call header(seasinit, noutdn, run)

      ! Write variable headers
      write (noutdn, '(A5, 1X, A3, 2(1X,A5), 4(1X,A7), 5(1X,A7), 6(1X,A7), 7(1X,A7))') &
         '@YEAR', 'DOY', 'DAS', 'DAP', &
         'RN%CD', 'VN%CD', 'RN%MD', 'VN%MD', &
         'RN%D', 'VN%D', 'SN%D', 'LN%D', 'GN%D', &
         'RNAD', 'VNAD', 'SNAD', 'LNAD', 'GNAD', 'CNAD', &
         'RNDEMD', 'TNDEMD', 'NDEMD', 'NSTD', 'NST1A', 'NST2A', 'SLN'

      ns1_av = 0.0
      ns2_av = 0.0
      count = 0

!----------------------------------------------------------------------
! DYNAMIC = OUTPUT
!----------------------------------------------------------------------
   else if (dynamic == output) then
      ! Don't print prior to planting date
      if (yrdoy < yrplt .OR. yrplt < 0) return

      ! Compute average stress factors since last printout
      ns1_av = ns1_av + (1.0 - nstres)
      ns2_av = ns2_av + (1.0 - agefac)
      count = count + 1

      ! Check frequency of output
      if ((mod(das, frop) == 0) &     !Daily output every FROP days,
          .OR. (yrdoy == yrplt) &     !on planting date, and
          .OR. (yrdoy == mdate)) then     !at harvest maturity
         ! Print
         dap = max(0, timdif(yrplt, yrdoy))
         call yr_doy(yrdoy, year, doy)

         ! Compute average stress factors since last printout
         if (count > 0) then
            ns1_av = ns1_av/count
            ns2_av = ns2_av/count
            count = 0
         end if

         write (noutdn, '(1X,I4,1X,I3.3,2(1X,I5), 9(1X,F7.2), 6(1X,F7.1), 3(1X,F7.2), 4(1X,F7.3))') &
            year, doy, das, dap, &
            rcnp*100., tcnp*100., rmnc*100., tmnc*100., &
            pcnrt, pcnveg, pcnst, pcnl, pcngrn, &
            wtnrt*10., wtnvg*10., wtnst*10., wtnlf*10., wtnsd*10., wtncan*10., &
            rndem*pltpop*10., tndem*pltpop*10., ndem*pltpop*10., (1.0 - nstres), ns1_av, ns2_av, sln

         ! Set average stress factors since last printout back to zero
         ns1_av = 0.0
         ns2_av = 0.0

      end if

!----------------------------------------------------------------------
! DYNAMIC = SEASEND
!----------------------------------------------------------------------
   else if (dynamic == seasend) then

      !Close daily output file.
      close (noutdn)

!-----------------------------------------------------------------------
! END OF DYNAMIC IF STRUCTURE
!-----------------------------------------------------------------------
   end if

!-----------------------------------------------------------------------
   return
end subroutine MZ_AG_OpNPlant
!=======================================================================

