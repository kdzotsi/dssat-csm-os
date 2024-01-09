!======================================================================================
! MZ_AG_ModuleDefs
! Constructed variables used by AgMaize subroutines
!----------------------------------------------------------------------
! Revision History
! 09/13/2012 KAD Created
!----------------------------------------------------------------------
!======================================================================================

!module MZ_AG_ModuleData
!implicit none
!public
![Parameters moved to cultivar and species files]
!contains
!The following methods
!end module MZ_AG_ModuleData

!*****************************************************************************************
module MZ_AG_ModuleDefs
implicit none
public

!-- Type definition for dssatinp variables (FILEIO)
type FileioType
     character(len=6)  :: varno, econo 
     character(len=12) :: filea, files, filee, filec
     character(len=16) :: vrname
     character(len=80) :: pathex, pathsr, pather, pathcr
     integer :: isens, trtnum, yremrg
     real :: pltpop, rowspc, azir, sdepth
     real :: crmat, rlamx, taint, gfpyr, photp, lftop, lalfx, ampli, asymp, sgfun, knpot, pgrmn, gasfn     
end type FileioType

!-- Type definition for variables read from species file
type SpeciesType     
     real :: swcg                                     !Seed germination
     real :: tceil, torla, tbrla, mldvs, bldvs        !Phenology
     real :: xtemp(5), ygdd(5)                        !APSIM thermal time parameters
     real :: tempCoef(3,3)                            !Leaf area
     real :: asmax, canh, xc                          !Photosynthesis
     real :: pormin, rwumx, rlwr, rwuep1              !Roots
     real :: res30c, r30c2                            !Respiration
     real :: pcarlf, pcarst, pcarrt, pcarea, pcarsd   !Carbohydrate fractions in leaves, stems, roots, reproductive organs and grain
     real :: pprolf, pprost, pprort, pproea, pprosd   !Protein fractions in leaves, stems, roots, reproductive organs and grain
     real :: pliplf, plipst, pliprt, plipea, plipsd   !Lipid fractions in leaves, stems, roots, reproductive organs and grain
     real :: pliglf, pligst, pligrt, pligea, pligsd   !Lignin fractions in leaves, stems, roots, reproductive organs and grain
     real :: poalf, poast, poart, poaea, poasd        !Organic acid fractions in leaves, stems, roots, reproductive organs and grain
     real :: pminlf, pminst, pminrt, pminea, pminsd   !Mineral fractions in leaves, stems, roots, reproductive organs and grain
     real :: tmnc, tance, rcnp, rance, ctcnp1, ctcnp2 !Plant and root nitrogen concentrations and parameters
end type SpeciesType

!-- Type definition for variables read from the ecotype file
type EcotypeType
     character(len=6) :: ecotyp
     character(len=16) :: econam
     integer :: cday
     real :: tbase, topt, ropt, tsen
end type EcotypeType

!-- Type definition for growth variables 
type GrowthType
     real :: gsh, grt, glv, glvsen, gst, egr, ggr                                     !Growth rates in kg[DM]/ha/day
     real :: gshplt, grtplt, glvplt, glvsenplt, gstplt, egrplt, ggrplt                !Growth rates in g/plt/day
     real :: tdrw, tadrw, stover, wrt, wlv, wlvsen, wst, we, grain                    !Weights of plant components in kg[DM]/ha
     real :: tdrwplt, tadrwplt, stoverplt, wrtplt, wlvplt, wlvsenplt, wstplt,  &      !Weights of plant components in g/plt
             weplt, grainplt
     real :: wlvtop, wlvbot, wsts, wstr, weolg, tadrwAnth, tadrwSilk                  !Special weights in kg[DM]/ha
     real :: kn, kw, seedno                                                           !Grain components
     real :: ptf, slw                                                                 !Other
end type GrowthType
      
end module MZ_AG_ModuleDefs
      
      