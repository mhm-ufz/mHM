!> \file mo_multi_param_reg.f90
!> \brief \copybrief mo_multi_param_reg
!> \details \copydetails mo_multi_param_reg

!> \brief Multiscale parameter regionalization (MPR).
!> \details This module provides the routines for multiscale parameter regionalization (MPR).
!> \authors Stephan Thober, Rohini Kumar
!> \date Dec 2012
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mpr
MODULE mo_multi_param_reg

  use mo_kind, only : i4, dp
  use mo_common_constants, only : nodata_dp, nodata_i4
  use mo_message, only : message, error_message

  implicit none

  private

  PUBLIC :: mpr                     ! calculates effective regionalised parameters
  PUBLIC :: canopy_intercept_param  ! estimate effective max. canopy interception

contains

  !> \brief Regionalizing and Upscaling process parameters
  !> \details calculating process parameters at L0 scale (Regionalization), like:
  !! - Baseflow recession parameter
  !! - Soil moisture parameters
  !! - PET correction for aspect
  !!
  !! and upscale these parameters to retrieve effective parameters at scale L1.
  !! Further parameter regionalizations are done for:
  !! - snow accumulation and melting parameters
  !! - threshold parameter for runoff generation on impervious layer
  !! - karstic percolation loss
  !! - setting up the Regionalized Routing Parameters
  !!
  !> \changelog
  !! - Stephan Thober           Jan 2013
  !!   - updated calling sequence for upscaling operators
  !! - Luis Samaniego           Feb 2013
  !!   - calling sequence, initial CHECK, call mpr_runoff
  !! - Stephan Thober           Feb 2013
  !!   - added subroutine for karstic percolation loss removed L1_, L0_ in variable names
  !! - Stephan Thober           Aug 2015
  !!   - moved regionalization of routing to mRM
  !! - Rohini Kumar             Mar 2016
  !!   - changes for handling multiple soil database options
  !! - Zink M. & Demirel M.C.   Mar 2017
  !!   - Added Jarvis soil water stress function at SM process(3)
  !! - Demirel M.C. & S. Stisen Apr 2017
  !!   - Added FC dependency on root fraction coefficient at SM process(3)
  !! - Robert Schweppe          Dec 2017
  !!   - added loop over LCscenes inside MPR, renamed variables rewrite
  !! - Robert Schweppe          Jun 2018
  !!   - refactoring and reformatting
  !! - Demirel M.C. & S. Stisen Jun 2020
  !!   - Added Feddes and global FC dependency on root fraction coefficient at SM process(3)=4
  !! - Rohini Kumar             Oct 2021
  !!   - Added Neutron count module to mHM integrate into develop branch (5.11.2)
  !! - Sebastian Müller         Mar 2023
  !!   - made L1_alpha, L1_kSlowFlow, L1_kBaseFlow and L1_kPerco land cover dependent
  !> \authors Stephan Thober, Rohini Kumar
  !> \date Dec 2012
subroutine mpr(mask0, geoUnit0, soilId0, Asp0, gridded_LAI0, LCover0, slope_emp0, y0, Id0, upper_bound1, lower_bound1, &
       left_bound1, right_bound1, n_subcells1, fSealed1, alpha1, degDayInc1, degDayMax1, degDayNoPre1, fAsp1, &
       HarSamCoeff1, PrieTayAlpha1, aeroResist1, surfResist1, fRoots1, kFastFlow1, kSlowFlow1, kBaseFlow1, &
       kPerco1, karstLoss1, soilMoistFC1, soilMoistSat1, soilMoistExp1, jarvis_thresh_c1, tempThresh1, &
       unsatThresh1, sealedThresh1, wiltingPoint1, maxInter1, petLAIcorFactor, &
       No_Count1, bulkDens1, latticeWater1, COSMICL31, &
       parameterset )

    use mo_common_variables, only : global_parameters, processMatrix
    use mo_mpr_SMhorizons, only : mpr_SMhorizons
    use mo_mpr_global_variables, only : HorizonDepth_mHM, fracSealed_CityArea, iFlag_soilDB, nSoilHorizons_mHM, &
         soilDB
    use mo_mpr_pet, only : bulksurface_resistance, pet_correctbyASP, pet_correctbyLAI, priestley_taylor_alpha
    use mo_mpr_runoff, only : mpr_runoff
    use mo_mpr_soilmoist, only : mpr_sm
    use mo_upscaling_operators, only : L0_fractionalCover_in_Lx, &
         upscale_arithmetic_mean
    use mo_mpr_neutrons,        only: mpr_neutrons
    implicit none

    !> mask at level 0 field
    logical, dimension(:, :), intent(in) :: mask0
    !> L0 geological units
    integer(i4), dimension(:), intent(in) :: geoUnit0
    !> soil Ids at level 0
    integer(i4), dimension(:, :), intent(in) :: soilId0
    !> [degree] Aspect at Level 0
    real(dp), dimension(:), intent(in) :: Asp0
    !> LAI grid at level 0, with dim2 = time
    real(dp), dimension(:, :), intent(in) :: gridded_LAI0
    !> land cover at level 0
    integer(i4), dimension(:, :), intent(in) :: LCOVER0
    !> Empirical quantiles of slope
    real(dp), dimension(:), intent(in) :: slope_emp0
    !> Cell ids at level 0
    integer(i4), dimension(:), intent(in) :: Id0
    !> Upper row of hi res block
    integer(i4), dimension(:), intent(in) :: upper_bound1
    !> Lower row of hi res block
    integer(i4), dimension(:), intent(in) :: lower_bound1
    !> Left column of hi res block
    integer(i4), dimension(:), intent(in) :: left_bound1
    !> Right column of hi res block
    integer(i4), dimension(:), intent(in) :: right_bound1
    !> Number of L0 cells within a L1 cell
    integer(i4), dimension(:), intent(in) :: n_subcells1
    !> y0 at level 0
    real(dp), dimension(:), intent(in) :: y0
    !> [1] fraction of sealed area
    real(dp), dimension(:, :, :), intent(inout) :: fSealed1
    !> Parameter that determines the rel. contribution to SM, upscal. Bulk den.
    real(dp), dimension(:, :, :), intent(inout) :: soilMoistExp1
    !> [1] jarvis critical value for norm SWC
    real(dp), dimension(:, :, :), intent(inout) :: jarvis_thresh_c1
    !> [10^-3 m] depth of saturated SM
    real(dp), dimension(:, :, :), intent(inout) :: soilMoistSat1
    !> [10^-3 m] field capacity
    real(dp), dimension(:, :, :), intent(inout) :: soilMoistFC1
    !> [10^-3 m] permanent wilting point
    real(dp), dimension(:, :, :), intent(inout) :: wiltingPoint1
    !> fraction of roots in soil horizon
    real(dp), dimension(:, :, :), intent(inout) :: fRoots1
    !> [degreeC] threshold temperature for snow rain
    real(dp), dimension(:, :, :), intent(inout) :: tempThresh1
    !> [mm-1 degreeC-1] Degree-day factor with no precipitation
    real(dp), dimension(:, :, :), intent(inout) :: degDayNoPre1
    !> [mm-1 degreeC-1] Maximum Degree-day factor
    real(dp), dimension(:, :, :), intent(inout) :: degDayMax1
    !> [d-1 degreeC-1]  Increase of the Degree-day factor per mm of increase in precipitation
    real(dp), dimension(:, :, :), intent(inout) :: degDayInc1
    !> [1]     PET correction for Aspect at level 1
    real(dp), dimension(:, :, :), intent(inout) :: fAsp1
    !> [1]     PET Hargreaves Samani coeff. at level 1
    real(dp), dimension(:, :, :), intent(inout) :: HarSamCoeff1
    !> [1]     PET Priestley Taylor coeff. at level 1
    real(dp), dimension(:, :, :), intent(inout) :: PrieTayAlpha1
    !> [s m-1] PET aerodynamical resitance at level 1
    real(dp), dimension(:, :, :), intent(inout) :: aeroResist1
    !> [s m-1] PET bulk surface resitance at level 1
    real(dp), dimension(:, :, :), intent(inout) :: surfResist1
    !> threshold parameter
    real(dp), dimension(:, :, :), intent(inout) :: sealedThresh1
    !> [10^-3 m] Threshhold water depth in upper reservoir (for Runoff contribution)
    real(dp), dimension(:, :, :), intent(inout) :: unsatThresh1
    !> [10^-3 m] Recession coefficient of the upper reservoir, upper outlet
    real(dp), dimension(:, :, :), intent(inout) :: kFastFlow1
    !> [10^-3 m] Recession coefficient of the upper reservoir, lower outlet
    real(dp), dimension(:, :, :), intent(inout) :: kSlowFlow1
    !> Level 1 baseflow recession
    real(dp), dimension(:, :, :), intent(inout) :: kBaseFlow1
    !> [1] Exponent for the upper reservoir
    real(dp), dimension(:, :, :), intent(inout) :: alpha1
    !> [d-1] percolation coefficient
    real(dp), dimension(:, :, :), intent(inout) :: kPerco1
    !> karstic percolation loss
    real(dp), dimension(:, :, :), intent(inout) :: karstLoss1
    !> max interception
    real(dp), dimension(:, :, :), intent(inout) :: maxInter1
    !> pet cor factor at level-1
    real(dp), dimension(:, :, :), intent(inout) :: petLAIcorFactor
    !> [-] inital neutron count
    real(dp), dimension(:, :, :), intent(inout) :: No_Count1
    !> [gcm-3] bulk density
    real(dp), dimension(:, :, :), intent(inout) :: bulkDens1
    !> [mm/mm] lattice water content
    real(dp), dimension(:, :, :), intent(inout) :: latticeWater1
    !> [-] cosmic L3 parameter
    real(dp), dimension(:, :, :), intent(inout) :: COSMICL31

    !> array of global parameters
    real(dp), dimension(:), intent(in), optional, target :: parameterset

    ! array of global parameters
    real(dp), dimension(:), pointer :: param

    real(dp), dimension(:, :, :), allocatable :: thetaS_till

    real(dp), dimension(:, :, :), allocatable :: thetaFC_till

    real(dp), dimension(:, :, :), allocatable :: thetaPW_till

    ! saturated hydraulic conductivity
    real(dp), dimension(:, :, :), allocatable :: Ks

    ! Bulk density
    real(dp), dimension(:, :, :), allocatable :: Db
    real(dp), dimension(:, :), allocatable :: thetaS
    real(dp), dimension(:, :), allocatable :: thetaFC
    real(dp), dimension(:, :), allocatable :: thetaPW

    ! neutron count
    real(dp), dimension(:,:,:), allocatable :: latWat_till
    real(dp), dimension(:,:,:), allocatable :: COSMIC_L3_till
    real(dp), dimension(:,:), allocatable   :: latWat         ! lattice water
    real(dp), dimension(:,:), allocatable   :: COSMIC_L3      ! COSMIC parameter L3


    ! relative variability of saturated
    ! hydraulic cound. for Horizantal flow
    real(dp), dimension(:), allocatable :: KsVar_H0

    ! relative variability of saturated
    ! hydraulic cound. for vertical flow
    real(dp), dimension(:), allocatable :: KsVar_V0

    ! soil mositure deficit from
    ! field cap. w.r.t to saturation
    real(dp), dimension(:), allocatable :: SMs_FC0

    ! L0 baseflow parameter
    real(dp), dimension(size(Id0, 1)) :: k2_0

    ! L1 baseflow parameter
    real(dp), dimension(:), allocatable :: k2_1

    ! L0 Aspect
    real(dp), dimension(size(Id0, 1)) :: fAsp0

    ! number of soil classes
    integer(i4) :: mSoil

    ! maximum of number of Tillage horizons
    integer(i4) :: mTill

    ! maximum number of horizons
    integer(i4) :: mHor

    ! number of Landcover classes
    integer(i4) :: mLC

    ! indexing of parameter vector - start
    integer(i4) :: iStart

    ! indexing of parameter vector - end
    integer(i4) :: iEnd

    ! 2nd indexing of parameter vector - start
    integer(i4) :: iStart2

    ! 2nd indexing of parameter vector - end
    integer(i4) :: iEnd2

    ! counter for looping over LCscenes
    integer(i4) :: iiLC

    ! [1]  Fraction of forest cover
    real(dp), dimension(size(fSealed1, dim = 1)) :: fForest1

    ! [1]  Fraction of permeable cover
    real(dp), dimension(size(fSealed1, dim = 1)) :: fPerm1


    if (present(parameterset)) then
       param => parameterset
    else
       param => global_parameters(:, 3)
    end if


    ! loop over all LCover scenes
    do iiLC = 1, size(LCover0, 2)

       ! estimate land cover fractions for dominant landcover class
       ! fSealed is intent inout, the rest only intent in
       fForest1(:) = L0_fractionalCover_in_Lx(LCover0(:, iiLC), 1, mask0, &
            upper_bound1, &
            lower_bound1, &
            left_bound1, &
            right_bound1, &
            n_subcells1)
       fSealed1(:, 1, iiLC) = L0_fractionalCover_in_Lx(LCover0(:, iiLC), 2, mask0, &
            upper_bound1, &
            lower_bound1, &
            left_bound1, &
            right_bound1, &
            n_subcells1)
       fPerm1(:) = L0_fractionalCover_in_Lx(LCover0(:, iiLC), 3, mask0, &
            upper_bound1, &
            lower_bound1, &
            left_bound1, &
            right_bound1, &
            n_subcells1)
       !---------------------------------------------------------
       ! Update fractions of sealed area fractions
       ! based on the sealing fraction[0-1] in cities
       !---------------------------------------------------------
       ! a factor is applied to the sealed area, effectively reducing it
       fSealed1(:, 1, iiLC) = fracSealed_CityArea * fSealed1(:, 1, iiLC)
       ! the forest area is kept constant, but the permeable area is increased so that the
       ! sum off all fractions equals 1 again
       fPerm1(:) = 1.0_dp - fSealed1(:, 1, iiLC) - fForest1(:)

       ! ------------------------------------------------------------------
       ! snow parameters
       ! ------------------------------------------------------------------
       select case(processMatrix(2,1))
       case(1)

          iStart = processMatrix(2, 3) - processMatrix(2, 2) + 1
          iEnd = processMatrix(2, 3)

          call snow_acc_melt_param(param(iStart : iEnd), & ! intent(in)
               fForest1, fSealed1(:, 1, iiLC), fPerm1, & ! intent(in)
               tempThresh1(:, 1, iiLC), degDayNoPre1(:, 1, iiLC), & ! intent(out)
               degDayInc1(:, 1, iiLC), degDayMax1(:, 1, iiLC) & ! intent(out)
               )
       case DEFAULT
          call error_message('***ERROR: Process description for process "snow pack" does not exist! mo_multi_param_reg')
       end select

       ! ------------------------------------------------------------------
       ! Soil moisture parametrization
       ! ------------------------------------------------------------------
       msoil = size(soilDB%is_present, 1)
       mLC = maxval(LCover0(:, iiLC), (LCover0(:, iiLC) .ne. nodata_i4))

       ! depending on which kind of soil database processing is to be performed
       if(iFlag_soilDB .eq. 0)then
          mtill = maxval(soilDB%nTillHorizons, (soilDB%nTillHorizons .ne. nodata_i4))
          mHor  = maxval(soilDB%nHorizons,     (soilDB%nHorizons     .ne. nodata_i4))
       else if(iFlag_soilDB .eq. 1) then
          ! here for each soil type both till and non-till soil hydraulic properties are to be estimated
          ! since a given soil type can lie in any horizon (till or non-till ones)
          ! adopt it in a way that it do not break the consistency of iFlag_soilDB = 0
          ! ** NOTE: SDB_nTillHorizons and SDB_nHorizons are also assigned in
          !          this flag option (see mo_soildatabase.f90 file - read_soil_LUT).
          !          But we are not using those variables here since in this case we have not
          !          varying number of soil horizons or either tillage horizons.
          !          So assigning them with a value = 1 is more than enough.
          mtill = 1
          mHor  = 1
       end if

       allocate(thetaS_till(msoil, mtill,  mLC))
       allocate(thetaFC_till(msoil, mtill, mLC))
       allocate(thetaPW_till(msoil, mtill, mLC))
       allocate(thetaS(msoil,  mHor))
       allocate(thetaFC(msoil, mHor))
       allocate(thetaPW(msoil, mHor))
       allocate(Ks(msoil, mHor, mLC))
       allocate(Db(msoil, mHor, mLC))

       ! neutron count related ones
       ! allocate and initalize here
       allocate(   latWat_till(msoil, mtill, mLC))
       allocate(COSMIC_L3_till(msoil, mtill, mLC))
       allocate(        latWat(msoil, mHor      ))
       allocate(     COSMIC_L3(msoil, mHor      ))
       latWat_till    = 0.000001_dp
       COSMIC_L3_till = 0.000001_dp
       COSMIC_L3      = 0.000001_dp
       latWat         = 0.000001_dp



       ! earlier these variables were allocated with  size(soilId0,1)
       ! in which the variable "soilId0" changes according to the iFlag_soilDB
       ! so better to use other variable which is common to both soilDB (0 AND 1) flags
       allocate(KsVar_H0(size(Id0, 1)))
       allocate(KsVar_V0(size(Id0, 1)))
       allocate( SMs_FC0(size(Id0, 1)))

       select case(processMatrix(3,1))
       case(1)
          ! first thirteen parameters go to this routine
          iStart = processMatrix(3, 3) - processMatrix(3, 2) + 1
          iEnd = processMatrix(3, 3) - 4

          ! next four parameters go here
          ! (the first three for the fRoots and the fourth one for the beta)
          iStart2 = processMatrix(3, 3) - 4 + 1
          iEnd2 = processMatrix(3, 3)

       case(2)
          ! first thirteen parameters go to this routine
          iStart = processMatrix(3, 3) - processMatrix(3, 2) + 1
          iEnd = processMatrix(3, 3) - 5

          ! next four parameters go here
          ! (the first three for the fRoots and the fourth one for the beta)
          iStart2 = processMatrix(3, 3) - 5 + 1
          iEnd2 = processMatrix(3, 3) - 1
          ! last parameter is jarvis parameter - no need to be regionalized
          jarvis_thresh_c1 = param(processMatrix(3, 3))
          !write(*,*) 'jarvis_thresh_c1 = ', jarvis_thresh_c1
          !write(*,*) 'iStart, iEnd, iStart2, iEnd2 = ', iStart, iEnd, iStart2, iEnd2

       case(3)
          ! first thirteen parameters go to this routine
          iStart = processMatrix(3, 3) - processMatrix(3, 2) + 1
          iEnd = processMatrix(3, 3) - 9

          ! next four parameters go here
          ! (the first three for the fRoots and the fourth one for the beta)
          iStart2 = processMatrix(3, 3) - 8
          iEnd2 = processMatrix(3, 3) - 1

          ! last parameter is jarvis parameter - no need to be regionalized
          jarvis_thresh_c1 = param(processMatrix(3, 3))

          !write(*,*) 'iStart, iEnd, iStart2, iEnd2 = ', iStart, iEnd, iStart2, iEnd2
          !write(*,*) 'jarvis_thresh_c1 = ', jarvis_thresh_c1

       case(4)
          ! first thirteen parameters go to this routine
          iStart = processMatrix(3, 3) - processMatrix(3, 2) + 1
          iEnd = processMatrix(3, 3) - 8

          ! next four parameters go here
          ! (the first three for the fRoots and the fourth one for the beta)
          iStart2 = processMatrix(3, 3) - 7
          iEnd2 = processMatrix(3, 3)
          !write(*,*) 'iStart, iEnd, iStart2, iEnd2 = ', iStart, iEnd, iStart2, iEnd2

       case DEFAULT
          call error_message('***ERROR: Process description for process "soil moisture parametrization"', &
               'does not exist! mo_multi_param_reg')
       end select

       call mpr_sm(param(iStart : iEnd), processMatrix, &
            soilDB%is_present, soilDB%nHorizons, soilDB%nTillHorizons, &
            soilDB%sand, soilDB%clay, soilDB%DbM, &
            Id0, soilId0, LCover0(:, iiLC), &
            thetaS_till, thetaFC_till, thetaPW_till, thetaS, &
            thetaFC, thetaPW, Ks, Db, KsVar_H0, KsVar_V0, SMs_FC0)

       ! >> neutron count related parameters
       if ( processMatrix(10,1) .GT. 0 ) &
            call mpr_neutrons( processMatrix(10,1), &  ! IN: processmatrix case
            param( processMatrix(10,3)-processMatrix(10,2)+1:processMatrix(10,3) ) , & ! IN:  global parameter set
            soilDB%is_present       , & ! IN:  flag indicating presence of soil
            soilDB%nHorizons        , & ! IN:  Number of Horizons of Soiltype
            soilDB%nTillHorizons    , & ! IN:  Number of tillage Horizons
            LCover0(:, iiLC)        , & ! IN:  land cover ids at level 0
            soilDB%clay             , & ! IN:  clay content
            soilDB%DbM              , & ! IN:  mineral Bulk density
            Db                      , & ! IN: Bulk density
            COSMIC_L3_till          , & ! OUT: COSMIC parameter L3 tillage layer
            latWat_till             , & ! OUT: COSMIC parameter Lattice Water tillage layer
            COSMIC_L3               , & ! OUT: COSMIC parameter L3
            latWat                    & ! OUT: COSMIC parameter Lattice Water
            )

       call mpr_SMhorizons(param(iStart2:iEnd2), processMatrix, &
            iFlag_soilDB, nSoilHorizons_mHM, HorizonDepth_mHM, &
            LCover0(:, iiLC), soilId0, &
            soilDB%nHorizons, soilDB%nTillHorizons, &
            thetaS_till, thetaFC_till, thetaPW_till, &
            thetaS, thetaFC, thetaPW, &
            soilDB%Wd, Db, soilDB%DbM, soilDB%RZdepth, &
            mask0, Id0, &
            upper_bound1, lower_bound1, left_bound1, right_bound1, n_subcells1, &
            soilMoistExp1(:, :, iiLC), soilMoistSat1(:, :, iiLC), soilMoistFC1(:, :, iiLC), &
            wiltingPoint1(:, :, iiLC), fRoots1(:, :, iiLC), &
            ! >>>>>> neutron count
            latWat_till, COSMIC_L3_till, latWat, COSMIC_L3, &
            bulkDens1(:,:,iiLC), latticeWater1(:,:,iiLC), COSMICL31(:,:,iiLC) &
            )

       deallocate(thetaS_till)
       deallocate(thetaFC_till)
       deallocate(thetaPW_till)
       deallocate(thetaS)
       deallocate(thetaFC)
       deallocate(thetaPW)
       deallocate(Ks)
       deallocate(Db)

       ! neutron count
       deallocate( latWat_till    )
       deallocate( COSMIC_L3_till )
       deallocate( latWat     )
       deallocate( COSMIC_L3  )

       ! ------------------------------------------------------------------
       ! potential evapotranspiration (PET)
       ! ------------------------------------------------------------------
       ! Penman-Monteith method is only method that is LCscene dependent
       if (processMatrix(5, 1) == 3) then
          iStart = processMatrix(5, 3) - processMatrix(5, 2) + 1
          iEnd = processMatrix(5, 3)
          call aerodynamical_resistance(gridded_LAI0, LCover0(:, iiLC), param(iStart : iEnd - 1), mask0, &
               Id0, n_subcells1, upper_bound1, lower_bound1, left_bound1, right_bound1, &
               aeroResist1(:, :, iiLC))
       else if (processMatrix(5, 1) == -1) then
          iStart = processMatrix(5, 3) - processMatrix(5, 2) + 1
          iEnd = processMatrix(5, 3)

          call pet_correctbyLAI(param(iStart : iEnd), nodata_dp, &
               LCover0(:, iiLC), gridded_LAI0, mask0, Id0, &
               upper_bound1, lower_bound1, left_bound1, &
               right_bound1, n_subcells1, petLAIcorFactor(:, :, iiLC))
       end if

       ! ------------------------------------------------------------------
       ! interflow
       ! ------------------------------------------------------------------
       select case(processMatrix(6, 1))
       case (1)
          !
          iStart = processMatrix(6, 3) - processMatrix(6, 2) + 1
          iEnd = processMatrix(6, 3)
          ! TODO: this subroutine should be split into each param (or at least extract kFastFlow1)
          ! because it is in the loop unnecessarily
          call mpr_runoff(LCover0(:, iiLC), mask0, SMs_FC0, slope_emp0, &
               KsVar_H0, param(iStart : iEnd), Id0, upper_bound1, lower_bound1, &
               left_bound1, right_bound1, n_subcells1, unsatThresh1(:, 1, 1), kFastFlow1(:, 1, iiLC), &
               kSlowFlow1(:, 1, iiLC), alpha1(:, 1, iiLC))
       case DEFAULT
          call error_message('***ERROR: Process description for process "interflow" does not exist! mo_multi_param_reg')
       END select

       ! ------------------------------------------------------------------
       ! percolation cofficient, karstic percolation loss
       ! ------------------------------------------------------------------
       select case(processMatrix(7, 1))
       case(1)

          iStart = processMatrix(7, 3) - processMatrix(7, 2) + 1
          iEnd = processMatrix(7, 3)
          call karstic_layer(& ! In
               param(iStart : iEnd), & ! In
               geoUnit0, mask0, & ! In
               SMs_FC0, KsVar_V0, Id0, & ! In
               n_subcells1, upper_bound1, lower_bound1, left_bound1, right_bound1, & ! In
               karstLoss1(:, 1, 1), kPerco1(:, 1, iiLC)                            & ! Out
               )

       case DEFAULT
          call error_message('***ERROR: Process description for process "percolation" does not exist! mo_multi_param_reg')
       end select

       deallocate(KsVar_H0)
       deallocate(KsVar_V0)
       deallocate(SMs_FC0)

    end do !! >>>>>>> LAND COVER SCENE LOOP


    ! ------------------------------------------------------------------
    ! sealed area threshold for runoff generation
    ! ------------------------------------------------------------------
    select case(processMatrix(4, 1))
    case (1)
       iStart = processMatrix(4, 3) - processMatrix(4, 2) + 1
       iEnd = processMatrix(4, 3)
       call iper_thres_runoff(param(iStart : iEnd), sealedThresh1)
    case DEFAULT
       call error_message('***ERROR: Process description for process "runoff_generation" does not exist! mo_multi_param_reg')
    end select

    ! ------------------------------------------------------------------
    ! potential evapotranspiration (PET)
    ! ------------------------------------------------------------------
    select case(processMatrix(5, 1))
    case(-1) ! LAI correction of input PET
       iEnd = -9999 ! dummy statement
    case(0) ! aspect correction of input PET
       iStart = processMatrix(5, 3) - processMatrix(5, 2) + 1
       iEnd = processMatrix(5, 3)
       call pet_correctbyASP(Id0, y0, Asp0, param(iStart : iEnd), nodata_dp, fAsp0)
       fAsp1(:, 1, 1) = upscale_arithmetic_mean(n_subcells1, upper_bound1, lower_bound1, &
            left_bound1, right_bound1, Id0, mask0, nodata_dp, fAsp0)
    case(1) ! Hargreaves-Samani method
       iStart = processMatrix(5, 3) - processMatrix(5, 2) + 1
       iEnd = processMatrix(5, 3)
       call pet_correctbyASP(Id0, y0, Asp0, param(iStart : iEnd - 1), nodata_dp, fAsp0)
       fAsp1(:, 1, 1) = upscale_arithmetic_mean(n_subcells1, upper_bound1, lower_bound1, &
            left_bound1, right_bound1, Id0, mask0, nodata_dp, fAsp0)
       HarSamCoeff1 = param(iEnd)
    case(2) ! Priestley-Taylor Method
       iStart = processMatrix(5, 3) - processMatrix(5, 2) + 1
       iEnd = processMatrix(5, 3)
       call priestley_taylor_alpha(gridded_LAI0, param(iStart : iEnd), &
            mask0, nodata_dp, Id0, n_subcells1, upper_bound1, lower_bound1, left_bound1, right_bound1, &
            PrieTayAlpha1(:, :, 1))
    case(3) ! Penman-Monteith method
       ! aerodynamic resistance is calculated inside LCscene loop
       iStart = processMatrix(5, 3) - processMatrix(5, 2) + 1
       iEnd = processMatrix(5, 3)
       call bulksurface_resistance(gridded_LAI0, param(iEnd), mask0, &
            nodata_dp, Id0, n_subcells1, upper_bound1, lower_bound1, left_bound1, right_bound1, &
            surfResist1(:, :, 1))
    case default
       call error_message('***ERROR: Process description for process "pet correction" does not exist! mo_multi_param_reg')
    end select
    ! ------------------------------------------------------------------
    ! baseflow recession parameter
    ! ------------------------------------------------------------------
    select case(processMatrix(9, 1))
      case(1)
        ! the number of process parameters, so the number in processMatrix(9,2) has
        ! to be equal to the size of geo_unit_list
        iStart = processMatrix(9, 3) - processMatrix(9, 2) + 1
        iEnd = processMatrix(9, 3)

        call baseflow_param(param(iStart : iEnd), geoUnit0, k2_0)

        ! Upscale by arithmetic mean
        allocate(k2_1(size(kBaseFlow1, 1)))
        k2_1 = upscale_arithmetic_mean(n_subcells1, upper_bound1, lower_bound1, &
            left_bound1, right_bound1, Id0, mask0, nodata_dp, k2_0)
        ! loop over all LCover scenes
        do iiLC = 1, size(LCover0, 2)
          kBaseFlow1(:, 1, iiLC) = k2_1
        end do
        deallocate(k2_1)

        ! correction and unit conversion
        ! if percolation is ON: correct K2 such that it is at least k1
        ! since kSlowFlow1 is LCover dependent, kBaseFlow1 is too
        if (processMatrix(7, 1) .gt. 0) kBaseFlow1 = merge(kSlowFlow1, kBaseFlow1, kBaseFlow1 .lt. kSlowFlow1)

      case DEFAULT
        call error_message('***ERROR: Process description for process "baseflow Recession" does not exist! mo_multi_param_reg')
    end select

    ! ------------------------------------------------------------------
    ! Neutron count related parameters
    ! >> only N0 parameter - others are defined above in soil parameters
    ! ------------------------------------------------------------------
    select case(processMatrix(10, 1))
    case(0)
       ! do nothing
    case(1)
       ! the number of process parameters, so the number in processMatrix(9,2) has
       iStart = processMatrix(10, 3) - processMatrix(10, 2) + 1
       iEnd = processMatrix(10, 3)
       No_Count1 = param(iStart)  ! >> 1st parameter --> N0 parameter
    case(2)
       ! the number of process parameters, so the number in processMatrix(9,2) has
       iStart = processMatrix(10, 3) - processMatrix(10, 2) + 1
       iEnd = processMatrix(10, 3)
       No_Count1 = param(iStart)  ! >> 1st parameter --> N0 parameter
    case DEFAULT
       call error_message('***ERROR: Process description for process "Neutron count" does not exist! mo_multi_param_reg')
    end select


    !-------------------------------------------------------------------
    ! call regionalization of parameters related to LAI
    ! it is now outside of mHM since LAI is now dynamic variable
    !-------------------------------------------------------------------
    call canopy_intercept_param(processMatrix, param(:), &
         gridded_LAI0, n_subcells1, upper_bound1, lower_bound1, left_bound1, right_bound1, Id0, mask0, &
         nodata_dp, maxInter1(:, :, 1))

  end subroutine mpr

  ! ----------------------------------------------------------------------------

  !    NAME
  !        baseflow_param

  !    PURPOSE
  !>       \brief baseflow recession parameter

  !>       \details This subroutine calculates the baseflow recession parameter
  !>       based on the geological units at the Level 0 scale. For each level 0
  !>       cell, it assigns the value specified in the parameter array param for the
  !>       geological unit in this cell.
  !>       Global parameters needed (see mhm_parameter.nml):
  !>       - param(1) = GeoParam(1,:)
  !>       - param(2) = GeoParam(2,:)
  !>       - ...

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: param"       list of required parameters
  !>       \param[in] "integer(i4), dimension(:) :: geoUnit0" ids of geological units at L0

  !    INTENT(OUT)
  !>       \param[out] "real(dp), dimension(:) :: k2_0" - baseflow recession parameter at Level 0

  !    HISTORY
  !>       \authors Stephan Thober, Rohini Kumar

  !>       \date Dec 2012

  ! Modifications:
  ! Stephan Thober  Dec 2013 - changed intent(inout) to intent(out)
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine baseflow_param(param, geoUnit0, k2_0)

    use mo_mpr_global_variables, only : GeoUnitList
    !$ use omp_lib

    implicit none

    ! list of required parameters
    real(dp), dimension(:), intent(in) :: param

    ! ids of geological units at L0
    integer(i4), dimension(:), intent(in) :: geoUnit0

    ! - baseflow recession parameter at Level 0
    real(dp), dimension(:), intent(out) :: k2_0

    ! loop variable
    integer(i4) :: ii

    ! geo unit
    integer(i4), dimension(1) :: gg


    if (size(param) .ne. size(geoUnitList)) &
          call error_message(' mo_multi_param_reg: baseflow_param: size mismatch, subroutine baseflow parameters ')

    k2_0 = nodata_dp

    !$OMP PARALLEL
    !$OMP DO PRIVATE(gg) SCHEDULE(STATIC)
    do ii = 1, size(k2_0)
       ! get parameter index in geoUnitList
       gg = minloc(abs(geoUnitList - geoUnit0(ii)))
       k2_0(ii) = param(gg(1))
    end do
    !$OMP END DO
    !$OMP END PARALLEL

  end subroutine baseflow_param

  ! ----------------------------------------------------------------------------

  !    NAME
  !        snow_acc_melt_param

  !    PURPOSE
  !>       \brief Calculates the snow parameters.

  !>       \details This subroutine calculates the snow parameters
  !>       threshold temperature (TT), degree-day factor without precipitation (DD)
  !>       and maximum degree-day factor (DDmax) as well as increase of degree-day
  !>       factor per mm of increase in precipitation (IDDP).

  !>       Global parameters needed (see mhm_parameter.nml):
  !>       - param(1) = snowTreshholdTemperature
  !>       - param(2) = degreeDayFactor_forest
  !>       - param(3) = degreeDayFactor_impervious
  !>       - param(4) = degreeDayFactor_pervious
  !>       - param(5) = increaseDegreeDayFactorByPrecip
  !>       - param(6) = maxDegreeDayFactor_forest
  !>       - param(7) = maxDegreeDayFactor_impervious
  !>       - param(8) = maxDegreeDayFactor_pervious
  !>       INTENT(IN)

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(8) :: param"    eight global parameters
  !>       \param[in] "real(dp), dimension(:) :: fForest1" [1] fraction of forest cover
  !>       \param[in] "real(dp), dimension(:) :: fIperm1"  [1] fraction of sealed area
  !>       \param[in] "real(dp), dimension(:) :: fPerm1"   [1] fraction of permeable area

  !    INTENT(OUT)
  !>       \param[out] "real(dp), dimension(:) :: tempThresh1"  [degreeC] threshold temperature for snow rain
  !>       \param[out] "real(dp), dimension(:) :: degDayNoPre1" [mm-1 degreeC-1] Degree-day factor with
  !>       \param[out] "real(dp), dimension(:) :: degDayInc1"   [d-1 degreeC-1]  Increase of the Degree-day
  !>       \param[out] "real(dp), dimension(:) :: degDayMax1"   [mm-1 degreeC-1] Maximum Degree-day factor

  !    HISTORY
  !>       \authors Stephan Thober, Rohini Kumar

  !>       \date Dec 2012

  ! Modifications:
  ! Juliane Mai    Oct 2013 - OLD parametrization --> param(1) = snowTreshholdTemperature
  !                                                --> param(2) = degreeDayFactor_forest
  !                                                --> param(3) = degreeDayFactor_impervious
  !                                                --> param(4) = degreeDayFactor_pervious
  !                                                --> param(5) = increaseDegreeDayFactorByPrecip
  !                                                --> param(6) = maxDegreeDayFactor_forest
  !                                                --> param(7) = maxDegreeDayFactor_impervious
  !                                                --> param(8) = maxDegreeDayFactor_pervious
  !                                             -------------------------------
  !                                             degreeDayFactor_impervious    = degreeDayFactor_forest + delta_1 + delta_2
  !                                             degreeDayFactor_pervious      = degreeDayFactor_forest + delta_1
  !                                             maxDegreeDayFactor_forest     = degreeDayFactor_forest + delta_3
  !                                             maxDegreeDayFactor_impervious = degreeDayFactor_impervious + delta_5
  !                                                                           = degreeDayFactor_forest + delta_1 + delta_2 + delta_5
  !                                             maxDegreeDayFactor_pervious   = degreeDayFactor_pervious + delta_4
  !                                                                           = degreeDayFactor_forest + delta_1 + delta_4
  !                                             -------------------------------
  !                                             NEW parametrization
  !                                                --> param(1) = snowTreshholdTemperature
  !                                                --> param(2) = degreeDayFactor_forest
  !                                                --> param(3) = delta_2
  !                                                --> param(4) = delta_1
  !                                                --> param(5) = increaseDegreeDayFactorByPrecip
  !                                                --> param(6) = delta_3
  !                                                --> param(7) = delta_5
  !                                                --> param(8) = delta_4
  ! Stephan Thober Dec 2013 - changed intent(inout) to intent(out)
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine snow_acc_melt_param(param, fForest1, fIperm1, fPerm1, tempThresh1, degDayNoPre1, degDayInc1, &
       degDayMax1)
    implicit none

    ! eight global parameters
    real(dp), dimension(8), intent(in) :: param

    ! [1] fraction of forest cover
    real(dp), dimension(:), intent(in) :: fForest1

    ! [1] fraction of sealed area
    real(dp), dimension(:), intent(in) :: fIperm1

    ! [1] fraction of permeable area
    real(dp), dimension(:), intent(in) :: fPerm1

    ! [degreeC] threshold temperature for snow rain
    real(dp), dimension(:), intent(out) :: tempThresh1

    ! [mm-1 degreeC-1] Degree-day factor with
    real(dp), dimension(:), intent(out) :: degDayNoPre1

    ! [mm-1 degreeC-1] Maximum Degree-day factor
    real(dp), dimension(:), intent(out) :: degDayMax1

    ! [d-1 degreeC-1]  Increase of the Degree-day
    real(dp), dimension(:), intent(out) :: degDayInc1

    real(dp) :: tmp_degreeDayFactor_forest, tmp_degreeDayFactor_impervious, tmp_degreeDayFactor_pervious

    real(dp) :: tmp_maxDegreeDayFactor_forest, tmp_maxDegreeDayFactor_impervious, tmp_maxDegreeDayFactor_pervious


    tmp_degreeDayFactor_forest = param(2)                                    ! OLD: param(2)
    tmp_degreeDayFactor_impervious = param(2) + param(4) + param(3)              ! OLD: param(3)
    tmp_degreeDayFactor_pervious = param(2) + param(4)                         ! OLD: param(4)
    tmp_maxDegreeDayFactor_forest = param(2) + param(6)   ! OLD: param(6)
    tmp_maxDegreeDayFactor_impervious = param(2) + param(4) + param(3) + param(7)   ! OLD: param(7)
    tmp_maxDegreeDayFactor_pervious = param(2) + param(4) + param(8)   ! OLD: param(8)

    tempThresh1 = param(1)
    degDayInc1 = param(5)

    degDayNoPre1 = (&
         tmp_degreeDayFactor_forest * fForest1 + &
         tmp_degreeDayFactor_impervious * fIperm1 + &
         tmp_degreeDayFactor_pervious * fPerm1)
    degDayMax1 = (&
         tmp_maxDegreeDayFactor_forest * fForest1 + &
         tmp_maxDegreeDayFactor_impervious * fIperm1 + &
         tmp_maxDegreeDayFactor_pervious * fPerm1)

  end subroutine snow_acc_melt_param

  ! ----------------------------------------------------------------------------

  !    NAME
  !        iper_thres_runoff

  !    PURPOSE
  !>       \brief sets the impervious layer threshold parameter for runoff generation

  !>       \details to be done by Kumar
  !>       ....

  !>       Global parameters needed (see mhm_parameter.nml):
  !>       - param(1) = imperviousStorageCapacity

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(1) :: param" - given threshold parameter

  !    INTENT(OUT)
  !>       \param[out] "real(dp), dimension(:, :, :) :: sealedThresh1"

  !    HISTORY
  !>       \authors Stephan Thober, Rohini Kumar

  !>       \date Dec 2012

  ! Modifications:
  ! Stephan Thober Dec 2013 - changed intent(inout) to intent(out)
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine iper_thres_runoff(param, sealedThresh1)
    implicit none

    ! - given threshold parameter
    real(dp), dimension(1), intent(in) :: param

    real(dp), dimension(:, :, :), intent(out) :: sealedThresh1


    sealedThresh1 = param(1)

  end subroutine iper_thres_runoff

  ! ----------------------------------------------------------------------------

  !    NAME
  !        karstic_layer

  !    PURPOSE
  !>       \brief calculates the Karstic percolation loss

  !>       \details This subroutine calls first the karstic_fraction upscaling
  !>       routine for determine the karstic fraction area for every Level 1
  !>       cell. Then, the karstic percolation loss is estimated given two
  !>       shape parameters by
  !>       \f[ karstLoss1 = 1 + ( fKarArea * param(1)) *( (-1)**INT(param(2),i4) ) \f]
  !>       where \f$ karstLoss1 \f$ is the karstic percolation loss and \f$ fKarArea \f$
  !>       is the fraction of karstic area at level 1
  !>       Global parameters needed (see mhm_parameter.nml):
  !>       - param(1) = rechargeCoefficient
  !>       - param(2) = rechargeFactor_karstic
  !>       - param(3) = gain_loss_GWreservoir_karstic

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(3) :: param"           parameters
  !>       \param[in] "integer(i4), dimension(:) :: geoUnit0"     id of the Karstic formation
  !>       \param[in] "logical, dimension(:, :) :: mask0"         mask at level 0
  !>       \param[in] "real(dp), dimension(:) :: SMs_FC0"         [-] soil mositure deficit from field
  !>       capacity w.r.t to saturation
  !>       \param[in] "real(dp), dimension(:) :: KsVar_V0"        [-] relative variability of saturated
  !>       \param[in] "integer(i4), dimension(:) :: Id0"          Cell ids of hi res field
  !>       \param[in] "integer(i4), dimension(:) :: n_subcells1"  number of l0 cells within a l1 cell
  !>       \param[in] "integer(i4), dimension(:) :: upper_bound1" upper row of a l1 cell in l0 grid
  !>       \param[in] "integer(i4), dimension(:) :: lower_bound1" lower row of a l1 cell in l0 grid
  !>       \param[in] "integer(i4), dimension(:) :: left_bound1"  left col of a l1 cell in l0 grid
  !>       \param[in] "integer(i4), dimension(:) :: right_bound1" right col of a l1 cell in l0 grid

  !    INTENT(OUT)
  !>       \param[out] "real(dp), dimension(:) :: karstLoss1" [-]    Karstic percolation loss
  !>       \param[out] "real(dp), dimension(:) :: L1_Kp"      [d-1] percolation coefficient

  !    HISTORY
  !>       \authors Rohini Kumar, Stephan Thober

  !>       \date Feb 2013

  ! Modifications:
  ! Stephan Thober Dec 2013 - changed intent(inout) to intent(out)
  ! Stephan Thober Dec 2013 - changed intent(inout) to intent(out)
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine karstic_layer(param, geoUnit0, mask0, SMs_FC0, KsVar_V0, Id0, n_subcells1, upper_bound1, lower_bound1, &
       left_bound1, right_bound1, karstLoss1, L1_Kp)

    use mo_mpr_global_variables, only : GeoUnitList, geoUnitKar
    use mo_upscaling_operators, only : L0_fractionalCover_in_Lx, upscale_arithmetic_mean
    !$ use omp_lib

    implicit none

    ! parameters
    real(dp), dimension(3), intent(in) :: param

    ! id of the Karstic formation
    integer(i4), dimension(:), intent(in) :: geoUnit0

    ! mask at level 0
    logical, dimension(:, :), intent(in) :: mask0

    ! [-] soil mositure deficit from field
    ! capacity w.r.t to saturation
    real(dp), dimension(:), intent(in) :: SMs_FC0

    ! [-] relative variability of saturated
    real(dp), dimension(:), intent(in) :: KsVar_V0

    ! Cell ids of hi res field
    integer(i4), dimension(:), intent(in) :: Id0

    ! number of l0 cells within a l1 cell
    integer(i4), dimension(:), intent(in) :: n_subcells1

    ! upper row of a l1 cell in l0 grid
    integer(i4), dimension(:), intent(in) :: upper_bound1

    ! lower row of a l1 cell in l0 grid
    integer(i4), dimension(:), intent(in) :: lower_bound1

    ! left col of a l1 cell in l0 grid
    integer(i4), dimension(:), intent(in) :: left_bound1

    ! right col of a l1 cell in l0 grid
    integer(i4), dimension(:), intent(in) :: right_bound1

    ! [-]    Karstic percolation loss
    real(dp), dimension(:), intent(out) :: karstLoss1

    ! [d-1] percolation coefficient
    real(dp), dimension(:), intent(out) :: L1_Kp

    ! fraction of karstic area
    real(dp), dimension(:), allocatable :: fKarArea

    ! temporal variable
    real(dp), dimension(size(SMs_FC0, 1)) :: tmp

    integer(i4) :: nGeoUnits

    integer(i4) :: i


    ! ------------------------------------------------------------------
    ! PERCOLATION; 1/Kp = f(Ks)
    ! ------------------------------------------------------------------
    ! Regionalise Kp with variability of last soil layer property
    !$OMP PARALLEL
    tmp = merge(param(1) * (1.0_dp + SMs_FC0) / (1.0_dp + KsVar_V0), &
         nodata_dp, Id0 .ne. nodata_i4)
    !$OMP END PARALLEL

    L1_Kp = upscale_arithmetic_mean(n_subcells1, upper_bound1, lower_bound1, &
         left_bound1, right_bound1, Id0, mask0, nodata_dp, tmp)

    ! minimum constrains
    L1_Kp = merge(2.0_dp, L1_Kp, L1_Kp .lt. 2.0_dp)

    nGeoUnits = size(geoUnitlist, 1)

    ! 1st calculate fraction of Karstic area
    allocate(fKarArea(size(karstLoss1, 1)))
    fKarArea = 0.0_dp

    do i = 1, nGeoUnits
       if(GeoUnitKar(i) .eq. 0) cycle
       fKarArea(:) = L0_fractionalCover_in_Lx(geoUnit0, geoUnitlist(i), mask0, &
            upper_bound1, lower_bound1, left_bound1, right_bound1, n_subcells1)
    end do

    ! 2nd calculate karstLoss1
    karstLoss1 = 1.0_dp - (fKarArea * param(2))

    deallocate(fKarArea)

  end subroutine karstic_layer

  ! ----------------------------------------------------------------------------

  !    NAME
  !        canopy_intercept_param

  !    PURPOSE
  !>       \brief estimate effective maximum interception capacity at L1

  !>       \details estimate effective maximum interception capacity at L1 for a given
  !>       Leaf Area Index field.
  !>       Global parameters needed (see mhm_parameter.nml):
  !>       Process Case 1:
  !>       - param(1) = canopyInterceptionFactor

  !    INTENT(IN)
  !>       \param[in] "integer(i4), dimension(:, :) :: processMatrix" indicate processes
  !>       \param[in] "real(dp), dimension(:) :: param"               array of global parameters
  !>       \param[in] "real(dp), dimension(:, :) :: LAI0"             LAI at level-0(nCells0, time)
  !>       \param[in] "integer(i4), dimension(:) :: n_subcells1"      Number of L0 cells within a L1 cell
  !>       \param[in] "integer(i4), dimension(:) :: upper_bound1"     Upper row of high resolution block
  !>       \param[in] "integer(i4), dimension(:) :: lower_bound1"     Lower row of high resolution block
  !>       \param[in] "integer(i4), dimension(:) :: left_bound1"      Left column of high resolution block
  !>       \param[in] "integer(i4), dimension(:) :: right_bound1"     Right column of high resolution block
  !>       \param[in] "integer(i4), dimension(:) :: Id0"              Cell ids at level 0
  !>       \param[in] "logical, dimension(:, :) :: mask0"             mask at level 0 field
  !>       \param[in] "real(dp) :: nodata"                            - nodata value

  !    INTENT(OUT)
  !>       \param[out] "real(dp), dimension(:, :) :: max_intercept1" max interception at level-1(nCells1, time)

  !    HISTORY
  !>       \authors Rohini Kumar

  !>       \date Aug. 2013

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine canopy_intercept_param(processMatrix, param, LAI0, n_subcells1, upper_bound1, lower_bound1, left_bound1, &
       right_bound1, Id0, mask0, nodata, max_intercept1)

    use mo_string_utils, only : num2str
    use mo_upscaling_operators, only : upscale_arithmetic_mean

    implicit none

    ! indicate processes
    integer(i4), dimension(:, :), intent(in) :: processMatrix

    ! array of global parameters
    real(dp), dimension(:), intent(in) :: param

    ! LAI at level-0(nCells0, time)
    real(dp), dimension(:, :), intent(in) :: LAI0

    ! Number of L0 cells within a L1 cell
    integer(i4), dimension(:), intent(in) :: n_subcells1

    ! Upper row of high resolution block
    integer(i4), dimension(:), intent(in) :: upper_bound1

    ! Lower row of high resolution block
    integer(i4), dimension(:), intent(in) :: lower_bound1

    ! Left column of high resolution block
    integer(i4), dimension(:), intent(in) :: left_bound1

    ! Right column of high resolution block
    integer(i4), dimension(:), intent(in) :: right_bound1

    ! Cell ids at level 0
    integer(i4), dimension(:), intent(in) :: Id0

    ! mask at level 0 field
    logical, dimension(:, :), intent(in) :: mask0

    ! - nodata value
    real(dp), intent(in) :: nodata

    ! max interception at level-1(nCells1, time)
    real(dp), dimension(:, :), intent(out) :: max_intercept1

    integer(i4) :: iStart, iEnd, it

    real(dp), dimension(:), allocatable :: max_intercept0

    real(dp), dimension(:), allocatable :: gamma_intercept


    ! ------------------------------------------------------------------
    ! Maximum interception parameter
    ! ------------------------------------------------------------------
    select case(processMatrix(1, 1))
    case(1)
       iStart = processMatrix(1, 3) - processMatrix(1, 2) + 1
       iEnd = processMatrix(1, 3)

       ! allocate space
       allocate(gamma_intercept(iEnd - iStart + 1))
       allocate(max_intercept0 (size(Id0, 1)))

       ! estimate max. intercept at Level-0
       gamma_intercept(:) = param(iStart : iEnd)

       do it = 1, size(LAI0, 2)
          !$OMP PARALLEL
          max_intercept0(:) = LAI0(:, it) * gamma_intercept(1)
          !$OMP END PARALLEL

          ! Upscale by arithmetic mean
          max_intercept1(:, it) = upscale_arithmetic_mean(n_subcells1, upper_bound1, lower_bound1, left_bound1, &
               right_bound1, Id0, mask0, nodata, max_intercept0(:))

       end do

       deallocate(gamma_intercept)
       deallocate(max_intercept0)
    CASE DEFAULT
       call error_message('mo_multi_param_reg: This processMatrix=', num2str(processMatrix(1, 1)), ' is not implemented!')
    end select

  end subroutine canopy_intercept_param


  ! ----------------------------------------------------------------------------

  !    NAME
  !        aerodynamical_resistance

  !    PURPOSE
  !>       \brief Regionalization of aerodynamic resistance

  !>       \details estimation of aerodynamical resistance
  !>       Global parameters needed (see mhm_parameter.nml):
  !>       - param(1) = canopyheigth_forest
  !>       - param(2) = canopyheigth_impervious
  !>       - param(3) = canopyheigth_pervious
  !>       - param(4) = displacementheight_coeff
  !>       - param(5) = roughnesslength_momentum_coeff
  !>       - param(6) = roughnesslength_heat_coeff

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:, :) :: LAI0"         LAI at level-0
  !>       \param[in] "integer(i4), dimension(:) :: LCover0"      land cover field
  !>       \param[in] "real(dp), dimension(6) :: param"           input parameter
  !>       \param[in] "logical, dimension(:, :) :: mask0"         mask at level 0
  !>       \param[in] "integer(i4), dimension(:) :: Id0"          Cell ids of hi res field
  !>       \param[in] "integer(i4), dimension(:) :: n_subcells1"  number of l0 cells within a l1 cell
  !>       \param[in] "integer(i4), dimension(:) :: upper_bound1" upper row of a l1 cell in l0 grid
  !>       \param[in] "integer(i4), dimension(:) :: lower_bound1" lower row of a l1 cell in l0 grid
  !>       \param[in] "integer(i4), dimension(:) :: left_bound1"  left col of a l1 cell in l0 grid
  !>       \param[in] "integer(i4), dimension(:) :: right_bound1" right col of a l1 cell in l0 grid

  !    INTENT(OUT)
  !>       \param[out] "real(dp), dimension(:, :) :: aerodyn_resistance1" aerodynmaical resistance

  !    HISTORY
  !>       \authors Matthias Zink

  !>       \date Apr 2013

  ! Modifications:
  ! Matthias Zink Jun 2017 - moved from mo_multi_scale_param_reg.f90 to mo_mpr_pet.f90
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine aerodynamical_resistance(LAI0, LCover0, param, mask0, Id0, n_subcells1, upper_bound1, lower_bound1, &
       left_bound1, right_bound1, aerodyn_resistance1)

    use mo_common_constants, only : eps_dp
    use mo_mpr_constants, only : WindMeasHeight, karman
    use mo_upscaling_operators, only : upscale_arithmetic_mean

    implicit none

    ! LAI at level-0
    real(dp), dimension(:, :), intent(in) :: LAI0

    ! land cover field
    integer(i4), dimension(:), intent(in) :: LCover0

    ! input parameter
    real(dp), dimension(6), intent(in) :: param

    ! mask at level 0
    logical, dimension(:, :), intent(in) :: mask0

    ! Cell ids of hi res field
    integer(i4), dimension(:), intent(in) :: Id0

    ! number of l0 cells within a l1 cell
    integer(i4), dimension(:), intent(in) :: n_subcells1

    ! upper row of a l1 cell in l0 grid
    integer(i4), dimension(:), intent(in) :: upper_bound1

    ! lower row of a l1 cell in l0 grid
    integer(i4), dimension(:), intent(in) :: lower_bound1

    ! left col of a l1 cell in l0 grid
    integer(i4), dimension(:), intent(in) :: left_bound1

    ! right col of a l1 cell in l0 grid
    integer(i4), dimension(:), intent(in) :: right_bound1

    ! aerodynmaical resistance
    real(dp), dimension(:, :), intent(out) :: aerodyn_resistance1

    integer(i4) :: tt

    real(dp), dimension(:), allocatable :: maxLAI

    real(dp), dimension(:), allocatable :: zm

    real(dp), dimension(:), allocatable :: canopy_height0

    real(dp), dimension(:), allocatable :: zm_zero, zh_zero, displace

    ! dim 1 = number of cells on level 0,
    ! dim2=month of year
    real(dp), dimension(:, :), allocatable :: aerodyn_resistance0


    ! initialize some things
    allocate(zm                  (size(LCover0, dim = 1))) ; zm = nodata_dp
    allocate(zm_zero             (size(LCover0, dim = 1))) ; zm_zero = nodata_dp
    allocate(zh_zero             (size(LCover0, dim = 1))) ; zh_zero = nodata_dp
    allocate(displace            (size(LCover0, dim = 1))) ; displace = nodata_dp
    allocate(canopy_height0      (size(LCover0, dim = 1))) ; canopy_height0 = nodata_dp
    allocate(aerodyn_resistance0 (size(LCover0, dim = 1), size(LAI0, 2))) ; aerodyn_resistance0 = nodata_dp
    allocate(maxLAI              (size(LCover0, dim = 1))) ; maxLAI = nodata_dp

    ! regionalization of canopy height
    ! substitute with canopy height
    canopy_height0 = merge(param(1), canopy_height0, LCover0 == 1)  ! forest
    canopy_height0 = merge(param(2), canopy_height0, LCover0 == 2)  ! impervious

    ! old implementation used values from LUT statically for all cells (Jan-Dec, 12 values):
    ! 7 Intensive-orchards: 2.0, 2.0,  2.0,  2.0,  3.0, 3.5,  4.0,  4.0,  4.0,  2.5,  2.0,  2.0
    maxLAI = MAXVAL(LAI0, dim=2)

    do tt = 1, size(LAI0, 2)

       ! pervious canopy height is scaled with LAI
       canopy_height0 = merge((param(3) * LAI0(:, tt) / maxLAI), canopy_height0, LCover0 == 3)  ! pervious

       ! estimation of the aerodynamic resistance on the lower level
       ! see FAO Irrigation and Drainage Paper No. 56 (p. 19 ff) for more information
       zm = WindMeasHeight
       ! correction: if wind measurement height is below canopy height loagarithm becomes negative
       zm = merge(canopy_height0 + zm, zm, ((abs(zm - nodata_dp) .GT. eps_dp) .AND. (zm .LT. canopy_height0)))

       ! zh       = zm
       displace = param(4) * canopy_height0
       zm_zero = param(5) * canopy_height0
       zh_zero = param(6) * zm_zero
       !
       ! calculate aerodynamic resistance (changes monthly)
       aerodyn_resistance0(:, tt) = log((zm - displace) / zm_zero) * log((zm - displace) / zh_zero) / (karman**2.0_dp)
       aerodyn_resistance1(:, tt) = upscale_arithmetic_mean(n_subcells1, upper_bound1, lower_bound1, &
            left_bound1, right_bound1, Id0, mask0, nodata_dp, aerodyn_resistance0(:, tt))

    end do

  end subroutine aerodynamical_resistance

END MODULE mo_multi_param_reg
