!> \file mo_mpr_smhorizons.f90
!> \brief \copybrief mo_mpr_smhorizons
!> \details \copydetails mo_mpr_smhorizons

!> \brief setting up the soil moisture horizons
!> \details This module sets up the soil moisture horizons
!> \authors Stephan Thober, Rohini Kumar
!> \date Dec 2012
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mpr
module mo_mpr_SMhorizons

  use mo_kind, only : i4, dp
  use mo_common_constants, only : nodata_dp

  implicit none

  public :: mpr_SMhorizons

  private

contains

  ! ----------------------------------------------------------------------------

  !    NAME
  !        mpr_SMhorizons

  !    PURPOSE
  !>       \brief upscale soil moisture horizons

  !>       \details calculate soil properties at the level 1.
  !>       Global parameters needed (see mhm_parameter.nml):
  !>       - param(1) = rootFractionCoefficient_forest
  !>       - param(2) = rootFractionCoefficient_impervious
  !>       - param(3) = rootFractionCoefficient_pervious
  !>       - param(4) = infiltrationShapeFactor

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: param"               parameters
  !>       \param[in] "integer(i4), dimension(:, :) :: processMatrix" - matrix specifying user defined processes
  !>       \param[in] "integer(i4) :: iFlag_soil"                     - flags for handling multiple soil databases
  !>       \param[in] "integer(i4) :: nHorizons_mHM"                  - number of horizons to model
  !>       \param[in] "real(dp), dimension(:) :: HorizonDepth"        [10^-3 m] horizon depth from
  !>       surface, postive downwards
  !>       \param[in] "integer(i4), dimension(:) :: LCOVER0"          Land cover at level 0
  !>       \param[in] "integer(i4), dimension(:, :) :: soilID0"       soil ID at level 0
  !>       \param[in] "integer(i4), dimension(:) :: nHorizons"        horizons per soil type
  !>       \param[in] "integer(i4), dimension(:) :: nTillHorizons"    Number of Tillage horizons
  !>       \param[in] "real(dp), dimension(:, :, :) :: thetaS_till"   saturated water content of soil
  !>       horizons upto tillage depth,
  !>       f(OM, management)
  !>       \param[in] "real(dp), dimension(:, :, :) :: thetaFC_till"  Field capacity of tillage
  !>       layers; LUC dependent,
  !>       f(OM, management)
  !>       \param[in] "real(dp), dimension(:, :, :) :: thetaPW_till"  Permament wilting point of
  !>       tillage layers; LUC dependent,
  !>       f(OM, management)
  !>       \param[in] "real(dp), dimension(:, :) :: thetaS"           saturated water content of soil
  !>       horizons after tillage depth
  !>       \param[in] "real(dp), dimension(:, :) :: thetaFC"          Field capacity of deeper layers
  !>       \param[in] "real(dp), dimension(:, :) :: thetaPW"          Permanent wilting point of
  !>       deeper layers
  !>       \param[in] "real(dp), dimension(:, :, :) :: Wd"            weights of mHM Horizons
  !>       according to horizons provided
  !>       in soil database
  !>       \param[in] "real(dp), dimension(:, :, :) :: Db"            Bulk density
  !>       \param[in] "real(dp), dimension(:, :) :: DbM"              mineral Bulk density
  !>       \param[in] "real(dp), dimension(:) :: RZdepth"             [mm]       Total soil depth
  !>       \param[in] "logical, dimension(:, :) :: mask0"             mask at L0
  !>       \param[in] "integer(i4), dimension(:) :: cell_id0"         Cell ids of hi res field
  !>       \param[in] "integer(i4), dimension(:) :: upp_row_L1"       Upper row of hi res block
  !>       \param[in] "integer(i4), dimension(:) :: low_row_L1"       Lower row of hi res block
  !>       \param[in] "integer(i4), dimension(:) :: lef_col_L1"       Left column of hi res block
  !>       \param[in] "integer(i4), dimension(:) :: rig_col_L1"       Right column of hi res block
  !>       \param[in] "integer(i4), dimension(:) :: nL0_in_L1"        Number of L0 cells within a L1 cel

  !    INTENT(INOUT)
  !>       \param[inout] "real(dp), dimension(:, :) :: L1_beta"   Parameter that determines the
  !>       relative contribution to SM, upscaled
  !>       Bulk density
  !>       \param[inout] "real(dp), dimension(:, :) :: L1_SMs"    [10^-3 m] depth of saturated SM cont
  !>       \param[inout] "real(dp), dimension(:, :) :: L1_FC"     [10^-3 m] field capacity
  !>       \param[inout] "real(dp), dimension(:, :) :: L1_PW"     [10^-3 m] permanent wilting point
  !>       \param[inout] "real(dp), dimension(:, :) :: L1_fRoots" fraction of roots in soil horizons

  !    HISTORY
  !>       \authors Luis Samaniego, Rohini Kumar, Stephan Thober

  !>       \date Dec 2012

  ! Modifications:
  ! Stephan Thober                  Jan 2013 - updated calling sequence for upscaling operators
  ! Juliane Mai                     Oct 2013 - OLD parametrization
  !                                                --> param(1) = rootFractionCoefficient_forest
  !                                                --> param(2) = rootFractionCoefficient_impervious
  !                                                --> param(3) = rootFractionCoefficient_pervious
  !                                                --> param(4) = infiltrationShapeFactor
  !                                             -------------------------------
  !                                             rootFractionCoeff_perv   = rootFractionCoeff_forest - delta_1
  !                                             -------------------------------
  !                                             NEW parametrization
  !                                                --> param(1) = rootFractionCoefficient_forest
  !                                                --> param(2) = rootFractionCoefficient_impervious
  !                                                --> param(3) = delta_1
  !                                                --> param(4) = infiltrationShapeFactor
  !                                                ! if processMatrix(3,1) = 3 and 4 additionally
  !                                                --> param(5) = rootFractionCoefficient_sand
  !                                                --> param(6) = rootFractionCoefficient_clay
  !                                                --> param(7) = FCmin_glob
  !                                                --> param(8) = FCdelta_glob
  ! Stephan Thober                  Mar 2014 - added omp parallelization
  ! Rohini Kumar                    Mar 2016 - changes for handling multiple soil database options
  ! M. Cuneyd Demirel, Simon Stisen Apr 2017 - added FC dependency on root fraction coefficient
  ! Robert Schweppe Jun 2018 - refactoring and reformatting
  ! M. Cuneyd Demirel, Simon Stisen Jun 2020 - added Feddes and FC dependency on root fraction coefficient processCase(3) = 4
  ! M. Cuneyd Demirel, Simon Stisen Feb 2021 - Bug fix normalization of FCnorm

  subroutine mpr_SMhorizons(param, processMatrix, iFlag_soil, nHorizons_mHM, HorizonDepth, LCOVER0, soilID0, nHorizons, &
                           nTillHorizons, thetaS_till, thetaFC_till, thetaPW_till, thetaS, thetaFC, thetaPW, Wd, Db, &
                           DbM, RZdepth, mask0, cell_id0, upp_row_L1, low_row_L1, lef_col_L1, rig_col_L1, nL0_in_L1, &
                           L1_beta, L1_SMs, L1_FC, L1_PW, L1_fRoots, &
                           ! neutron count
                           latWat_till   , & ! lattice water upto tillage depth
                           COSMIC_L3_till, & ! COSMIC parameter L3 upto tillage depth
                           latWat        , & ! lattice water
                           COSMIC_L3     , & ! COSMIC paramter L3
                           L1_bulkDens   , & ! L bulk density
                           L1_latticeWater,& ! L1 lattice water content
                           L1_COSMICL3   )   ! L1 COSMIC L3 parameter from neutron module

    use mo_message, only : message, error_message
    use mo_string_utils, only : num2str
    use mo_upscaling_operators, only : upscale_harmonic_mean
    !$ use omp_lib

    implicit none

    ! parameters
    real(dp), dimension(:), intent(in) :: param

    ! - matrix specifying user defined processes
    integer(i4), dimension(:, :), intent(in) :: processMatrix

    ! - flags for handling multiple soil databases
    integer(i4), intent(in) :: iFlag_soil

    ! - number of horizons to model
    integer(i4), intent(in) :: nHorizons_mHM

    ! [10^-3 m] horizon depth from
    ! surface, postive downwards
    real(dp), dimension(:), intent(in) :: HorizonDepth

    ! Land cover at level 0
    integer(i4), dimension(:), intent(in) :: LCOVER0

    ! soil ID at level 0
    integer(i4), dimension(:, :), intent(in) :: soilID0

    ! horizons per soil type
    integer(i4), dimension(:), intent(in) :: nHorizons

    ! Number of Tillage horizons
    integer(i4), dimension(:), intent(in) :: nTillHorizons

    ! saturated water content of soil
    ! horizons upto tillage depth,
    ! f(OM, management)
    real(dp), dimension(:, :, :), intent(in) :: thetaS_till

    ! Field capacity of tillage
    ! layers; LUC dependent,
    ! f(OM, management)
    real(dp), dimension(:, :, :), intent(in) :: thetaFC_till

    ! Permament wilting point of
    ! tillage layers; LUC dependent,
    ! f(OM, management)
    real(dp), dimension(:, :, :), intent(in) :: thetaPW_till

    ! saturated water content of soil
    ! horizons after tillage depth
    real(dp), dimension(:, :), intent(in) :: thetaS

    ! Field capacity of deeper layers
    real(dp), dimension(:, :), intent(in) :: thetaFC

    ! Permanent wilting point of
    ! deeper layers
    real(dp), dimension(:, :), intent(in) :: thetaPW

    ! weights of mHM Horizons
    ! according to horizons provided
    ! in soil database
    real(dp), dimension(:, :, :), intent(in) :: Wd

    ! Bulk density
    real(dp), dimension(:, :, :), intent(in) :: Db

    ! mineral Bulk density
    real(dp), dimension(:, :), intent(in) :: DbM

    ! [mm]       Total soil depth
    real(dp), dimension(:), intent(in) :: RZdepth

    ! mask at L0
    logical, dimension(:, :), intent(in) :: mask0

    ! Cell ids of hi res field
    integer(i4), dimension(:), intent(in) :: cell_id0

    ! Upper row of hi res block
    integer(i4), dimension(:), intent(in) :: upp_row_L1

    ! Lower row of hi res block
    integer(i4), dimension(:), intent(in) :: low_row_L1

    ! Left column of hi res block
    integer(i4), dimension(:), intent(in) :: lef_col_L1

    ! Right column of hi res block
    integer(i4), dimension(:), intent(in) :: rig_col_L1

    ! Number of L0 cells within a L1 cel
    integer(i4), dimension(:), intent(in) :: nL0_in_L1

    ! Parameter that determines the
    ! relative contribution to SM, upscaled
    ! Bulk density
    real(dp), dimension(:, :), intent(inout) :: L1_beta

    ! [10^-3 m] depth of saturated SM cont
    real(dp), dimension(:, :), intent(inout) :: L1_SMs

    ! [10^-3 m] field capacity
    real(dp), dimension(:, :), intent(inout) :: L1_FC

    ! [10^-3 m] permanent wilting point
    real(dp), dimension(:, :), intent(inout) :: L1_PW

    ! fraction of roots in soil horizons
    real(dp), dimension(:, :), intent(inout) :: L1_fRoots

    ! neutron count
    real(dp), dimension(:,:,:), intent(in) :: latWat_till   ! lattice water
    real(dp), dimension(:,:,:), intent(in) :: COSMIC_L3_till! COSMIC parameter L3
    real(dp), dimension(:,:),   intent(in) :: latWat        ! lattice water
    real(dp), dimension(:,:),   intent(in) :: COSMIC_L3     ! COSMIC paramter L3
    ! out
    real(dp), dimension(:,:), intent(inout) :: L1_bulkDens
    real(dp), dimension(:,:), intent(inout) :: L1_latticeWater
    real(dp), dimension(:,:), intent(inout) :: L1_COSMICL3


    ! loop index
    integer(i4) :: h

    ! loop index
    integer(i4) :: k

    ! loop index
    integer(i4) :: l

    ! loop index
    integer(i4) :: s

    real(dp) :: dpth_f

    real(dp) :: dpth_t

    real(dp) :: fTotRoots

    ! beta 0
    real(dp), dimension(size(LCOVER0, 1)) :: beta0

    ! [10^3 kg/m3] Bulk density
    real(dp), dimension(size(LCOVER0, 1)) :: Bd0

    ! [10^-3 m] depth of saturated SM cont
    real(dp), dimension(size(LCOVER0, 1)) :: SMs0

    ! [10^-3 m] field capacity
    real(dp), dimension(size(LCOVER0, 1)) :: FC0

    ! [10^-3 m] permanent wilting point
    real(dp), dimension(size(LCOVER0, 1)) :: PW0

    ! neutron count
    real(dp), dimension(size(LCOVER0,1))    :: LW0     ! lattice water
    real(dp), dimension(size(LCOVER0,1))    :: L30     ! COSMIC parameter L3


    ! fraction of roots in soil horizons
    real(dp), dimension(size(LCOVER0, 1)) :: fRoots0

    real(dp) :: tmp_rootFractionCoefficient_forest

    real(dp) :: tmp_rootFractionCoefficient_impervious

    real(dp) :: tmp_rootFractionCoefficient_pervious

    ! Field capacity dependent
    ! root frac coeffiecient
    real(dp) :: tmp_rootFractionCoefficient_perviousFC

    ! Model parameter describing the threshold for
    ! actual ET reduction for sand
    real(dp) :: tmp_rootFractionCoefficient_sand

    ! Model parameter describing the threshold for actual
    ! ET reduction for clay
    real(dp) :: tmp_rootFractionCoefficient_clay

    real(dp) :: FCmin_glob

    ! real(dp) :: FCdelta_glob

    real(dp) :: FCmax_glob

    real(dp) :: FCnorm
    ! the minimum number of till horizons
    integer(i4) :: min_nTH


    min_nTH = minval(nTillHorizons(:))
    tmp_rootFractionCoefficient_forest     = param(1)     ! min(1.0_dp, param(2) + param(3) + param(1))
    tmp_rootFractionCoefficient_impervious = param(2)

    ! decide which parameterization should be used for route fraction:
    select case (processMatrix(3, 1))
    case(1,2)
       ! 1 and 2 - dependent on land cover
       tmp_rootFractionCoefficient_pervious = param(1) - param(3) ! min(1.0_dp, param(2) + param(3))
       !write(*,*) 'tmp_rootFractionCoefficient_forest     = ', tmp_rootFractionCoefficient_forest
       !write(*,*) 'tmp_rootFractionCoefficient_impervious = ', tmp_rootFractionCoefficient_impervious
       !write(*,*) 'tmp_rootFractionCoefficient_pervious   = ', tmp_rootFractionCoefficient_pervious
    case(3,4)
       ! 3 and 4 - dependent on land cover and additionally soil texture
       tmp_rootFractionCoefficient_pervious = param(3) ! min(1.0_dp, param(2) + param(3))
       !delta approach is used as in tmp_rootFractionCoefficient_pervious
       tmp_rootFractionCoefficient_sand = param(6) - param(5)
       !the value in parameter namelist is before substraction i.e. param(5)
       tmp_rootFractionCoefficient_clay = param(6)
       FCmin_glob=param(7)
       FCmax_glob=param(7)+param(8)
       !write(*,*) 'tmp_rootFractionCoefficient_forest     = ', tmp_rootFractionCoefficient_forest
       !write(*,*) 'tmp_rootFractionCoefficient_impervious = ', tmp_rootFractionCoefficient_impervious
       !write(*,*) 'tmp_rootFractionCoefficient_pervious   = ', tmp_rootFractionCoefficient_pervious
       !write(*,*) 'tmp_rootFractionCoefficient_sand       = ', tmp_rootFractionCoefficient_sand
       !write(*,*) 'tmp_rootFractionCoefficient_clay       = ', tmp_rootFractionCoefficient_clay
       !write(*,*) 'FCmin_glob = ', FCmin_glob
       !write(*,*) 'FCmax_glob = ', FCmax_glob
    end select


    ! select case according to a given soil database flag
    SELECT CASE(iFlag_soil)
      ! classical mHM soil database format
    CASE(0)
       do h = 1, nHorizons_mHM

          Bd0  = nodata_dp
          SMs0 = nodata_dp
          FC0  = nodata_dp
          PW0  = nodata_dp
          fRoots0 = nodata_dp
          tmp_rootFractionCoefficient_perviousFC = nodata_dp

          ! neutron count
          LW0 = nodata_dp
          L30 = nodata_dp

          ! Initalise mHM horizon depth
          ! Last layer depth is soil type dependent, and hence it assigned within the inner loop
          ! by default for the first soil layer
          dpth_f = 0.0_dp
          dpth_t = HorizonDepth(H)
          ! check for the layer (2, ... n-1 layers) update depth
          if( (H .gt. 1) .and. (H .lt. nHorizons_mHM) ) then
             dpth_f = HorizonDepth(H-1)
             dpth_t = HorizonDepth(H)
          end if

          !$OMP PARALLEL
          !$OMP DO PRIVATE( l, s ) SCHEDULE( STATIC )
          cellloop0 : do k = 1, size(LCOVER0, 1)
             l = LCOVER0(k)
             s = soilID0(k, 1)  ! >> in this case the second dimension of soilId0 = 1
             ! depth weightage bulk density
             Bd0(k) = sum(Db(s, : nTillHorizons(s), L) * Wd(S, H, 1 : nTillHorizons(S)), &
                  Wd(S, H, 1 : nTillHorizons(S)) > 0.0_dp) &
                  + sum(dbM(S, nTillHorizons(S) + 1 : nHorizons(S)) &
                  * Wd(S, H, nTillHorizons(S) + 1 : nHorizons(S)), &
                  Wd(S, H, nTillHorizons(S) + 1 : nHorizons(S)) >= 0.0_dp)
             ! depth weightage thetaS
             SMs0(k) = sum(thetaS_till(S, : nTillHorizons(s), L) &
                  * Wd(S, H, 1 : nTillHorizons(S)), &
                  Wd(S, H, 1 : nTillHorizons(S)) > 0.0_dp) &
                  + sum(thetaS(S, nTillHorizons(S) + 1 - min_nTH : nHorizons(s) - min_nTH) &
                  * Wd(S, H, nTillHorizons(S) + 1 : nHorizons(S)), &
                  Wd(S, H, nTillHorizons(S) + 1 : nHorizons(S)) > 0.0_dp)
             ! depth weightage FC
             FC0(k) = sum(thetaFC_till(S, : nTillHorizons(s), L) &
                  * Wd(S, H, 1 : nTillHorizons(S)), &
                  Wd(S, H, 1 : nTillHorizons(S)) > 0.0_dp) &
                  + sum(thetaFC(S, nTillHorizons(S) + 1 - min_nTH : nHorizons(s) - min_nTH) &
                  * Wd(S, H, nTillHorizons(S) + 1 : nHorizons(S)), &
                  Wd(S, H, nTillHorizons(S) + 1 : nHorizons(S)) > 0.0_dp)
             ! depth weightage PWP
             PW0(k) = sum(thetaPW_till(S, : nTillHorizons(s), L) &
                  * Wd(S, H, 1 : nTillHorizons(S)), &
                  Wd(S, H, 1 : nTillHorizons(S)) > 0.0_dp) &
                  + sum(thetaPW(S, nTillHorizons(S) + 1 - min_nTH : nHorizons(s) - min_nTH) &
                  * Wd(S, H, nTillHorizons(S) + 1 : nHorizons(S)), &
                  Wd(S, H, nTillHorizons(S) + 1 : nHorizons(S)) > 0.0_dp)

             ! neutron count --> depth weightage LW and L30
             LW0(k) = sum( latWat_till(S, : nTillHorizons(s), L) &
                  * Wd(S, H, 1 : nTillHorizons(S) ), &
                  Wd(S, H, 1 : nTillHorizons(S)) > 0.0_dp ) &
                  + sum( latWat(S,nTillHorizons(S) + 1 - min_nTH : nHorizons(s) - min_nTH) &
                  * Wd(S, H, nTillHorizons(S) + 1 : nHorizons(S)), &
                  Wd(S, H, nTillHorizons(S) + 1 : nHorizons(S)) > 0.0_dp )

             L30(k) = sum( COSMIC_L3_till(S, : nTillHorizons(s), L) &
                  * Wd(S, H, 1 : nTillHorizons(S) ), &
                  Wd(S, H, 1 : nTillHorizons(S)) > 0.0_dp ) &
                  + sum( COSMIC_L3(S,nTillHorizons(S) + 1 - min_nTH : nHorizons(s) - min_nTH) &
                  * Wd(S, H, nTillHorizons(S) + 1 : nHorizons(S)), &
                  Wd(S, H, nTillHorizons(S) + 1 : nHorizons(S)) > 0.0_dp )


             ! Horizon depths: last soil horizon is varying, and thus the depth
             ! of the horizon too...
             if(H .eq. nHorizons_mHM) then
                dpth_f = HorizonDepth(nHorizons_mHM - 1)
                dpth_t = RZdepth(S)
             end if
             ! other soil properties [SMs, FC, PWP in mm]
             SMs0(k) = SMs0(k) * (dpth_t - dpth_f)
             FC0(k)  = FC0(k)  * (dpth_t - dpth_f)
             PW0(k)  = PW0(k)  * (dpth_t - dpth_f)
             LW0(k)  = LW0(k)  * (dpth_t - dpth_f)
          end do cellloop0
          !$OMP END DO
          !$OMP END PARALLEL


          !$OMP PARALLEL
          !$OMP DO PRIVATE( l, tmp_rootFractionCoefficient_perviousFC, FCnorm ) SCHEDULE( STATIC )
          celllloop0 : do k = 1, size(LCOVER0, 1)
             l = LCOVER0(k)
             !---------------------------------------------------------------------
             ! Effective root fractions in soil horizon...
             !  as weightage sum (according to LC fraction)
             !---------------------------------------------------------------------
             ! vertical root distribution = f(LC), following asymptotic equation
             ! [for refrence see, Jackson et. al., Oecologia, 1996. 108(389-411)]

             ! Roots(H) = 1 - beta^d
             !  where,
             !   Roots(H) = cumulative root fraction [-], range: 0-1
             !   beta     = fitted extinction cofficient parameter [-], as a f(LC)
             !   d        = soil surface to depth [cm]

             ! NOTES **
             !  sum(fRoots) for soil horions = 1.0

             !  if [sum(fRoots) for soil horions < 1.0], then
             !  normalise fRoot distribution such that all roots end up
             !  in soil horizon and thus satisfying the constrain that
             !  sum(fRoots) = 1

             !  The above constrains means that there are not roots below the soil horizon.
             !  This may or may not be realistic but it has been coded here to satisfy the
             !  conditions of the EVT vales, otherwise which the EVT values would be lesser
             !  than the acutal EVT from whole soil layers.

             !  Code could be modified in a way that a portion of EVT comes from the soil layer
             !  which is in between unsaturated and saturated zone or if necessary the saturated
             !  layer (i.e. Groundwater layer) can also contribute to EVT. Note that the above
             !  modification should be done only if and only if [sum(fRoots) for soil horions < 1.0].
             !  In such cases, you have to judiciously decide which layers (either soil layer between
             !  unsaturated and saturated zone or saturated zone) will contribute to EVT and in which
             !  proportions. Also note that there are no obervations on the depth avialable ata a
             !  moment on these layers.
             !------------------------------------------------------------------------
             select case(L)
             case(1)
                ! forest
                fRoots0(k) = (1.0_dp - tmp_rootFractionCoefficient_forest**(dpth_t * 0.1_dp)) &
                     - (1.0_dp - tmp_rootFractionCoefficient_forest**(dpth_f * 0.1_dp))
             case(2)
                ! impervious
                fRoots0(k) = (1.0_dp - tmp_rootFractionCoefficient_impervious**(dpth_t * 0.1_dp)) &
                     - (1.0_dp - tmp_rootFractionCoefficient_impervious**(dpth_f * 0.1_dp))
             case(3)
                ! permeable
                select case (processMatrix(3, 1))
                case(1,2)
                   ! permeable
                   fRoots0(k) = (1.0_dp - tmp_rootFractionCoefficient_pervious**(dpth_t * 0.1_dp)) &
                        - (1.0_dp - tmp_rootFractionCoefficient_pervious**(dpth_f * 0.1_dp))
                case(3,4)
                   ! permeable
                   ! introducing global FC dependency on root frac. coef. by Simon Stisen and M. Cuneyd Demirel from GEUS.dk
                   ! The normalization is based on Demirel et al 2018 (doi: 10.5194/hess-22-1299-2018)
                   ! Case 3 is based on Jarvis (doi: 10.1016/0022-1694(89)90050-4)
                   ! Case 4 is based on Feddes (doi: 10.1016/0022-1694(76)90017-2)

                   FCnorm = (((FC0(k) / (dpth_t - dpth_f)) - FCmin_glob) / (FCmax_glob - FCmin_glob))

                   if(FCnorm .lt. 0.0_dp) then
                      ! print*, "FCnorm is below 0, will become 0", FCnorm
                      FCnorm=0.0_dp
                   else if(FCnorm .gt. 1.0_dp) then
                      ! print*, "FCnorm is above 1, will become 1", FCnorm
                      FCnorm=1.0_dp
                   end if

                   tmp_rootFractionCoefficient_perviousFC = (FCnorm * tmp_rootFractionCoefficient_clay) &
                        + ((1 - FCnorm) * tmp_rootFractionCoefficient_sand)

                   fRoots0(k) = (1.0_dp - tmp_rootFractionCoefficient_perviousFC**(dpth_t * 0.1_dp)) &
                        - (1.0_dp - tmp_rootFractionCoefficient_perviousFC**(dpth_f * 0.1_dp))

                end select

                if((fRoots0(k) .lt. 0.0_dp) .OR. (fRoots0(k) .gt. 1.0_dp)) then
                  ! why is this not stopping here?
                  call message('***ERROR: Fraction of roots out of range [0,1]. Cell', &
                        num2str(k), ' has value ', num2str(fRoots0(k)))
                  ! stop
                end if
             end select

          end do celllloop0
          !$OMP END DO
          !$OMP END PARALLEL

          beta0 = Bd0 * param(4)

          !---------------------------------------------
          ! Upscale the soil related parameters
          !---------------------------------------------
          L1_SMs(:, h) = upscale_harmonic_mean(nL0_in_L1, Upp_row_L1, Low_row_L1, &
               Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata_dp, SMs0)
          L1_beta(:, h) = upscale_harmonic_mean(nL0_in_L1, Upp_row_L1, Low_row_L1, &
               Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata_dp, beta0)
          L1_PW(:, h) = upscale_harmonic_mean(nL0_in_L1, Upp_row_L1, Low_row_L1, &
               Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata_dp, PW0)
          L1_FC(:, h) = upscale_harmonic_mean(nL0_in_L1, Upp_row_L1, Low_row_L1, &
               Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata_dp, FC0)
          L1_fRoots(:, h) = upscale_harmonic_mean(nL0_in_L1, Upp_row_L1, Low_row_L1, &
               Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata_dp, fRoots0)

          ! !> neutron count
          L1_bulkDens(:,h) = upscale_harmonic_mean( nL0_in_L1, Upp_row_L1, Low_row_L1, &
               Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata_dp, Bd0 )
          L1_latticeWater(:,h) = upscale_harmonic_mean( nL0_in_L1, Upp_row_L1, Low_row_L1, &
               Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata_dp, LW0 )
          L1_COSMICL3(:,h) = upscale_harmonic_mean( nL0_in_L1, Upp_row_L1, Low_row_L1, &
               Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata_dp, L30 )

       end do

      ! to handle multiple soil horizons with unique soil class
    CASE(1)
      ! horizon wise calculation
      do h = 1, nHorizons_mHM
        Bd0 = nodata_dp
        SMs0 = nodata_dp
        FC0 = nodata_dp
        PW0 = nodata_dp
        fRoots0 = nodata_dp
        tmp_rootFractionCoefficient_perviousFC = nodata_dp

        ! neutron count
        LW0     = nodata_dp
        L30     = nodata_dp

        ! initalise mHM horizon depth
        if (h .eq. 1) then
          dpth_f = 0.0_dp
          dpth_t = HorizonDepth(h)
          ! check for the layer (2, ... n-1 layers) update depth
        else
          dpth_f = HorizonDepth(h - 1)
          dpth_t = HorizonDepth(h)
        end if
        ! need to be done for every layer to get fRoots
        !$OMP PARALLEL
        !$OMP DO PRIVATE( l, s ) SCHEDULE( STATIC )
        cellloop1 : do k = 1, size(LCOVER0, 1)
          L = LCOVER0(k)
          s = soilID0(k, h)
          if (h .le. nTillHorizons(1)) then
            Bd0(k) = Db(s, 1, L)
            SMs0(k)= thetaS_till (s, 1, L) * (dpth_t - dpth_f)  ! in mm
            FC0(k) = thetaFC_till(s, 1, L) * (dpth_t - dpth_f)  ! in mm
            PW0(k) = thetaPW_till(s, 1, L) * (dpth_t - dpth_f)  ! in mm
            LW0(k) = latWat_till(s, 1, L)  * (dpth_t - dpth_f)  ! in mm  ! >> neutron count
          else
            Bd0(k) = DbM(s, 1)
            SMs0(k)= thetaS (s, 1) * (dpth_t - dpth_f)  ! in mm
            FC0(k) = thetaFC(s, 1) * (dpth_t - dpth_f)  ! in mm
            PW0(k) = thetaPW(s, 1) * (dpth_t - dpth_f)  ! in mm
            LW0(k) = latWat(s, 1)  * (dpth_t - dpth_f)  ! in mm  ! >> neutron count
          end if
        end do cellloop1
        !$OMP END DO
        !$OMP END PARALLEL


        !$OMP PARALLEL
        !$OMP DO PRIVATE( l, tmp_rootFractionCoefficient_perviousFC, FCnorm ) SCHEDULE( STATIC )

        celllloop1 : do k = 1, size(LCOVER0, 1)
          l = LCOVER0(k)
          !================================================================================
          ! fRoots = f[LC] --> (fRoots(H) = 1 - beta^d)
          ! see below for comments and references for the use of this simple equation
          ! NOTE that in this equation the unit of soil depth is in cm
          !================================================================================

          select case(L)
          case(1)
            ! forest
            fRoots0(k) = (1.0_dp - tmp_rootFractionCoefficient_forest**(dpth_t * 0.1_dp)) &
                    - (1.0_dp - tmp_rootFractionCoefficient_forest**(dpth_f * 0.1_dp))
          case(2)
            ! impervious
            fRoots0(k) = (1.0_dp - tmp_rootFractionCoefficient_impervious**(dpth_t * 0.1_dp)) &
                    - (1.0_dp - tmp_rootFractionCoefficient_impervious**(dpth_f * 0.1_dp))
          case(3)

            select case (processMatrix(3, 1))
            case(1,2)
              ! permeable
              fRoots0(k) = (1.0_dp - tmp_rootFractionCoefficient_pervious**(dpth_t * 0.1_dp)) &
                      - (1.0_dp - tmp_rootFractionCoefficient_pervious**(dpth_f * 0.1_dp))
            case(3,4)
              ! permeable
              ! introducing global FC dependency on root frac. coef. by Simon Stisen and M. Cuneyd Demirel from GEUS.dk
              ! The normalization is based on Demirel et al 2018 (doi: 10.5194/hess-22-1299-2018)
              ! Case 3 is based on Jarvis (doi: 10.1016/0022-1694(89)90050-4)
              ! Case 4 is based on Feddes (doi: 10.1016/0022-1694(76)90017-2)
              FCnorm = (((FC0(k) / (dpth_t - dpth_f)) - FCmin_glob) / (FCmax_glob - FCmin_glob))

              if(FCnorm .lt. 0.0_dp) then
                ! print*, "FCnorm is below 0, will become 0", FCnorm
                FCnorm=0.0_dp
              else if(FCnorm .gt. 1.0_dp) then
                ! print*, "FCnorm is above 1, will become 1", FCnorm
                FCnorm=1.0_dp
              end if

              tmp_rootFractionCoefficient_perviousFC = (FCnorm * tmp_rootFractionCoefficient_clay) &
                      + ((1 - FCnorm) * tmp_rootFractionCoefficient_sand)


              fRoots0(k) = (1.0_dp - tmp_rootFractionCoefficient_perviousFC**(dpth_t * 0.1_dp)) &
                      - (1.0_dp - tmp_rootFractionCoefficient_perviousFC**(dpth_f * 0.1_dp))

            end select

            if((fRoots0(k) .lt. 0.0_dp) .OR. (fRoots0(k) .gt. 1.0_dp)) then
              ! why is this not stopping here?
              call message('***ERROR: Fraction of roots out of range [0,1]. Cell', &
                    num2str(k), ' has value ', num2str(fRoots0(k)))
              ! stop
            end if
          end select

        end do celllloop1
        !$OMP END DO
        !$OMP END PARALLEL

        ! beta parameter
        beta0 = Bd0 * param(4)

        ! Upscale the soil related parameters
        L1_SMs(:, h) = upscale_harmonic_mean(nL0_in_L1, Upp_row_L1, Low_row_L1, &
                Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata_dp, SMs0)
        L1_beta(:, h) = upscale_harmonic_mean(nL0_in_L1, Upp_row_L1, Low_row_L1, &
                Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata_dp, beta0)
        L1_PW(:, h) = upscale_harmonic_mean(nL0_in_L1, Upp_row_L1, Low_row_L1, &
                Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata_dp, PW0)
        L1_FC(:, h) = upscale_harmonic_mean(nL0_in_L1, Upp_row_L1, Low_row_L1, &
                Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata_dp, FC0)
        L1_fRoots(:, h) = upscale_harmonic_mean(nL0_in_L1, Upp_row_L1, Low_row_L1, &
             Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata_dp, fRoots0)

        ! neutron count
        L1_bulkDens(:,h) = upscale_harmonic_mean( nL0_in_L1, Upp_row_L1, Low_row_L1, &
             Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata_dp, Bd0 )
       L1_latticeWater(:,h)   = upscale_harmonic_mean( nL0_in_L1, Upp_row_L1, Low_row_L1, &
             Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata_dp, LW0  )
        L1_COSMICL3(:,h)   = upscale_harmonic_mean( nL0_in_L1, Upp_row_L1, Low_row_L1, &
             Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata_dp, L30  )

      end do
      ! anything else
    CASE DEFAULT
      call error_message('***ERROR: iFlag_soilDB option given does not exist. Only 0 and 1 is taken at the moment.')
    END SELECT


    ! below operations are common to all soil databases flags
    !$OMP PARALLEL
    !------------------------------------------------------------------------
    ! CHECK LIMITS OF PARAMETERS
    !   [PW <= FC <= ThetaS]
    !  units of all variables are now in [mm]
    ! If by any means voilation of this rule appears (e.g. numerical errors)
    ! than correct it --> threshold limit = 1% of the upper ones
    !------------------------------------------------------------------------
    L1_FC = merge(L1_SMs - 0.01_dp * L1_SMs, L1_FC, L1_FC .gt. L1_SMs)
    L1_PW = merge(L1_FC - 0.01_dp  * L1_FC,  L1_PW, L1_PW .gt. L1_FC )
    ! check the physical limit
    L1_SMs = merge(0.0001_dp, L1_SMs, L1_SMs .lt. 0.0_dp)
    L1_FC  = merge(0.0001_dp, L1_FC,  L1_FC  .lt. 0.0_dp)
    L1_PW  = merge(0.0001_dp, L1_PW,  L1_PW  .lt. 0.0_dp)
    ! Normalise the vertical root distribution profile such that [sum(fRoots) = 1.0]
    !$OMP DO PRIVATE( fTotRoots ) SCHEDULE( STATIC )
    do k = 1, size(L1_fRoots, 1)
      fTotRoots = sum(L1_fRoots(k, :), L1_fRoots(k, :) .gt. 0.0_dp)
      ! This if clause is necessary for test program but may be redundant in actual program
      if (fTotRoots .gt. 0.0_dp) then
        L1_fRoots(k, :) = L1_fRoots(k, :) / fTotRoots
      else
        L1_fRoots(k, :) = 0.0_dp
      end If
    end do

    !$OMP END DO
    !$OMP END PARALLEL
!close(1)
  end subroutine mpr_SMhorizons

end module mo_mpr_SMhorizons
