!> \file mo_mpr_soilmoist.f90
!> \brief \copybrief mo_mpr_soilmoist
!> \details \copydetails mo_mpr_soilmoist

!> \brief Multiscale parameter regionalization (MPR) for soil moisture
!> \details This module contains all routines required for parametrizing soil moisture processes.
!> \authors Stephan Thober, Rohini Kumar
!> \date Dec 2012
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mpr
module mo_mpr_soilmoist

  use mo_kind, only : i4, dp
  use mo_message, only : message, error_message

  implicit none

  public :: mpr_sm

  private

contains
  ! ----------------------------------------------------------------------------

  !    NAME
  !        mpr_sm

  !    PURPOSE
  !>       \brief multiscale parameter regionalization for soil moisture

  !>       \details This subroutine is a wrapper around all soil moisture
  !>       parameter routines. This subroutine requires 13 parameters. These
  !>       parameters have to correspond to the parameters in the original
  !>       parameter array at the following locations: 10-12, 13-18, 27-30.
  !>       Global parameters needed (see mhm_parameter.nml):
  !>       - param( 1) = orgMatterContent_forest
  !>       - param( 2) = orgMatterContent_impervious
  !>       - param( 3) = orgMatterContent_pervious
  !>       - param( 4) = PTF_lower66_5_constant
  !>       - param( 5) = PTF_lower66_5_clay
  !>       - param( 6) = PTF_lower66_5_Db
  !>       - param( 7) = PTF_higher66_5_constant
  !>       - param( 8) = PTF_higher66_5_clay
  !>       - param( 9) = PTF_higher66_5_Db
  !>       - param(10) = PTF_Ks_constant
  !>       - param(11) = PTF_Ks_sand
  !>       - param(12) = PTF_Ks_clay
  !>       - param(13) = PTF_Ks_curveSlope

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(13) :: param"           global parameters
  !>       \param[in] "integer(i4), dimension(:) :: is_present"    indicates whether soiltype is present
  !>       \param[in] "integer(i4), dimension(:) :: nHorizons"     Number of Horizons per soiltype
  !>       \param[in] "integer(i4), dimension(:) :: nTillHorizons" Number of Tillage Horizons
  !>       \param[in] "real(dp), dimension(:, :) :: sand"          sand content
  !>       \param[in] "real(dp), dimension(:, :) :: clay"          clay content
  !>       \param[in] "real(dp), dimension(:, :) :: DbM"           mineral Bulk density
  !>       \param[in] "integer(i4), dimension(:) :: ID0"           cell ids at level 0
  !>       \param[in] "integer(i4), dimension(:, :) :: soilId0"    soil ids at level 0
  !>       \param[in] "integer(i4), dimension(:) :: LCOVER0"       land cover ids at level 0

  !    INTENT(OUT)
  !>       \param[out] "real(dp), dimension(:, :, :) :: thetaS_till"  saturated soil moisture tillage layer
  !>       \param[out] "real(dp), dimension(:, :, :) :: thetaFC_till" field capacity tillage layer
  !>       \param[out] "real(dp), dimension(:, :, :) :: thetaPW_till" permanent wilting point tillage layer
  !>       \param[out] "real(dp), dimension(:, :) :: thetaS"          saturated soil moisture
  !>       \param[out] "real(dp), dimension(:, :) :: thetaFC"         field capacity
  !>       \param[out] "real(dp), dimension(:, :) :: thetaPW"         permanent wilting point
  !>       \param[out] "real(dp), dimension(:, :, :) :: Ks"           saturated hydraulic conductivity
  !>       \param[out] "real(dp), dimension(:, :, :) :: Db"           Bulk density
  !>       \param[out] "real(dp), dimension(:) :: KsVar_H0"           rel. var. of Ks for horizontal flow
  !>       \param[out] "real(dp), dimension(:) :: KsVar_V0"           rel. var. of Ks for vertical flow
  !>       \param[out] "real(dp), dimension(:) :: SMs_FC0"            soil mositure deficit from
  !>       field cap. w.r.t to saturation

  !    HISTORY
  !>       \authors Stephan Thober, Rohini Kumar

  !>       \date Dec 2012

  ! Modifications:
  ! Juliane Mai    Oct 2013 - OLD parametrization --> param(1) = orgMatterContent_forest
  !                                               --> param(2) = orgMatterContent_impervious
  !                                               --> param(3) = orgMatterContent_pervious
  !                                               --> param(4:13) = ... orgMatterContent_forest = orgMatterContent_perv
  !                                                                 + delta_1 NEW parametrization
  !                                               --> param(1) = delta_1 --> param(2) = orgMatterContent_impervious
  !                                               --> param(3) = orgMatterContent_pervious --> param(4:13) = ...
  ! Matthias Zink  Nov 2013 - documentation, inouts --> out moved constants to mhm_constants
  ! Stephan Thober Mar 2014 - separated cell loop from soil loop for better scaling in parallelisation
  ! David Schaefer Mar 2015 - Added dummy variable to avoid redundant computations
  !                           -> Total number of instruction is reduced by ~25%
  !                           (tested on packaged example/gnu48/{release,debug})
  ! Rohini Kumar   Mar 2016 - changes for handling multiple soil database options
  ! Robert Schweppe Jun 2018 - refactoring and reformatting
  ! M. Cuneyd Demirel, Simon Stisen Jun 2020 - added Feddes and FC dependency on root fraction coefficient processCase(3) = 4
  ! Rohini Kumar                    Oct 2021 - Neutron count module to mHM integrate into develop branch (5.11.2)

  subroutine mpr_sm(param, processMatrix, is_present, nHorizons, nTillHorizons, sand, clay, DbM, ID0, soilId0, LCover0, &
                   thetaS_till, thetaFC_till, thetaPW_till, thetaS, thetaFC, thetaPW, Ks, Db, KsVar_H0, KsVar_V0, SMs_FC0)

    use mo_common_constants, only : nodata_dp, nodata_i4
    use mo_mpr_constants, only : BulkDens_OrgMatter
    use mo_mpr_global_variables, only : iFlag_soilDB
    !$ use omp_lib

    implicit none

    ! global parameters
    real(dp), dimension(13), intent(in) :: param

    ! - matrix specifying user defined processes
    integer(i4), dimension(:, :), intent(in) :: processMatrix


    ! indicates whether soiltype is present
    integer(i4), dimension(:), intent(in) :: is_present

    ! Number of Horizons per soiltype
    integer(i4), dimension(:), intent(in) :: nHorizons

    ! Number of Tillage Horizons
    integer(i4), dimension(:), intent(in) :: nTillHorizons

    ! sand content
    real(dp), dimension(:, :), intent(in) :: sand

    ! clay content
    real(dp), dimension(:, :), intent(in) :: clay

    ! mineral Bulk density
    real(dp), dimension(:, :), intent(in) :: DbM

    ! cell ids at level 0
    integer(i4), dimension(:), intent(in) :: ID0

    ! soil ids at level 0
    integer(i4), dimension(:, :), intent(in) :: soilId0

    ! land cover ids at level 0
    integer(i4), dimension(:), intent(in) :: LCOVER0

    ! saturated soil moisture tillage layer
    real(dp), dimension(:, :, :), intent(out) :: thetaS_till

    ! field capacity tillage layer
    real(dp), dimension(:, :, :), intent(out) :: thetaFC_till

    ! permanent wilting point tillage layer
    real(dp), dimension(:, :, :), intent(out) :: thetaPW_till

    ! saturated soil moisture
    real(dp), dimension(:, :), intent(out) :: thetaS

    ! field capacity
    real(dp), dimension(:, :), intent(out) :: thetaFC

    ! permanent wilting point
    real(dp), dimension(:, :), intent(out) :: thetaPW

    ! saturated hydraulic conductivity
    real(dp), dimension(:, :, :), intent(out) :: Ks

    ! Bulk density
    real(dp), dimension(:, :, :), intent(out) :: Db

    ! rel. var. of Ks for horizontal flow
    real(dp), dimension(:), intent(out) :: KsVar_H0

    ! rel. var. of Ks for vertical flow
    real(dp), dimension(:), intent(out) :: KsVar_V0

    ! soil mositure deficit from
    ! field cap. w.r.t to saturation
    real(dp), dimension(:), intent(out) :: SMs_FC0

    ! loop index
    integer(i4) :: i

    ! loop index
    integer(i4) :: j

    ! loop index
    integer(i4) :: l

    ! dummy variable for storing soil class
    integer(i4) :: s

    integer(i4) :: tmp_minSoilHorizon

    real(dp) :: pM

    real(dp) :: pOM

    ! temporal saturated hydr. cond
    real(dp) :: Ks_tmp

    ! van Genuchten shape param
    real(dp) :: Genu_Mual_n

    ! van Genuchten shape param
    real(dp) :: Genu_Mual_alpha

    real(dp) :: tmp_orgMatterContent_forest

    real(dp) :: tmp_orgMatterContent_pervious

    real(dp) :: tmp_orgMatterContent_impervious

    ! total saturated soil mositure content
    real(dp), dimension(:), allocatable :: SMs_tot0

    ! maximum LCover class in L0
    integer(i4) :: max_LCover

    ! saturated hydraulic conductivity
    real(dp), dimension(:, :), allocatable :: Ks_non_till

    ! Case 1 and 4 is based on Jarvis https://doi.org/10.1016/0022-1694(89)90050-4
    ! Case 2 and 3 is based on Feddes https://doi.org/10.1016/0022-1694(76)90017-2
    select case (processMatrix(3, 1))
    case(1,2)
       tmp_orgMatterContent_forest = param(3) + param(1)
    case(3,4)
       tmp_orgMatterContent_forest = param(1)
    end select

    tmp_orgMatterContent_impervious = param(2)
    tmp_orgMatterContent_pervious   = param(3)
    !write(*,*) 'tmp_orgMatterContent_forest = ', tmp_orgMatterContent_forest
    !write(*,*) 'tmp_orgMatterContent_impervious = ', tmp_orgMatterContent_impervious
    !write(*,*) 'tmp_orgMatterContent_pervious = ', tmp_orgMatterContent_pervious

    tmp_minSoilHorizon = minval(nTillHorizons(:))

    ! allocatable local variables
    ! total saturated soil moisture content
    allocate(SMs_tot0(size(ID0, 1)))

    ! some additional variables for iFlag_soil = 1
    if(iFlag_soilDB .eq. 1) then
      s = size(is_present, 1)
      ! note: although second dimension is not required
      ! we assign value 1 to be comparable to other assigned variables
      allocate(Ks_non_till(s, 1))
    end if

    ! initializing soil hydraulic properties related params
    KsVar_H0 = merge(0.0_dp, nodata_dp, ID0 .ne. nodata_i4)
    KsVar_V0 = merge(0.0_dp, nodata_dp, ID0 .ne. nodata_i4)
    SMs_tot0 = merge(0.0_dp, nodata_dp, ID0 .ne. nodata_i4)
    SMs_FC0 = merge(0.0_dp, nodata_dp, ID0 .ne. nodata_i4)

    ! initialization
    thetaS_till = 0.0_dp
    thetaFC_till = 0.0_dp
    thetaPW_till = 0.0_dp
    thetaS = 0.0_dp
    thetaFC = 0.0_dp
    thetaPW = 0.0_dp
    Ks = 0.0_dp
    Db = 0.0_dp
    if(allocated(Ks_non_till)) Ks_non_till = 0.0_dp
    max_LCover = maxval(LCOVER0)
    ! select case according to a given soil database flag
    SELECT CASE(iFlag_soilDB)
      ! classical mHM soil database format
    CASE(0)
      !$OMP PARALLEL default(shared)
      !$OMP DO &
      !$OMP PRIVATE( i, j, L, pOM, pM, Ks_tmp, Genu_Mual_alpha, Genu_Mual_n ) &
      !$OMP SCHEDULE( STATIC )
      do i = 1, size(is_present)
        if (is_present(i) .lt. 1) cycle
        horizon : do j = 1, nHorizons(i)
          ! calculating vertical hydraulic conductivity
          call hydro_cond(Ks_tmp, param(10 : 13), sand(i, j), clay(i, j))
          Ks(i, j, :) = Ks_tmp
          ! calculating other soil hydraulic properties
          ! tillage horizons
          if (j .le. nTillHorizons(i)) then
            ! LC class
            do L = 1, max_LCover
              select case (L)
              case(1)               ! forest
                pOM = tmp_orgMatterContent_forest
              case(2)               ! impervious
                pOM = tmp_orgMatterContent_impervious !param(2)
              case(3)               ! permeable
                pOM = tmp_orgMatterContent_pervious
              case default
                 call error_message('Error mpr_sm: pOM used uninitialized.')
              end select
              pM = 100.0_dp - pOM
              ! bulk density acording to Rawl's (1982) paper
              Db(i, j, L) = 100.0_dp / ((pOM / BulkDens_OrgMatter) + (pM / DbM(i, j)))
              ! Effect of organic matter content
              ! This is taken into account in a simplified form by using
              ! the ratio of(Bd / BdOM)
              Ks(i, j, L) = Ks(i, j, L) * (DbM(i, j) / Db(i, j, L))
              ! estimated SMs_till & van Genuchten's shape parameter (n)
              call Genuchten(thetaS_till(i, j, L), Genu_Mual_n, Genu_Mual_alpha, &
                      param(4 : 9), sand(i, j), clay(i, j), Db(i, j, L))
              ! estimating field capacity
              call field_cap(thetaFC_till(i, j, L), Ks(i, j, L), thetaS_till(i, j, L), Genu_Mual_n)
              ! estimating permanent wilting point
              call PWP(Genu_Mual_n, Genu_Mual_alpha, thetaS_till(i, j, L), thetaPW_till(i, j, L))
            end do
            ! deeper layers
          else
            ! estimate SMs & van Genuchten's shape parameter (n)
            call Genuchten(thetaS(i, j - tmp_minSoilHorizon), Genu_Mual_n, Genu_Mual_alpha, &
                    param(4 : 9), sand(i, j), clay(i, j), DbM(i, j))
            ! estimate field capacity
            call field_cap(thetaFC(i, j - tmp_minSoilHorizon), &
                    Ks(i, j, 1), thetaS(i, j - tmp_minSoilHorizon), Genu_Mual_n)
            ! estimate permanent wilting point
            call PWP(Genu_Mual_n, Genu_Mual_alpha, thetaS(i, j - tmp_minSoilHorizon), &
                    thetaPW(i, j - tmp_minSoilHorizon))
          end if
        end do horizon
      end do
      !$OMP END DO

      ! calculate other soil properties at each location [L0] for regionalising model parameters
      !$OMP DO PRIVATE( s, j ) SCHEDULE( STATIC )
      cellloop : do i = 1, size(soilId0, 1) ! >> here = ncells0
        s = soilId0(i, 1)                    ! >> in this case the second dimension of soilId0 = 1
        do j = 1, nHorizons(s)
          if (j .le. nTillHorizons(s)) then
            ! Soil properties over the whole soil coloum depth
            KsVar_H0(i) = KsVar_H0(i) + thetaS_till(s, j, LCover0(i)) * Ks(s, j, LCover0(i))
            KsVar_V0(i) = KsVar_V0(i) + thetaS_till(s, j, LCover0(i)) / Ks(s, j, LCover0(i))
            SMs_FC0(i) = SMs_FC0(i) + thetaFC_till(s, j, LCover0(i))
            SMs_tot0(i) = SMs_tot0(i) + thetaS_till (s, j, LCover0(i))
          else
            ! soil_properties over the whole soil column
            KsVar_H0(i) = KsVar_H0(i) + thetaS(s, j - tmp_minSoilHorizon) * Ks(s, j, 1)
            KsVar_V0(i) = KsVar_V0(i) + thetaS(s, j - tmp_minSoilHorizon) / Ks(s, j, 1)
            SMs_FC0(i) = SMs_FC0(i) + thetaFC(s, j - tmp_minSoilHorizon)
            SMs_tot0(i) = SMs_tot0(i) + thetaS (s, j - tmp_minSoilHorizon)
          end if
        end do
        ! ------------------------------------------------------------------
        ! DETERMINE RELATIVE VARIABILITIES OF
        !   Ks FOR HORIZONTAL FLOW (KsVar_H)
        !               &
        !   Ks FOR VERTICAL FLOW (KsVar_V)
        ! ------------------------------------------------------------------
        ! soil moisture saturation deficit relative to the field capacity soil moisture
        SMs_FC0(i) = (SMs_tot0(i) - SMs_FC0(i)) / SMs_tot0(i)
        ! Ks variability over the whole soil coloum depth for
        ! both horizontal and vertical flows including relative variabilities
        KsVar_H0(i) = KsVar_H0(i) / SMs_tot0(i) / param(13)
        KsVar_V0(i) = SMs_tot0(i) / KsVar_V0(i) / param(13)
      end do cellloop
      !$OMP END DO
      !$OMP END PARALLEL

      ! to handle multiple soil horizons with unique soil class
    CASE(1)
      do i = 1, size(is_present)
        if (is_present(i) .lt. 1) cycle
        ! **** FOR THE TILLAGE TYPE OF SOIL *****
        ! there is actually no soil horizons/soil type in this case
        ! but we assign of j = 1 to use variables as defined in the classical option (iFlag_soil = 0)
        do j = 1, 1
          ! calculating vertical hydraulic conductivity
          call hydro_cond(Ks_tmp, param(10 : 13), sand(i, j), clay(i, j))
          Ks_non_till(i, j) = Ks_tmp  ! >> non-till
          Ks(i, j, :) = Ks_tmp  ! >> till layers
          ! calculating other soil hydraulic properties
          ! tillage horizons properties depending on the LC class
          do L = 1, max_LCover
            select case (L)
            case(1)               ! forest
              pOM = tmp_orgMatterContent_forest
            case(2)               ! impervious
              pOM = tmp_orgMatterContent_impervious !param(2)
            case(3)               ! permeable
              pOM = tmp_orgMatterContent_pervious
            case default
               call error_message('Error mpr_sm: pOM used is not initialized.')
            end select
            pM = 100.0_dp - pOM
            ! bulk density acording to Rawl's (1982) paper
            Db(i, j, L) = 100.0_dp / ((pOM / BulkDens_OrgMatter) + (pM / DbM(i, j)))
            ! Effect of organic matter content on Ks estimates
            ! This is taken into account in a simplified form by using
            ! the ratio of (Bd/BdOM)
            Ks_tmp = Ks_tmp * (DbM(i, j) / Db(i, j, L))
            Ks(i, j, L) = Ks_tmp
            ! estimated SMs_till & van Genuchten's shape parameter (n)
            call Genuchten(thetaS_till(i, j, L), Genu_Mual_n, Genu_Mual_alpha, &
                    param(4 : 9), sand(i, j), clay(i, j), Db(i, j, L))
            ! estimating field capacity
            call field_cap(thetaFC_till(i, j, L), Ks_tmp, thetaS_till(i, j, L), Genu_Mual_n)
            ! estimating permanent wilting point
            call PWP(Genu_Mual_n, Genu_Mual_alpha, thetaS_till(i, j, L), thetaPW_till(i, j, L))
          end do

          ! *** FOR NON-TILLAGE TYPE OF SOILS ***
          ! note j = 1
          ! since Ks_tmp has changed earlier ... get the original Ks once again
          Ks_tmp = Ks_non_till(i, j)
          ! estimate SMs & van Genuchten's shape parameter (n)
          call Genuchten(thetaS(i, j), Genu_Mual_n, Genu_Mual_alpha, param(4 : 9), sand(i, j), clay(i, j), DbM(i, j))
          ! estimate field capacity
          call field_cap(thetaFC(i, j), Ks_tmp, thetaS(i, j), Genu_Mual_n)
          ! estimate permanent wilting point
          call PWP(Genu_Mual_n, Genu_Mual_alpha, thetaS(i, j), thetaPW(i, j))

        end do  ! >> HORIZON
      end do   ! >> SOIL TYPE

      ! calculate other soil properties at each location [L0] for regionalising model parameters
      do i = 1, size(soilId0, 1)     !! over all cells
        do j = 1, size(soilId0, 2)  !! over horizons
          s = soilId0(i, j)
          if (j .le. nTillHorizons(1)) then
            ! soil properties over the whole soil coloum depth
            KsVar_H0(i) = KsVar_H0(i) + thetaS_till (s, 1, LCover0(i)) * Ks(s, 1, LCover0(i))
            KsVar_V0(i) = KsVar_V0(i) + thetaS_till (s, 1, LCover0(i)) / Ks(s, 1, LCover0(i))
            SMs_FC0(i) = SMs_FC0 (i) + thetaFC_till(s, 1, LCover0(i))
            SMs_tot0(i) = SMs_tot0(i) + thetaS_till (s, 1, LCover0(i))
          else
            ! soil_properties over the whole soil column
            KsVar_H0(i) = KsVar_H0(i) + thetaS (s, 1) * Ks_non_till(s, 1)
            KsVar_V0(i) = KsVar_V0(i) + thetaS (s, 1) / Ks_non_till(s, 1)
            SMs_FC0(i) = SMs_FC0 (i) + thetaFC(s, 1)
            SMs_tot0(i) = SMs_tot0(i) + thetaS (s, 1)
          end if
        end do
        ! ------------------------------------------------------------------
        ! DETERMINE RELATIVE VARIABILITIES OF
        ! Ks FOR HORIZONTAL FLOW (KsVar_H) & Ks FOR VERTICAL FLOW (KsVar_V)
        ! ------------------------------------------------------------------
        ! soil moisture saturation deficit relative to the field capacity soil moisture
        SMs_FC0(i) = (SMs_tot0(i) - SMs_FC0(i)) / SMs_tot0(i)
        ! Ks variability over the whole soil coloum depth for
        ! both horizontal and vertical flows including relative variabilities
        KsVar_H0(i) = KsVar_H0(i) / SMs_tot0(i) / param(13)
        KsVar_V0(i) = SMs_tot0(i) / KsVar_V0(i) / param(13)
      end do

    CASE DEFAULT
      call error_message('***ERROR: iFlag_soilDB option given does not exist. Only 0 and 1 is taken at the moment.')
    END SELECT

    ! free space **
    deallocate(SMs_tot0)
    if(allocated(Ks_non_till)) deallocate(Ks_non_till)

  end subroutine mpr_sm

  ! ------------------------------------------------------------------

  !    NAME
  !        PWP

  !    PURPOSE
  !>       \brief Permanent Wilting point

  !>       \details This subroutine calculates the permanent wilting
  !>       point according to Zacharias et al. (2007, Soil Phy.) and
  !>       using van Genuchten 1980's equation. For the water retention curve at
  !>       a matrix potential of -1500 kPa, it is assumed that thetaR = 0.
  !>       ADDITIONAL INFORMATION
  !>       Zacharias et al. 2007, Soil Phy.

  !    INTENT(IN)
  !>       \param[in] "real(dp) :: Genu_Mual_n"     - Genuchten shape parameter
  !>       \param[in] "real(dp) :: Genu_Mual_alpha" - Genuchten shape parameter
  !>       \param[in] "real(dp) :: thetaS"          - saturated water content

  !    INTENT(OUT)
  !>       \param[out] "real(dp) :: thetaPWP" - Permanent Wilting point

  !    HISTORY
  !>       \authors Stephan Thober, Rohini Kumar

  !>       \date Dec, 2012

  ! Modifications:
  ! Matthias Zink Nov 2013 - documentation, moved constants to mhm_constants
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  elemental pure subroutine PWP(Genu_Mual_n, Genu_Mual_alpha, thetaS, thetaPWP)

    use mo_mpr_constants, only : PWP_c, PWP_matPot_ThetaR

    implicit none

    ! - Genuchten shape parameter
    real(dp), intent(in) :: Genu_Mual_n

    ! - Genuchten shape parameter
    real(dp), intent(in) :: Genu_Mual_alpha

    ! - saturated water content
    real(dp), intent(in) :: thetaS

    ! - Permanent Wilting point
    real(dp), intent(out) :: thetaPWP

    real(dp) :: x

    ! Genuchten shape parameter
    real(dp) :: Genu_Mual_m


    Genu_Mual_m = PWP_c - (PWP_c / Genu_Mual_n)
    x = PWP_c + exp(Genu_Mual_n * log(Genu_Mual_alpha * PWP_matPot_ThetaR))
    x = exp(Genu_Mual_m * log(x))
    ! constrain
    if (x < 1.0_dp) x = 1.0_dp
    thetaPWP = thetaS / x

  end subroutine PWP

  ! ----------------------------------------------------------------------------

  !    NAME
  !        field_cap

  !    PURPOSE
  !>       \brief calculates the field capacity

  !>       \details estimate Field capacity; FC -- Flux based
  !>       approach (Twarakavi, et. al. 2009, WRR)
  !>       According to the
  !>       above reference FC is defined as the soil water content at
  !>       which the drainage from a profile ceases under natural
  !>       conditions. Since drainage from a soil profile in a simulation
  !>       never becomes zero, we assume that drainage ceases when the
  !>       bottom flux from the soil reaches a value that is equivalent to
  !>       the minimum amount of precipitation that could be recorded
  !>       (i.e. 0.01 cm/d == 1 mm/d). It is assumed that ThetaR = 0.0_dp
  !>       ADDITIONAL INFORMATION
  !>       Twarakavi, et. al. 2009, WRR

  !    INTENT(OUT)
  !>       \param[out] "real(dp) :: thetaFC" - Field capacity

  !    INTENT(IN)
  !>       \param[in] "real(dp) :: Ks"          - saturated hydraulic conductivity
  !>       \param[in] "real(dp) :: thetaS"      - saturated water content
  !>       \param[in] "real(dp) :: Genu_Mual_n" - Genuchten shape parameter

  !    HISTORY
  !>       \authors Stephan Thober, Rohini Kumar

  !>       \date Dec 2012

  ! Modifications:
  ! Matthias Zink Nov 2013 - documentation, moved constants to mhm_constants
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  elemental pure subroutine field_cap(thetaFC, Ks, thetaS, Genu_Mual_n)

    use mo_mpr_constants, only : field_cap_c1, field_cap_c2

    implicit none

    ! - saturated hydraulic conductivity
    real(dp), intent(in) :: Ks

    ! - saturated water content
    real(dp), intent(in) :: thetaS

    ! - Genuchten shape parameter
    real(dp), intent(in) :: Genu_Mual_n

    ! - Field capacity
    real(dp), intent(out) :: thetaFC

    real(dp) :: x


    x = (field_cap_c1) * (field_cap_c2 + log10(Ks))
    thetaFC = thetaS * exp(x * log(Genu_Mual_n))

  end subroutine field_cap

  ! ----------------------------------------------------------------------------

  !    NAME
  !        Genuchten

  !    PURPOSE
  !>       \brief calculates the Genuchten shape parameter

  !>       \details estimate SMs_till & van Genuchten's shape parameter (n)
  !>       (Zacharias et al, 2007, soil Phy.)
  !>       Global parameters needed (see mhm_parameter.nml):
  !>       - param( 1) = PTF_lower66_5_constant
  !>       - param( 2) = PTF_lower66_5_clay
  !>       - param( 3) = PTF_lower66_5_Db
  !>       - param( 4) = PTF_higher66_5_constant
  !>       - param( 5) = PTF_higher66_5_clay
  !>       - param( 6) = PTF_higher66_5_Db
  !>       ADDITIONAL INFORMATION
  !>       Zacharias et al, 2007, soil Phy.

  !    INTENT(OUT)
  !>       \param[out] "real(dp) :: thetaS"          - saturated water content
  !>       \param[out] "real(dp) :: Genu_Mual_n"     - van Genuchten shape parameter
  !>       \param[out] "real(dp) :: Genu_Mual_alpha" - van Genuchten shape parameter

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(6) :: param" parameters
  !>       \param[in] "real(dp) :: sand"                - [%] sand content
  !>       \param[in] "real(dp) :: clay"                - [%] clay content
  !>       \param[in] "real(dp) :: Db"                  - [10^3 kg/m3] bulk density

  !    HISTORY
  !>       \authors Stephan Thober, Rohini Kumar

  !>       \date Dec 2012

  ! Modifications:
  ! Rohini Kumar Mar 2014 - ThetaS limit changed from 0 to 0.001
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine Genuchten(thetaS, Genu_Mual_n, Genu_Mual_alpha, param, sand, clay, Db)

    use mo_mpr_constants, only : vGenuchtenN_c1, vGenuchtenN_c10, vGenuchtenN_c11, vGenuchtenN_c12, vGenuchtenN_c13, &
                                 vGenuchtenN_c14, vGenuchtenN_c15, vGenuchtenN_c16, vGenuchtenN_c17, vGenuchtenN_c18, &
                                 vGenuchtenN_c2, vGenuchtenN_c3, vGenuchtenN_c4, vGenuchtenN_c5, vGenuchtenN_c6, &
                                 vGenuchtenN_c7, vGenuchtenN_c8, vGenuchtenN_c9, vGenuchten_sandtresh

    implicit none

    ! parameters
    real(dp), dimension(6), intent(in) :: param

    ! - [%] sand content
    real(dp), intent(in) :: sand

    ! - [%] clay content
    real(dp), intent(in) :: clay

    ! - [10^3 kg/m3] bulk density
    real(dp), intent(in) :: Db

    ! - saturated water content
    real(dp), intent(out) :: thetaS

    ! - van Genuchten shape parameter
    real(dp), intent(out) :: Genu_Mual_n

    ! - van Genuchten shape parameter
    real(dp), intent(out) :: Genu_Mual_alpha

    ! temporal variable
    real(dp) :: x


    ! estimate SMs_till & van Genuchten's parameters (alpha and n)
    if (sand < vGenuchten_sandtresh) then
      thetaS = param(1) + param(2) * clay + param(3) * Db
      Genu_Mual_n = vGenuchtenN_c1 - vGenuchtenN_c2 * (sand**(vGenuchtenN_c3)) + &
              vGenuchtenN_c4 * (clay**(vGenuchtenN_c5))
      x = vGenuchtenN_c6 + vGenuchtenN_c7 * sand + vGenuchtenN_c8 * clay - &
              vGenuchtenN_c9 * Db
    else
      thetaS = param(4) + param(5) * clay + param(6) * Db
      Genu_Mual_n = vGenuchtenN_c10 + vGenuchtenN_c11 * (sand**(vGenuchtenN_c12)) + &
              vGenuchtenN_c13 * (clay**(vGenuchtenN_c14))
      x = vGenuchtenN_c15 + vGenuchtenN_c16 * sand + vGenuchtenN_c17 * clay - &
              vGenuchtenN_c18 * Db
    end if

    ! Mualem alpha
    Genu_Mual_alpha = exp(x)

    ! hard coded limits, according to (Zacharias et al, 2007, soil Phy.)
    if (thetaS < 0.01_dp) then
      call message('thetaS below threshold limit 1e-2, reset.')
      ! Put constrains on theta_S
      thetaS = 0.01_dp
    end if
    if (thetaS > 1.0_dp) then
      call message('thetaS above 1, reset.')
      ! Put constrains on theta_S
      thetaS = 1.0_dp
    end if
    if (Genu_Mual_n < 1.01000_dp) then
      call message('Genu_Mual_n below threshold limit 1.01, reset.')
      Genu_Mual_n = 1.01000_dp
    end if
    if (Genu_Mual_alpha < 0.00001_dp) then
      call message('Genu_Mual_alpha below threshold limit 1e-5, reset.')
      Genu_Mual_alpha = 0.00001_dp
    end if

  end subroutine Genuchten

  ! ----------------------------------------------------------------------------

  !    NAME
  !        hydro_cond

  !    PURPOSE
  !>       \brief calculates the hydraulic conductivity Ks

  !>       \details By default save this value of Ks, particularly for the
  !>       deeper layers where OM content plays relatively low or no role
  !>       Global parameters needed (see mhm_parameter.nml):
  !>       - param(1) = PTF_Ks_constant
  !>       - param(2) = PTF_Ks_sand
  !>       - param(3) = PTF_Ks_clay
  !>       - param(4) = PTF_Ks_curveSlope
  !>       ADDITIONAL INFORMATION
  !>       Written,  Stephan Thober, Dec 2012

  !    INTENT(OUT)
  !>       \param[out] "real(dp) :: KS"

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(4) :: param"
  !>       \param[in] "real(dp) :: sand"                - [%] sand content
  !>       \param[in] "real(dp) :: clay"                - [%] clay content

  !    HISTORY
  !>       \authors Stephan Thober, Rohini Kumar

  !>       \date Dec 2012

  ! Modifications:
  ! Matthias Zink  Nov 2013 - documentation, moved constants to mhm_constants
  ! Matthias Cuntz Jun 2014 - suggested to fix param(4)
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine hydro_cond(KS, param, sand, clay)

    use mo_mpr_constants, only : Ks_c

    implicit none

    real(dp), dimension(4), intent(in) :: param

    ! - [%] sand content
    real(dp), intent(in) :: sand

    ! - [%] clay content
    real(dp), intent(in) :: clay

    real(dp), intent(out) :: KS

    ! temporal variable
    real(dp) :: x


    ! saturated vertical hydraulic conductivity, Ks (cm/d)
    !   from Cosby et. al. (WRR 1984) Table 4
    ! param(4) is the unit conversion from inch/h to cm/d and should be a constant.
    ! Fix it in the namelist, i.e. in
    ! mhm_parameter.nml set the 4th value (=FLAG) to 0 and the third value to 60.96
    !   PTF_Ks_curveSlope = 60.96, 60.96, 60.96, 0, 1
    x = param(1) + param(2) * sand - param(3) * clay
    Ks = param(4) * exp(X * log(Ks_c))

    if (Ks < 1.10_dp) then
      call message('JMJMJM-Ks-BAD')
    end if

    ! minimum value of Ks = 1.1cm/d
    if (Ks < 1.10_dp) Ks = 1.10_dp

  end subroutine hydro_cond

end module mo_mpr_soilmoist
