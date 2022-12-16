!> \file mo_mpr_pet.f90
!> \brief \copybrief mo_mpr_pet
!> \details \copydetails mo_mpr_pet

!> \brief MPR routine for PET.
!> \details This module sets up pet correction factor at level-1 based on LAI
!> \authors Mehmet Cuneyd Demirel, Simon Stisen
!> \date May 2017
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mpr
module mo_mpr_pet

  use mo_kind, only : i4, dp

  implicit none

  PUBLIC :: pet_correctbyLAI          ! estimate PET correction factor with distributed LAI
  PUBLIC :: pet_correctbyASP          ! estimate PET correction factor with distributed Aspect
  PUBLIC :: priestley_taylor_alpha    ! factor (alpha) for Presley-Taylor ET estimation
  !PUBLIC :: aerodynamical_resistance  ! aerodynamical resistance (ra) for Penman-Monteith ET estimation
  PUBLIC :: bulksurface_resistance    ! bulk surface (stomatal) resistance (rs) for Penman-Monteith ET estimation

  private

contains


  ! ----------------------------------------------------------------------------

  !    NAME
  !        pet_correctbyLAI

  !    PURPOSE
  !>       \brief estimate PET correction factor based on LAI at L1

  !>       \details estimate PET correction factor based on LAI at L1 for a given
  !>       Leaf Area Index field.

  !>       Global parameters needed (see mhm_parameter.nml):
  !>       Process Case 5:
  !>       - param(1) = PET_a_forest
  !>       - param(2) = PET_a_impervious
  !>       - param(3) = PET_a_pervious
  !>       - param(4) = PET_b
  !>       - param(5) = PET_c

  !>       Example DSF=PET_a+PET_b*(1-exp(PET_c*LAI))

  !>       Similar to the crop coefficient concept Kc=a+b*(1-exp(c*LAI)) by Allen, R. G., L. S. Pereira,
  !>       D. Raes, and M. Smith (1998), Crop evapotranspiration - Guidelines for computing crop water requirements.,
  !>       FAO
  !>       Irrigation and drainage paper 56. See Chapter 9, Equation 97  <http://www.fao.org/docrep/X0490E/x0490e0f.htm>
  !>       Date: 17/5/2017

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(5) :: param"         parameters
  !>       \param[in] "real(dp) :: nodata"                      - nodata value
  !>       \param[in] "integer(i4), dimension(:) :: LCOVER0"    Land cover at level 0
  !>       \param[in] "real(dp), dimension(:, :) :: LAI0"       LAI at level-0
  !>       \param[in] "logical, dimension(:, :) :: mask0"       mask at L0
  !>       \param[in] "integer(i4), dimension(:) :: cell_id0"   Cell ids of hi res field
  !>       \param[in] "integer(i4), dimension(:) :: upp_row_L1" Upper row of hi res block
  !>       \param[in] "integer(i4), dimension(:) :: low_row_L1" Lower row of hi res block
  !>       \param[in] "integer(i4), dimension(:) :: lef_col_L1" Left column of hi res block
  !>       \param[in] "integer(i4), dimension(:) :: rig_col_L1" Right column of hi res block
  !>       \param[in] "integer(i4), dimension(:) :: nL0_in_L1"  Number of L0 cells within a L1 cel

  !    INTENT(INOUT)
  !>       \param[inout] "real(dp), dimension(:, :) :: L1_petLAIcorFactor" pet cor factor at level-1

  !    HISTORY
  !>       \authors M. Cuneyd Demirel and Simon Stisen from GEUS.dk

  !>       \date May. 2017

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine pet_correctbyLAI(param, nodata, LCOVER0, LAI0, mask0, cell_id0, upp_row_L1, low_row_L1, lef_col_L1, &
                             rig_col_L1, nL0_in_L1, L1_petLAIcorFactor)

    use mo_upscaling_operators, only : upscale_harmonic_mean
    !$ use omp_lib

    implicit none

    ! parameters
    real(dp), dimension(5), intent(in) :: param

    ! - nodata value
    real(dp), intent(in) :: nodata

    ! Land cover at level 0
    integer(i4), dimension(:), intent(in) :: LCOVER0

    ! LAI at level-0
    real(dp), dimension(:, :), intent(in) :: LAI0

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

    ! pet cor factor at level-1
    real(dp), dimension(:, :), intent(inout) :: L1_petLAIcorFactor

    ! pet cor factor at level-0
    real(dp), dimension(size(LCOVER0, 1), size(LAI0, 2)) :: petLAIcorFactor_0

    ! loop index
    integer(i4) :: kk, tt

    ! loop index
    integer(i4) :: LL


    ! ------------------------------------------------------------------
    ! Estimate DSF=PET_a+PET_b*(1-exp(PET_c*LAI)) to correct PET as PET=DSF*PET
    ! ------------------------------------------------------------------
    !$OMP PARALLEL
    !$OMP DO PRIVATE( LL ) SCHEDULE( STATIC )
    ! need to be done for every landcover to get DSF
    do kk = 1, size(LCOVER0, 1)

      LL = LCOVER0(kk)

      ! TODO: memory order of arrays is not optimal, how to improve?
      select case(LL)
      case(1) ! forest
        petLAIcorFactor_0(kk, :) = param(1) + (param(4) * (1.0_dp - exp(param(5) * LAI0(kk, :))))
      case(2) ! impervious
        petLAIcorFactor_0(kk, :) = param(2) + (param(4) * (1.0_dp - exp(param(5) * LAI0(kk, :))))
      case(3) ! permeable
        petLAIcorFactor_0(kk, :) = param(3) + (param(4) * (1.0_dp - exp(param(5) * LAI0(kk, :))))
      end select

    end do
    !$OMP END DO
    !$OMP END PARALLEL

    do tt = 1, size(LAI0, 2)
      L1_petLAIcorFactor(:, tt) = upscale_harmonic_mean(nL0_in_L1, Upp_row_L1, Low_row_L1, &
              Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata, petLAIcorFactor_0(:, tt))
    end do


  end subroutine pet_correctbyLAI

  ! ----------------------------------------------------------------------------

  !    NAME
  !        pet_correctbyASP

  !    PURPOSE
  !>       \brief correction of PET

  !>       \details Correction of PET based on L0 aspect data.


  !>       Global parameters needed (see mhm_parameter.nml):
  !>       - param(1) = minCorrectionFactorPET
  !>       - param(2) = maxCorrectionFactorPET
  !>       - param(3) = aspectTresholdPET

  !    INTENT(IN)
  !>       \param[in] "integer(i4), dimension(:) :: id0"      Level 0 cell id
  !>       \param[in] "real(dp), dimension(:) :: latitude_l0" latitude on l0
  !>       \param[in] "real(dp), dimension(:) :: Asp0"        [degree] Aspect at Level 0
  !>       \param[in] "real(dp), dimension(3) :: param"       process parameters
  !>       \param[in] "real(dp) :: nodata"                    - no data value

  !    INTENT(OUT)
  !>       \param[out] "real(dp), dimension(:) :: fAsp0" PET correction for Aspect

  !    HISTORY
  !>       \authors Stephan Thober, Rohini Kumar

  !>       \date Dec 2012

  ! Modifications:
  ! Juliane Mai    Oct 2013 - OLD parametrization  --> param(1) = minCorrectionFactorPET
  !                                                --> param(2) = maxCorrectionFactorPET
  !                                                --> param(3) = aspectTresholdPET
  !                                             -------------------------------
  !                                             maxCorrectionFactorPET = minCorrectionFactorPET + delta
  !                                             -------------------------------
  !                                             NEW parametrization
  !                                                --> param(1) = minCorrectionFactorPET
  !                                                --> param(2) = delta
  !                                                --> param(3) = aspectTresholdPET
  ! Stephan Thober Dec 2013 - changed intent(inout) to intent(out)
  ! Stephan Thober Sep 2015 - Mapping L1 to Lo, latitude on L0
  ! Luis Samaniego Sep 2015 - PET correction on the southern hemisphere
  ! Matthias Zink  Jun 2017 - renamed and moved from mo_multi_scale_param_reg.f90 to mo_mpr_pet.f90
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine pet_correctbyASP(Id0, latitude_l0, Asp0, param, nodata, fAsp0)
    !$ use omp_lib

    implicit none

    ! Level 0 cell id
    integer(i4), dimension(:), intent(in) :: id0

    ! latitude on l0
    real(dp), dimension(:), intent(in) :: latitude_l0

    ! - no data value
    real(dp), intent(in) :: nodata

    ! process parameters
    real(dp), dimension(3), intent(in) :: param

    ! [degree] Aspect at Level 0
    real(dp), dimension(:), intent(in) :: Asp0

    ! PET correction for Aspect
    real(dp), dimension(:), intent(out) :: fAsp0

    ! PET correction for Aspect, south
    real(dp), dimension(size(id0, 1)) :: fAsp0S

    logical, dimension(size(id0, 1)) :: mask_north_hemisphere_l0

    real(dp) :: tmp_maxCorrectionFactorPET


    mask_north_hemisphere_l0 = merge(.TRUE., .FALSE., latitude_l0 .gt. 0.0_dp)

    tmp_maxCorrectionFactorPET = param(1) + param(2)

    ! for cells on the northern hemisphere
    !$OMP PARALLEL
    fAsp0 = merge(&
            param(1) + (tmp_maxCorrectionFactorPET - param(1)) / param(3) * asp0, &
            param(1) + (tmp_maxCorrectionFactorPET - param(1)) / (360._dp - param(3)) * (360._dp - Asp0), &
            !         ( asp0 < param(3) ) .and. mask_north_hemisphere_l0  )
            asp0 < param(3))
    fAsp0 = merge(fAsp0, nodata, Id0 /= int(nodata, i4))
    !$OMP END PARALLEL

    ! for cells on the southern hemisphere
    !$OMP PARALLEL
    fAsp0S = merge(&
            param(1) + (tmp_maxCorrectionFactorPET - param(1)) / (360._dp - param(3)) * (360._dp - Asp0), &
            param(1) + (tmp_maxCorrectionFactorPET - param(1)) / param(3) * asp0, &
            asp0 < param(3))
    fAsp0S = merge(fAsp0S, nodata, Id0 /= int(nodata, i4))
    !$OMP END PARALLEL

    !$OMP PARALLEL
    fAsp0 = merge(fAsp0, fAsp0S, mask_north_hemisphere_l0)
    !$OMP END PARALLEL

  end subroutine pet_correctbyASP

  ! ----------------------------------------------------------------------------

  !    NAME
  !        priestley_taylor_alpha

  !    PURPOSE
  !>       \brief Regionalization of priestley taylor alpha

  !>       \details estimation of priestley taylor alpha
  !>       Global parameters needed (see mhm_parameter.nml):
  !>       - param(1) = PriestleyTaylorCoeff
  !>       - param(2) = PriestleyTaylorLAIcorr

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:, :) :: LAI0"       LAI at level-0
  !>       \param[in] "real(dp), dimension(:) :: param"         input parameter
  !>       \param[in] "logical, dimension(:, :) :: mask0"       mask at level 0
  !>       \param[in] "real(dp) :: nodata"                      - nodata value
  !>       \param[in] "integer(i4), dimension(:) :: cell_id0"   Cell ids of hi res field
  !>       \param[in] "integer(i4), dimension(:) :: nL0_in_L1"  number of l0 cells within a l1 cell
  !>       \param[in] "integer(i4), dimension(:) :: Upp_row_L1" upper row of a l1 cell in l0 grid
  !>       \param[in] "integer(i4), dimension(:) :: Low_row_L1" lower row of a l1 cell in l0 grid
  !>       \param[in] "integer(i4), dimension(:) :: Lef_col_L1" left col of a l1 cell in l0 grid
  !>       \param[in] "integer(i4), dimension(:) :: Rig_col_L1" right col of a l1 cell in l0 grid

  !    INTENT(OUT)
  !>       \param[out] "real(dp), dimension(:, :) :: priestley_taylor_alpha1" bulk surface resistance

  !    HISTORY
  !>       \authors Matthias Zink

  !>       \date Apr 2013

  ! Modifications:
  ! Matthias Zink Jun 2017 - moved from mo_multi_scale_param_reg.f90 to mo_mpr_pet.f90
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine priestley_taylor_alpha(LAI0, param, mask0, nodata, cell_id0, nL0_in_L1, Upp_row_L1, Low_row_L1, Lef_col_L1, &
                                   Rig_col_L1, priestley_taylor_alpha1)

    use mo_upscaling_operators, only : upscale_arithmetic_mean

    implicit none

    ! LAI at level-0
    real(dp), dimension(:, :), intent(in) :: LAI0

    ! input parameter
    real(dp), dimension(:), intent(in) :: param

    ! mask at level 0
    logical, dimension(:, :), intent(in) :: mask0

    ! - nodata value
    real(dp), intent(in) :: nodata

    ! Cell ids of hi res field
    integer(i4), dimension(:), intent(in) :: cell_id0

    ! number of l0 cells within a l1 cell
    integer(i4), dimension(:), intent(in) :: nL0_in_L1

    ! upper row of a l1 cell in l0 grid
    integer(i4), dimension(:), intent(in) :: Upp_row_L1

    ! lower row of a l1 cell in l0 grid
    integer(i4), dimension(:), intent(in) :: Low_row_L1

    ! left col of a l1 cell in l0 grid
    integer(i4), dimension(:), intent(in) :: Lef_col_L1

    ! right col of a l1 cell in l0 grid
    integer(i4), dimension(:), intent(in) :: Rig_col_L1

    ! bulk surface resistance
    real(dp), dimension(:, :), intent(out) :: priestley_taylor_alpha1

    integer(i4) :: tt

    ! dim 1 = number of cells on level 0, time
    real(dp), dimension(:, :), allocatable :: priestley_taylor_alpha0


    ! initialize some things
    allocate(priestley_taylor_alpha0 (size(LAI0, 1), size(LAI0, 2))) ; priestley_taylor_alpha0 = nodata
    priestley_taylor_alpha1 = nodata
    !
    do tt = 1, size(LAI0, 2)
      priestley_taylor_alpha0(:, tt) = param(1) + param(2) * LAI0(:, tt)

      priestley_taylor_alpha1(:, tt) = upscale_arithmetic_mean(nL0_in_L1, Upp_row_L1, Low_row_L1, &
              Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata, priestley_taylor_alpha0(:, tt))

    end do

  end subroutine priestley_taylor_alpha

  ! ----------------------------------------------------------------------------

  !    NAME
  !        bulksurface_resistance

  !    PURPOSE
  !>       \brief Regionalization of bulk surface resistance

  !>       \details estimation of bulk surface resistance
  !>       Global parameters needed (see mhm_parameter.nml):
  !>       - param(1) = stomatal_resistance

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:, :) :: LAI0"       LAI at level-0
  !>       \param[in] "real(dp) :: param"                       - global parameter
  !>       \param[in] "logical, dimension(:, :) :: mask0"       mask at level 0
  !>       \param[in] "real(dp) :: nodata"                      - nodata value
  !>       \param[in] "integer(i4), dimension(:) :: cell_id0"   Cell ids of hi res field
  !>       \param[in] "integer(i4), dimension(:) :: nL0_in_L1"  number of l0 cells within a l1 cell
  !>       \param[in] "integer(i4), dimension(:) :: Upp_row_L1" upper row of a l1 cell in l0 grid
  !>       \param[in] "integer(i4), dimension(:) :: Low_row_L1" lower row of a l1 cell in l0 grid
  !>       \param[in] "integer(i4), dimension(:) :: Lef_col_L1" left col of a l1 cell in l0 grid
  !>       \param[in] "integer(i4), dimension(:) :: Rig_col_L1" right col of a l1 cell in l0 grid

  !    INTENT(OUT)
  !>       \param[out] "real(dp), dimension(:, :) :: bulksurface_resistance1" bulk surface resistance

  !    HISTORY
  !>       \authors Matthias Zink

  !>       \date Apr 2013

  ! Modifications:
  ! Matthias Zink Jun 2017 - moved from mo_multi_scale_param_reg.f90 to mo_mpr_pet.f90
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine bulksurface_resistance(LAI0, param, mask0, nodata, cell_id0, nL0_in_L1, Upp_row_L1, Low_row_L1, Lef_col_L1, &
                                   Rig_col_L1, bulksurface_resistance1)

    use mo_mpr_constants, only : LAI_factor_surfResi, LAI_offset_surfResi, max_surfResist
    use mo_upscaling_operators, only : upscale_arithmetic_mean

    implicit none

    ! LAI at level-0
    real(dp), dimension(:, :), intent(in) :: LAI0

    ! - global parameter
    real(dp), intent(in) :: param

    ! mask at level 0
    logical, dimension(:, :), intent(in) :: mask0

    ! - nodata value
    real(dp), intent(in) :: nodata

    ! Cell ids of hi res field
    integer(i4), dimension(:), intent(in) :: cell_id0

    ! number of l0 cells within a l1 cell
    integer(i4), dimension(:), intent(in) :: nL0_in_L1

    ! upper row of a l1 cell in l0 grid
    integer(i4), dimension(:), intent(in) :: Upp_row_L1

    ! lower row of a l1 cell in l0 grid
    integer(i4), dimension(:), intent(in) :: Low_row_L1

    ! left col of a l1 cell in l0 grid
    integer(i4), dimension(:), intent(in) :: Lef_col_L1

    ! right col of a l1 cell in l0 grid
    integer(i4), dimension(:), intent(in) :: Rig_col_L1

    ! bulk surface resistance
    real(dp), dimension(:, :), intent(out) :: bulksurface_resistance1

    integer(i4) :: tt

    ! dim 1 = number of cells on level 0,
    ! dim 2 = number of months in year (12)
    real(dp), dimension(:, :), allocatable :: bulksurface_resistance0


    ! initialize some things
    allocate(bulksurface_resistance0 (size(LAI0, 1), size(LAI0, 2))) ; bulksurface_resistance0 = nodata
    bulksurface_resistance1 = nodata
    !

    do tt = 1, size(LAI0, 2)
      bulksurface_resistance0(:, tt) = param / (LAI0(:, tt) / &
              (LAI_factor_surfResi * LAI0(:, tt) + LAI_offset_surfResi))
      ! efeective LAI from McMahon et al ,2013 , HESS supplements

      ! since LAI may be very low, rs becomes very high
      ! thus the values are restricted to maximum literaure values (i.e. McMahon et al ,2013 , HESS)
      bulksurface_resistance0(:, tt) = merge(max_surfResist, bulksurface_resistance0(:, tt), &
              bulksurface_resistance0(:, tt) .GT. max_surfResist)

      bulksurface_resistance1(:, tt) = upscale_arithmetic_mean(nL0_in_L1, Upp_row_L1, Low_row_L1, &
              Lef_col_L1, Rig_col_L1, cell_id0, mask0, nodata, bulksurface_resistance0(:, tt))

    end do

  end subroutine bulksurface_resistance

end module mo_mpr_pet
