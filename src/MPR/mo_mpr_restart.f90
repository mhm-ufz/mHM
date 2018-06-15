!> \file mo_restart.f90

!> \brief reading and writing states, fluxes and configuration for restart of mHM.

!> \details routines are seperated for reading and writing variables for:\n
!>          - states and fluxes, and \n
!>          - configuration.\n
!>          Reading of L11 configuration is also seperated from the rest, 
!>          since it is only required when routing is activated.

!> \authors Stephan Thober
!> \date Jul 2013

MODULE mo_mpr_restart

  ! This module is a restart for the UFZ CHS mesoscale hydrologic model mHM.

  ! Written  Stephan Thober, Apr 2011

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: write_eff_params        ! read restart files for configuration from a given path
  PUBLIC :: write_mpr_restart_files     ! write restart files for configuration to a given path

  INTERFACE unpack_field_and_write
    MODULE PROCEDURE unpack_field_and_write_1d_i4, &
            unpack_field_and_write_1d_dp, &
            unpack_field_and_write_2d_dp, &
            unpack_field_and_write_3d_dp
  end interface unpack_field_and_write


CONTAINS
  ! ------------------------------------------------------------------

  !      NAME
  !         write_restart

  !     PURPOSE
  !>        \brief write restart files for each basin

  !>        \details write restart files for each basin. For each basin
  !>        three restart files are written. These are xxx_states.nc, 
  !>        xxx_L11_config.nc, and xxx_config.nc (xxx being the three digit
  !>        basin index). If a variable is added here, it should also be added
  !>        in the read restart routines below.

  !     INTENT(IN)
  !>        \param[in] "character(256), dimension(:) :: OutPath"     Output Path for each basin

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN

  !     RESTRICTIONS 
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author   Stephan Thober
  !>        \date     Jun 2014
  !         Modified  Matthias Zink   Nov. 2014  - added PET related parameter writing
  !                   Stephan Thober  Aug  2015  - moved write of routing states to mRM
  !                   David Schaefer  Nov  2015  - mo_netcdf
  !                   Stephan Thober  Nov  2016  - moved processMatrix to common variables
  !                   Zink M. Demirel C.,Mar 2017 - Added Jarvis soil water stress function at SM process(3)  

  ! ------------------------------------------------------------------ 
  subroutine write_mpr_restart_files(OutPath)

    use mo_kind, only : i4
    use mo_message, only : message
    use mo_string_utils, only : num2str
    use mo_netcdf, only : NcDataset, NcDimension
    use mo_mpr_global_variables, only : &
            nSoilHorizons_mHM, &  ! number of soil horizons
            nLAI
    use mo_common_variables, only : &
            level1, &
            nLCoverScene ! number of land cover scenes
    use mo_common_restart, only : write_grid_info

    implicit none

    character(256) :: Fname
    character(256), dimension(:), intent(in) :: OutPath ! list of Output paths per Basin
    integer(i4) :: iBasin
    integer(i4) :: s1       ! start index at level 1
    integer(i4) :: e1       ! end index at level 1
    logical, dimension(:, :), allocatable :: mask1    ! mask at level 1

    type(NcDataset) :: nc
    type(NcDimension) :: rows1, cols1, soil1, lcscenes, lais

    basin_loop : do iBasin = 1, size(OutPath)

      ! write restart file for iBasin
      Fname = trim(OutPath(iBasin)) // "mHM_restart_" // trim(num2str(iBasin, "(i3.3)")) // ".nc"
      ! print a message
      call message("    Writing Restart-file: ", trim(adjustl(Fname)), " ...")

      nc = NcDataset(fname, "w")

      call write_grid_info(level1(iBasin), "1", nc)

      rows1 = nc%getDimension("nrows1")
      cols1 = nc%getDimension("ncols1")

      soil1 = nc%setDimension("L1_soilhorizons", nSoilHorizons_mHM)
      lcscenes = nc%setDimension("LCoverScenes", nLCoverScene)
      lais = nc%setDimension("LAI_timesteps", nLAI)

      ! for appending and intialization
      allocate(mask1(rows1%getLength(), cols1%getLength()))
      s1 = level1(iBasin)%iStart
      e1 = level1(iBasin)%iEnd
      mask1 = level1(iBasin)%mask

      call write_eff_params(mask1, s1, e1, rows1, cols1, soil1, lcscenes, lais, nc)
      deallocate(mask1)
      call nc%close()

    end do basin_loop

  end subroutine write_mpr_restart_files

  subroutine write_eff_params(mask1, s1, e1, rows1, cols1, soil1, lcscenes, lais, nc)

    use mo_kind, only : i4
    use mo_netcdf, only : NcDataset, NcDimension, NcVariable
    use mo_common_constants, only : nodata_dp, nodata_i4
    use mo_mpr_global_variables, only : &
            L1_fSealed, &
            L1_alpha, &
            L1_degDayInc, &
            L1_degDayMax, &
            L1_degDayNoPre, &
            L1_degDay, &
            L1_karstLoss, &
            L1_fAsp, &
            L1_HarSamCoeff, &
            L1_PrieTayAlpha, &
            L1_aeroResist, &
            L1_surfResist, &
            L1_fRoots, &
            L1_maxInter, &
            L1_kfastFlow, &
            L1_kSlowFlow, &
            L1_kBaseFlow, &
            L1_kPerco, &
            L1_soilMoistFC, &
            L1_soilMoistSat, &
            L1_soilMoistExp, &
            L1_jarvis_thresh_c1, &
            L1_petLAIcorFactor, &
            L1_tempThresh, &
            L1_unsatThresh, &
            L1_sealedThresh, &
            L1_wiltingPoint
    use mo_common_variables, only : &
            LC_year_start, LC_year_end, & ! LCscenes details
            processMatrix ! process configuration

    implicit none

    logical, dimension(:, :), allocatable, intent(in) :: mask1    ! mask at level 1
    integer(i4), intent(in) :: s1       ! start index at level 1
    integer(i4), intent(in) :: e1       ! end index at level 1
    type(NcDimension), intent(in) :: rows1, cols1, soil1, lcscenes, lais
    type(NcDataset), intent(inout) :: nc
    type(NcVariable) :: var

    ! set variable
    var = nc%setVariable("LC_year_start", "i32", (/lcscenes/))
    call var%setFillValue(nodata_i4)
    call var%setData(LC_year_start)
    call var%setAttribute("long_name", "start year of land cover scene")

    var = nc%setVariable("LC_year_end", "i32", (/lcscenes/))
    call var%setFillValue(nodata_i4)
    call var%setData(LC_year_end)
    call var%setAttribute("long_name", "end year of land cover scene")

    !-------------------------------------------
    ! EFFECTIVE PARAMETERS
    !-------------------------------------------
    call unpack_field_and_write(nc, "L1_fSealed", &
            (/rows1, cols1, lcscenes/), nodata_dp, L1_fSealed(s1 : e1, 1, :), mask1, &
            "fraction of Sealed area at level 1")

    call unpack_field_and_write(nc, "L1_alpha", &
            (/rows1, cols1/), nodata_dp, L1_alpha(s1 : e1, 1, 1), mask1, &
            "exponent for the upper reservoir at level 1")

    call unpack_field_and_write(nc, "L1_degDayInc", &
            (/rows1, cols1, lcscenes/), nodata_dp, L1_degDayInc(s1 : e1, 1, :), mask1, &
            "increase of the Degree-day factor per mm of increase in precipitation at level 1")

    call unpack_field_and_write(nc, "L1_degDayMax", &
            (/rows1, cols1, lcscenes/), nodata_dp, L1_degDayMax(s1 : e1, 1, :), mask1, &
            "maximum degree-day factor at level 1")

    call unpack_field_and_write(nc, "L1_degDayNoPre", &
            (/rows1, cols1, lcscenes/), nodata_dp, L1_degDayNoPre(s1 : e1, 1, :), mask1, &
            "degree-day factor with no precipitation at level 1")

    call unpack_field_and_write(nc, "L1_degDay", &
            (/rows1, cols1, lcscenes/), nodata_dp, L1_degDay(s1 : e1, 1, :), mask1, &
            "degree-day factor with no precipitation at level 1")

    call unpack_field_and_write(nc, "L1_karstLoss", &
            (/rows1, cols1/), nodata_dp, L1_karstLoss(s1 : e1, 1, 1), mask1, &
            "Karstic percolation loss at level 1")

    call unpack_field_and_write(nc, "L1_fRoots", &
            (/rows1, cols1, soil1, lcscenes/), nodata_dp, L1_fRoots(s1 : e1, :, :), mask1, &
            "Fraction of roots in soil horizons at level 1")

    call unpack_field_and_write(nc, "L1_maxInter", &
            (/rows1, cols1, lais/), nodata_dp, L1_maxInter(s1 : e1, :, 1), mask1, &
            "Maximum interception at level 1")

    call unpack_field_and_write(nc, "L1_kfastFlow", &
            (/rows1, cols1, lcscenes/), nodata_dp, L1_kfastFlow(s1 : e1, 1, :), mask1, &
            "fast interflow recession coefficient at level 1")

    call unpack_field_and_write(nc, "L1_kSlowFlow", &
            (/rows1, cols1/), nodata_dp, L1_kSlowFlow(s1 : e1, 1, 1), mask1, &
            "slow interflow recession coefficient at level 1")

    call unpack_field_and_write(nc, "L1_kBaseFlow", &
            (/rows1, cols1/), nodata_dp, L1_kBaseFlow(s1 : e1, 1, 1), mask1, &
            "baseflow recession coefficient at level 1")

    call unpack_field_and_write(nc, "L1_kPerco", &
            (/rows1, cols1/), nodata_dp, L1_kPerco(s1 : e1, 1, 1), mask1, &
            "percolation coefficient at level 1")

    call unpack_field_and_write(nc, "L1_soilMoistFC", &
            (/rows1, cols1, soil1, lcscenes/), nodata_dp, L1_soilMoistFC(s1 : e1, :, :), mask1, &
            "SM below which actual ET is reduced linearly till PWP at level 1 for processCase(3)=1")

    call unpack_field_and_write(nc, "L1_soilMoistSat", &
            (/rows1, cols1, soil1, lcscenes/), nodata_dp, L1_soilMoistSat(s1 : e1, :, :), mask1, &
            "Saturation soil moisture for each horizon [mm] at level 1")

    call unpack_field_and_write(nc, "L1_soilMoistExp", &
            (/rows1, cols1, soil1, lcscenes/), nodata_dp, L1_soilMoistExp(s1 : e1, :, :), mask1, &
            "Exponential parameter to how non-linear is the soil water retention at level 1")

    if (processMatrix(3, 1) == 2) then
      call unpack_field_and_write(nc, "L1_jarvis_thresh_c1", &
              (/rows1, cols1/), nodata_dp, L1_jarvis_thresh_c1(s1 : e1, 1, 1), mask1, &
              "jarvis critical value for normalized soil water content")
    end if

    if (processMatrix(5, 1) == -1) then
      call unpack_field_and_write(nc, "L1_petLAIcorFactor", &
              (/rows1, cols1, lais, lcscenes/), nodata_dp, L1_petLAIcorFactor(s1 : e1, :, :), mask1, &
              "PET correction factor based on LAI")
    end if

    call unpack_field_and_write(nc, "L1_tempThresh", &
            (/rows1, cols1, lcscenes/), nodata_dp, L1_tempThresh(s1 : e1, 1, :), mask1, &
            "Threshold temperature for snow/rain at level 1")

    call unpack_field_and_write(nc, "L1_unsatThresh", &
            (/rows1, cols1/), nodata_dp, L1_unsatThresh(s1 : e1, 1, 1), mask1, &
            "Threshold water depth controlling fast interflow at level 1")

    call unpack_field_and_write(nc, "L1_sealedThresh", &
            (/rows1, cols1/), nodata_dp, L1_sealedThresh(s1 : e1, 1, 1), mask1, &
            "Threshold water depth for surface runoff in sealed surfaces at level 1")

    call unpack_field_and_write(nc, "L1_wiltingPoint", &
            (/rows1, cols1, soil1, lcscenes/), nodata_dp, L1_wiltingPoint(s1 : e1, :, :), mask1, &
            "Permanent wilting point at level 1")

    select case (processMatrix(5, 1))
    case(-1 : 0) ! PET is input
      call unpack_field_and_write(nc, "L1_fAsp", &
              (/rows1, cols1/), nodata_dp, L1_fAsp(s1 : e1, 1, 1), mask1, &
              "PET correction factor due to terrain aspect at level 1")

    case(1) ! Hargreaves-Samani
      call unpack_field_and_write(nc, "L1_fAsp", &
              (/rows1, cols1/), nodata_dp, L1_fAsp(s1 : e1, 1, 1), mask1, &
              "PET correction factor due to terrain aspect at level 1")

      call unpack_field_and_write(nc, "L1_HarSamCoeff", &
              (/rows1, cols1/), nodata_dp, L1_HarSamCoeff(s1 : e1, 1, 1), mask1, &
              "Hargreaves-Samani coefficient")

    case(2) ! Priestley-Taylor
      call unpack_field_and_write(nc, "L1_PrieTayAlpha", &
              (/rows1, cols1, lais/), nodata_dp, L1_PrieTayAlpha(s1 : e1, :, 1), mask1, &
              "Priestley Taylor coeffiecient (alpha)")

    case(3) ! Penman-Monteith
      call unpack_field_and_write(nc, "L1_aeroResist", &
              (/rows1, cols1, lais, lcscenes/), nodata_dp, L1_aeroResist(s1 : e1, :, :), mask1, &
              "aerodynamical resitance")

      call unpack_field_and_write(nc, "L1_surfResist", &
              (/rows1, cols1, lais/), nodata_dp, L1_surfResist(s1 : e1, :, 1), mask1, &
              "bulk surface resitance")

    end select

  end subroutine write_eff_params

  subroutine unpack_field_and_write_1d_i4 (nc, var_name, var_dims, fill_value, data, mask, var_long_name)

    use mo_kind, only : i4
    use mo_netcdf, only : NcDataset, NcDimension, NcVariable

    implicit none
    ! intent inout
    type(NcDataset), intent(inout) :: nc             ! NcDataset to add variable to
    ! intent in
    character(*), intent(in) :: var_name       ! variable name
    type(NcDimension), dimension(:), intent(in) :: var_dims       ! vector of Variable dimensions
    integer(i4), intent(in) :: fill_value     ! fill value used for missing values
    integer(i4), dimension(:), intent(in) :: data           ! packed data to be set to variable
    logical, dimension(:, :), intent(in) :: mask           ! mask used for unpacking
    character(*), optional, intent(in) :: var_long_name  ! variable long name attribute
    ! local variables
    type(NcVariable) :: var

    ! set variable
    var = nc%setVariable(var_name, "i32", var_dims)
    call var%setFillValue(fill_value)

    ! set the unpacked data
    call var%setData(unpack(data, mask, fill_value))

    ! optionally set attributes
    if (present(var_long_name)) then
      call var%setAttribute("long_name", trim(var_long_name))
    end if

  end subroutine

  subroutine unpack_field_and_write_1d_dp (nc, var_name, var_dims, fill_value, data, mask, var_long_name)

    use mo_kind, only : dp
    use mo_netcdf, only : NcDataset, NcDimension, NcVariable

    implicit none
    ! intent inout
    type(NcDataset), intent(inout) :: nc             ! NcDataset to add variable to
    ! intent in
    character(*), intent(in) :: var_name       ! variable name
    type(NcDimension), dimension(:), intent(in) :: var_dims       ! vector of Variable dimensions
    real(dp), intent(in) :: fill_value     ! fill value used for missing values
    real(dp), dimension(:), intent(in) :: data           ! packed data to be set to variable
    logical, dimension(:, :), intent(in) :: mask           ! mask used for unpacking
    character(*), optional, intent(in) :: var_long_name  ! variable long name attribute
    ! local variables
    type(NcVariable) :: var

    ! set variable
    var = nc%setVariable(var_name, "f64", var_dims)
    call var%setFillValue(fill_value)

    ! set the unpacked data
    call var%setData(unpack(data, mask, fill_value))

    ! optionally set attributes
    if (present(var_long_name)) then
      call var%setAttribute("long_name", trim(var_long_name))
    end if

  end subroutine

  subroutine unpack_field_and_write_2d_dp (nc, var_name, var_dims, fill_value, data, mask, var_long_name)

    use mo_kind, only : dp, i4
    use mo_netcdf, only : NcDataset, NcDimension, NcVariable

    implicit none
    ! intent inout
    type(NcDataset), intent(inout) :: nc             ! NcDataset to add variable to
    ! intent in
    character(*), intent(in) :: var_name       ! variable name
    type(NcDimension), dimension(:), intent(in) :: var_dims       ! vector of Variable dimensions
    real(dp), intent(in) :: fill_value     ! fill value used for missing values
    real(dp), dimension(:, :), intent(in) :: data           ! packed data to be set to variable
    logical, dimension(:, :), intent(in) :: mask           ! mask used for unpacking
    character(*), optional, intent(in) :: var_long_name  ! variable long name attribute
    ! local variables
    type(NcVariable) :: var
    real(dp), dimension(:, :, :), allocatable :: dummy_arr
    integer(i4), dimension(3) :: dim_length
    integer(i4) :: ii

    ! set variable
    var = nc%setVariable(var_name, "f64", var_dims)
    call var%setFillValue(fill_value)

    dim_length = var%getShape()
    allocate(dummy_arr(dim_length(1), dim_length(2), dim_length(3)))
    do ii = 1, size(data, 2)
      dummy_arr(:, :, ii) = unpack(data(:, ii), mask, fill_value)
    end do

    ! set the unpacked data
    call var%setData(dummy_arr)

    ! optionally set attributes
    if (present(var_long_name)) then
      call var%setAttribute("long_name", trim(var_long_name))
    end if

  end subroutine

  subroutine unpack_field_and_write_3d_dp (nc, var_name, var_dims, fill_value, data, mask, var_long_name)

    use mo_kind, only : dp, i4
    use mo_netcdf, only : NcDataset, NcDimension, NcVariable

    implicit none
    ! intent inout
    type(NcDataset), intent(inout) :: nc             ! NcDataset to add variable to
    ! intent in
    character(*), intent(in) :: var_name       ! variable name
    type(NcDimension), dimension(:), intent(in) :: var_dims       ! vector of Variable dimensions
    real(dp), intent(in) :: fill_value     ! fill value used for missing values
    real(dp), dimension(:, :, :), intent(in) :: data           ! packed data to be set to variable
    logical, dimension(:, :), intent(in) :: mask           ! mask used for unpacking
    character(*), optional, intent(in) :: var_long_name  ! variable long name attribute
    ! local variables
    type(NcVariable) :: var
    real(dp), dimension(:, :, :, :), allocatable :: dummy_arr
    integer(i4), dimension(4) :: dim_length
    integer(i4) :: ii, jj

    ! set variable
    var = nc%setVariable(var_name, "f64", var_dims)
    call var%setFillValue(fill_value)

    dim_length = var%getShape()
    allocate(dummy_arr(dim_length(1), dim_length(2), dim_length(3), dim_length(4)))
    do ii = 1, size(data, 2)
      do jj = 1, size(data, 3)
        dummy_arr(:, :, ii, jj) = unpack(data(:, ii, jj), mask, fill_value)
      end do
    end do

    ! set the unpacked data
    call var%setData(dummy_arr)

    ! optionally set attributes
    if (present(var_long_name)) then
      call var%setAttribute("long_name", trim(var_long_name))
    end if

  end subroutine

END MODULE mo_mpr_restart
