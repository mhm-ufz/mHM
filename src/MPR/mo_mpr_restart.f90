!> \file mo_mpr_restart.f90
!> \brief \copybrief mo_mpr_restart
!> \details \copydetails mo_mpr_restart

!> \brief reading and writing states, fluxes and configuration for restart of mHM.
!> \details routines are seperated for reading and writing variables for:
!!  - states and fluxes, and
!!  - configuration.
!!
!! Reading of L11 configuration is also seperated from the rest, since it is only required when routing is activated.
!> \authors Stephan Thober
!> \date Jul 2013
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mpr
MODULE mo_mpr_restart

  ! This module is a restart for the UFZ CHS mesoscale hydrologic model mHM.

  ! Written  Stephan Thober, Apr 2011

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: write_eff_params            ! read restart files for configuration from a given path
  PUBLIC :: write_mpr_restart_files     ! write restart files for configuration to a given path


  !> \brief unpack parameter fields and write them to file
  !> \param[inout] "type(NcDataset) :: nc"                    NcDataset to add variable to
  !> \param[in] "character(*) :: var_name"                    variable name
  !> \param[in] "type(NcDimension), dimension(:) :: var_dims" vector of Variable dimensions
  !> \param[in] "integer(i4) :: fill_value"                   fill value used for missing values
  !> \param[in] "integer(i4/dp), dimension(...) :: data"      packed data to be set to variable
  !> \param[in] "logical, dimension(:, :) :: mask"            mask used for unpacking
  !> \param[in] "character(*), optional :: var_long_name"     variable long name attribute
  !> \authors Robert Schweppe
  !> \date Jun 2018
  INTERFACE unpack_field_and_write
    MODULE PROCEDURE unpack_field_and_write_1d_i4, &
            unpack_field_and_write_1d_dp, &
            unpack_field_and_write_2d_dp, &
            unpack_field_and_write_3d_dp
  end interface unpack_field_and_write


CONTAINS

  !> \brief write restart files for each domain
  !> \details write restart files for each domain. For each domain
  !! three restart files are written. These are xxx_states.nc,
  !! xxx_L11_config.nc, and xxx_config.nc (xxx being the three digit
  !! domain index). If a variable is added here, it should also be added
  !! in the read restart routines below.
  !! ADDITIONAL INFORMATION
  !! write_restart
  !> \changelog
  !! - Stephan Thober Aug 2015
  !!   - moved write of routing states to mRM
  !! - David Schaefer Nov 2015
  !!   - mo_netcdf
  !! - Stephan Thober Nov 2016
  !!   - moved processMatrix to common variables
  !! - Zink M. Demirel C. Mar 2017
  !!   - Added Jarvis soil water stress function at SM process(3)
  !> \authors Stephan Thober
  !> \date Jun 2014
  subroutine write_mpr_restart_files(OutFile)

    use mo_common_restart, only : write_grid_info
    use mo_common_variables, only : level1, nLCoverScene, domainMeta, LC_year_start, LC_year_end
    use mo_kind, only : i4, dp
    use mo_message, only : message
    use mo_mpr_global_variables, only : nLAI, nSoilHorizons_mHM, HorizonDepth_mHM
    use mo_netcdf, only : NcDataset, NcDimension
    use mo_string_utils, only : num2str
    use mo_common_constants, only : soilHorizonsVarName, landCoverPeriodsVarName, LAIVarName

    implicit none

    character(256) :: Fname

    !> Output Path for each domain
    character(256), dimension(:), intent(in) :: OutFile

    integer(i4) :: iDomain, domainID

    ! start index at level 1
    integer(i4) :: s1

    ! end index at level 1
    integer(i4) :: e1

    ! mask at level 1
    logical, dimension(:, :), allocatable :: mask1

    type(NcDataset) :: nc

    type(NcDimension) :: rows1, cols1, soil1, lcscenes, lais

    real(dp), dimension(:), allocatable :: dummy_1D


    domain_loop : do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)

      ! write restart file for iDomain
      Fname = trim(OutFile(iDomain))
      ! print a message
      call message("    Writing Restart-file: ", trim(adjustl(Fname)), " ...")

      nc = NcDataset(fname, "w")

      call write_grid_info(level1(iDomain), "1", nc)

      rows1 = nc%getDimension("nrows1")
      cols1 = nc%getDimension("ncols1")

      ! write the dimension to the file and also save bounds
      allocate(dummy_1D(nSoilHorizons_mHM+1))
      dummy_1D(1) = 0.0_dp
      dummy_1D(2:nSoilHorizons_mHM+1) = HorizonDepth_mHM(:)
      soil1 = nc%setCoordinate(trim(soilHorizonsVarName), nSoilHorizons_mHM, dummy_1D, 2_i4)
      deallocate(dummy_1D)
      allocate(dummy_1D(nLCoverScene+1))
      dummy_1D(1:nLCoverScene) = LC_year_start(:)
      ! this is done because bounds are always stored as real so e.g.
      ! 1981-1990,1991-2000 is thus saved as 1981.0-1991.0,1991.0-2001.0
      ! it is translated back into ints correctly during reading
      dummy_1D(nLCoverScene+1) = LC_year_end(nLCoverScene) + 1
      lcscenes = nc%setCoordinate(trim(landCoverPeriodsVarName), nLCoverScene, dummy_1D, 0_i4)
      deallocate(dummy_1D)
      ! write the dimension to the file
      lais = nc%setDimension(trim(LAIVarName), nLAI)

     ! for appending and intialization
      allocate(mask1(rows1%getLength(), cols1%getLength()))
      s1 = level1(iDomain)%iStart
      e1 = level1(iDomain)%iEnd
      mask1 = level1(iDomain)%mask

      call write_eff_params(mask1, s1, e1, rows1, cols1, soil1, lcscenes, lais, nc)
      deallocate(mask1)
      call nc%close()

    end do domain_loop

  end subroutine write_mpr_restart_files


  !> \brief write effective parameter fields to given restart file
  !> \changelog
  !!  - Rohini Kumar Oct 2021
  !!    - Added Neutron count module to mHM integrate into develop branch (5.11.2)
  !!  - Sebastian Müller Mar 2023
  !!    - made L1_alpha, L1_kSlowFlow, L1_kBaseFlow and L1_kPerco land cover dependent
  !> \authors Robert Schweppe
  !> \date Jun 2018
  subroutine write_eff_params(mask1, s1, e1, rows1, cols1, soil1, lcscenes, lais, nc)

    use mo_common_constants, only : nodata_dp, nodata_i4
    use mo_common_variables, only : LC_year_end, LC_year_start, processMatrix
    use mo_kind, only : i4
    use mo_mpr_global_variables, only : L1_HarSamCoeff, L1_PrieTayAlpha, L1_aeroResist, &
                                        L1_alpha, L1_degDay, L1_degDayInc, L1_degDayMax, L1_degDayNoPre, L1_fAsp, &
                                        L1_fRoots, L1_fSealed, L1_jarvis_thresh_c1, L1_kBaseFlow, L1_kPerco, &
                                        L1_kSlowFlow, L1_karstLoss, L1_kfastFlow, L1_maxInter, L1_petLAIcorFactor, &
                                        L1_sealedThresh, L1_soilMoistExp, L1_soilMoistFC, L1_soilMoistSat, L1_surfResist, &
                                        L1_tempThresh, L1_unsatThresh, L1_wiltingPoint, &
                                        ! neutron count
                                        L1_No_Count, L1_bulkDens, L1_latticeWater, L1_COSMICL3

    use mo_netcdf, only : NcDataset, NcDimension, NcVariable

    implicit none

    logical, dimension(:, :), allocatable, intent(in) :: mask1 !< mask at level 1
    integer(i4), intent(in) :: s1 !< start index at level 1
    integer(i4), intent(in) :: e1 !< end index at level 1
    type(NcDimension), intent(in) :: rows1 !< y dimension
    type(NcDimension), intent(in) :: cols1 !< x dimension
    type(NcDimension), intent(in) :: soil1 !< soil dimension
    type(NcDimension), intent(in) :: lcscenes !< land conver scenes dimension
    type(NcDimension), intent(in) :: lais !< LAI dimension
    type(NcDataset), intent(inout) :: nc !< NetCDF file to write to

    type(NcVariable) :: var

    !-------------------------------------------
    ! EFFECTIVE PARAMETERS
    !-------------------------------------------
    call unpack_field_and_write(nc, "L1_fSealed", &
            (/rows1, cols1, lcscenes/), nodata_dp, L1_fSealed(s1 : e1, 1, :), mask1, &
            "fraction of Sealed area at level 1")

    call unpack_field_and_write(nc, "L1_alpha", &
            (/rows1, cols1, lcscenes/), nodata_dp, L1_alpha(s1 : e1, 1, :), mask1, &
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
            (/rows1, cols1, lcscenes/), nodata_dp, L1_kSlowFlow(s1 : e1, 1, :), mask1, &
            "slow interflow recession coefficient at level 1")

    call unpack_field_and_write(nc, "L1_kBaseFlow", &
            (/rows1, cols1, lcscenes/), nodata_dp, L1_kBaseFlow(s1 : e1, 1, :), mask1, &
            "baseflow recession coefficient at level 1")

    call unpack_field_and_write(nc, "L1_kPerco", &
            (/rows1, cols1, lcscenes/), nodata_dp, L1_kPerco(s1 : e1, 1, :), mask1, &
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

    if (any(processMatrix(3, 1) == (/2, 3/))) then
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

   ! neutron count
   select case (processMatrix(10, 1))
   case(1) ! deslet
      call unpack_field_and_write(nc, "L1_No_Count", &
           (/rows1, cols1/), nodata_dp, L1_No_Count(s1:e1, 1, 1), mask1, &
           "N0 count at level 1")
      call unpack_field_and_write(nc, "L1_bulkDens", &
           (/rows1, cols1, soil1, lcscenes/), nodata_dp, L1_bulkDens(s1:e1, :, :), mask1, &
           "Bulk density at level 1 for processCase(10)")
      call unpack_field_and_write(nc, "L1_latticeWater", &
           (/rows1, cols1, soil1, lcscenes/), nodata_dp, L1_latticeWater(s1:e1, :, :), mask1, &
           "Lattice water content at level 1 for processCase(10)")

   case(2) ! COSMIC
      call unpack_field_and_write(nc, "L1_No_Count", &
           (/rows1, cols1/), nodata_dp, L1_No_Count(s1 : e1, 1, 1), mask1, &
           "N0 count at level 1")
      call unpack_field_and_write(nc, "L1_bulkDens", &
           (/rows1, cols1, soil1, lcscenes/), nodata_dp, L1_bulkDens(s1 : e1, :, :), mask1, &
           "Bulk density at level 1 for processCase(10)")
      call unpack_field_and_write(nc, "L1_latticeWater", &
           (/rows1, cols1, soil1, lcscenes/), nodata_dp, L1_latticeWater(s1 : e1, :, :), mask1, &
           "Lattice water content at level 1 for processCase(10)")
      call unpack_field_and_write(nc, "L1_COSMICL3", &
           (/rows1, cols1, soil1, lcscenes/), nodata_dp, L1_COSMICL3(s1 : e1, :, :), mask1, &
           "COSMIC L3 parameter at level 1 for processCase(10)")
   end select

  end subroutine write_eff_params

  !> \copydoc unpack_field_and_write
  subroutine unpack_field_and_write_1d_i4(nc, var_name, var_dims, fill_value, data, mask, var_long_name)

    use mo_kind, only : i4
    use mo_netcdf, only : NcDataset, NcDimension, NcVariable

    implicit none

    ! NcDataset to add variable to
    type(NcDataset), intent(inout) :: nc

    ! variable name
    character(*), intent(in) :: var_name

    ! vector of Variable dimensions
    type(NcDimension), dimension(:), intent(in) :: var_dims

    ! fill value used for missing values
    integer(i4), intent(in) :: fill_value

    ! packed data to be set to variable
    integer(i4), dimension(:), intent(in) :: data

    ! mask used for unpacking
    logical, dimension(:, :), intent(in) :: mask

    ! variable long name attribute
    character(*), optional, intent(in) :: var_long_name

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

  !> \copydoc unpack_field_and_write
  subroutine unpack_field_and_write_1d_dp(nc, var_name, var_dims, fill_value, data, mask, var_long_name)

    use mo_kind, only : dp
    use mo_netcdf, only : NcDataset, NcDimension, NcVariable

    implicit none

    ! NcDataset to add variable to
    type(NcDataset), intent(inout) :: nc

    ! variable name
    character(*), intent(in) :: var_name

    ! vector of Variable dimensions
    type(NcDimension), dimension(:), intent(in) :: var_dims

    ! fill value used for missing values
    real(dp), intent(in) :: fill_value

    ! packed data to be set to variable
    real(dp), dimension(:), intent(in) :: data

    ! mask used for unpacking
    logical, dimension(:, :), intent(in) :: mask

    ! variable long name attribute
    character(*), optional, intent(in) :: var_long_name

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

  !> \copydoc unpack_field_and_write
  subroutine unpack_field_and_write_2d_dp(nc, var_name, var_dims, fill_value, data, mask, var_long_name)

    use mo_kind, only : dp, i4
    use mo_netcdf, only : NcDataset, NcDimension, NcVariable

    implicit none

    ! NcDataset to add variable to
    type(NcDataset), intent(inout) :: nc

    ! variable name
    character(*), intent(in) :: var_name

    ! vector of Variable dimensions
    type(NcDimension), dimension(:), intent(in) :: var_dims

    ! fill value used for missing values
    real(dp), intent(in) :: fill_value

    ! packed data to be set to variable
    real(dp), dimension(:, :), intent(in) :: data

    ! mask used for unpacking
    logical, dimension(:, :), intent(in) :: mask

    ! variable long name attribute
    character(*), optional, intent(in) :: var_long_name

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

  !> \copydoc unpack_field_and_write
  subroutine unpack_field_and_write_3d_dp(nc, var_name, var_dims, fill_value, data, mask, var_long_name)

    use mo_kind, only : dp, i4
    use mo_netcdf, only : NcDataset, NcDimension, NcVariable

    implicit none

    ! NcDataset to add variable to
    type(NcDataset), intent(inout) :: nc

    ! variable name
    character(*), intent(in) :: var_name

    ! vector of Variable dimensions
    type(NcDimension), dimension(:), intent(in) :: var_dims

    ! fill value used for missing values
    real(dp), intent(in) :: fill_value

    ! packed data to be set to variable
    real(dp), dimension(:, :, :), intent(in) :: data

    ! mask used for unpacking
    logical, dimension(:, :), intent(in) :: mask

    ! variable long name attribute
    character(*), optional, intent(in) :: var_long_name

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
