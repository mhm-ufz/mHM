!>       \file mo_restart.f90

!>       \brief reading and writing states, fluxes and configuration for restart of mHM.

!>       \details routines are seperated for reading and writing variables for:
!>       - states and fluxes, and
!>       - configuration.
!>       Reading of L11 configuration is also seperated from the rest,
!>       since it is only required when routing is activated.

!>       \authors Stephan Thober

!>       \date Jul 2013

! Modifications:

MODULE mo_restart

  ! This module is a restart for the UFZ CHS mesoscale hydrologic model mHM.

  ! Written  Stephan Thober, Apr 2011
  use mo_common_constants, only : soilHorizonsVarName, landCoverPeriodsVarName, LAIVarName

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_restart_states     ! read restart files for state variables from a given path
  PUBLIC :: write_restart_files     ! write restart files for configuration to a given path

  !    NAME
  !        unpack_field_and_write

  !    PURPOSE
  !>       \brief TODO: add description

  !>       \details TODO: add description

  !    INTENT(INOUT)
  !>       \param[inout] "type(NcDataset) :: nc" NcDataset to add variable to

  !    INTENT(IN)
  !>       \param[in] "character(*) :: var_name"                    variable name
  !>       \param[in] "type(NcDimension), dimension(:) :: var_dims" vector of Variable dimensions
  !>       \param[in] "integer(i4) :: fill_value"                   fill value used for missing values
  !>       \param[in] "integer(i4), dimension(:) :: data"           packed data to be set to variable
  !>       \param[in] "logical, dimension(:, :) :: mask"            mask used for unpacking

  !    INTENT(IN), OPTIONAL
  !>       \param[in] "character(*), optional :: var_long_name" variable long name attribute

  !    HISTORY
  !>       \authors Robert Schweppe

  !>       \date Jun 2018

  ! Modifications:


  INTERFACE unpack_field_and_write
    MODULE PROCEDURE unpack_field_and_write_1d_i4, &
            unpack_field_and_write_1d_dp, &
            unpack_field_and_write_2d_dp, &
            unpack_field_and_write_3d_dp
  end interface unpack_field_and_write


CONTAINS
  ! ------------------------------------------------------------------

  !    NAME
  !        write_restart_files

  !    PURPOSE
  !>       \brief write restart files for each domain

  !>       \details write restart files for each domain. For each domain
  !>       three restart files are written. These are xxx_states.nc,
  !>       xxx_L11_config.nc, and xxx_config.nc (xxx being the three digit
  !>       domain index). If a variable is added here, it should also be added
  !>       in the read restart routines below.

  !    INTENT(IN)
  !>       \param[in] "character(256), dimension(:) :: OutFile" Output Path for each domain

  !    HISTORY
  !>       \authors Stephan Thober

  !>       \date Jun 2014

  ! Modifications:
  ! Stephan Thober     Aug  2015 - moved write of routing states to mRM
  ! David Schaefer     Nov  2015 - mo_netcdf
  ! Stephan Thober     Nov  2016 - moved processMatrix to common variables
  ! Zink M. Demirel C. Mar 2017 - Added Jarvis soil water stress function at SM process(3)
  ! Robert Schweppe    Feb 2018 - Removed all L0 references
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine write_restart_files(OutFile)

    use mo_common_constants, only : nodata_dp, nodata_i4
    use mo_grid, only : write_grid_info
    use mo_common_variables, only : level1, nLandCoverPeriods, domainMeta, landCoverPeriodBoundaries
    use mo_common_datetime_type, only: simPer, LCyearId
    use mo_global_variables, only : L1_Inter, L1_Throughfall, L1_aETCanopy, L1_aETSealed, L1_aETSoil, L1_baseflow, &
                                    L1_fastRunoff, L1_infilSoil, L1_melt, L1_percol, L1_preEffect, L1_rain, &
                                    L1_runoffSeal, L1_satSTW, L1_sealSTW, L1_slowRunoff, L1_snow, L1_snowPack, &
                                    L1_soilMoist, L1_total_runoff, L1_unsatSTW, L1_degDay, &
                                    nLAIs, nSoilHorizons, soilHorizonBoundaries, LAIBoundaries
    use mo_kind, only : dp, i4
    use mo_message, only : message
    use mo_netcdf, only : NcDataset, NcDimension, NcVariable
    use mo_string_utils, only : num2str

    implicit none

    character(256) :: Fname

    ! Output Path for each domain
    character(256), dimension(:), intent(in) :: OutFile

    integer(i4) :: iDomain, domainID, iYear, iBoundary

    integer(i4) :: ii

    ! start index at level 1
    integer(i4) :: s1

    ! end index at level 1
    integer(i4) :: e1

    ! number of landcoverperiods for current domain
    integer(i4) :: iDomainNLandCoverPeriods

    ! mask at level 1
    logical, dimension(:, :), allocatable :: mask1
    logical, dimension(:, :, :), allocatable :: mask_soil1

    ! dummy variable
    real(dp), dimension(:, :, :), allocatable :: dummy_3D
    real(dp), dimension(:), allocatable :: dummy_1D

    integer(i4) :: max_extent

    type(NcDataset) :: nc

    type(NcDimension) :: rows1, cols1, soil1, lcscenes, lais

    type(NcVariable) :: var


    ! get maximum extent of one dimension 2 or 3
    max_extent = max(nSoilHorizons, nLandCoverPeriods, nLAIs)

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
      soil1 = nc%setCoordinate(trim(soilHorizonsVarName), nSoilHorizons, soilHorizonBoundaries, 2_i4)
      ! write the dimension to the file
      lais = nc%setCoordinate(trim(LAIVarName), nLAIs, LAIBoundaries, 0_i4)

      iDomainNLandCoverPeriods = maxval(LCyearId(:, iDomain), mask=LCyearId(:, iDomain) /= nodata_i4)
      allocate(dummy_1D(size(landCoverPeriodBoundaries, dim=1)))
      dummy_1D = real(landCoverPeriodBoundaries(:, iDomain), dp)
      lcscenes = nc%setCoordinate(trim(landCoverPeriodsVarName), iDomainNLandCoverPeriods, &
              dummy_1D, 0_i4)
      deallocate(dummy_1D)

      ! for appending and intialization
      allocate(mask1(rows1%getLength(), cols1%getLength()))
      s1 = level1(iDomain)%iStart
      e1 = level1(iDomain)%iEnd

      mask1 = level1(iDomain)%mask
      allocate(mask_soil1 (rows1%getLength(), cols1%getLength(), nSoilHorizons))
      mask_soil1 = spread(mask1, 3, nSoilHorizons)

      var = nc%setVariable("L1_Inter", "f64", [rows1, cols1])
      call var%setFillValue(nodata_dp)
      call var%setData(unpack(L1_inter(s1 : e1), mask1, nodata_dp))
      call var%setAttribute("long_name", "Interception storage at level 1")

      var = nc%setVariable("L1_snowPack", "f64", [rows1, cols1])
      call var%setFillValue(nodata_dp)
      call var%setData(unpack(L1_snowPack(s1 : e1), mask1, nodata_dp))
      call var%setAttribute("long_name", "Snowpack at level 1")

      var = nc%setVariable("L1_sealSTW", "f64", [rows1, cols1])
      call var%setFillValue(nodata_dp)
      call var%setData(unpack(L1_sealSTW(s1 : e1), mask1, nodata_dp))
      call var%setAttribute("long_name", "Retention storage of impervious areas at level 1")

      var = nc%setVariable("L1_soilMoist", "f64", [rows1, cols1, soil1])
      call var%setFillValue(nodata_dp)
      call var%setData(unpack(pack(L1_soilMoist(s1 : e1, :), .true.), mask_soil1, nodata_dp))
      call var%setAttribute("long_name", "soil moisture at level 1")

      var = nc%setVariable("L1_unsatSTW", "f64", [rows1, cols1])
      call var%setFillValue(nodata_dp)
      call var%setData(unpack(L1_unsatSTW(s1 : e1), mask1, nodata_dp))
      call var%setAttribute("long_name", "upper soil storage at level 1")

      var = nc%setVariable("L1_satSTW", "f64", [rows1, cols1])
      call var%setFillValue(nodata_dp)
      call var%setData(unpack(L1_satSTW(s1 : e1), mask1, nodata_dp))
      call var%setAttribute("long_name", "groundwater storage at level 1")

      var = nc%setVariable("L1_aETSoil", "f64", [rows1, cols1, soil1])
      call var%setFillValue(nodata_dp)
      call var%setData(unpack(pack(L1_aETSoil(s1 : e1, :), .true.), mask_soil1, nodata_dp))
      call var%setAttribute("long_name", "soil actual ET at level 1")

      var = nc%setVariable("L1_aETCanopy", "f64", [rows1, cols1])
      call var%setFillValue(nodata_dp)
      call var%setData(unpack(L1_aETCanopy(s1 : e1), mask1, nodata_dp))
      call var%setAttribute("long_name", "canopy actual ET at level 1")

      var = nc%setVariable("L1_aETSealed", "f64", [rows1, cols1])
      call var%setFillValue(nodata_dp)
      call var%setData(unpack(L1_aETSealed(s1 : e1), mask1, nodata_dp))
      call var%setAttribute("long_name", "sealed actual ET at level 1")

      var = nc%setVariable("L1_baseflow", "f64", [rows1, cols1])
      call var%setFillValue(nodata_dp)
      call var%setData(unpack(L1_baseflow(s1 : e1), mask1, nodata_dp))
      call var%setAttribute("long_name", "baseflow at level 1")

      var = nc%setVariable("L1_infilSoil", "f64", [rows1, cols1, soil1])
      call var%setFillValue(nodata_dp)
      call var%setData(unpack(pack(L1_infilSoil(s1 : e1, :), .true.), mask_soil1, nodata_dp))
      call var%setAttribute("long_name", "soil in-exfiltration at level 1")

      var = nc%setVariable("L1_fastRunoff", "f64", [rows1, cols1])
      call var%setFillValue(nodata_dp)
      call var%setData(unpack(L1_fastRunoff(s1 : e1), mask1, nodata_dp))
      call var%setAttribute("long_name", "fast runoff")

      var = nc%setVariable("L1_percol", "f64", [rows1, cols1])
      call var%setFillValue(nodata_dp)
      call var%setData(unpack(L1_percol(s1 : e1), mask1, nodata_dp))
      call var%setAttribute("long_name", "percolation at level 1")

      var = nc%setVariable("L1_melt", "f64", [rows1, cols1])
      call var%setFillValue(nodata_dp)
      call var%setData(unpack(L1_melt(s1 : e1), mask1, nodata_dp))
      call var%setAttribute("long_name", "snow melt at level 1")

      var = nc%setVariable("L1_preEffect", "f64", [rows1, cols1])
      call var%setFillValue(nodata_dp)
      call var%setData(unpack(L1_preEffect(s1 : e1), mask1, nodata_dp))
      call var%setAttribute("long_name", "effective precip. depth (snow melt + rain) at level 1")

      var = nc%setVariable("L1_rain", "f64", [rows1, cols1])
      call var%setFillValue(nodata_dp)
      call var%setData(unpack(L1_rain(s1 : e1), mask1, nodata_dp))
      call var%setAttribute("long_name", "rain (liquid water) at level 1")

      var = nc%setVariable("L1_runoffSeal", "f64", [rows1, cols1])
      call var%setFillValue(nodata_dp)
      call var%setData(unpack(L1_runoffSeal(s1 : e1), mask1, nodata_dp))
      call var%setAttribute("long_name", "runoff from impervious area at level 1")

      var = nc%setVariable("L1_slowRunoff", "f64", [rows1, cols1])
      call var%setFillValue(nodata_dp)
      call var%setData(unpack(L1_slowRunoff(s1 : e1), mask1, nodata_dp))
      call var%setAttribute("long_name", "slow runoff at level 1")

      var = nc%setVariable("L1_snow", "f64", [rows1, cols1])
      call var%setFillValue(nodata_dp)
      call var%setData(unpack(L1_snow(s1 : e1), mask1, nodata_dp))
      call var%setAttribute("long_name", "snow (solid water) at level 1")

      var = nc%setVariable("L1_Throughfall", "f64", [rows1, cols1])
      call var%setFillValue(nodata_dp)
      call var%setData(unpack(L1_Throughfall(s1 : e1), mask1, nodata_dp))
      call var%setAttribute("long_name", "throughfall at level 1")

      var = nc%setVariable("L1_total_runoff", "f64", [rows1, cols1])
      call var%setFillValue(nodata_dp)
      call var%setData(unpack(L1_total_runoff(s1 : e1), mask1, nodata_dp))
      call var%setAttribute("long_name", "total runoff at level 1")

       var = nc%setVariable("L1_degDay", "f64", [rows1, cols1])
      call var%setFillValue(nodata_dp)
      call var%setData(unpack(L1_degDay(s1 : e1), mask1, nodata_dp))
      call var%setAttribute("long_name", "degree-day factor with no precipitation at level 1")

      call write_eff_params(mask1, s1, e1, rows1, cols1, soil1, lcscenes, lais, nc, iDomainNLandCoverPeriods)

      call nc%close()

      deallocate(mask1, mask_soil1)
    end do domain_loop

  end subroutine write_restart_files

  ! ------------------------------------------------------------------

  !    NAME
  !        read_restart_states

  !    PURPOSE
  !>       \brief reads fluxes and state variables from file

  !>       \details read fluxes and state variables from given
  !>       restart directory and initialises all state variables
  !>       that are initialized in the subroutine initialise,
  !>       contained in module mo_startup.

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iDomain"    number of domains
  !>       \param[in] "character(256) :: InFile" Input Path including trailing slash

  !    HISTORY
  !>       \authors Stephan Thober

  !>       \date Apr 2013

  ! Modifications:
  ! Stephan Thober Aug  2015 - moved read of routing states to mRM
  ! David Schaefer Nov  2015 - mo_netcdf
  ! Stephan Thober Nov  2016 - moved processMatrix to common variables
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine read_restart_states(iDomain, InFile, do_read_states_arg, do_read_dims_arg)

    use mo_common_variables, only : level1, nLandCoverPeriods, processMatrix
    use mo_kind, only : dp, i4
    use mo_global_variables, only : L1_Inter, L1_Throughfall, L1_aETCanopy, &
                                    L1_aETSealed, L1_aETSoil, L1_baseflow, L1_fastRunoff, L1_infilSoil, L1_melt, &
                                    L1_percol, L1_preEffect, L1_rain, L1_runoffSeal, L1_satSTW, L1_sealSTW, &
                                    L1_slowRunoff, L1_snow, L1_snowPack, L1_soilMoist, L1_total_runoff, L1_unsatSTW, &
            nSoilHorizons, are_parameter_initialized, nLAIs, &
            L1_HarSamCoeff, &
                                        L1_PrieTayAlpha, L1_aeroResist, L1_alpha, L1_degDay, L1_degDayInc, L1_degDayMax, &
                                        L1_degDayNoPre, L1_fAsp, L1_fRoots, L1_fSealed, L1_jarvis_thresh_c1, &
                                    L1_kBaseFlow, L1_kPerco, L1_kSlowFlow, L1_karstLoss, L1_kFastFlow, L1_maxInter, &
                                        L1_petLAIcorFactor, L1_sealedThresh, L1_soilMoistExp, L1_soilMoistFC, &
                                        L1_soilMoistSat, L1_surfResist, L1_tempThresh, L1_unsatThresh, L1_wiltingPoint, &
                                        L1_latitude
    use mo_netcdf, only : NcDataset, NcDimension, NcVariable
    use mo_string_utils, only : num2str
    use mo_read_nc, only: check_dimension_consistency

    implicit none

    ! number of domain
    integer(i4), intent(in) :: iDomain

    ! Input Path including trailing slash
    character(256), intent(in) :: InFile

    logical, intent(in), optional :: do_read_states_arg, do_read_dims_arg

    character(256) :: Fname

    ! loop index
    integer(i4) :: ii, jj

    ! start index at level 1
    integer(i4) :: s1

    ! end index at level 1
    integer(i4) :: e1

    logical :: do_read_states, do_read_dims

    ! mask at level 1
    logical, dimension(:, :), allocatable :: mask1
    logical, dimension(:, :, :), allocatable :: mask_soil1

    ! dummy, 2 dimension
    real(dp), dimension(:, :), allocatable :: dummyD2

    ! dummy, 3 dimension
    real(dp), dimension(:, :, :), allocatable :: dummyD3

    ! dummy, 3 dimension
    real(dp), dimension(:, :, :, :), allocatable :: dummyD4

    type(NcDataset) :: nc

    type(NcVariable) :: var

    type(NcDimension) :: nc_dim

    integer(i4) :: nSoilHorizons_temp, nLAIs_temp, nLandCoverPeriods_temp
    real(dp), dimension(:), allocatable :: landCoverPeriodBoundaries_temp, soilHorizonBoundaries_temp, &
            LAIBoundaries_temp
    integer(i4), dimension(:), allocatable :: landCoverPeriodSelect
    integer(i4), dimension(:), allocatable :: varShape


    do_read_states = .true.
    if (present(do_read_states_arg)) do_read_states = do_read_states_arg
    do_read_dims = .true.
    if (present(do_read_dims_arg)) do_read_dims = do_read_dims_arg

    Fname = trim(InFile)
    ! call message('    Reading states from ', trim(adjustl(Fname)),' ...')

    ! get domain information at level 1
    allocate(mask1 (level1(iDomain)%nrows, level1(iDomain)%ncols))
    mask1 = level1(iDomain)%mask
    s1 = level1(iDomain)%iStart
    e1 = level1(iDomain)%iEnd

    nc = NcDataset(fname, "r")

    ! get the dimensions
    var = nc%getVariable(trim(soilHorizonsVarName)//'_bnds')
    call var%getData(dummyD2)
    nSoilHorizons_temp = size(dummyD2, 2)
    allocate(soilHorizonBoundaries_temp(nSoilHorizons_temp+1))
    soilHorizonBoundaries_temp(1:nSoilHorizons_temp) = dummyD2(1,:)
    soilHorizonBoundaries_temp(nSoilHorizons_temp+1) = dummyD2(2, nSoilHorizons_temp)

    ! get the landcover dimension
    var = nc%getVariable(trim(landCoverPeriodsVarName)//'_bnds')
    call var%getData(dummyD2)
    nLandCoverPeriods_temp = size(dummyD2, 2)
    allocate(landCoverPeriodBoundaries_temp(nLandCoverPeriods_temp+1))
    landCoverPeriodBoundaries_temp(1:nLandCoverPeriods_temp) = dummyD2(1,:)
    landCoverPeriodBoundaries_temp(nLandCoverPeriods_temp+1) = dummyD2(2, nLandCoverPeriods_temp)

    ! get the LAI dimension
    if (nc%hasVariable(trim(LAIVarName)//'_bnds')) then
      var = nc%getVariable(trim(LAIVarName)//'_bnds')
      call var%getData(dummyD2)
      nLAIs_temp = size(dummyD2, 2)
      allocate(LAIBoundaries_temp(nLAIs_temp+1))
      LAIBoundaries_temp(1:nLAIs_temp) = dummyD2(1,:)
      LAIBoundaries_temp(nLAIs_temp+1) = dummyD2(2,nLAIs_temp)
    else if (nc%hasDimension('L1_LAITimesteps')) then
      nc_dim = nc%getDimension('L1_LAITimesteps')
      nLAIs_temp = nc_dim%getLength()
      allocate(LAIBoundaries_temp(nLAIs_temp+1))
      LAIBoundaries_temp = [(ii, ii=1, nLAIs_temp+1)]
    end if

    call check_dimension_consistency(iDomain, nSoilHorizons_temp, soilHorizonBoundaries_temp, &
          nLAIs_temp, LAIBoundaries_temp, nLandCoverPeriods_temp, landCoverPeriodBoundaries_temp, &
            landCoverPeriodSelect, do_read_dims)

    if (do_read_states) then
      allocate(mask_soil1 (level1(iDomain)%nrows, level1(iDomain)%ncols, nSoilHorizons))
      mask_soil1 = spread(mask1, 3, nSoilHorizons)

      if (nc%hasVariable('L1_Inter')) then
        !-------------------------------------------
        ! STATE VARIABLES (optionally)
        !-------------------------------------------

        ! Interception
        var = nc%getVariable("L1_Inter")
        call var%getData(dummyD2)
        L1_inter(s1 : e1) = pack(dummyD2, mask1)

        ! Snowpack
        var = nc%getVariable("L1_snowPack")
        call var%getData(dummyD2)
        L1_snowPack(s1 : e1) = pack(dummyD2, mask1)

        ! Retention storage of impervious areas
        var = nc%getVariable("L1_sealSTW")
        call var%getData(dummyD2)
        L1_sealSTW(s1 : e1) = pack(dummyD2, mask1)

        ! upper soil storage
        var = nc%getVariable("L1_unsatSTW")
        call var%getData(dummyD2)
        L1_unsatSTW(s1 : e1) = pack(dummyD2, mask1)

        ! groundwater storage
        var = nc%getVariable("L1_satSTW")
        call var%getData(dummyD2)
        L1_satSTW(s1 : e1) = pack(dummyD2, mask1)

        ! degree-day factor
        var = nc%getVariable("L1_degDay")
        call var%getData(dummyD2)
        L1_degDay(s1 : e1) = pack(dummyD2, mask1)

        ! Soil moisture of each horizon
        var = nc%getVariable("L1_soilMoist")
        call var%getData(dummyD3)
        L1_soilMoist(s1 : e1, :) = reshape(pack(dummyD3(:, :, :), mask_soil1), [e1-s1+1, nSoilHorizons])

        !-------------------------------------------
        ! FLUXES (optionally)
        !-------------------------------------------

        !  soil actual ET
        var = nc%getVariable("L1_aETSoil")
        call var%getData(dummyD3)
        L1_aETSoil(s1 : e1, : ) = reshape(pack(dummyD3(:, :, :), mask_soil1), [e1-s1+1, nSoilHorizons])

        ! canopy actual ET
        var = nc%getVariable("L1_aETCanopy")
        call var%getData(dummyD2)
        L1_aETCanopy(s1 : e1) = pack(dummyD2, mask1)

        ! sealed area actual ET
        var = nc%getVariable("L1_aETSealed")
        call var%getData(dummyD2)
        L1_aETSealed(s1 : e1) = pack(dummyD2, mask1)

        ! baseflow
        var = nc%getVariable("L1_baseflow")
        call var%getData(dummyD2)
        L1_baseflow(s1 : e1) = pack(dummyD2, mask1)

        ! soil in-exfiltration
        var = nc%getVariable("L1_infilSoil")
        call var%getData(dummyD3)
        L1_infilSoil(s1 : e1, :) = reshape(pack(dummyD3(:, :, :), mask_soil1), [e1-s1+1, nSoilHorizons])

        ! fast runoff
        var = nc%getVariable("L1_fastRunoff")
        call var%getData(dummyD2)
        L1_fastRunoff(s1 : e1) = pack(dummyD2, mask1)

        ! snow melt
        var = nc%getVariable("L1_melt")
        call var%getData(dummyD2)
        L1_melt(s1 : e1) = pack(dummyD2, mask1)

        ! percolation
        var = nc%getVariable("L1_percol")
        call var%getData(dummyD2)
        L1_percol(s1 : e1) = pack(dummyD2, mask1)

        ! effective precip. depth (snow melt + rain)
        var = nc%getVariable("L1_preEffect")
        call var%getData(dummyD2)
        L1_preEffect(s1 : e1) = pack(dummyD2, mask1)

        ! rain (liquid water)
        var = nc%getVariable("L1_rain")
        call var%getData(dummyD2)
        L1_rain(s1 : e1) = pack(dummyD2, mask1)

        ! runoff from impervious area
        var = nc%getVariable("L1_runoffSeal")
        call var%getData(dummyD2)
        L1_runoffSeal(s1 : e1) = pack(dummyD2, mask1)

        ! slow runoff
        var = nc%getVariable("L1_slowRunoff")
        call var%getData(dummyD2)
        L1_slowRunoff(s1 : e1) = pack(dummyD2, mask1)

        ! snow (solid water)
        var = nc%getVariable("L1_snow")
        call var%getData(dummyD2)
        L1_snow(s1 : e1) = pack(dummyD2, mask1)

        ! throughfall
        var = nc%getVariable("L1_Throughfall")
        call var%getData(dummyD2)
        L1_Throughfall(s1 : e1) = pack(dummyD2, mask1)

        ! total runoff
        var = nc%getVariable("L1_total_runoff")
        call var%getData(dummyD2)
        L1_total_runoff(s1 : e1) = pack(dummyD2, mask1)
      end if

      if (nc%hasVariable('L1_fSealed')) then
        ! Parameter fields have to be allocated in any case
        call init_eff_params(level1(iDomain)%nCells)

        !-------------------------------------------
        ! EFFECTIVE PARAMETERS
        !-------------------------------------------
        var = nc%getVariable("L1_fSealed")
        call var%getData(dummyD3)
        dummyD3 = dummyD3(:, :, landCoverPeriodSelect)
        do ii = 1, size(dummyD3, 3)
          L1_fSealed(s1 : e1, ii) = pack(dummyD3(:, :, ii), mask1)
        end do

        ! exponent for the upper reservoir
        var = nc%getVariable("L1_alpha")
        call var%getData(dummyD2)
        L1_alpha(s1 : e1, 1) = pack(dummyD2, mask1)
        ! call var%getData(dummyD3)
        ! dummyD3 = dummyD3(:, :, landCoverPeriodSelect)
        ! do ii = 1, size(dummyD3, 3)
        !   L1_alpha(s1 : e1, ii) = pack(dummyD3(:, :, ii), mask1)
        ! end do

        ! increase of the Degree-day factor per mm of increase in precipitation
        var = nc%getVariable("L1_degDayInc")
        call var%getData(dummyD3)
        dummyD3 = dummyD3(:, :, landCoverPeriodSelect)
        do ii = 1, size(dummyD3, 3)
          L1_degDayInc(s1 : e1, ii) = pack(dummyD3(:, :, ii), mask1)
        end do

        ! maximum degree-day factor
        var = nc%getVariable("L1_degDayMax")
        call var%getData(dummyD3)
        dummyD3 = dummyD3(:, :, landCoverPeriodSelect)
        do ii = 1, size(dummyD3, 3)
          L1_degDayMax(s1 : e1, ii) = pack(dummyD3(:, :, ii), mask1)
        end do

        ! degree-day factor with no precipitation
        var = nc%getVariable("L1_degDayNoPre")
        call var%getData(dummyD3)
        dummyD3 = dummyD3(:, :, landCoverPeriodSelect)
        do ii = 1, size(dummyD3, 3)
          L1_degDayNoPre(s1 : e1, ii) = pack(dummyD3(:, :, ii), mask1)
        end do

        ! Karstic percolation loss
        var = nc%getVariable("L1_karstLoss")
        call var%getData(dummyD2)
        L1_karstLoss(s1 : e1) = pack(dummyD2, mask1)

        ! Fraction of roots in soil horizons
        var = nc%getVariable("L1_fRoots")
        call var%getData(dummyD4)
        dummyD4 = dummyD4(:, :, :, landCoverPeriodSelect)
        do jj = 1, size(dummyD4, 4)
          L1_fRoots(s1 : e1, :, jj) = reshape(pack(dummyD4(:, :, :, jj), mask_soil1), [e1-s1+1, nSoilHorizons])
        end do

        ! Maximum interception
        var = nc%getVariable("L1_maxInter")
        call var%getData(dummyD3)
        do ii = 1, nLAIs
          L1_maxInter(s1 : e1, ii) = pack(dummyD3(:, :, ii), mask1)
        end do

        ! fast interflow recession coefficient
        var = nc%getVariable("L1_kFastFlow")
        call var%getData(dummyD3)
        dummyD3 = dummyD3(:, :, landCoverPeriodSelect)
        do ii = 1, size(dummyD3, 3)
          L1_kFastFlow(s1 : e1, ii) = pack(dummyD3(:, :, ii), mask1)
        end do

        ! slow interflow recession coefficient
        var = nc%getVariable("L1_kSlowFlow")
        call var%getData(dummyD2)
        L1_kSlowFlow(s1 : e1, 1) = pack(dummyD2, mask1)
        !call var%getData(dummyD3)
        !dummyD3 = dummyD3(:, :, landCoverPeriodSelect)
        ! do ii = 1, size(dummyD3, 3)
        !   L1_kSlowFlow(s1 : e1, ii) = pack(dummyD3(:, :, ii), mask1)
        ! end do

        ! baseflow recession coefficient
        var = nc%getVariable("L1_kBaseFlow")
        call var%getData(dummyD2)
        L1_kBaseFlow(s1 : e1, 1) = pack(dummyD2, mask1)
        ! call var%getData(dummyD3)
        ! dummyD3 = dummyD3(:, :, landCoverPeriodSelect)
        ! do ii = 1, size(dummyD3, 3)
        !   L1_kBaseFlow(s1 : e1, ii) = pack(dummyD3(:, :, ii), mask1)
        ! end do

        ! percolation coefficient
        var = nc%getVariable("L1_kPerco")
        call var%getData(dummyD2)
        L1_kPerco(s1 : e1, 1) = pack(dummyD2, mask1)
        ! call var%getData(dummyD3)
        ! dummyD3 = dummyD3(:, :, landCoverPeriodSelect)
        ! do ii = 1, size(dummyD3, 3)
        !   L1_kPerco(s1 : e1, ii) = pack(dummyD3(:, :, ii), mask1)
        ! end do

        ! Soil moisture below which actual ET is reduced linearly till PWP
        ! for processCase(3) = 1
        var = nc%getVariable("L1_soilMoistFC")
        call var%getData(dummyD4)
        dummyD4 = dummyD4(:, :, :, landCoverPeriodSelect)
        do jj = 1, size(dummyD4, 4)
          L1_soilMoistFC(s1 : e1, :, jj) = reshape(pack(dummyD4(:, :, :, jj), mask_soil1), [e1-s1+1, nSoilHorizons])
        end do

        ! Saturation soil moisture for each horizon [mm]
        var = nc%getVariable("L1_soilMoistSat")
        call var%getData(dummyD4)
        dummyD4 = dummyD4(:, :, :, landCoverPeriodSelect)
        do jj = 1, size(dummyD4, 4)
          L1_soilMoistSat(s1 : e1, :, jj) = reshape(pack(dummyD4(:, :, :, jj), mask_soil1), [e1-s1+1, nSoilHorizons])
        end do

        ! Exponential parameter to how non-linear is the soil water retention
        var = nc%getVariable("L1_soilMoistExp")
        call var%getData(dummyD4)
        dummyD4 = dummyD4(:, :, :, landCoverPeriodSelect)
        do jj = 1, size(dummyD4, 4)
          L1_soilMoistExp(s1 : e1, :, jj) = reshape(pack(dummyD4(:, :, :, jj), mask_soil1), [e1-s1+1, nSoilHorizons])
        end do

        if (any(processMatrix(3, 1) == [2, 3])) then
          ! jarvis critical value for normalized soil water content
          var = nc%getVariable("L1_jarvis_thresh_c1")
          call var%getData(dummyD2)
          L1_jarvis_thresh_c1(s1 : e1) = pack(dummyD2, mask1)
        end if

        ! Threshold temperature for snow/rain
        var = nc%getVariable("L1_tempThresh")
        call var%getData(dummyD3)
        dummyD3 = dummyD3(:, :, landCoverPeriodSelect)
        do ii = 1, size(dummyD3, 3)
          L1_tempThresh(s1 : e1, ii) = pack(dummyD3(:, :, ii), mask1)
        end do

        ! Threshold water depth controlling fast interflow
        var = nc%getVariable("L1_unsatThresh")
        call var%getData(dummyD2)
        L1_unsatThresh(s1 : e1, 1) = pack(dummyD2, mask1)
        ! call var%getData(dummyD3)
        ! dummyD3 = dummyD3(:, :, landCoverPeriodSelect)
        ! do ii = 1, size(dummyD3, 3)
        !   L1_unsatThresh(s1 : e1, ii) = pack(dummyD3(:, :, ii), mask1)
        ! end do

        ! Threshold water depth for surface runoff in sealed surfaces
        var = nc%getVariable("L1_sealedThresh")
        call var%getData(dummyD2)
        L1_sealedThresh(s1 : e1) = pack(dummyD2, mask1)

        ! Permanent wilting point
        var = nc%getVariable("L1_wiltingPoint")
        call var%getData(dummyD4)
        dummyD4 = dummyD4(:, :, :, landCoverPeriodSelect)
        do jj = 1, size(dummyD4, 4)
          L1_wiltingPoint(s1 : e1, :, jj) = reshape(pack(dummyD4(:, :, :, jj), mask_soil1), [e1-s1+1, nSoilHorizons])
        end do

        ! different parameters dependent on PET formulation
        select case (processMatrix(5, 1))
        case(-1) ! PET is input

          ! PET correction factor due to LAI
          var = nc%getVariable("L1_petLAIcorFactor")
          call var%getData(dummyD4)
          dummyD4 = dummyD4(:, :, :, landCoverPeriodSelect)
          do jj = 1, size(dummyD4, 4)
            do ii = 1, nLAIs
              L1_petLAIcorFactor(s1 : e1, ii, jj) = pack(dummyD4(:, :, ii, jj), mask1)
            end do
          end do

        case(0) ! PET is input

          ! PET correction factor due to terrain aspect
          var = nc%getVariable("L1_fAsp")
          call var%getData(dummyD2)
          L1_fAsp(s1 : e1) = pack(dummyD2, mask1)

        case(1) ! Hargreaves-Samani

          ! PET correction factor due to terrain aspect
          var = nc%getVariable("L1_fAsp")
          call var%getData(dummyD2)
          L1_fAsp(s1 : e1) = pack(dummyD2, mask1)

          ! Hargreaves Samani coeffiecient
          var = nc%getVariable("L1_HarSamCoeff")
          call var%getData(dummyD2)
          L1_HarSamCoeff(s1 : e1) = pack(dummyD2, mask1)

        case(2) ! Priestely-Taylor

          ! Priestley Taylor coeffiecient (alpha)
          var = nc%getVariable("L1_PrieTayAlpha")
          call var%getData(dummyD3)
          do ii = 1, nLAIs
            L1_PrieTayAlpha(s1 : e1, ii) = pack(dummyD3(:, :, ii), mask1)
          end do

        case(3) ! Penman-Monteith

          ! aerodynamical resitance
          var = nc%getVariable("L1_aeroResist")
          call var%getData(dummyD4)
          dummyD4 = dummyD4(:, :, :, landCoverPeriodSelect)
          do jj = 1, size(dummyD4, 4)
            do ii = 1, nLAIs
              L1_aeroResist(s1 : e1, ii, jj) = pack(dummyD4(:, :, ii, jj), mask1)
            end do
          end do

          ! bulk surface resitance
          var = nc%getVariable("L1_surfResist")
          call var%getData(dummyD3)
          do ii = 1, nLAIs
            L1_surfResist(s1 : e1, ii) = pack(dummyD3(:, :, ii), mask1)
          end do

        end select
        are_parameter_initialized = .true.
      end if

      deallocate(mask_soil1, mask1)
    end if

    call nc%close()

  end subroutine read_restart_states

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
    deallocate(dummy_arr)

  end subroutine

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
    do jj = 1, size(data, 3)
      do ii = 1, size(data, 2)
        dummy_arr(:, :, ii, jj) = unpack(data(:, ii, jj), mask, fill_value)
      end do
    end do

    ! set the unpacked data
    call var%setData(dummy_arr)

    ! optionally set attributes
    if (present(var_long_name)) then
      call var%setAttribute("long_name", trim(var_long_name))
    end if
    deallocate(dummy_arr)

  end subroutine

  subroutine unpack_field_and_write_soil(nc, var_name, var_dims, fill_value, data, mask, var_long_name)

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
    logical, dimension(:, :, :), allocatable :: mask_soil1

    integer(i4), dimension(4) :: dim_length

    integer(i4) :: jj


    ! set variable
    var = nc%setVariable(var_name, "f64", var_dims)
    call var%setFillValue(fill_value)
    dim_length = var%getShape()

    allocate(mask_soil1 (dim_length(1), dim_length(2), dim_length(3)))
    mask_soil1 = spread(mask, 3, dim_length(3))

    allocate(dummy_arr(dim_length(1), dim_length(2), dim_length(3), dim_length(4)))
    do jj = 1, size(data, 3)
      dummy_arr(:, :, :, jj) = unpack(pack(data(:, :, jj), .true.), mask_soil1, fill_value)
    end do

    ! set the unpacked data
    call var%setData(dummy_arr)

    ! optionally set attributes
    if (present(var_long_name)) then
      call var%setAttribute("long_name", trim(var_long_name))
    end if
    deallocate(dummy_arr, mask_soil1)

  end subroutine

  !    NAME
  !        write_eff_params

  !    PURPOSE
  !>       \brief TODO: add description

  !>       \details TODO: add description

  !    INTENT(IN)
  !>       \param[in] "logical, dimension(:, :) :: mask1"                        mask at level 1
  !>       \param[in] "integer(i4) :: s1"                                        start index at level 1
  !>       \param[in] "integer(i4) :: e1"                                        end index at level 1
  !>       \param[in] "type(NcDimension) :: rows1, cols1, soil1, lcscenes, lais"
  !>       \param[in] "type(NcDimension) :: rows1, cols1, soil1, lcscenes, lais"
  !>       \param[in] "type(NcDimension) :: rows1, cols1, soil1, lcscenes, lais"
  !>       \param[in] "type(NcDimension) :: rows1, cols1, soil1, lcscenes, lais"
  !>       \param[in] "type(NcDimension) :: rows1, cols1, soil1, lcscenes, lais"

  !    INTENT(INOUT)
  !>       \param[inout] "type(NcDataset) :: nc"

  !    HISTORY
  !>       \authors Robert Schweppe

  !>       \date Jun 2018

  ! Modifications:

  subroutine write_eff_params(mask1, s1, e1, rows1, cols1, soil1, lcscenes, lais, nc, iDomainNLandCoverPeriods)

    use mo_constants, only : nodata_dp
    use mo_common_variables, only : processMatrix
    use mo_kind, only : i4
    use mo_global_variables, only : L1_HarSamCoeff, L1_PrieTayAlpha, L1_aeroResist, &
                                        L1_alpha, L1_degDayInc, L1_degDayMax, L1_degDayNoPre, L1_fAsp, L1_latitude, &
                                        L1_fRoots, L1_fSealed, L1_jarvis_thresh_c1, L1_kBaseFlow, L1_kPerco, &
                                        L1_kSlowFlow, L1_karstLoss, L1_kFastFlow, L1_maxInter, L1_petLAIcorFactor, &
                                        L1_sealedThresh, L1_soilMoistExp, L1_soilMoistFC, L1_soilMoistSat, L1_surfResist, &
                                        L1_tempThresh, L1_unsatThresh, L1_wiltingPoint, nLAIs, nSoilHorizons
    use mo_netcdf, only : NcDataset, NcDimension

    implicit none

    ! mask at level 1
    logical, dimension(:, :), intent(in) :: mask1

    ! start and end index at level 1, number of land cover periods
    integer(i4), intent(in) :: s1, e1, iDomainNLandCoverPeriods

    type(NcDimension), intent(in) :: rows1, cols1, soil1, lcscenes, lais

    type(NcDataset), intent(inout) :: nc

    !-------------------------------------------
    ! EFFECTIVE PARAMETERS
    !-------------------------------------------
    ! TODO: MPR add dims
    call unpack_field_and_write(nc, "L1_fSealed", &
            [rows1, cols1, lcscenes], nodata_dp, L1_fSealed(s1 : e1, 1:iDomainNLandCoverPeriods), mask1, &
            "fraction of Sealed area at level 1")

    call unpack_field_and_write(nc, "L1_alpha", &
            ![rows1, cols1, lcscenes], nodata_dp, L1_alpha(s1 : e1, 1:iDomainNLandCoverPeriods), mask1, &
            [rows1, cols1], nodata_dp, L1_alpha(s1 : e1, 1), mask1, &
            "exponent for the upper reservoir at level 1")

    call unpack_field_and_write(nc, "L1_degDayInc", &
            [rows1, cols1, lcscenes], nodata_dp, L1_degDayInc(s1 : e1, 1:iDomainNLandCoverPeriods), mask1, &
            "increase of the Degree-day factor per mm of increase in precipitation at level 1")

    call unpack_field_and_write(nc, "L1_degDayMax", &
            [rows1, cols1, lcscenes], nodata_dp, L1_degDayMax(s1 : e1, 1:iDomainNLandCoverPeriods), mask1, &
            "maximum degree-day factor at level 1")

    call unpack_field_and_write(nc, "L1_degDayNoPre", &
            [rows1, cols1, lcscenes], nodata_dp, L1_degDayNoPre(s1 : e1, 1:iDomainNLandCoverPeriods), mask1, &
            "degree-day factor with no precipitation at level 1")

    call unpack_field_and_write(nc, "L1_karstLoss", &
            [rows1, cols1], nodata_dp, L1_karstLoss(s1 : e1), mask1, &
            "Karstic percolation loss at level 1")

    call unpack_field_and_write_soil(nc, "L1_fRoots", &
            [rows1, cols1, soil1, lcscenes], nodata_dp, L1_fRoots(s1 : e1, 1:nSoilHorizons, 1:iDomainNLandCoverPeriods), mask1, &
            "Fraction of roots in soil horizons at level 1")

    call unpack_field_and_write(nc, "L1_maxInter", &
            [rows1, cols1, lais], nodata_dp, L1_maxInter(s1 : e1, 1:nLAIs), mask1, &
            "Maximum interception at level 1")

    call unpack_field_and_write(nc, "L1_kFastFlow", &
            [rows1, cols1, lcscenes], nodata_dp, L1_kFastFlow(s1 : e1, 1:iDomainNLandCoverPeriods), mask1, &
            "fast interflow recession coefficient at level 1")

    call unpack_field_and_write(nc, "L1_kSlowFlow", &
            ! [rows1, cols1, lcscenes], nodata_dp, L1_kSlowFlow(s1 : e1, 1:iDomainNLandCoverPeriods), mask1, &
            [rows1, cols1], nodata_dp, L1_kSlowFlow(s1 : e1, 1), mask1, &
            "slow interflow recession coefficient at level 1")

    call unpack_field_and_write(nc, "L1_kBaseFlow", &
            ! [rows1, cols1, lcscenes], nodata_dp, L1_kBaseFlow(s1 : e1, 1:iDomainNLandCoverPeriods), mask1, &
            [rows1, cols1], nodata_dp, L1_kBaseFlow(s1 : e1, 1), mask1, &
            "baseflow recession coefficient at level 1")

    call unpack_field_and_write(nc, "L1_kPerco", &
            ! [rows1, cols1, lcscenes], nodata_dp, L1_kPerco(s1 : e1, 1:iDomainNLandCoverPeriods), mask1, &
            [rows1, cols1], nodata_dp, L1_kPerco(s1 : e1, 1), mask1, &
            "percolation coefficient at level 1")

    call unpack_field_and_write_soil(nc, "L1_soilMoistFC", &
            [rows1, cols1, soil1, lcscenes], nodata_dp, &
            L1_soilMoistFC(s1 : e1, 1:nSoilHorizons, 1:iDomainNLandCoverPeriods), mask1, &
            "SM below which actual ET is reduced linearly till PWP at level 1 for processCase(3)=1")

    call unpack_field_and_write_soil(nc, "L1_soilMoistSat", &
            [rows1, cols1, soil1, lcscenes], nodata_dp, &
            L1_soilMoistSat(s1 : e1, 1:nSoilHorizons, 1:iDomainNLandCoverPeriods), mask1, &
            "Saturation soil moisture for each horizon [mm] at level 1")

    call unpack_field_and_write_soil(nc, "L1_soilMoistExp", &
            [rows1, cols1, soil1, lcscenes], nodata_dp, &
            L1_soilMoistExp(s1 : e1, 1:nSoilHorizons, 1:iDomainNLandCoverPeriods), mask1, &
            "Exponential parameter to how non-linear is the soil water retention at level 1")

    if (processMatrix(3, 1) == 2 .or. processMatrix(3, 1) == 3) then
      call unpack_field_and_write(nc, "L1_jarvis_thresh_c1", &
              [rows1, cols1], nodata_dp, L1_jarvis_thresh_c1(s1 : e1), mask1, &
              "jarvis critical value for normalized soil water content")
    end if

    if (processMatrix(5, 1) == -1) then
      call unpack_field_and_write(nc, "L1_petLAIcorFactor", &
              [rows1, cols1, lais, lcscenes], nodata_dp, &
              L1_petLAIcorFactor(s1 : e1, 1:nLAIs, 1:iDomainNLandCoverPeriods), mask1, &
              "PET correction factor based on LAI")
    end if

    call unpack_field_and_write(nc, "L1_tempThresh", &
            [rows1, cols1, lcscenes], nodata_dp, L1_tempThresh(s1 : e1, 1:iDomainNLandCoverPeriods), mask1, &
            "Threshold temperature for snow/rain at level 1")

    call unpack_field_and_write(nc, "L1_unsatThresh", &
            ! [rows1, cols1, lcscenes], nodata_dp, L1_unsatThresh(s1 : e1, 1:iDomainNLandCoverPeriods), mask1, &
            [rows1, cols1], nodata_dp, L1_unsatThresh(s1 : e1, 1), mask1, &
            "Threshold water depth controlling fast interflow at level 1")

    call unpack_field_and_write(nc, "L1_sealedThresh", &
            [rows1, cols1], nodata_dp, L1_sealedThresh(s1 : e1), mask1, &
            "Threshold water depth for surface runoff in sealed surfaces at level 1")

    call unpack_field_and_write_soil(nc, "L1_wiltingPoint", &
            [rows1, cols1, soil1, lcscenes], nodata_dp, &
            L1_wiltingPoint(s1 : e1, 1:nSoilHorizons, 1:iDomainNLandCoverPeriods), mask1, &
            "Permanent wilting point at level 1")

    !call unpack_field_and_write(nc, "L1_latitude", &
    !        [rows1, cols1], nodata_dp, L1_latitude(s1 : e1), mask1, &
    !        "Latitude at level 1")

    select case (processMatrix(5, 1))
    case(-1 : 0) ! PET is input
      call unpack_field_and_write(nc, "L1_fAsp", &
              [rows1, cols1], nodata_dp, L1_fAsp(s1 : e1), mask1, &
              "PET correction factor due to terrain aspect at level 1")

    case(1) ! Hargreaves-Samani
      call unpack_field_and_write(nc, "L1_fAsp", &
              [rows1, cols1], nodata_dp, L1_fAsp(s1 : e1), mask1, &
              "PET correction factor due to terrain aspect at level 1")

      call unpack_field_and_write(nc, "L1_HarSamCoeff", &
              [rows1, cols1], nodata_dp, L1_HarSamCoeff(s1 : e1), mask1, &
              "Hargreaves-Samani coefficient")

    case(2) ! Priestley-Taylor
      call unpack_field_and_write(nc, "L1_PrieTayAlpha", &
              [rows1, cols1, lais], nodata_dp, L1_PrieTayAlpha(s1 : e1, 1:nLAIs), mask1, &
              "Priestley Taylor coeffiecient (alpha)")

    case(3) ! Penman-Monteith
      call unpack_field_and_write(nc, "L1_aeroResist", &
              [rows1, cols1, lais, lcscenes], nodata_dp, &
              L1_aeroResist(s1 : e1, 1:nLAIs, 1:iDomainNLandCoverPeriods), mask1, &
              "aerodynamical resitance")

      call unpack_field_and_write(nc, "L1_surfResist", &
              [rows1, cols1, lais], nodata_dp, L1_surfResist(s1 : e1, 1:nLAIs), mask1, &
              "bulk surface resitance")

    end select

  end subroutine write_eff_params

  subroutine init_eff_params(ncells1)

    use mo_append, only : append
    use mo_kind, only: i4, dp
    use mo_common_constants, only : P1_InitStateFluxes
    use mo_common_variables, only: nLandCoverPeriods
    use mo_global_variables, only : L1_HarSamCoeff, L1_PrieTayAlpha, L1_aeroResist, L1_alpha, &
                                    L1_degDayInc, L1_degDayMax, L1_degDayNoPre, L1_fAsp, L1_fRoots, L1_fSealed, &
                                    L1_jarvis_thresh_c1, L1_kBaseFlow, L1_kPerco, L1_kSlowFlow, L1_karstLoss, &
                                    L1_kFastFlow, L1_maxInter, L1_petLAIcorFactor, L1_sealedThresh, L1_soilMoistExp, &
                                    L1_soilMoistFC, L1_soilMoistSat, L1_surfResist, L1_tempThresh, L1_unsatThresh, &
                                    L1_wiltingPoint, L1_latitude, nLAIs, nSoilHorizons

    implicit none

    integer(i4), intent(in) :: ncells1

    real(dp), dimension(:, :, :), allocatable :: dummy_3D
    real(dp), dimension(:, :), allocatable :: dummy_2D
    real(dp), dimension(:), allocatable :: dummy_1D

    !-------------------------------------------
    ! EFFECTIVE PARAMETERS
    !-------------------------------------------
    ! for appending and intialization
    allocate(dummy_3D(nCells1, nSoilHorizons, nLandCoverPeriods))
    dummy_3D = P1_InitStateFluxes
    ! Fraction of roots in soil horizons
    call append(L1_fRoots, dummy_3D)
    ! Soil moisture below which actual ET is reduced linearly till PWP
    call append(L1_soilMoistFC, dummy_3D)
    ! Saturation soil moisture for each horizon [mm]
    call append(L1_soilMoistSat, dummy_3D)
    ! Exponential parameter to how non-linear is the soil water retention
    call append(L1_soilMoistExp, dummy_3D)
    ! Permanent wilting point
    call append(L1_wiltingPoint, dummy_3D)
    deallocate(dummy_3D)

    allocate(dummy_3D(nCells1, nLAIs, nLandCoverPeriods))
    dummy_3D = P1_InitStateFluxes
    ! PET correction factor due to LAI
    call append(L1_petLAIcorFactor, dummy_3D)
    ! PET aerodynamical resistance
    call append(L1_aeroResist, dummy_3D)
    deallocate(dummy_3D)

    allocate(dummy_2D(nCells1, nLandCoverPeriods))
    dummy_2D = P1_InitStateFluxes
    call append(L1_fSealed, dummy_2D)
    ! increase of the Degree-day factor per mm of increase in precipitation
    call append(L1_degDayInc, dummy_2D)
    ! maximum degree-day factor
    call append(L1_degDayMax, dummy_2D)
    ! degree-day factor with no precipitation
    call append(L1_degDayNoPre, dummy_2D)
    ! fast interflow recession coefficient
    call append(L1_kFastFlow, dummy_2D)
    ! Threshold temperature for snow/rain
    call append(L1_tempThresh, dummy_2D)
    ! Threshold water depth controlling fast interflow
    call append(L1_unsatThresh, dummy_2D)
    ! slow interflow recession coefficient
    call append(L1_kSlowFlow, dummy_2D)
    ! percolation coefficient
    call append(L1_kPerco, dummy_2D)
    ! exponent for the upper reservoir
    call append(L1_alpha, dummy_2D)
    ! baseflow recession coefficient
    call append(L1_kBaseFlow, dummy_2D)
    deallocate(dummy_2D)

    allocate(dummy_2D(nCells1, nLAIs))
    dummy_2D = P1_InitStateFluxes
    ! PET Prietley Taylor coefficient
    call append(L1_PrieTayAlpha, dummy_2D)
    ! PET bulk surface resistance
    call append(L1_surfResist, dummy_2D)
    ! Maximum interception
    call append(L1_maxInter, dummy_2D)
    deallocate(dummy_2D)

    allocate(dummy_1D(nCells1))
    dummy_1D = P1_InitStateFluxes
    ! Karstic percolation loss
    call append(L1_karstLoss, dummy_1D)
    ! PET correction factor due to terrain aspect
    call append(L1_fAsp, dummy_1D)
    ! latitude
    call append(L1_latitude, dummy_1D)
    ! PET Hargreaves Samani coefficient
    call append(L1_HarSamCoeff, dummy_1D)
    ! jarvis critical value for normalized soil water content
    call append(L1_jarvis_thresh_c1, dummy_1D)
    ! Threshhold water depth for surface runoff in sealed surfaces
    call append(L1_sealedThresh, dummy_1D)
    deallocate(dummy_1D)


  end subroutine init_eff_params

END MODULE mo_restart
