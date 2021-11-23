!>       \file mo_mrm_read_data.f90

!>       \brief This module contains all routines to read mRM data from file.
!>       \details

!>       \details TODO: add description

!>       \authors Stephan Thober

!>       \date Aug 2015

! Modifications:

module mo_mrm_read_data
  use mo_kind, only : i4, dp
  use mo_netcdf, only : NcDataset, NcVariable, NcDimension
  use mo_read_nc, only: check_sort_order

  implicit none
  public :: mrm_read_L0_data
  public :: mrm_read_discharge
  public :: mrm_read_total_runoff
  public :: mrm_read_bankfull_runoff
  private
contains
  ! ------------------------------------------------------------------

  !    NAME
  !        mrm_read_L0_data

  !    PURPOSE
  !>       \brief read L0 data from file

  !>       \details With the exception of L0_mask, L0_elev, and L0_LCover, all
  !>       L0 variables are read from file. The former three are only read if they
  !>       are not provided as variables.

  !    INTENT(IN)
  !>       \param[in] "logical :: do_reinit"
  !>       \param[in] "logical :: do_readlatlon"
  !>       \param[in] "logical :: do_readlcover"

  !    HISTORY
  !>       \authors Juliane Mai, Matthias Zink, and Stephan Thober

  !>       \date Aug 2015

  ! Modifications:
  ! Stephan Thober Sep 2015 - added L0_mask, L0_elev, and L0_LCover
  ! Robert Schweppe Jun 2018 - refactoring and reformatting
  ! Stephan Thober Jun 2018 - including varying celerity functionality

  subroutine mrm_read_L0_data(do_readlcover)

    use mo_append, only : append, add_nodata_slice
    use mo_common_constants, only : nodata_i4
    use mo_common_read_data, only : read_dem, read_lcover
    use mo_common_variables, only : L0_elev, L0_LCover, level0, domainMeta, processMatrix
    use mo_grid, only: Grid
    use mo_message, only : error_message, message
    use mo_mrm_file, only : file_facc, file_fdir, file_gaugeloc
    use mo_mrm_global_variables, only : L0_InflowGaugeLoc, L0_fAcc, L0_fDir, L0_gaugeLoc, domain_mrm, dirGauges
    use mo_string_utils, only : num2str
    use mo_common_datetime_type, only: get_land_cover_period_indices, simPer, nLandCoverPeriods

    implicit none

    logical, intent(in) :: do_readlcover
    integer(i4) ::domainID, iDomain
    integer(i4) :: iVar
    integer(i4) :: iGauge
    character(256) :: fname, varName
    integer(i4) :: nunit
    integer(i4) :: nCells
    integer(i4) :: nLandCoverPeriods_temp
    real(dp), dimension(:), allocatable :: landCoverPeriodBoundaries_temp
    integer(i4), dimension(:), allocatable :: landCoverSelect


    integer(i4), dimension(:, :), allocatable :: data_i4_2d
    real(dp), dimension(:, :), allocatable :: data_dp_2d
    logical, dimension(:, :), allocatable :: mask_2d

    type(Grid), pointer :: level0_iDomain => null()
    type(NcDataset)                        :: nc           ! netcdf file
    type(NcVariable)                       :: ncVar        ! variables for data form netcdf
    integer(i4)                            :: nodata_value ! data nodata value


    do iDomain = 1, domainMeta%nDomains
      ! ************************************************
      ! READ SPATIAL DATA FOR EACH DOMAIN
      ! ************************************************
      domainID = domainMeta%indices(iDomain)

      level0_iDomain => level0(domainMeta%L0DataFrom(iDomain))

      ! check whether L0 data is shared
      if (domainMeta%L0DataFrom(iDomain) < iDomain) then
        call message('      Using data of domain ', &
                trim(adjustl(num2str(domainMeta%indices(domainMeta%L0DataFrom(iDomain))))), ' for domain: ',&
                trim(adjustl(num2str(domainID))), '...')
        cycle
      end if

      call message('      Reading data for domain: ', trim(adjustl(num2str(domainID))), ' ...')
      call read_dem(iDomain, level0_iDomain, data_dp_2d)
      ! put data in variable
      call append(L0_elev, pack(data_dp_2d, level0_iDomain%mask))

      if (do_readlcover) then
        ! if case 8==1, then we need lcover for parameter estimation
        if (processMatrix(8, 1) .eq. 1) then
          ! read the land cover file, all periods
          call read_lcover(iDomain, data_i4_2d, nLandCoverPeriods_temp, landCoverPeriodBoundaries_temp)
          ! compare the simulation period and the land cover periods in the file,
          ! get a boolean vector with periods to select
          call get_land_cover_period_indices(simPer(iDomain), landCoverPeriodBoundaries_temp, &
                  selectIndices=landCoverSelect)
          ! select the needed periods and fill remaining slices so appending works
          ! background: all domains can have different number of land cover periods but data are in one big
          ! pre-allocated array for all domains
          data_i4_2d = data_i4_2d(:, landCoverSelect)
          call add_nodata_slice(data_i4_2d, nLandCoverPeriods - size(landCoverSelect), nodata_i4)
        else if ((processMatrix(8, 1) .eq. 2) .or. (processMatrix(8, 1) .eq. 3)) then
          allocate(data_i4_2d(level0_iDomain%nCells, nLandCoverPeriods))
          data_i4_2d = nodata_i4
        end if
        call append(L0_LCover, data_i4_2d)
        ! free memory
        deallocate(data_i4_2d)
      end if

      ! read fAcc, fDir, gaugeLoc
      do iVar = 1, 3
        select case (iVar)
        case(1) ! flow accumulation
          fName = trim(adjustl(dirGauges(iDomain))) // trim(adjustl(file_facc))
          varName = 'facc'
        case(2) ! flow direction
          fName = trim(adjustl(dirGauges(iDomain))) // trim(adjustl(file_fdir))
          varName = 'fdir'
        case(3) ! location of gauging stations
          fName = trim(adjustl(dirGauges(iDomain))) // trim(adjustl(file_gaugeloc))
          varName = 'idgauges'
        end select
        if (iVar == 3 .and. domain_mrm(iDomain)%nGauges < 1_i4) then
          ! set to nodata, but not omit because data association between arrays and domains might break
          allocate(data_i4_2d(level0_iDomain%nrows, level0_iDomain%ncols))
          data_i4_2d = nodata_i4
        else
          ! read the Dataset
          nc = NcDataset(fname, "r")
          ! get the variable
          ncVar = nc%getVariable(trim(varName))

          call ncVar%getData(data_i4_2d, mask=mask_2d)
          if ( size(data_i4_2d, 1) /= level0_iDomain%nrows .or. size(data_i4_2d, 2) /= level0_iDomain%ncols) then
            call error_message('***ERROR: read_forcing_nc: mHM generated x and y: ', &
                    num2str(level0_iDomain%nrows), num2str(level0_iDomain%ncols) , &
                    'are not matching NetCDF dimensions: ', num2str(size(data_i4_2d, 1)), num2str(size(data_i4_2d, 2)))
          end if

          ! flip the data if any dimension is not sorted correctly
          call check_sort_order(data_i4_2d, ncVar)

          ! put global nodata value into array (probably not all grid cells have values)
          data_i4_2d = merge(data_i4_2d, nodata_i4, mask=mask_2d)

          call nc%close()
        end if
        select case (iVar)
        case(1) ! flow accumulation
          call append(L0_fAcc, pack(data_i4_2d, level0_iDomain%mask))
        case(2) ! flow direction
          ! TODO: if a flipping occurs, we have to rotate the fdir variable
          call append(L0_fDir, pack(data_i4_2d, level0_iDomain%mask))
        case(3) ! location of evaluation and inflow gauging stations
          ! evaluation gauges
          ! Input data check
          do iGauge = 1, domain_mrm(iDomain)%nGauges
            ! If gaugeId is found in gauging location file?
            if (.not. any(data_i4_2d .EQ. domain_mrm(iDomain)%gaugeIdList(iGauge))) then
              call error_message(&
                      '***ERROR: Gauge ID "', trim(adjustl(num2str(domain_mrm(iDomain)%gaugeIdList(iGauge)))), &
                      '" not found in gauge location input file: ', &
                      trim(adjustl(dirGauges(iDomain))) // trim(adjustl(file_gaugeloc)))
            end if
          end do

          call append(L0_gaugeLoc, pack(data_i4_2d, level0_iDomain%mask))

          ! inflow gauges
          ! if no inflow gauge for this subdomain exists still matirx with dim of subdomain has to be paded
          if (domain_mrm(iDomain)%nInflowGauges .GT. 0_i4) then
            ! Input data check
            do iGauge = 1, domain_mrm(iDomain)%nInflowGauges
              ! If InflowGaugeId is found in gauging location file?
              if (.not. any(data_i4_2d .EQ. domain_mrm(iDomain)%InflowGaugeIdList(iGauge))) then
                call error_message('***ERROR: Inflow Gauge ID "', &
                        trim(adjustl(num2str(domain_mrm(iDomain)%InflowGaugeIdList(iGauge)))), &
                        '" not found in gauge location input file: ', &
                        trim(adjustl(dirGauges(iDomain))) // trim(adjustl(file_gaugeloc)))
              end if
            end do
          end if

          call append(L0_InflowGaugeLoc, pack(data_i4_2d, level0_iDomain%mask))

        end select
        !
        ! deallocate arrays
        deallocate(data_i4_2d)
        !
      end do
    end do

    ! some sanity checks
    if ( any (L0_fAcc == nodata_i4)) then
      call error_message(' Error: flow accumulation field has missing values within the valid masked area')
    end if
    if ( any (L0_fDir == nodata_i4)) then
      call error_message(' Error: flow direction field has missing values within the valid masked area')
    end if


  end subroutine mrm_read_L0_data
  ! ---------------------------------------------------------------------------

  !    NAME
  !        mrm_read_discharge

  !    PURPOSE
  !>       \brief Read discharge timeseries from file

  !>       \details Read Observed discharge at the outlet of a catchment
  !>       and at the inflow of a catchment. Allocate global runoff
  !>       variable that contains the simulated runoff after the simulation.

  !    HISTORY
  !>       \authors Matthias Zink & Stephan Thober

  !>       \date Aug 2015

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine mrm_read_discharge

    use mo_append, only : paste
    use mo_common_constants, only : nodata_dp
    use mo_common_variables, only : domainMeta, evalPer, opti_function, optimize
    use mo_common_datetime_type, only: nTstepDay, simPer
    use mo_message, only : error_message
    use mo_mrm_file, only : udischarge
    use mo_mrm_global_variables, only : InflowGauge, gauge, mRM_runoff, nGaugesLocal, &
                                        nInflowGaugesTotal, nMeasPerDay, &
                                        riv_temp_pcs
    use mo_read_timeseries, only : read_timeseries
    use mo_string_utils, only : num2str

    implicit none

    integer(i4) :: iGauge

    integer(i4) :: iDomain

    integer(i4) :: maxTimeSteps

    ! file name of file to read
    character(256) :: fName

    integer(i4), dimension(3) :: start_tmp, end_tmp

    real(dp), dimension(:), allocatable :: data_dp_1d

    logical, dimension(:), allocatable :: mask_1d


    !----------------------------------------------------------
    ! INITIALIZE RUNOFF
    !----------------------------------------------------------
    maxTimeSteps = maxval(simPer(1 : domainMeta%nDomains)%julEnd - simPer(1 : domainMeta%nDomains)%julStart + 1) * nTstepDay
    allocate(mRM_runoff(maxTimeSteps, nGaugesLocal))
    mRM_runoff = nodata_dp

    ! READ GAUGE DATA
    do iGauge = 1, nGaugesLocal
      ! get domain id
      iDomain = gauge%domainId(iGauge)
      ! get start and end dates
      start_tmp = (/evalPer(iDomain)%yStart, evalPer(iDomain)%mStart, evalPer(iDomain)%dStart/)
      end_tmp = (/evalPer(iDomain)%yEnd, evalPer(iDomain)%mEnd, evalPer(iDomain)%dEnd  /)
      ! evaluation gauge
      fName = trim(adjustl(gauge%fname(iGauge)))
      call read_timeseries(trim(fName), udischarge, &
              start_tmp, end_tmp, optimize, opti_function, &
              data_dp_1d, mask = mask_1d, nMeasPerDay = nMeasPerDay)
      data_dp_1d = merge(data_dp_1d, nodata_dp, mask_1d)
      call paste(gauge%Q, data_dp_1d, nodata_dp)
      deallocate (data_dp_1d)
      ! TODO-RIV-TEMP: read temperature at gauge
    end do
    !
    ! inflow gauge
    !
    ! in mhm call InflowGauge%Q has to be initialized -- dummy allocation with period of domain 1 and initialization
    if (nInflowGaugesTotal .EQ. 0) then
      allocate(data_dp_1d(maxval(simPer(:)%julEnd - simPer(:)%julStart + 1)))
      data_dp_1d = nodata_dp
      call paste(InflowGauge%Q, data_dp_1d, nodata_dp)
    else
      do iGauge = 1, nInflowGaugesTotal
        ! get domain id
        iDomain = InflowGauge%domainId(iGauge)
        ! get start and end dates
        start_tmp = (/simPer(iDomain)%yStart, simPer(iDomain)%mStart, simPer(iDomain)%dStart/)
        end_tmp = (/simPer(iDomain)%yEnd, simPer(iDomain)%mEnd, simPer(iDomain)%dEnd  /)
        ! inflow gauge
        fName = trim(adjustl(InflowGauge%fname(iGauge)))
        call read_timeseries(trim(fName), udischarge, &
                start_tmp, end_tmp, optimize, opti_function, &
                data_dp_1d, mask = mask_1d, nMeasPerDay = nMeasPerDay)
        if (.NOT. (all(mask_1d))) then
          call error_message('***ERROR: Nodata values in inflow gauge time series. File: ', trim(fName), &
                  new_line('a'), '          During simulation period from ', num2str(simPer(iDomain)%yStart) &
                  , ' to ', num2str(simPer(iDomain)%yEnd))
        end if
        data_dp_1d = merge(data_dp_1d, nodata_dp, mask_1d)
        call paste(InflowGauge%Q, data_dp_1d, nodata_dp)
        deallocate (data_dp_1d)
      end do
    end if

  end subroutine mrm_read_discharge

  ! ---------------------------------------------------------------------------

  !    NAME
  !        mrm_read_total_runoff

  !    PURPOSE
  !>       \brief read simulated runoff that is to be routed

  !>       \details read spatio-temporal field of total runoff that has been
  !>       simulated by a hydrologic model or land surface model. This
  !>       total runoff will then be aggregated to the level 11 resolution
  !>       and then routed through the stream network.

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iDomain" domain id

  !    HISTORY
  !>       \authors Stephan Thober

  !>       \date Sep 2015

  ! Modifications:
  ! Stephan Thober  Feb 2016 - refactored deallocate statements
  ! Stephan Thober  Sep 2016 - added ALMA convention
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine mrm_read_total_runoff(iDomain)

    use mo_append, only : append
    use mo_constants, only : HourSecs
    use mo_common_constants, only : nodata_dp
    use mo_common_variables, only : ALMA_convention, level1
    use mo_common_datetime_type, only: simPer, timestep
    use mo_mrm_global_variables, only : L1_total_runoff_in, dirTotalRunoff, filenameTotalRunoff, &
                                        varnameTotalRunoff
    use mo_read_nc, only : read_nc

    implicit none

    ! domain id
    integer(i4), intent(in) :: iDomain

    integer(i4) :: tt

    integer(i4) :: nTimeSteps

    ! tell nc file to read daily or hourly values
    integer(i4) :: nctimestep

    ! read data from file
    real(dp), dimension(:, :, :), allocatable :: L1_data

    real(dp), dimension(:, :), allocatable :: L1_data_packed


    if (timestep .eq. 1) nctimestep = -4 ! hourly input
    if (timestep .eq. 24) nctimestep = -1 ! daily input
    call read_nc(trim(dirTotalRunoff(iDomain)), level1(iDomain)%nrows, level1(iDomain)%ncols, &
            varnameTotalRunoff, level1(iDomain)%mask, L1_data, target_period = simPer(iDomain), &
            nctimestep = nctimestep, filename = filenameTotalRunoff)
    ! pack variables
    nTimeSteps = size(L1_data, 3)
    allocate(L1_data_packed(level1(iDomain)%nCells, nTimeSteps))
    do tt = 1, nTimeSteps
      L1_data_packed(:, tt) = pack(L1_data(:, :, tt), mask = level1(iDomain)%mask)
    end do
    ! free space immediately
    deallocate(L1_data)

    ! convert if ALMA conventions have been given
    if (ALMA_convention) then
      ! convert from kg m-2 s-1 to mm TS-1
      ! 1 kg m-2 -> 1 mm depth
      ! multiply with time to obtain per timestep
      L1_data_packed = L1_data_packed * timestep * HourSecs
    end if

    ! append
    call append(L1_total_runoff_in, L1_data_packed(:, :), nodata_dp)

    !free space
    deallocate(L1_data_packed)

  end subroutine mrm_read_total_runoff

  subroutine mrm_read_bankfull_runoff(iDomain)
  ! ---------------------------------------------------------------------------

  !      NAME
  !          mrm_read_bankfull_runoff

  !>         \brief reads the bankfull runoff for approximating the channel widths

  !>         \details reads the bankfull runoff, which can be calculated with
  !>         the script in mhm/post_proc/bankfull_discharge.py

  !     INTENT(IN)
  !>        \param[in] "integer(i4)               :: iDomain"  domain id

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
  !         None

  !     RESTRICTIONS
  !>        \note The file read in must contain a double precision float variable with the name
  !>        "Q_bkfl".

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !         \author Lennart Schueler
  !         \date    May 2018

    use mo_mrm_global_variables, only: level11
    use mo_read_nc, only: read_const_nc
    use mo_mrm_global_variables, only: &
         dirBankfullRunoff, &   ! directory of bankfull_runoff file for each domain
         L11_bankfull_runoff_in ! bankfull runoff at L1

    implicit none

    ! input variables
    integer(i4), intent(in) :: iDomain

    ! local variables
    real(dp), dimension(:,:), allocatable :: L11_data ! read data from file
    real(dp), dimension(:), allocatable :: L11_data_packed

    call read_const_nc(trim(dirBankfullRunoff(iDomain)), &
                               "Q_bkfl", L11_data, &
                               nRows=level11(iDomain)%nrows, &
                               nCols=level11(iDomain)%ncols &
            )

    allocate(L11_data_packed(level11(iDomain)%nCells))
    L11_data_packed(:) = pack(L11_data(:,:), mask=level11(iDomain)%mask)

    ! append
    if (allocated(L11_bankfull_runoff_in)) then
        L11_bankfull_runoff_in = [L11_bankfull_runoff_in, L11_data_packed]
    else
        allocate(L11_bankfull_runoff_in(size(L11_data_packed)))
        L11_bankfull_runoff_in = L11_data_packed
    end if

    deallocate(L11_data)
    deallocate(L11_data_packed)

  end subroutine mrm_read_bankfull_runoff

  !    NAME
  !        rotate_fdir_variable

  !    PURPOSE
  !>       \brief TODO: add description

  !>       \details TODO: add description

  !    INTENT(INOUT)
  !>       \param[inout] "integer(i4), dimension(:, :) :: x"

  !    HISTORY
  !>       \authors L. Samaniego & R. Kumar

  !>       \date Jun 2018

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine rotate_fdir_variable(x)

    use mo_common_constants, only : nodata_i4

    implicit none

    integer(i4), dimension(:, :), intent(INOUT) :: x


    ! 0  -1   0 |
    integer(i4) :: i, j


    !-------------------------------------------------------------------
    ! NOTE:
    !
    !     Since the DEM was transposed from (lat,lon) to (lon,lat), i.e.
    !     new DEM = transpose(DEM_old), then
    !
    !     the flow direction X (which was read) for every i, j cell needs
    !     to be rotated as follows
    !
    !                 X(i,j) = [R] * {uVector}
    !
    !     with
    !     {uVector} denoting the fDir_old (read) in vector form, and
    !               e.g. dir 8 => {-1, -1, 0 }
    !     [R]       denting a full rotation matrix to transform the flow
    !               direction into the new coordinate system (lon,lat).
    !
    !     [R] = [rx][rz]
    !
    !     with
    !           |      1       0      0 |
    !     [rx] =|      0   cos T  sin T | = elemental rotation along x axis
    !           |      0  -sin T  cos T |
    !
    !           |  cos F   sin F      0 |
    !     [rz] =| -sin F   cos F      0 | = elemental rotation along z axis
    !           |      0       0      1 |
    !
    !     and T = pi, F = - pi/2
    !     thus
    !          |  0  -1   0 |
    !     [R] =| -1   0   0 |
    !          |  0   0  -1 |
    !     making all 8 directions the following transformation were
    !     obtained.
    ! Addon, RS (2018-12-12): inverting latitude (now ascending order)
    ! made another conversion necessary
    ! non-inverted replacements are commented!
    !-------------------------------------------------------------------
    do i = 1, size(x, 1)
      do j = 1, size(x, 2)
        if (x(i, j)  == nodata_i4) cycle
        select case (x(i, j))
        case(1)
          !x(i, j) = 4
          x(i, j) = 4
        case(2)
          !x(i, j) = 2
          x(i, j) = 8
        case(4)
          !x(i, j) = 1
          x(i, j) = 16
        case(8)
          !x(i, j) = 128
          x(i, j) = 32
        case(16)
          !x(i, j) = 64
          x(i, j) = 64
        case(32)
          !x(i, j) = 32
          x(i, j) = 128
        case(64)
          !x(i, j) = 16
          x(i, j) = 1
        case(128)
          !x(i, j) = 8
          x(i, j) = 2
        end select
      end do
    end do

  end subroutine rotate_fdir_variable

end module mo_mrm_read_data
