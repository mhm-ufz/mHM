module mo_mrm_river_head
  use mo_common_variables,     only : level0, domainMeta
  use mo_mrm_global_variables, only : L0_L11_remap, L11_bankfull_runoff_in, &
                                      L0_slope, L0_channel_elevation
  use mo_kind,                 only : i4, dp
  use mo_common_constants,     only : nodata_dp
  use mo_append,               only : append
  use mo_netcdf,               only : NcVariable


  implicit none

  private

  public ::  init_masked_zeros_l0 ! allocates memory
  ! calculates the channel elevation from the bankfull river discharge
  public :: calc_channel_elevation
  public :: calc_river_head ! calculates the river head
  ! calculates the monthly averages and writes the river heads to nc file
  public :: avg_and_write_timestep
  public :: create_output ! creates the NetCDF file
  ! resets the summation variable for the river head averages
  private :: reset_sum


  ! All dims = nDomains
  ! NC variables kept in memory for writing timesteps
  type(NcVariable), dimension(:), allocatable :: nc_time, nc_riverhead
  ! counter for the written timesteps
  integer(i4), dimension(:), allocatable :: time_counter
  ! counter for the monthly averages of the river head
  integer(i4), dimension(:), allocatable :: sum_counter


  contains


  subroutine init_masked_zeros_l0(iDomain, data)
    integer(i4), intent(in) :: iDomain
    real(dp), dimension(:), allocatable, intent(inout) :: data
    real(dp), dimension(:), allocatable :: dummy_1D

    allocate(dummy_1D(level0(iDomain)%nCells))
    dummy_1D = nodata_dp
    call append(data, dummy_1D)
    deallocate(dummy_1D)
    call reset_sum(iDomain, data)
  end subroutine init_masked_zeros_l0


  subroutine reset_sum(iDomain, data)
    ! Sets all elements of data to 0.0 where mask0 is true
    integer(i4), intent(in) :: iDomain
    real(dp), dimension(:), intent(inout) :: data
    integer(i4) :: i, j, k

    do k = 1, level0(iDomain)%nCells
      i = level0(iDomain)%CellCoor(k, 1)
      j = level0(iDomain)%CellCoor(k, 2)
      if (level0(iDomain)%mask(i,j)) then
        data(k) = 0.0
      end if
    end do
  end subroutine reset_sum


  subroutine calc_channel_elevation()
    use mo_common_constants, only : nodata_i4
    use mo_common_variables, only : domainMeta, L0_elev
    use mo_mrm_global_variables, only : L0_fDir, L0_fAcc, L0_channel_depth
    real(dp), dimension(:,:), allocatable :: channel_dpth
    real(dp), dimension(:,:), allocatable :: channel_elev
    real(dp), dimension(:,:), allocatable :: slope
    real(dp), dimension(:,:), allocatable :: elev0
    integer(i4), dimension(:,:), allocatable :: fDir0
    integer(i4), dimension(:,:), allocatable :: fAcc0
    real(dp) n ! Manning's roughness coefficient
    integer(i4) :: nrows0, ncols0
    integer(i4) :: s0, e0
    integer(i4) i, j, k
    integer(i4) iDomain

    n = .045_dp ! m^-1/3 s from Sutanudjaja et al. 2011

    do iDomain = 1, domainMeta%nDomains
      nrows0 = level0(iDomain)%nrows
      ncols0 = level0(iDomain)%ncols
      s0 = level0(iDomain)%iStart
      e0 = level0(iDomain)%iEnd
      allocate(channel_dpth(nrows0, ncols0))
      allocate(channel_elev(nrows0, ncols0))
      allocate(elev0(nrows0, ncols0))
      allocate(fDir0(nrows0, ncols0))
      allocate(fAcc0(nrows0, ncols0))
      allocate(slope(nrows0, ncols0))
      channel_dpth(:,:) = nodata_dp
      channel_elev(:,:) = nodata_dp
      slope(:,:) = nodata_dp

      elev0(:,:) = unpack(L0_elev(s0:e0), level0(iDomain)%mask, nodata_dp)
      fDir0(:,:) = unpack(L0_fDir(s0:e0), level0(iDomain)%mask, nodata_i4)
      fAcc0(:,:) = unpack(L0_fAcc(s0:e0), level0(iDomain)%mask, nodata_i4)

      do k = 1, level0(iDomain)%nCells
        i = level0(iDomain)%CellCoor(k, 1)
        j = level0(iDomain)%CellCoor(k, 2)
        if (fAcc0(i,j) > 1) then
          slope(i,j) = calc_slope(iDomain, elev0, fDir0, i, j)
          channel_dpth(i,j) = ((n * &
  L11_bankfull_runoff_in(L0_L11_remap(iDomain)%lowres_id_on_highres(i,j))) / &
  (4.8 * slope(i,j)**.5))**.6
          channel_elev(i,j) = elev0(i,j) - channel_dpth(i,j)
        end if
      end do

      call append(L0_channel_depth, pack(channel_dpth(:,:), &
                  level0(iDomain)%mask))
      call append(L0_channel_elevation, pack(channel_elev(:,:), &
                  level0(iDomain)%mask))
      call append(L0_slope, pack(slope(:,:), &
                  level0(iDomain)%mask))

      deallocate(channel_dpth)
      deallocate(channel_elev)
      deallocate(elev0)
      deallocate(fDir0)
      deallocate(fAcc0)
      deallocate(slope)
    end do
  end subroutine calc_channel_elevation


  subroutine calc_river_head(iDomain, L11_Qmod, river_head)
    integer(i4), intent(in) :: iDomain
    real(dp), dimension(:), intent(in) :: L11_Qmod
    real(dp), dimension(:), allocatable, intent(inout) :: river_head
    real(dp) :: n ! Manning's roughness coefficient
    integer(i4) :: s0, e0
    integer(i4) i, j, k, L11_ind

    n = .045_dp ! m^-1/3 s from Sutanudjaja et al. 2011

    s0 = level0(iDomain)%iStart
    e0 = level0(iDomain)%iEnd

    do k = s0, e0
      i = level0(iDomain)%CellCoor(k, 1)
      j = level0(iDomain)%CellCoor(k, 2)
      if (i >= 0 .and. i < 99999 .and. j  >= 0 .and. j < 99999) then
        L11_ind = L0_L11_remap(iDomain)%lowres_id_on_highres(i,j)
        ! TODO L11_Qmid(L11_ind) causes IEEE_UNDERFLOW_FLAG IEEE_DENORMAL
        river_head(k) = river_head(k) + L0_channel_elevation(k) + &
            (n * L11_Qmod(L11_ind) / L11_bankfull_runoff_in(L11_ind) / &
            L0_slope(k)**.5)**.6
      end if
    end do

    sum_counter(iDomain) = sum_counter(iDomain) + 1
  end subroutine calc_river_head


  function calc_slope(iDomain, elev0, fDir0, i, j) result(slope)
    use mo_common_variables, only: iFlag_cordinate_sys
    use mo_mrm_net_startup,      only: cellLength, moveDownOneCell
    integer(i4), intent(in) :: iDomain
    integer(i4), intent(in) :: i, j
    integer(i4), intent(in), dimension(:,:), allocatable :: fDir0
    real(dp), intent(in), dimension(:,:), allocatable :: elev0
    real(dp) :: slope, length
    integer(i4) :: i_down, j_down

    call cellLength(iDomain, fDir0(i, j), i, j, &
                    iFlag_cordinate_sys, length)
    i_down = i
    j_down = j
    call moveDownOneCell(fDir0(i, j), i_down, j_down)

    slope = (elev0(i,j) - elev0(i_down, j_down)) / length
    ! TODO: as soon as current gfortran compiler is available on EVE,
    ! use ieee_isnan from ieee_arithmetic module, instead of
    ! slope /= slope
    if(slope < 0.0001_dp .OR. slope > 20._dp .OR. slope /= slope) then
      slope = 0.0001_dp
    end if
  end function calc_slope


  subroutine avg_and_write_timestep(iDomain, timestep, data)
    integer(i4), intent(in) :: iDomain
    integer(i4), intent(in) :: timestep
    real(dp), intent(inout), dimension(:) :: data

    call nc_time(iDomain)%setData(timestep, [time_counter(iDomain)])
    ! -> index along the time dimension of the netcdf variable
    call nc_riverhead(iDomain)%setData(unpack(data / sum_counter(iDomain), &
        level0(iDomain)%mask, nodata_dp), [1, 1, time_counter(iDomain)])
    time_counter(iDomain) = time_counter(iDomain) + 1
    call reset_sum(iDomain, data)
    sum_counter(iDomain) = 0
  end subroutine avg_and_write_timestep


  subroutine create_output(iDomain, OutPath)
    use mo_common_mhm_mrm_variables, only : evalPer
    use mo_netcdf, only : NcDataset, NcDimension
    use mo_string_utils, only : num2str
    use mo_julian, only : dec2date
    use mo_grid, only : geoCoordinates, mapCoordinates
    use mo_file, only : version
    use mo_common_variables, only : project_details, setup_description, &
        simulation_type, Conventions, contact, mHM_details, history
    use mo_message, only : message

    ! number of domains
    integer(i4), intent(in) :: iDomain
    ! list of Output paths per Domain
    character(256), intent(in) :: OutPath
    character(256) :: filename
    real(dp), allocatable, dimension(:) :: easting, northing
    real(dp), allocatable, dimension(:, :) :: lat, lon
    character(128) :: unit, date, time, datetime
    integer(i4) :: day, month, year
    type(NcDataset) :: nc
    type(NcDimension), dimension(3) :: dimids0
    type(NcVariable) :: var

    filename = trim(OutPath) // 'mRM_riverhead_' // trim(num2str(iDomain, '(i3.3)')) // '.nc'

    call message('    Writing mRM groundwater coupling file to ' // trim(filename) // ' ...')

    call mapCoordinates( level0(iDomain), northing, easting )
    call geoCoordinates( level0(iDomain), lat, lon )

    nc = NcDataset(filename, "w")

    dimids0 = [ nc%setDimension("easting", size(easting)), &
                nc%setDimension("northing", size(northing)), &
                nc%setDimension("time", 0) &
              ]

    if (.not. allocated(sum_counter)) then
        allocate(sum_counter(domainMeta%nDomains))
    end if
    sum_counter(iDomain) = 0

    ! time units
    call dec2date( real(evalPer(iDomain)%julStart, dp), dd = day, mm = month, yy = year )
    write(unit, "('hours since ', i4, '-' ,i2.2, '-', i2.2, 1x, '00:00:00')") year, month, day

    ! time
    if (.not. allocated(time_counter)) then
        allocate(time_counter(domainMeta%nDomains))
    end if
    time_counter(iDomain) = 1
    if (.not. allocated(nc_time)) then
        allocate(nc_time(domainMeta%nDomains))
    end if
    var = nc%setVariable( "time", "i32", [ dimids0(3) ] )
    call var%setAttribute( "units", unit )
    call var%setAttribute( "long_name", "time" )
    nc_time(iDomain) = var

    ! northing
    var = nc%setVariable( "northing", "f64", [ dimids0(2) ] )
    call var%setFillValue(nodata_dp)
    call var%setData( northing )
    call var%setAttribute( "units", "m or deg. dec." )
    call var%setAttribute( "long_name", "y-coordinate in the given coordinate system" )
    call var%setAttribute( "missing_value", nodata_dp )

    ! easting
    var = nc%setVariable( "easting", "f64", [ dimids0(1) ] )
    call var%setFillValue(nodata_dp)
    call var%setData( easting )
    call var%setAttribute( "units", "m or deg. dec." )
    call var%setAttribute( "long_name", "x-coordinate in the given coordinate system" )
    call var%setAttribute( "missing_value", nodata_dp )

    ! lon
    var = nc%setVariable( "lon", "f64", dimids0(1 : 2) )
    call var%setFillValue(nodata_dp)
    call var%setData( lon )
    call var%setAttribute( "units", "deg. dec." )
    call var%setAttribute( "long_name", "longitude" )
    call var%setAttribute( "missing_value", nodata_dp )

    ! lat
    var = nc%setVariable( "lat", "f64", dimids0(1 : 2) )
    call var%setFillValue(nodata_dp)
    call var%setData( lat )
    call var%setAttribute( "units", "deg. dec." )
    call var%setAttribute( "long_name", "latitude" )
    call var%setAttribute( "missing_value", nodata_dp )

    if (.not. allocated(nc_riverhead)) then
        allocate(nc_riverhead(domainMeta%nDomains))
    end if
    var = nc%setVariable( "riverhead", "f64", dimids0(:) )
    call var%setAttribute( "units", "m" )
    call var%setAttribute( "long_name", "simulated riverhead at each node at level 0" )
    call var%setAttribute( "missing_value", nodata_dp )
    call var%setFillValue(nodata_dp)
    nc_riverhead(iDomain) = var

    ! global attributes
    call date_and_time( date = date, time = time )
    write(datetime, "(a4,'-',a2,'-',a2,1x,a2,':',a2,':',a2)") date(1 : 4), &
            date(5 : 6), date(7 : 8), time(1 : 2), time(3 : 4), time(5 : 6)

    call nc%setAttribute( "project", project_details )
    call nc%setAttribute( "setup_description", setup_description )
    call nc%setAttribute( "simulation_type", simulation_type )
    call nc%setAttribute( "Conventions", Conventions )
    call nc%setAttribute( "contact", contact )
    call nc%setAttribute( "mHM_details", trim(mHM_details) // ", release mHMv" // trim(version) )
    call nc%setAttribute( "history", trim(datetime) // ", " // history )
    call nc%setAttribute( "title", "mHMv"//trim(version)//" "//trim(simulation_type)//" outputs" )
    call nc%setAttribute( "creation_date", datetime )
  end subroutine create_output


end module mo_mrm_river_head
