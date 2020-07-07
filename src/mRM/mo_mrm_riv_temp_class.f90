!>       \file mo_mrm_riv_temp_class.f90

!>       \brief Class for the river temperature calculations
!>       \details

!>       \authors Sebastian Mueller

!>       \date Apr 2020

module mo_mrm_riv_temp_class

  use mo_kind, only: dp, i4

  implicit none

  private
  public :: riv_temp_type

  ! This is a container to define the river temperature routing in the current time step
  type riv_temp_type
    ! config settings
    character(256) :: nml_name = 'config_riv_temp' ! file name for river widths
    ! riv geometry
    character(256), dimension(:), allocatable :: dir_riv_widths ! Directory where river widths are stored
    character(256) :: riv_widths_file ! file name for river widths
    character(256) :: riv_widths_name ! variable name for river widths
    real(dp), public, dimension(:), allocatable :: L11_riv_widths !river widths in L11
    real(dp), public, dimension(:), allocatable :: L11_riv_areas !river area in L11
    ! PET related vars
    real(dp) :: albedo_water ! albedo of open water
    real(dp) :: pt_a_water ! priestley taylor alpha parameter for PET on open water
    ! \f$ E_L \f$ Generated lateral temperature energy flux [m3 s-1 K] on L1
    ! accumulated later fluxes
    real(dp), dimension(:), allocatable :: L1_runoff_E
    real(dp), dimension(:), allocatable :: L1_acc_srad
    real(dp), dimension(:), allocatable :: L1_acc_lrad
    real(dp), dimension(:), allocatable :: L1_acc_temp
    real(dp), dimension(:), allocatable :: L1_acc_tann
    integer(i4) :: ts_cnt ! time-step counter for accumulation of meteo
    ! vars for routing
    integer(i4) :: s11 ! starting index for current L11 domain
    integer(i4) :: e11 ! ending index for current L11 domain
    real(dp), dimension(:,:), allocatable :: netNode_E_IN ! Total energy inputs at t-1 and t
    real(dp), dimension(:,:), allocatable :: netNode_E_R ! energy leaving at t-1 and t
    real(dp), dimension(:), allocatable :: netNode_E_mod ! Simulated routed energy
    real(dp), dimension(:), allocatable :: netNode_E_out ! total energy source from cell in L11
    ! results
    real(dp), dimension(:), allocatable :: river_temp ! variable containing river temp for each domain and gauge
  contains
    procedure :: config => meth_config
    procedure :: init => meth_init
    procedure :: init_area => meth_init_area
    procedure :: acc_source_E => meth_acc_source_E
    procedure :: alloc_lateral => meth_alloc_lateral
    procedure :: dealloc_lateral => meth_dealloc_lateral
  end type riv_temp_type

contains

  subroutine meth_config( &
    self, &
    file_namelist, &
    unamelist, &
    file_namelist_param, &
    unamelist_param &
  )

    use mo_common_constants, only : maxNoDomains, nodata_i4
    use mo_common_variables, only : domainMeta
    use mo_nml, only : close_nml, open_nml, position_nml

    implicit none

    class(riv_temp_type), intent(inout) :: self
    character(*), intent(in) :: file_namelist, file_namelist_param
    integer, intent(in) :: unamelist, unamelist_param

    character(256), dimension(maxNoDomains) :: dir_riv_widths
    ! parameter to read in
    real(dp) :: albedo_water ! albedo of open water
    real(dp) :: pt_a_water ! priestley taylor alpha parameter for PET on open water
    character(256) :: riv_widths_file ! file name for river widths
    character(256) :: riv_widths_name ! variable name for river widths

    ! namelist for river temperature configuration
    namelist /config_riv_temp/ &
      albedo_water, &
      pt_a_water, &
      riv_widths_file, &
      riv_widths_name, &
      dir_riv_widths

    ! TODO-RIV-TEMP:
    print *, '    read config: river temperature routing'

    ! allocate the directory arrays
    allocate(self%dir_riv_widths(domainMeta%nDomains))

    ! open the namelist file
    call open_nml(file_namelist, unamelist, quiet=.true.)
    ! find the river-temp config namelist
    call position_nml(self%nml_name, unamelist)
    ! read the river-temp config namelist
    read(unamelist, nml=config_riv_temp)

    self%albedo_water = albedo_water
    self%pt_a_water = pt_a_water
    self%riv_widths_file = riv_widths_file
    self%riv_widths_name = riv_widths_name
    self%dir_riv_widths = dir_riv_widths

    ! closing the namelist file
    call close_nml(unamelist)

  end subroutine meth_config


  subroutine meth_init( &
    self, &
    nCells &
  )
    use mo_append, only : append
    use mo_mrm_constants, only : nRoutingStates
    ! use mo_common_variables, only : level0, domainMeta

    implicit none

    class(riv_temp_type), intent(inout) :: self
    integer(i4), intent(in) :: nCells

    real(dp), dimension(:), allocatable :: dummy_Vector11
    real(dp), dimension(:, :), allocatable :: dummy_Matrix11_IT

    print *, 'init river temp'
    ! TODO-RIV-TEMP:
    ! - init_condition riv_temp (0)
    ! - see: variables_alloc_routing

    ! set sub time-step counter to 0
    self%ts_cnt = 0_i4

    ! dummy vector and matrix
    allocate(dummy_Vector11   (nCells))
    allocate(dummy_Matrix11_IT(nCells, nRoutingStates))

    ! simulated energy flux at each node
    dummy_Vector11(:) = 0.0_dp
    call append(self%netNode_E_mod, dummy_Vector11)
    ! simulated river temperature at each node
    dummy_Vector11(:) = 0.0_dp
    call append(self%river_temp, dummy_Vector11)
    ! Total outflow from cells L11 at time tt
    dummy_Vector11(:) = 0.0_dp
    call append(self%netNode_E_out, dummy_Vector11)
    ! Total discharge inputs at t-1 and t
    dummy_Matrix11_IT(:, :) = 0.0_dp
    call append(self%netNode_E_IN, dummy_Matrix11_IT)
    !  Routed outflow leaving a node
    dummy_Matrix11_IT(:, :) = 0.0_dp
    call append(self%netNode_E_R, dummy_Matrix11_IT)

    ! free space
    if (allocated(dummy_Vector11)) deallocate(dummy_Vector11)
    if (allocated(dummy_Matrix11_IT)) deallocate(dummy_Matrix11_IT)

  end subroutine meth_init

  subroutine meth_init_area( &
    self, &
    iDomain, &
    L11_netPerm, &
    L11_fromN, &
    L11_length, &
    nLinks, &
    nCells, &
    nrows, &
    ncols, &
    L11_mask &
  )
    use mo_append, only : append
    use mo_read_forcing_nc, only: read_const_forcing_nc

    implicit none

    class(riv_temp_type), intent(inout) :: self
    ! Domain counter
    integer(i4), intent(in) :: iDomain
    ! L11 routing order
    integer(i4), dimension(:), intent(in) :: L11_netPerm
    ! L11 source grid cell order
    integer(i4), dimension(:), intent(in) :: L11_fromN
    ! L11 link length
    real(dp), dimension(:), intent(in) :: L11_length
    ! number of L11 links in the current domain
    integer(i4), intent(in) :: nLinks
    ! number of L11 cells of the current domain
    integer(i4), intent(in) :: nCells
    integer(i4), intent(in) :: ncols     ! Number of columns
    integer(i4), intent(in) :: nrows     ! Number of rows
    ! the mask for valid cells in the original grid (nrows*ncols)
    logical, dimension(:, :), intent(in) :: L11_mask

    logical, dimension(:,:), allocatable :: mask
    real(dp), dimension(:,:), allocatable :: L11_data ! read data from file
    real(dp), dimension(:), allocatable :: L11_riv_widths, L11_riv_areas

    integer(i4) :: i, k, iNode

    print *, 'init river area'
    call read_const_forcing_nc(&
      trim(self%dir_riv_widths(iDomain)), &
      nrows, &
      ncols, &
      self%riv_widths_name, &
      mask, &
      L11_data, &
      self%riv_widths_file &
    )

    allocate(L11_riv_widths(nCells))
    L11_riv_widths(:) = pack(L11_data(:,:), mask=L11_mask)
    call append(self%L11_riv_widths, L11_riv_widths)

    allocate(L11_riv_areas(nCells))
    ! at area at Outlet is 0 (correct?)
    L11_riv_areas = 0.0_dp
    do k = 1, nLinks
      ! get LINK routing order -> i
      i = L11_netPerm(k)
      iNode = L11_fromN(i)
      L11_riv_areas(iNode) = L11_length(i) * L11_riv_widths(iNode)
    end do
    call append(self%L11_riv_areas, L11_riv_areas)

    deallocate(L11_data)
    deallocate(L11_riv_widths)
    deallocate(L11_riv_areas)

  end subroutine meth_init_area


  subroutine meth_reset_time(self)
    implicit none

    class(riv_temp_type), intent(inout) :: self

    self%L1_acc_inter = 0.0_dp
    self%L1_acc_base = 0.0_dp
    self%L1_acc_direct = 0.0_dp
    self%L1_runoff_E = 0.0_dp

  end subroutine meth_reset_time

  subroutine meth_alloc_lateral( &
    self, &
    nCells &
  )
    implicit none

    class(riv_temp_type), intent(inout) :: self
    integer(i4), intent(in) :: nCells

    print *, 'riv-temp: allocate later runoff components from L1'
    allocate(self%L1_acc_inter(nCells))
    allocate(self%L1_acc_base(nCells))
    allocate(self%L1_acc_direct(nCells))
    allocate(self%L1_runoff_E(nCells))
    ! do we need different arrays for these?!
    self%L1_acc_inter = 0.0_dp
    self%L1_acc_base = 0.0_dp
    self%L1_acc_direct = 0.0_dp
    self%L1_runoff_E = 0.0_dp

  end subroutine meth_alloc_lateral

  subroutine meth_dealloc_lateral( &
    self &
  )
    implicit none

    class(riv_temp_type), intent(inout) :: self

    print *, 'riv-temp: deallocate later runoff components from L1'
    deallocate(self%L1_acc_inter)
    deallocate(self%L1_acc_base)
    deallocate(self%L1_acc_direct)
    deallocate(self%L1_runoff_E)

  end subroutine meth_dealloc_lateral

  subroutine meth_acc_source_E( &
    self, &
    time, &
    rout_case, &
    tsRoutFactor, &
    fSealed_area_fraction, &
    fast_interflow, &
    slow_interflow, &
    baseflow, &
    direct_runoff, &
    temp_air, &
    mean_temp_air, &
    ssrd_day, &
    strd_day, &
    fday_ssrd, &
    fnight_ssrd, &
    fday_strd, &
    fnight_strd, &
    do_rout &
  )

    use mo_julian, only : dec2date

    implicit none

    class(riv_temp_type), intent(inout) :: self
    ! current decimal Julian day
    real(dp), intent(in) :: time
    ! processMatrix(8, 1) -> routing process case
    integer(i4), intent(in) :: rout_case
    ! factor between routing and hydrological modelling resolution
    real(dp), intent(in) :: tsRoutFactor
    ! sealed area fraction [1]
    real(dp), dimension(:), intent(in) :: fSealed_area_fraction
    ! \f$ q_0 \f$ Fast runoff component [mm tst-1]
    real(dp), dimension(:), intent(in) :: fast_interflow
    ! \f$ q_1 \f$ Slow runoff component [mm tst-1]
    real(dp), dimension(:), intent(in) :: slow_interflow
    ! \f$ q_2 \f$ Baseflow [mm tsts-1]
    real(dp), dimension(:), intent(in) :: baseflow
    ! \f$ q_D \f$ Direct runoff from impervious areas  [mm tst-1]
    real(dp), dimension(:), intent(in) :: direct_runoff
    ! air temperature [K]
    real(dp), dimension(:), intent(in) :: temp_air
    ! annual mean air temperature [K]
    real(dp), dimension(:), intent(in) :: mean_temp_air
    ! Daily mean short radiation
    real(dp), dimension(:), intent(in) :: ssrd_day
    ! Daily mean longwave radiation
    real(dp), dimension(:), intent(in) :: strd_day
    ! Daytime fraction of ssrd
    real(dp), dimension(:), intent(in) :: fday_ssrd
    ! Nighttime fraction of ssrd
    real(dp), dimension(:), intent(in) :: fnight_ssrd
    ! Daytime fraction of strd
    real(dp), dimension(:), intent(in) :: fday_strd
    ! Nighttime fraction of strd
    real(dp), dimension(:), intent(in) :: fnight_strd
    ! flag for performing routing
    logical, intent(in) :: do_rout

    ! is day or night
    logical :: isday
    ! current hour of a given day
    integer(i4) :: hour

    call dec2date(time, hh=hour)
    ! flag for day or night depending on hours of the day
    isday = (hour .gt. 6) .AND. (hour .le. 18)

    ! switch to turn on adding energy if routing TS is larger then mhm TS

    ! if ( do_add ) then
    !     self%L1_lateral_E = self%L1_lateral_E + lateral_E
    ! else
    !     self%L1_lateral_E = lateral_E
    ! end if

    ! if (rout_case .eq. 1) then
    !   RunToRout = L1_total_runoff(s1 : e1) ! runoff [mm TST-1] mm per timestep
    !   InflowDischarge = InflowGauge%Q(iDischargeTS, :) ! inflow discharge in [m3 s-1]
    ! else if ((rout_case .eq. 2) .or. (rout_case .eq. 3)) then
    !   if (tsRoutFactor .lt. 1._dp) then
    !     ! routing timesteps are shorter than hydrologic time steps
    !     RunToRout = L1_total_runoff(s1 : e1) ! runoff [mm TST-1] mm per timestep
    !   else
    !     ! routing timesteps are longer than hydrologic time steps
    !     RunToRout = RunToRout + L1_total_runoff(s1 : e1)
    !   end if
    ! end if

  end subroutine meth_acc_source_E

end module mo_mrm_riv_temp_class

