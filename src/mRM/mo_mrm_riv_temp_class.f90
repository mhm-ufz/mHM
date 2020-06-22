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
    real(dp), public, dimension(:), allocatable :: L11_riv_area !river area in L11
    ! PET related vars
    real(dp) :: albedo_water ! albedo of open water
    real(dp) :: pt_a_water ! priestley taylor alpha parameter for PET on open water
    ! \f$ E_L \f$ Generated lateral temperature energy flux [m3 s-1 K] on L1
    real(dp), dimension(:), allocatable :: L1_runoff_E
    ! accumulated later fluxes
    real(dp), dimension(:), allocatable :: L1_acc_inter
    real(dp), dimension(:), allocatable :: L1_acc_direct
    real(dp), dimension(:), allocatable :: L1_acc_base
    ! vars for routing
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
    procedure :: init_cond_from_L1 => meth_init_cond_from_L1
    procedure :: calc_source_E => meth_calc_source_E
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
    use mo_common_variables, only : level0, domainMeta

    implicit none

    class(riv_temp_type), intent(inout) :: self
    integer(i4), intent(in) :: nCells

    real(dp), dimension(:), allocatable :: dummy_Vector11
    real(dp), dimension(:, :), allocatable :: dummy_Matrix11_IT

    print *, 'init river temp'
    ! TODO-RIV-TEMP:
    ! - init_condition riv_temp (0)
    ! - see: variables_alloc_routing

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
    self &
  )
    implicit none

    class(riv_temp_type), intent(inout) :: self

    print *, 'init river area'
    ! TODO-RIV-TEMP:
    ! - needs riv-widths
    ! - walk the riv-network in order to get link-length

  end subroutine meth_init_area


  subroutine meth_prepare_meteo( &
    self &
  )
    implicit none

    class(riv_temp_type), intent(inout) :: self

    print *, 'riv-temp: prepare_meteo'
    ! TODO-RIV-TEMP:
    ! - see: prepare_meteo_forcings_data

  end subroutine meth_prepare_meteo


  subroutine meth_init_cond_from_L1( &
    self &
  )
    implicit none

    class(riv_temp_type), intent(inout) :: self

    print *, 'init river temperature with L1 air temp'
    ! TODO-RIV-TEMP:
    ! - dis-agg L1 air temp for tt=1
    ! - init with 0 for now

  end subroutine meth_init_cond_from_L1

  subroutine meth_alloc_lateral( &
    self, &
    size &
  )
    implicit none

    class(riv_temp_type), intent(inout) :: self
    integer(i4), intent(in) :: size

    print *, 'riv-temp: allocate later runoff components from L1'
    allocate(self%L1_acc_inter(size))
    allocate(self%L1_acc_base(size))
    allocate(self%L1_acc_direct(size))
    allocate(self%L1_runoff_E(size))
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

  subroutine meth_calc_source_E( &
    self, &
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
    fday_strd &
  )

    implicit none

    class(riv_temp_type), intent(inout) :: self
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
    real(dp), intent(in) :: ssrd_day
    ! Daily mean longwave radiation
    real(dp), intent(in) :: strd_day
    ! Daytime fraction of ssrd
    real(dp), intent(in) :: fday_ssrd
    ! Daytime fraction of strd
    real(dp), intent(in) :: fday_strd
    ! Nighttime fraction of ssrd
    ! switch to turn on adding energy if routing TS is larger then mhm TS
    ! logical, intent(in), optional :: add

    ! if ( do_add ) then
    !     self%L1_lateral_E = self%L1_lateral_E + lateral_E
    ! else
    !     self%L1_lateral_E = lateral_E
    ! end if

  end subroutine meth_calc_source_E

  SUBROUTINE calc_L1_runoff_E( &
    fSealed_area_fraction, &
    fast_interflow, &
    slow_interflow, &
    baseflow, &
    direct_runoff, &
    temp_air, &
    mean_temp_air, &
    lateral_E &
  )

    implicit none

    ! sealed area fraction [1]
    REAL(dp), dimension(:), INTENT(IN) :: fSealed_area_fraction
    ! \f$ q_0 \f$ Fast runoff component [mm tst-1]
    REAL(dp), dimension(:), INTENT(IN) :: fast_interflow
    ! \f$ q_1 \f$ Slow runoff component [mm tst-1]
    REAL(dp), dimension(:), INTENT(IN) :: slow_interflow
    ! \f$ q_2 \f$ Baseflow [mm tsts-1]
    REAL(dp), dimension(:), INTENT(IN) :: baseflow
    ! \f$ q_D \f$ Direct runoff from impervious areas  [mm tst-1]
    REAL(dp), dimension(:), INTENT(IN) :: direct_runoff
    ! air temperature [K]
    real(dp), dimension(:), intent(in) :: temp_air
    ! annual mean air temperature [K]
    real(dp), dimension(:), intent(in) :: mean_temp_air
    ! \f$ q_T \f$ Generated lateral Energy [K mm tst-1]
    REAL(dp), dimension(:), INTENT(OUT) :: lateral_E

    lateral_E = ( &
      (baseflow * mean_temp_air + (slow_interflow + fast_interflow) * temp_air) * (1.0_dp - fSealed_area_fraction) &
      + direct_runoff * temp_air * fSealed_area_fraction &
    )

  END SUBROUTINE calc_L1_runoff_E

  elemental pure subroutine temporal_disagg_radiation( &
      isday, ntimesteps_day, &
      ssrd_day, strd_day, &
      fday_ssrd, fday_strd, &
      fnight_ssrd, fnight_strd, &
      ssrd, strd &
  )
    implicit none

    ! is day or night
    logical, intent(in) :: isday
    ! # of time steps per day
    real(dp), intent(in) :: ntimesteps_day
    ! Daily mean shortwave radiation
    real(dp), intent(in) :: ssrd_day
    ! Daily mean longwave radiation
    real(dp), intent(in) :: strd_day
    ! Daytime fraction of ssrd
    real(dp), intent(in) :: fday_ssrd
    ! Daytime fraction of strd
    real(dp), intent(in) :: fday_strd
    ! Nighttime fraction of ssrd
    real(dp), intent(in) :: fnight_ssrd
    ! Nighttime fraction of strd
    real(dp), intent(in) :: fnight_strd
    ! Actual ssrd
    real(dp), intent(out) :: ssrd
    ! Actual strd
    real(dp), intent(out) :: strd

    ! default vaule used if ntimesteps_day = 1 (i.e., e.g. daily values)
    ssrd = ssrd_day
    strd = strd_day

    ! Distribute into time steps night/day
    if(ntimesteps_day .gt. 1.0_dp) then
      if (isday) then ! DAY-TIME
        ssrd = 2.0_dp * ssrd_day * fday_ssrd / ntimesteps_day
        strd = 2.0_dp * strd_day * fday_strd / ntimesteps_day
      else            ! NIGHT-TIME
        ssrd = 2.0_dp * ssrd_day * fnight_ssrd / ntimesteps_day
        strd = 2.0_dp * strd_day * fnight_strd / ntimesteps_day
      end if
    end if

  end subroutine temporal_disagg_radiation

  SUBROUTINE L11_E_acc(E_L1, efecArea, L1_L11_Id, L11_areaCell, L11_L1_Id, TS, map_flag, E_acc)

    use mo_common_constants, only : HourSecs, nodata_dp

    implicit none

    ! total runoff L1 [mm tst-1]
    real(dp), intent(in), dimension(:) :: E_L1
    ! effective area in [km2] at Level 1
    real(dp), intent(in), dimension(:) :: efecarea
    ! L11 Ids mapped on L1
    integer(i4), intent(in), dimension(:) :: L1_L11_Id
    ! effective area in [km2] at Level 11
    real(dp), intent(in), dimension(:) :: L11_areacell
    ! L1 Ids mapped on L11
    integer(i4), intent(in), dimension(:) :: L11_L1_Id
    ! time step in [s]
    integer(i4), intent(in) :: TS
    ! Flag indicating whether routing resolution is higher than hydrologic one
    logical, intent(in) :: map_flag
    ! aggregated runoff at L11 [m3 s-1]
    real(dp), intent(out), dimension(:) :: E_acc
    integer(i4) :: k
    ! [s] time step
    real(dp) :: TST

    ! ------------------------------------------------------------------
    ! ACCUMULATION OF DISCHARGE TO A ROUTING CELL
    ! ------------------------------------------------------------------
    ! Hydrologic timestep in seconds
    TST = HourSecs * TS

    if (map_flag) then
      E_acc = 0._dp
      ! loop over high-resolution cells (L1) and add discharge to
      ! corresponding low-resolution cells (L11)
      do k = 1, size(E_L1, 1)
        E_acc(L1_L11_Id(k)) = E_acc(L1_L11_Id(k)) + E_L1(k) * efecArea(k)
      end do
      E_acc = E_acc * 1000.0_dp / TST
      !
    else
      ! initialize qout
      E_acc = nodata_dp
      do k = 1, size(E_acc, 1)
        ! map energy flux from coarse L1 resolution to fine L11 resolution
        E_acc(k) = E_L1(L11_L1_Id(k))
      end do
      ! adjust energy flux by area cell
      E_acc(:) = E_acc(:) * L11_areaCell(:) * 1000.0_dp / TST
    end if

  END SUBROUTINE L11_E_acc

end module mo_mrm_riv_temp_class

