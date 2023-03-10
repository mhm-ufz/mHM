!> \file    mo_mrm_riv_temp_class.f90
!> \brief   \copybrief mo_mrm_riv_temp_class
!> \details \copydetails mo_mrm_riv_temp_class

!> \brief   Class for the river temperature calculations
!> \details River temperature routing on top of mRM.
!> \warning This feature is still experimental! River freezing is still missing.
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    Sep 2020
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mrm
module mo_mrm_riv_temp_class

  use mo_kind, only: dp, i4

  implicit none

  private

  !> \class   riv_temp_type
  !> \brief   This is a container to define the river temperature routing in the current time step
  !> \details This class provides all procedures to rout river tmperature through the river network.
  !> \warning This feature is still experimental!
  !!          This first version doesn't provide ice covering, which means, that the river temperature
  !!          can drop below 0 [deg C] in winter.
  type, public :: riv_temp_type
    ! config settings
    logical :: active = .false. !< state if this process is active
    integer(i4) :: case = 0_i4 !< the selected process-case option
    character(256) :: nml_name = 'config_riv_temp' !< namelist name in mhm namelist
    ! riv geometry
    character(256), dimension(:), allocatable :: dir_riv_widths !< Directory where river widths are stored
    character(256) :: riv_widths_file !< file name for river widths
    character(256) :: riv_widths_name !< variable name for river widths
    real(dp), public, dimension(:), allocatable :: L11_riv_widths !< river widths in L11
    real(dp), public, dimension(:), allocatable :: L11_riv_areas !< river area in L11
    ! PET related vars
    real(dp) :: albedo_water !< albedo of open water
    real(dp) :: pt_a_water   !< priestley taylor alpha parameter for PET on open water
    ! sensible heat flux related
    real(dp) :: emissivity_water   !< emissivity of water
    real(dp) :: turb_heat_ex_coeff !< lateral heat exchange coefficient water <-> air
    !> cutoff value for temperature
    real(dp) :: delta_T = 0.1_dp
    ! controlling variables for iterative solver
    integer(i4) :: max_iter  !< input: maximal number of iterations done
    real(dp) :: delta_iter   !< input: convergence criteria for iterative solver
    real(dp) :: step_iter    !< input: step-size for linear search
    logical :: first_iter    !< whether it is at the first iteration (to determine search direction)
    logical :: up_iter       !< whether the search direction is upwards
    logical :: bisect_iter   !< whether to do the bisection search part (after the interval is found)
    real(dp) :: up_bnd_iter  !< upper bound for the current bisection step
    real(dp) :: low_bnd_iter !< lower bound for the current bisection step
    ! accumulated later fluxes (in current time-step)
    real(dp), dimension(:), allocatable :: L1_runoff_E !< runoff energy at L1 level
    real(dp), dimension(:), allocatable :: L1_acc_ssrd !< accumulated shortwave radiation at L1 level
    real(dp), dimension(:), allocatable :: L1_acc_strd !< accumulated longwave radiation at L1 level
    real(dp), dimension(:), allocatable :: L1_acc_temp !< accumulated air temperature at L1 level
    real(dp), dimension(:), allocatable :: L1_ssrd_calc !< current shortwave radiation at L1 level
    real(dp), dimension(:), allocatable :: L1_strd_calc !< current longwave radiation at L1 level
    real(dp), dimension(:), allocatable :: L1_tann_calc !< current mean air temperature at L1 level
    integer(i4) :: ts_cnt !< sub time-step counter for accumulation of meteo
    ! vars for routing
    integer(i4) :: s11 !< starting index for current L11 domain
    integer(i4) :: e11 !< ending index for current L11 domain
    real(dp), dimension(:,:), allocatable :: netNode_E_IN !< Total energy inputs at t-1 and t
    real(dp), dimension(:,:), allocatable :: netNode_E_R  !< energy leaving at t-1 and t
    real(dp), dimension(:), allocatable :: netNode_E_mod  !< Simulated routed energy
    real(dp), dimension(:), allocatable :: netNode_E_out  !< total energy source from cell in L11
    real(dp), dimension(:), allocatable :: L11_srad_net   !< net short wave radiation at L11
    real(dp), dimension(:), allocatable :: L11_lrad_in    !< incoming long wave radiation at L11
    real(dp), dimension(:), allocatable :: L11_air_temp   !< air temp at L11
    ! variable containing river temp for each domain at L11 level
    real(dp), dimension(:), allocatable :: river_temp !< resulting river temp at L11 in [deg C]
  contains
    ! config and inits
    !> \copydoc mo_mrm_riv_temp_class::config
    procedure :: config !< \see mo_mrm_riv_temp_class::config
    !> \copydoc mo_mrm_riv_temp_class::init
    procedure :: init !< \see mo_mrm_riv_temp_class::init
    !> \copydoc mo_mrm_riv_temp_class::init_area
    procedure :: init_area !< \see mo_mrm_riv_temp_class::init_area
    !> \copydoc mo_mrm_riv_temp_class::init_riv_temp
    procedure :: init_riv_temp !< \see mo_mrm_riv_temp_class::init_riv_temp

    ! source accumulations
    !> \copydoc mo_mrm_riv_temp_class::acc_source_e
    procedure :: acc_source_E !< \see mo_mrm_riv_temp_class::acc_source_e
    !> \copydoc mo_mrm_riv_temp_class::finalize_source_e
    procedure :: finalize_source_E !< \see mo_mrm_riv_temp_class::finalize_source_e

    ! temp-energy routing routines
    !> \copydoc mo_mrm_riv_temp_class::get_lrad_out
    procedure :: get_lrad_out !< \see mo_mrm_riv_temp_class::get_lrad_out
    !> \copydoc mo_mrm_riv_temp_class::get_lat_heat
    procedure :: get_lat_heat !< \see mo_mrm_riv_temp_class::get_lat_heat
    !> \copydoc mo_mrm_riv_temp_class::get_sens_heat
    procedure :: get_sens_heat !< \see mo_mrm_riv_temp_class::get_sens_heat
    !> \copydoc mo_mrm_riv_temp_class::get_e_io
    procedure :: get_E_IO !< \see mo_mrm_riv_temp_class::get_e_io
    !> \copydoc mo_mrm_riv_temp_class::l11_routing_e
    procedure :: L11_routing_E !< \see mo_mrm_riv_temp_class::l11_routing_e

    ! helper for iterative solver
    !> \copydoc mo_mrm_riv_temp_class::init_iter
    procedure :: init_iter !< \see mo_mrm_riv_temp_class::init_iter
    !> \copydoc mo_mrm_riv_temp_class::next_iter
    procedure :: next_iter !< \see mo_mrm_riv_temp_class::next_iter

    ! care taker
    !> \copydoc mo_mrm_riv_temp_class::reset_timestep
    procedure :: reset_timestep !< \see mo_mrm_riv_temp_class::reset_timestep
    !> \copydoc mo_mrm_riv_temp_class::alloc_lateral
    procedure :: alloc_lateral !< \see mo_mrm_riv_temp_class::alloc_lateral
    !> \copydoc mo_mrm_riv_temp_class::dealloc_lateral
    procedure :: dealloc_lateral !< \see mo_mrm_riv_temp_class::dealloc_lateral
    !> \copydoc mo_mrm_riv_temp_class::clean_up
    procedure :: clean_up !< \see mo_mrm_riv_temp_class::clean_up

  end type riv_temp_type

contains


  !> \brief clean up
  subroutine clean_up( &
    self &
  )
    implicit none

    class(riv_temp_type), intent(inout) :: self

    if ( allocated(self%L1_runoff_E) ) deallocate(self%L1_runoff_E)
    if ( allocated(self%L1_acc_strd) ) deallocate(self%L1_acc_strd)
    if ( allocated(self%L1_acc_ssrd) ) deallocate(self%L1_acc_ssrd)
    if ( allocated(self%L1_acc_temp) ) deallocate(self%L1_acc_temp)
    if ( allocated(self%dir_riv_widths) ) deallocate(self%dir_riv_widths)
    if ( allocated(self%L11_riv_widths) ) deallocate(self%L11_riv_widths)
    if ( allocated(self%L11_riv_areas) ) deallocate(self%L11_riv_areas)
    if ( allocated(self%netNode_E_IN) ) deallocate(self%netNode_E_IN)
    if ( allocated(self%netNode_E_R) ) deallocate(self%netNode_E_R)
    if ( allocated(self%netNode_E_mod) ) deallocate(self%netNode_E_mod)
    if ( allocated(self%netNode_E_out) ) deallocate(self%netNode_E_out)
    if ( allocated(self%L11_srad_net) ) deallocate(self%L11_srad_net)
    if ( allocated(self%L11_lrad_in) ) deallocate(self%L11_lrad_in)
    if ( allocated(self%L11_air_temp) ) deallocate(self%L11_air_temp)
    if ( allocated(self%river_temp) ) deallocate(self%river_temp)
    ! meteo arrays
    if ( allocated(self%L1_ssrd_calc) ) deallocate(self%L1_ssrd_calc)
    if ( allocated(self%L1_strd_calc) ) deallocate(self%L1_strd_calc)
    if ( allocated(self%L1_tann_calc) ) deallocate(self%L1_tann_calc)

  end subroutine clean_up

  !> \brief configure the \ref riv_temp_type class from the mhm namelist
  subroutine config( &
    self, &
    file_namelist, &
    unamelist, &
    file_namelist_param, &
    unamelist_param &
  )

    use mo_common_constants, only : maxNoDomains, nodata_i4
    use mo_common_variables, only : domainMeta
    use mo_nml, only : close_nml, open_nml, position_nml
    use mo_check, only : check_dir
    USE mo_string_utils, ONLY : num2str

    implicit none

    class(riv_temp_type), intent(inout) :: self
    character(*), intent(in) :: file_namelist !< mhm namelist file
    character(*), intent(in) :: file_namelist_param !< mhm parameter namelist file
    integer, intent(in) :: unamelist
    integer, intent(in) :: unamelist_param

    ! parameter to read in
    real(dp) :: albedo_water ! albedo of open water
    real(dp) :: pt_a_water ! priestley taylor alpha parameter for PET on open water
    real(dp) :: emissivity_water ! emissivity of water
    real(dp) :: turb_heat_ex_coeff ! lateral heat exchange coefficient water <-> air
    ! controlling variables for iterative solver
    integer(i4) :: max_iter
    real(dp) :: delta_iter
    real(dp) :: step_iter
    ! files for river widths
    character(256), dimension(maxNoDomains) :: dir_riv_widths
    character(256) :: riv_widths_file ! file name for river widths
    character(256) :: riv_widths_name ! variable name for river widths

    integer(i4) :: iDomain, domainID

    ! namelist for river temperature configuration
    namelist /config_riv_temp/ &
      albedo_water, &
      pt_a_water, &
      emissivity_water, &
      turb_heat_ex_coeff, &
      max_iter, &
      delta_iter, &
      step_iter, &
      riv_widths_file, &
      riv_widths_name, &
      dir_riv_widths

    ! allocate the directory arrays
    allocate(self%dir_riv_widths(domainMeta%nDomains))

    ! default values
    albedo_water = 0.15_dp
    pt_a_water = 1.26_dp
    emissivity_water = 0.96_dp
    turb_heat_ex_coeff = 20.0_dp
    max_iter = 20_i4
    delta_iter = 1.0e-02_dp
    step_iter = 5.0_dp

    ! open the namelist file
    call open_nml(file_namelist, unamelist, quiet=.true.)
    ! find the river-temp config namelist
    call position_nml(self%nml_name, unamelist)
    ! read the river-temp config namelist
    read(unamelist, nml=config_riv_temp)

    self%albedo_water = albedo_water
    self%pt_a_water = pt_a_water
    self%emissivity_water = emissivity_water
    self%turb_heat_ex_coeff = turb_heat_ex_coeff
    self%delta_iter = delta_iter
    self%max_iter = max_iter
    self%step_iter = step_iter
    self%riv_widths_file = riv_widths_file
    self%riv_widths_name = riv_widths_name
    self%dir_riv_widths = dir_riv_widths

    ! closing the namelist file
    call close_nml(unamelist)

    do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)
      call check_dir( &
        path=self%dir_riv_widths(iDomain), &
        text="(domain "//trim(num2str(domainID,'(I3)'))//") River widths directory:", &
        raise=.true., &
        tab=4, &
        text_length=40 &
      )
    end do

  end subroutine config

  !> \brief initalize the \ref riv_temp_type class for the current domain
  subroutine init( &
    self, &
    nCells &
  )
    use mo_append, only : append
    use mo_mrm_constants, only : nRoutingStates
    ! use mo_common_variables, only : level0, domainMeta

    implicit none

    class(riv_temp_type), intent(inout) :: self
    integer(i4), intent(in) :: nCells !< number of level-11 cells for the current domain

    real(dp), dimension(:), allocatable :: dummy_Vector11
    real(dp), dimension(:, :), allocatable :: dummy_Matrix11_IT

    ! dummy vector and matrix
    allocate(dummy_Vector11   (nCells))
    allocate(dummy_Matrix11_IT(nCells, nRoutingStates))
    dummy_Vector11(:) = 0.0_dp
    dummy_Matrix11_IT(:, :) = 0.0_dp

    ! simulated energy flux at each node
    call append(self%netNode_E_mod, dummy_Vector11)
    ! simulated river temperature at each node
    call append(self%river_temp, dummy_Vector11)
    ! Total outflow from cells L11 at time tt
    call append(self%netNode_E_out, dummy_Vector11)
    ! Total discharge inputs at t-1 and t
    call append(self%netNode_E_IN, dummy_Matrix11_IT)
    !  Routed outflow leaving a node
    call append(self%netNode_E_R, dummy_Matrix11_IT)

    ! free space
    deallocate(dummy_Vector11)
    deallocate(dummy_Matrix11_IT)

  end subroutine init

  !> \brief initialize the river area of \ref riv_temp_type class for the current domain
  subroutine init_area( &
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
    use mo_read_nc, only: read_const_nc

    implicit none

    class(riv_temp_type), intent(inout) :: self
    integer(i4), intent(in) :: iDomain !< Domain ID
    integer(i4), dimension(:), intent(in) :: L11_netPerm !< L11 routing order
    integer(i4), dimension(:), intent(in) :: L11_fromN !< L11 source grid cell order
    real(dp), dimension(:), intent(in) :: L11_length !< L11 link length
    integer(i4), intent(in) :: nLinks !< number of L11 links in the current domain
    integer(i4), intent(in) :: nCells !< number of L11 cells of the current domain
    integer(i4), intent(in) :: ncols  !< Number of columns
    integer(i4), intent(in) :: nrows  !< Number of rows
    logical, dimension(:, :), intent(in) :: L11_mask !< the mask for valid cells in the original grid (nrows, ncols)

    real(dp), dimension(:,:), allocatable :: L11_data ! read data from file
    real(dp), dimension(:), allocatable :: L11_riv_widths, L11_riv_areas

    integer(i4) :: i, k, iNode

    call read_const_nc(&
      trim(self%dir_riv_widths(iDomain)), &
      nrows, &
      ncols, &
      self%riv_widths_name, &
      L11_data, &
      self%riv_widths_file &
    )

    allocate(L11_riv_widths(nCells))
    L11_riv_widths(:) = pack(L11_data(:,:), mask=L11_mask)
    call append(self%L11_riv_widths, L11_riv_widths)

    allocate(L11_riv_areas(nCells))
    ! at area at Outlet is 0
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

  end subroutine init_area

  !> \brief initialize the river temperature of \ref riv_temp_type class for the current domain
  subroutine init_riv_temp( &
    self, &
    temp_air, &
    efecarea, &
    L1_L11_Id, &
    L11_areacell, &
    L11_L1_Id, &
    map_flag &
  )

    use mo_mrm_pre_routing, only : L11_meteo_acc

    implicit none

    class(riv_temp_type), intent(inout) :: self
    real(dp), dimension(:), intent(in) :: temp_air !< air temperature [degC] for current timestep
    real(dp), intent(in), dimension(:) :: efecarea !< effective area in [km2] at Level 1
    integer(i4), intent(in), dimension(:) :: L1_L11_Id !< L11 Ids mapped on L1
    real(dp), intent(in), dimension(:) :: L11_areacell !< effective area in [km2] at Level 11
    integer(i4), intent(in), dimension(:) :: L11_L1_Id !< L1 Ids mapped on L11
    logical, intent(in) :: map_flag !< Flag indicating whether routing resolution is higher than hydrologic one

    ! map temperature from L1 to L11
    call L11_meteo_acc( &
      temp_air, efecarea, L1_L11_Id, L11_areacell, L11_L1_Id, map_flag, self%river_temp(self%s11 : self%e11))

    ! assure positive temperature
    self%river_temp(self%s11 : self%e11) = max(self%delta_T, self%river_temp(self%s11 : self%e11))

  end subroutine init_riv_temp

  !> \brief reset \ref riv_temp_type class for next timestep
  subroutine reset_timestep(self)
    implicit none

    class(riv_temp_type), intent(inout) :: self

    ! set sub time-step counter to 0
    self%ts_cnt = 0_i4

    self%L1_runoff_E = 0.0_dp
    self%L1_acc_strd = 0.0_dp
    self%L1_acc_ssrd = 0.0_dp
    self%L1_acc_temp = 0.0_dp

  end subroutine reset_timestep

  !> \brief allocate lateral temp components of \ref riv_temp_type class for current domain
  subroutine alloc_lateral( &
    self, &
    nCells &
  )
    implicit none

    class(riv_temp_type), intent(inout) :: self
    integer(i4), intent(in) :: nCells !< number of level-1 cells for the current domain

    allocate(self%L1_runoff_E(nCells))
    allocate(self%L1_acc_strd(nCells))
    allocate(self%L1_acc_ssrd(nCells))
    allocate(self%L1_acc_temp(nCells))
    ! meteo arrays
    allocate(self%L1_ssrd_calc(nCells))
    allocate(self%L1_strd_calc(nCells))
    allocate(self%L1_tann_calc(nCells))
    ! init these arrays to 0
    call self%reset_timestep()

  end subroutine alloc_lateral

  !> \brief deallocate lateral temp components of \ref riv_temp_type
  subroutine dealloc_lateral( &
    self &
  )
    implicit none

    class(riv_temp_type), intent(inout) :: self

    deallocate(self%L1_runoff_E)
    deallocate(self%L1_acc_strd)
    deallocate(self%L1_acc_ssrd)
    deallocate(self%L1_acc_temp)
    deallocate(self%L1_ssrd_calc)
    deallocate(self%L1_strd_calc)
    deallocate(self%L1_tann_calc)

  end subroutine dealloc_lateral

  !> \brief accumulate energy sources of \ref riv_temp_type
  subroutine acc_source_e( &
    self, &
    fSealed_area_fraction, &
    fast_interflow, &
    slow_interflow, &
    baseflow, &
    direct_runoff, &
    temp_air &
  )

    use mo_mrm_pre_routing, only : calc_L1_runoff_E

    implicit none

    class(riv_temp_type), intent(inout) :: self
    real(dp), dimension(:), intent(in) :: fSealed_area_fraction !< sealed area fraction [1]
    real(dp), dimension(:), intent(in) :: fast_interflow !< \f$ q_0 \f$ Fast runoff component [mm TS-1]
    real(dp), dimension(:), intent(in) :: slow_interflow !< \f$ q_1 \f$ Slow runoff component [mm TS-1]
    real(dp), dimension(:), intent(in) :: baseflow !< \f$ q_2 \f$ Baseflow [mm TS-1]
    real(dp), dimension(:), intent(in) :: direct_runoff !< \f$ q_D \f$ Direct runoff from impervious areas  [mm TS-1]
    real(dp), dimension(:), intent(in) :: temp_air !< air temperature [K]

    ! increase the sub time-step counter
    self%ts_cnt = self%ts_cnt + 1_i4

    ! caclucate the temperature energy of the runoffs at L1 in [K mm]
    ! automatically accumulate them
    call calc_L1_runoff_E( &
      fSealed_area_fraction, &
      fast_interflow, slow_interflow, baseflow, direct_runoff, &
      temp_air, self%L1_tann_calc, &
      self%L1_runoff_E & ! will be added here
    )
    ! accumulate meteo forcings (will be averaged with sub time-step counter later)
    self%L1_acc_ssrd = self%L1_acc_ssrd + self%L1_ssrd_calc
    self%L1_acc_strd = self%L1_acc_strd + self%L1_strd_calc
    self%L1_acc_temp = self%L1_acc_temp + temp_air

  end subroutine acc_source_e

  !> \brief finalize energy sources of \ref riv_temp_type
  subroutine finalize_source_E( &
    self, &
    efecarea, &
    L1_L11_Id, &
    L11_areacell, &
    L11_L1_Id, &
    timestep, &
    map_flag &
  )

    use mo_mrm_pre_routing, only : L11_meteo_acc, L11_runoff_acc

    implicit none

    class(riv_temp_type), intent(inout) :: self
    !> effective area in [km2] at Level 1
    real(dp), intent(in), dimension(:) :: efecarea
    !> L11 Ids mapped on L1
    integer(i4), intent(in), dimension(:) :: L1_L11_Id
    !> effective area in [km2] at Level 11
    real(dp), intent(in), dimension(:) :: L11_areacell
    !> L1 Ids mapped on L11
    integer(i4), intent(in), dimension(:) :: L11_L1_Id
    !> simulation timestep in [h]
    integer(i4), intent(in) :: timestep
    !> Flag indicating whether routing resolution is higher than hydrologic one
    logical, intent(in) :: map_flag

    ! prepare temporal variables for temp-routing
    if ( allocated(self%L11_srad_net) ) deallocate(self%L11_srad_net)
    if ( allocated(self%L11_lrad_in) ) deallocate(self%L11_lrad_in)
    if ( allocated(self%L11_air_temp) ) deallocate(self%L11_air_temp)
    allocate(self%L11_srad_net(size(L11_L1_Id, dim = 1)))
    allocate(self%L11_lrad_in(size(L11_L1_Id, dim = 1)))
    allocate(self%L11_air_temp(size(L11_L1_Id, dim = 1)))

    ! convert the L1 runoff temperature energy from [K mm] to [K m3 s-1] on L11
    call L11_runoff_acc( &
      self%L1_runoff_E, &
      efecarea, &
      L1_L11_Id, &
      L11_areacell, &
      L11_L1_Id, &
      timestep, &
      map_flag, &
      self%netNode_E_out(self%s11 : self%e11) &
    )
    ! all meteo forcings are accumulated over sub time-steps and need to be divided by the sub time-step counter
    ! short wave radiation on L11
    self%L1_acc_ssrd = self%L1_acc_ssrd / real(self%ts_cnt, dp)
    call L11_meteo_acc(self%L1_acc_ssrd, efecarea, L1_L11_Id, L11_areacell, L11_L1_Id, map_flag, self%L11_srad_net)
    ! calculate the net-radiation from incoming/outging short radiation
    self%L11_srad_net = self%L11_srad_net * (1._dp - self%albedo_water)
    ! long wave radiation on L11
    self%L1_acc_strd = self%L1_acc_strd / real(self%ts_cnt, dp)
    call L11_meteo_acc(self%L1_acc_strd, efecarea, L1_L11_Id, L11_areacell, L11_L1_Id, map_flag, self%L11_lrad_in)
    ! air temperature on L11
    self%L1_acc_temp = self%L1_acc_temp / real(self%ts_cnt, dp)
    call L11_meteo_acc(self%L1_acc_temp, efecarea, L1_L11_Id, L11_areacell, L11_L1_Id, map_flag, self%L11_air_temp)

  end subroutine finalize_source_E

  !> \brief get outgoing longwave radiation of \ref riv_temp_type
  !> \return outgoing longwave radiation
  real(dp) function get_lrad_out(self, riv_temp) result(lrad_out)

    use mo_constants, only: sigma_dp

    implicit none

    class(riv_temp_type), intent(in) :: self
    !> river temperature in K
    real(dp), intent(in) :: riv_temp

    ! outgoing longwave radiation from Boltzmann equation
    lrad_out = self%emissivity_water * sigma_dp * riv_temp ** 4_i4

  end function get_lrad_out

  !> \brief latent heat flux of \ref riv_temp_type
  !> \return latent heat flux
  real(dp) function get_lat_heat(self, air_temp, netrad) result(lat_heat)

    use mo_constants, only : Psychro_dp
    use mo_pet, only : slope_satpressure

    implicit none

    class(riv_temp_type), intent(in) :: self
    !> air temperature in deg C
    real(dp), intent(in) :: air_temp
    !> net radiation in W * m-2
    real(dp), intent(in) :: netrad

    ! save slope of saturation vapor pressure curve
    real(dp) :: delta

    delta = slope_satpressure(air_temp) ! slope of saturation vapor pressure curve
    ! latent heat for priestley taylor PET formula
    lat_heat = self%pt_a_water * delta / (Psychro_dp + delta) * max(netrad, 0._dp)

  end function get_lat_heat

  !> \brief sensible heat flux of \ref riv_temp_type
  !> \return sensible heat flux
  real(dp) function get_sens_heat(self, air_temp, riv_temp) result(sens_heat)

    implicit none

    class(riv_temp_type), intent(in) :: self
    real(dp), intent(in) :: air_temp !< air temperature in [deg C]
    real(dp), intent(in) :: riv_temp !< river temperature in [deg C]

    ! sensible heat flux resulting for temp. diff.: river <-> air
    sens_heat = self%turb_heat_ex_coeff * (riv_temp - air_temp)

  end function get_sens_heat

  !> \brief get complete energy source of \ref riv_temp_type at given cell
  !> \return energy IO
  real(dp) function get_E_IO(self, riv_temp, cell) result(E_IO)

    use mo_constants, only : T0_dp, cp_w_dp
    use mo_mhm_constants, only : H2Odens

    implicit none

    class(riv_temp_type), intent(in) :: self
    !> given river temperature in K to calculate heat fluxes
    real(dp), intent(in) :: riv_temp
    !> cell index in the current domain
    integer(i4), intent(in) :: cell
    !> net radiation calc from short and longwave radiation in/out
    real(dp) :: netrad, sens_heat, lat_heat

    ! net radiation
    netrad = self%L11_srad_net(cell) + self%L11_lrad_in(cell) - self%get_lrad_out(riv_temp)
    ! sensible heat flux
    sens_heat = self%get_sens_heat(self%L11_air_temp(cell), riv_temp - T0_dp)
    ! latent heat flux
    lat_heat = self%get_lat_heat(self%L11_air_temp(cell), netrad)

    ! convert energy flux [W m-2] to [K m s-1]
    ! needs to be multiplied with river-area to result in temp-energy flux [K m3 s-1]
    E_IO = (netrad - sens_heat - lat_heat) / (H2Odens * cp_w_dp)

  end function get_E_IO

  !> \brief execute the temperature routing of \ref riv_temp_type
  subroutine L11_routing_E( &
    self, &
    nLinks, &
    netPerm, &
    netLink_fromN, &
    netLink_toN, &
    netLink_C1, &
    netLink_C2, &
    nInflowGauges, &
    InflowHeadwater, &
    InflowNodeList, &
    L11_qTR, &
    L11_Qmod &
  )
    use mo_constants, only : T0_dp

    implicit none

    class(riv_temp_type), intent(inout) :: self
    !> number of stream segment (reaches)
    integer(i4), intent(in) :: nLinks
    !> routing order of a given domain (permutation)
    integer(i4), dimension(:), intent(in) :: netPerm
    !> from node
    integer(i4), dimension(:), intent(in) :: netLink_fromN
    !> to node
    integer(i4), dimension(:), intent(in) :: netLink_toN
    !> routing parameter  C1 (\cite CMM1988 p. 25-41)
    real(dp), dimension(:), intent(in) :: netLink_C1
    !> routing parameters C2 (id)
    real(dp), dimension(:), intent(in) :: netLink_C2
    !> [-]      number of inflow points
    integer(i4), intent(in) :: nInflowGauges
    !> [-]      if to consider headwater cells of inflow gauge
    logical, dimension(:), intent(in) :: InflowHeadwater
    !> [-]      L11 ID of inflow points
    integer(i4), dimension(:), intent(in) :: InflowNodeList
    !> [m3 s-1] Transformed outflow leaving node I at current timestep(Muskingum)
    real(dp), intent(in), dimension(:) :: L11_qTR
    !> [m3 s-1] Simulated routed discharge
    real(dp), intent(in), dimension(:) :: L11_Qmod

    integer(i4) :: i, k, m, iNode, tNode, L11in, L11to
    real(dp) :: E_IO_in, E_IO, T_est, T_rout

    ! initialize total input at point time (t+1) in all nodes
    self%netNode_E_IN(self%s11 : self%e11, 2) = 0.0_dp

    riv_loop: do k = 1, nLinks
      ! get LINK routing order -> i
      i = netPerm(k)
      iNode = netLink_fromN(i)
      tNode = netLink_toN(i)
      ! indices for concaternated domains
      L11in = iNode + self%s11 - 1_i4
      L11to = tNode + self%s11 - 1_i4

      ! accumulate all inputs in iNode
      self%netNode_E_IN(L11in, 2) = self%netNode_E_IN(L11in, 2) + self%netNode_E_out(L11in)
      ! calculate the river temp at the upstream node
      self%river_temp(L11in) = self%netNode_E_IN(L11in, 2) / L11_Qmod(iNode) - T0_dp

      ! first guess for the routed temperature is the upstream temp
      T_est = self%river_temp(L11in) + T0_dp
      ! calculate the meteo energy source at the upstream node
      E_IO_in = self%get_E_IO(T_est, iNode)

      ! initialize iteration
      call self%init_iter()
      ! perform iteration
      iterloop: do m=1, self%max_iter
        ! meteo energy flux depending on the estimated routed temperature (avg. IN- and OUT-node)
        E_IO = (E_IO_in + self%get_E_IO(T_est, tNode)) * 0.5_dp * self%L11_riv_areas(L11in)
        ! routing iNode
        self%netNode_E_R(L11in, 2) = self%netNode_E_R(L11in, 1) &
          + netLink_C1(i) * (self%netNode_E_IN(L11in, 1) - self%netNode_E_R(L11in, 1)) &
          + netLink_C2(i) * (self%netNode_E_IN(L11in, 2) + E_IO - self%netNode_E_IN(L11in, 1))
        ! calculate the routed temperature
        T_rout = self%netNode_E_R(L11in, 2) / L11_qTR(iNode)
        ! check for convergence
        if ( abs(T_rout - T_est) < self%delta_iter ) exit iterloop
        ! get new estimate for next iteration
        call self%next_iter(T_est, T_rout)
      end do iterloop

      ! add the meteo energy flux to the accumulated incoming energy at the IN-node
      self%netNode_E_IN(L11in, 2) = self%netNode_E_IN(L11in, 2) + E_IO

      ! check if the inflow from upstream cells should be deactivated
      if (nInflowGauges .GT. 0) then
        do i = 1, nInflowGauges
          ! check if downstream Node (tNode) is inflow gauge and headwaters should be ignored
          if ((tNode == InflowNodeList(i)) .AND. (.NOT. InflowHeadwater(i))) &
            self%netNode_E_R(L11in, 2) = 0.0_dp
        end do
      end if

      ! add routed temp energy to downstream node
      self%netNode_E_IN(L11to, 2) = self%netNode_E_IN(L11to, 2) + self%netNode_E_R(L11in, 2)

    end do riv_loop

    ! Accumulate all inputs in tNode (netNode_E_out) ONLY for last link
    tNode = netLink_toN(netPerm(nLinks))
    L11to = tNode + self%s11 - 1_i4
    self%netNode_E_IN(L11to, 2) = self%netNode_E_IN(L11to, 2) + self%netNode_E_out(L11to)
    ! calculate riv temp at last link
    self%river_temp(L11to) = self%netNode_E_IN(L11to, 2) / L11_Qmod(tNode) - T0_dp

    ! backflow t-> t-1
    self%netNode_E_R(self%s11 : self%e11, 1) = self%netNode_E_R(self%s11 : self%e11, 2)
    self%netNode_E_IN(self%s11 : self%e11, 1) = self%netNode_E_IN(self%s11 : self%e11, 2)

  end subroutine L11_routing_E

  !> \brief initialize iterative solver of \ref riv_temp_type
  subroutine init_iter(self)

    implicit none

    class(riv_temp_type), intent(inout) :: self

    ! first iteration is used to determine the search direction for the interval
    self%first_iter = .true.
    ! after the interval is found, we will perform a biseciton to nail down the temperature
    self%bisect_iter = .false.

  end subroutine init_iter

  !> \brief execute next iteration with iterative solver of \ref riv_temp_type
  subroutine next_iter(self, T_est, T_rout)

    implicit none

    class(riv_temp_type), intent(inout) :: self
    real(dp), intent(inout) :: T_est !< estimated river temperature
    real(dp), intent(in) :: T_rout !< calculated (routed) river temperature

    ! before performing a bisection we need to search for the interval (with given step-size)
    if ( .not. self%bisect_iter ) then
      ! to determine the direction of interval search, we have to wait for the first iteration
      if ( self%first_iter ) then
        self%first_iter = .false.
        ! determine search direction
        self%up_iter = ( T_rout > T_est )
      end if
      ! check if we need to take another step
      if ( self%up_iter .eqv. ( T_rout > T_est ) ) then
        if ( self%up_iter ) then
          T_est = T_est + self%step_iter
        else
          T_est = T_est - self%step_iter
        end if
      ! if direction changes, we have found the interval for bisection
      else
        ! start bisection in next iteration
        self%bisect_iter = .true.
        ! set interval bounds for bisection
        if ( self%up_iter ) then
          self%up_bnd_iter = T_est
          self%low_bnd_iter = T_est - self%step_iter
        else
          self%up_bnd_iter = T_est + self%step_iter
          self%low_bnd_iter = T_est
        end if
        ! set estimation to interval center
        T_est = (self%up_bnd_iter + self%low_bnd_iter) / 2.0_dp
      end if
    ! perform bisection
    else
      ! if routed temp is lower than estimated, use lower interval
      if ( T_rout < T_est ) then
        self%up_bnd_iter = T_est
      ! if routed temp is higher than estimated, use upper interval
      else
        self%low_bnd_iter = T_est
      end if
      ! set estimation to interval center
      T_est = (self%up_bnd_iter + self%low_bnd_iter) / 2.0_dp
    end if

  end subroutine next_iter

end module mo_mrm_riv_temp_class
