!>       \file mo_mrm_riv_temp_class.f90

!>       \brief Class for the river temperature calculations
!>       \details

!>       \authors Sebastian Mueller

!>       \date Apr 2020

module mo_mrm_riv_temp_class

  use mo_kind, only: dp

  implicit none

  private
  public :: riv_temp_type

  ! This is a container to define the river temperature routing in the current time step
  type riv_temp_type
    ! riv geometry
    character(256), dimension(:), allocatable :: dirWidths ! Directory where river widths are stored
    character(256) :: riv_widths_file ! file name for river widths
    character(256) :: riv_widths_name ! variable name for river widths
    real(dp), public, dimension(:), allocatable :: L11_riv_widths !river widths in L11
    real(dp), public, dimension(:), allocatable :: L11_riv_area !river area in L11
    ! PET related vars
    real(dp) :: albedo_water ! albedo of open water
    real(dp) :: pt_a_water ! priestley taylor alpha parameter for PET on open water
    ! \f$ E_L \f$ Generated lateral temperature energy flux [m3 s-1 K] on L1
    real(dp), dimension(:), allocatable :: L1_lateral_E
    ! vars for temp routing
    real(dp), dimension(:), allocatable :: L11_lateral_E
    real(dp), dimension(:,:), allocatable :: netNode_E_IN ! Total energy inputs at t-1 and t
    real(dp), dimension(:,:), allocatable :: netNode_E_R ! energy leaving at t-1 and t
    real(dp), dimension(:), allocatable :: netNode_E_mod ! Simulated routed energy
    real(dp), dimension(:), allocatable :: netNode_E_out ! total energy source from cell in L11
    ! results
    real(dp), dimension(:, :), allocatable :: mRM_river_temp ! variable containing river temp for each domain and gauge
  contains
    procedure :: config => meth_config
    procedure :: init => meth_init
    procedure :: init_area => meth_init_area
    procedure :: init_cond_from_L1 => meth_init_cond_from_L1
    procedure :: calc_lateral_E => meth_calc_lateral_E
  end type riv_temp_type

contains

  subroutine meth_config( &
    self &
  )

    implicit none

    class(riv_temp_type), intent(inout) :: self

    print *, 'config river temp'
    ! TODO-RIV-TEMP:
    ! - set variables from nml

  end subroutine meth_config


  subroutine meth_init( &
    self &
  )
    use mo_append, only : append
    use mo_kind, only : dp, i4
    use mo_mrm_constants, only : nRoutingStates
    use mo_common_variables, only : level0, domainMeta

    implicit none

    class(riv_temp_type), intent(inout) :: self

    real(dp), dimension(:), allocatable :: dummy_Vector11
    real(dp), dimension(:, :), allocatable :: dummy_Matrix11_IT

    print *, 'init river temp'
    ! TODO-RIV-TEMP:
    ! - riv-width
    ! - init_condition riv_temp (air temp)
    ! - see: variables_alloc_routing

    ! dummy vector and matrix
    ! allocate(dummy_Vector11   (level11(iDomain)%nCells))
    ! allocate(dummy_Matrix11_IT(level11(iDomain)%nCells, nRoutingStates))
    ! dummy_Vector11(:) = 0.0_dp

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


  subroutine meth_init_cond_from_L1( &
    self &
  )
    implicit none

    class(riv_temp_type), intent(inout) :: self

    print *, 'init river temperature with L1 air temp'
    ! TODO-RIV-TEMP:
    ! - dis-agg L1 air temp for tt=1

  end subroutine meth_init_cond_from_L1

  subroutine meth_calc_lateral_E( &
    self, &
    fSealed_area_fraction, &
    fast_interflow, &
    slow_interflow, &
    baseflow, &
    direct_runoff, &
    temp_air, &
    mean_temp_air, &
    add &
  )

    implicit none

    class(riv_temp_type), intent(inout) :: self
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
    ! switch to turn on adding energy if routing TS is larger then mhm TS
    logical, intent(in), optional :: add

    ! internal switch for adding energy
    logical :: do_add
    ! lateral energy
    real(dp), dimension(size(baseflow)) :: lateral_E

    ! check optional argument add
    if (present(add)) then
        do_add = add
    else
        do_add = .false.
    end if

    lateral_E = ( &
        (baseflow * mean_temp_air + (slow_interflow + fast_interflow) * temp_air) * (1.0_dp - fSealed_area_fraction) &
        + direct_runoff * temp_air * fSealed_area_fraction &
    )

    if ( do_add ) then
        self%L1_lateral_E = self%L1_lateral_E + lateral_E
    else
        self%L1_lateral_E = lateral_E
    end if

  end subroutine meth_calc_lateral_E

end module mo_mrm_riv_temp_class

