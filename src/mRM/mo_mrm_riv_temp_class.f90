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

  type riv_temp_type
    ! This is a container to define the river temperature routing in the current time step
    character(256), dimension(:), allocatable :: dirWidths ! Directory where river widths are stored
    character(256) :: riv_widths_file ! file name for river widths
    character(256) :: riv_widths_name ! variable name for river widths
    real(dp), public, dimension(:), allocatable :: L11_riv_widths !river widths in L11
    real(dp) :: albedo_water ! albedo of open water
    real(dp) :: pt_a_water ! priestley taylor alpha parameter for PET on open water
    ! \f$ E_L \f$ Generated lateral temperature energy flux [m3 s-1 K]
    real(dp), dimension(:), allocatable :: L1_lateral_E
    real(dp), dimension(:, :), allocatable :: mRM_river_temp ! variable containing river temp for each domain and gauge
  contains
    procedure :: config => config_riv_temp
    procedure :: init => init_riv_temp
    procedure :: calc_lateral_E => pro_lateral_E
  end type riv_temp_type

contains

  subroutine config_riv_temp( &
    self &
  )

    implicit none

    class(riv_temp_type), intent(inout) :: self

    print *, 'config river temp'
    ! TODO-RIV-TEMP:
    ! - set variables from nml

  end subroutine config_riv_temp


  subroutine init_riv_temp ( &
    self &
  )

    implicit none

    class(riv_temp_type), intent(inout) :: self

    print *, 'init river temp'
    ! TODO-RIV-TEMP:
    ! - riv-width
    ! - init_condition riv_temp (air temp)
    ! -

  end subroutine init_riv_temp


  subroutine pro_lateral_E( &
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

  end subroutine pro_lateral_E

end module mo_mrm_riv_temp_class

