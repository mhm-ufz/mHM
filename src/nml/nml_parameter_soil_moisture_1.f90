!> \file nml_parameter_soil_moisture_1.f90
!> \copydoc nml_soil_moisture_1

!> \brief Soil moisture - Case 1
!> \details Feddes ET reduction with multi-layer infiltration and Brooks-Corey parameterization.
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_soil_moisture_1
  use nml_helper, only: &
    nml_file_t, &
    nml_line_buffer, &
    NML_OK, &
    NML_ERR_FILE_NOT_FOUND, &
    NML_ERR_OPEN, &
    NML_ERR_NOT_OPEN, &
    NML_ERR_NML_NOT_FOUND, &
    NML_ERR_READ, &
    NML_ERR_CLOSE, &
    NML_ERR_REQUIRED, &
    NML_ERR_ENUM, &
    NML_ERR_BOUNDS, &
    NML_ERR_NOT_SET, &
    NML_ERR_INVALID_NAME, &
    NML_ERR_INVALID_INDEX, &
    idx_check, &
    to_lower
  use ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
  ! kind specifiers listed in the nml-tools configuration file
  use mo_kind, only: &
    dp
  use mo_parameter_types, only: parameter_t

  implicit none

  !> \class nml_soil_moisture_1_t
  !> \brief Soil moisture - Case 1
  !> \details Feddes ET reduction with multi-layer infiltration and Brooks-Corey parameterization.
  type, public :: nml_soil_moisture_1_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    type(parameter_t) :: organic_matter_forest !< Organic matter content for forest
    type(parameter_t) :: organic_matter_impervious !< Organic matter content for impervious areas
    type(parameter_t) :: organic_matter_pervious !< Organic matter content for pervious areas
    type(parameter_t) :: ptf_lower_66_5_constant !< PTF constant below 66.5 percent sand
    type(parameter_t) :: ptf_lower_66_5_clay !< PTF clay multiplier below 66.5 percent sand
    type(parameter_t) :: ptf_lower_66_5_bulk_density !< PTF bulk-density multiplier below 66.5 percent sand
    type(parameter_t) :: ptf_upper_66_5_constant !< PTF constant above 66.5 percent sand
    type(parameter_t) :: ptf_upper_66_5_clay !< PTF clay multiplier above 66.5 percent sand
    type(parameter_t) :: ptf_upper_66_5_bulk_density !< PTF bulk-density multiplier above 66.5 percent sand
    type(parameter_t) :: ptf_ks_constant !< PTF constant for saturated hydraulic conductivity
    type(parameter_t) :: ptf_ks_sand !< PTF sand multiplier for saturated hydraulic conductivity
    type(parameter_t) :: ptf_ks_clay !< PTF clay multiplier for saturated hydraulic conductivity
    type(parameter_t) :: root_fraction_forest !< Root-fraction coefficient for forest
    type(parameter_t) :: root_fraction_impervious !< Root-fraction coefficient for impervious areas
    type(parameter_t) :: root_fraction_pervious !< Root-fraction coefficient for pervious areas
    type(parameter_t) :: infiltration_shape_factor !< Shape factor partitioning effective precipitation into runoff and infiltration
  contains
    procedure :: init => nml_soil_moisture_1_init
    procedure :: init_type => nml_soil_moisture_1_init_type
    procedure :: from_file => nml_soil_moisture_1_from_file
    procedure :: set => nml_soil_moisture_1_set
    procedure :: is_set => nml_soil_moisture_1_is_set
    procedure :: is_valid => nml_soil_moisture_1_is_valid
  end type nml_soil_moisture_1_t

contains

  !> \brief Initialize defaults and sentinels for soil_moisture_1
  integer function nml_soil_moisture_1_init(this, errmsg) result(status)
    class(nml_soil_moisture_1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! derived values
    status = this%init_type( &
      organic_matter_forest=this%organic_matter_forest, &
      organic_matter_impervious=this%organic_matter_impervious, &
      organic_matter_pervious=this%organic_matter_pervious, &
      ptf_lower_66_5_constant=this%ptf_lower_66_5_constant, &
      ptf_lower_66_5_clay=this%ptf_lower_66_5_clay, &
      ptf_lower_66_5_bulk_density=this%ptf_lower_66_5_bulk_density, &
      ptf_upper_66_5_constant=this%ptf_upper_66_5_constant, &
      ptf_upper_66_5_clay=this%ptf_upper_66_5_clay, &
      ptf_upper_66_5_bulk_density=this%ptf_upper_66_5_bulk_density, &
      ptf_ks_constant=this%ptf_ks_constant, &
      ptf_ks_sand=this%ptf_ks_sand, &
      ptf_ks_clay=this%ptf_ks_clay, &
      root_fraction_forest=this%root_fraction_forest, &
      root_fraction_impervious=this%root_fraction_impervious, &
      root_fraction_pervious=this%root_fraction_pervious, &
      infiltration_shape_factor=this%infiltration_shape_factor, &
      errmsg=errmsg)
    if (status /= NML_OK) return
  end function nml_soil_moisture_1_init

  !> \brief Initialize derived values with their field-specific defaults
  integer function nml_soil_moisture_1_init_type(this, &
    organic_matter_forest, &
    organic_matter_impervious, &
    organic_matter_pervious, &
    ptf_lower_66_5_constant, &
    ptf_lower_66_5_clay, &
    ptf_lower_66_5_bulk_density, &
    ptf_upper_66_5_constant, &
    ptf_upper_66_5_clay, &
    ptf_upper_66_5_bulk_density, &
    ptf_ks_constant, &
    ptf_ks_sand, &
    ptf_ks_clay, &
    root_fraction_forest, &
    root_fraction_impervious, &
    root_fraction_pervious, &
    infiltration_shape_factor, &
    errmsg) result(status)
    class(nml_soil_moisture_1_t), intent(in) :: this !< parent namelist instance
    type(parameter_t), intent(inout), optional :: organic_matter_forest !< Organic matter content for forest
    type(parameter_t), intent(inout), optional :: organic_matter_impervious !< Organic matter content for impervious areas
    type(parameter_t), intent(inout), optional :: organic_matter_pervious !< Organic matter content for pervious areas
    type(parameter_t), intent(inout), optional :: ptf_lower_66_5_constant !< PTF constant below 66.5 percent sand
    type(parameter_t), intent(inout), optional :: ptf_lower_66_5_clay !< PTF clay multiplier below 66.5 percent sand
    type(parameter_t), intent(inout), optional :: ptf_lower_66_5_bulk_density !< PTF bulk-density multiplier below 66.5 percent sand
    type(parameter_t), intent(inout), optional :: ptf_upper_66_5_constant !< PTF constant above 66.5 percent sand
    type(parameter_t), intent(inout), optional :: ptf_upper_66_5_clay !< PTF clay multiplier above 66.5 percent sand
    type(parameter_t), intent(inout), optional :: ptf_upper_66_5_bulk_density !< PTF bulk-density multiplier above 66.5 percent sand
    type(parameter_t), intent(inout), optional :: ptf_ks_constant !< PTF constant for saturated hydraulic conductivity
    type(parameter_t), intent(inout), optional :: ptf_ks_sand !< PTF sand multiplier for saturated hydraulic conductivity
    type(parameter_t), intent(inout), optional :: ptf_ks_clay !< PTF clay multiplier for saturated hydraulic conductivity
    type(parameter_t), intent(inout), optional :: root_fraction_forest !< Root-fraction coefficient for forest
    type(parameter_t), intent(inout), optional :: root_fraction_impervious !< Root-fraction coefficient for impervious areas
    type(parameter_t), intent(inout), optional :: root_fraction_pervious !< Root-fraction coefficient for pervious areas
    type(parameter_t), intent(inout), optional :: infiltration_shape_factor !< Shape factor partitioning effective precipitation into runoff and infiltration
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (present(organic_matter_forest)) then
      organic_matter_forest%value = ieee_value(organic_matter_forest%value, ieee_quiet_nan) ! sentinel for derived component value
      organic_matter_forest%optimize = .false.
      organic_matter_forest%min = ieee_value(organic_matter_forest%min, ieee_quiet_nan) ! sentinel for derived component min
      organic_matter_forest%max = ieee_value(organic_matter_forest%max, ieee_quiet_nan) ! sentinel for derived component max
      organic_matter_forest%min = 5.0_dp
      organic_matter_forest%max = 10.0_dp
    end if
    if (present(organic_matter_impervious)) then
      organic_matter_impervious%value = ieee_value(organic_matter_impervious%value, ieee_quiet_nan) ! sentinel for derived component value
      organic_matter_impervious%optimize = .false.
      organic_matter_impervious%min = ieee_value(organic_matter_impervious%min, ieee_quiet_nan) ! sentinel for derived component min
      organic_matter_impervious%max = ieee_value(organic_matter_impervious%max, ieee_quiet_nan) ! sentinel for derived component max
      organic_matter_impervious%min = 0.0_dp
      organic_matter_impervious%max = 1.0_dp
    end if
    if (present(organic_matter_pervious)) then
      organic_matter_pervious%value = ieee_value(organic_matter_pervious%value, ieee_quiet_nan) ! sentinel for derived component value
      organic_matter_pervious%optimize = .false.
      organic_matter_pervious%min = ieee_value(organic_matter_pervious%min, ieee_quiet_nan) ! sentinel for derived component min
      organic_matter_pervious%max = ieee_value(organic_matter_pervious%max, ieee_quiet_nan) ! sentinel for derived component max
      organic_matter_pervious%min = 1.0_dp
      organic_matter_pervious%max = 5.0_dp
    end if
    if (present(ptf_lower_66_5_constant)) then
      ptf_lower_66_5_constant%value = ieee_value(ptf_lower_66_5_constant%value, ieee_quiet_nan) ! sentinel for derived component value
      ptf_lower_66_5_constant%optimize = .false.
      ptf_lower_66_5_constant%min = ieee_value(ptf_lower_66_5_constant%min, ieee_quiet_nan) ! sentinel for derived component min
      ptf_lower_66_5_constant%max = ieee_value(ptf_lower_66_5_constant%max, ieee_quiet_nan) ! sentinel for derived component max
      ptf_lower_66_5_constant%min = 0.75_dp
      ptf_lower_66_5_constant%max = 0.8_dp
    end if
    if (present(ptf_lower_66_5_clay)) then
      ptf_lower_66_5_clay%value = ieee_value(ptf_lower_66_5_clay%value, ieee_quiet_nan) ! sentinel for derived component value
      ptf_lower_66_5_clay%optimize = .false.
      ptf_lower_66_5_clay%min = ieee_value(ptf_lower_66_5_clay%min, ieee_quiet_nan) ! sentinel for derived component min
      ptf_lower_66_5_clay%max = ieee_value(ptf_lower_66_5_clay%max, ieee_quiet_nan) ! sentinel for derived component max
      ptf_lower_66_5_clay%min = 0.0008_dp
      ptf_lower_66_5_clay%max = 0.0012_dp
    end if
    if (present(ptf_lower_66_5_bulk_density)) then
      ptf_lower_66_5_bulk_density%value = ieee_value(ptf_lower_66_5_bulk_density%value, ieee_quiet_nan) ! sentinel for derived component value
      ptf_lower_66_5_bulk_density%optimize = .false.
      ptf_lower_66_5_bulk_density%min = ieee_value(ptf_lower_66_5_bulk_density%min, ieee_quiet_nan) ! sentinel for derived component min
      ptf_lower_66_5_bulk_density%max = ieee_value(ptf_lower_66_5_bulk_density%max, ieee_quiet_nan) ! sentinel for derived component max
      ptf_lower_66_5_bulk_density%min = -0.27_dp
      ptf_lower_66_5_bulk_density%max = -0.25_dp
    end if
    if (present(ptf_upper_66_5_constant)) then
      ptf_upper_66_5_constant%value = ieee_value(ptf_upper_66_5_constant%value, ieee_quiet_nan) ! sentinel for derived component value
      ptf_upper_66_5_constant%optimize = .false.
      ptf_upper_66_5_constant%min = ieee_value(ptf_upper_66_5_constant%min, ieee_quiet_nan) ! sentinel for derived component min
      ptf_upper_66_5_constant%max = ieee_value(ptf_upper_66_5_constant%max, ieee_quiet_nan) ! sentinel for derived component max
      ptf_upper_66_5_constant%min = 0.8_dp
      ptf_upper_66_5_constant%max = 0.9_dp
    end if
    if (present(ptf_upper_66_5_clay)) then
      ptf_upper_66_5_clay%value = ieee_value(ptf_upper_66_5_clay%value, ieee_quiet_nan) ! sentinel for derived component value
      ptf_upper_66_5_clay%optimize = .false.
      ptf_upper_66_5_clay%min = ieee_value(ptf_upper_66_5_clay%min, ieee_quiet_nan) ! sentinel for derived component min
      ptf_upper_66_5_clay%max = ieee_value(ptf_upper_66_5_clay%max, ieee_quiet_nan) ! sentinel for derived component max
      ptf_upper_66_5_clay%min = -0.0012_dp
      ptf_upper_66_5_clay%max = -0.0008_dp
    end if
    if (present(ptf_upper_66_5_bulk_density)) then
      ptf_upper_66_5_bulk_density%value = ieee_value(ptf_upper_66_5_bulk_density%value, ieee_quiet_nan) ! sentinel for derived component value
      ptf_upper_66_5_bulk_density%optimize = .false.
      ptf_upper_66_5_bulk_density%min = ieee_value(ptf_upper_66_5_bulk_density%min, ieee_quiet_nan) ! sentinel for derived component min
      ptf_upper_66_5_bulk_density%max = ieee_value(ptf_upper_66_5_bulk_density%max, ieee_quiet_nan) ! sentinel for derived component max
      ptf_upper_66_5_bulk_density%min = -0.35_dp
      ptf_upper_66_5_bulk_density%max = -0.3_dp
    end if
    if (present(ptf_ks_constant)) then
      ptf_ks_constant%value = ieee_value(ptf_ks_constant%value, ieee_quiet_nan) ! sentinel for derived component value
      ptf_ks_constant%optimize = .false.
      ptf_ks_constant%min = ieee_value(ptf_ks_constant%min, ieee_quiet_nan) ! sentinel for derived component min
      ptf_ks_constant%max = ieee_value(ptf_ks_constant%max, ieee_quiet_nan) ! sentinel for derived component max
      ptf_ks_constant%min = -1.2_dp
      ptf_ks_constant%max = -0.285_dp
    end if
    if (present(ptf_ks_sand)) then
      ptf_ks_sand%value = ieee_value(ptf_ks_sand%value, ieee_quiet_nan) ! sentinel for derived component value
      ptf_ks_sand%optimize = .false.
      ptf_ks_sand%min = ieee_value(ptf_ks_sand%min, ieee_quiet_nan) ! sentinel for derived component min
      ptf_ks_sand%max = ieee_value(ptf_ks_sand%max, ieee_quiet_nan) ! sentinel for derived component max
      ptf_ks_sand%min = 0.006_dp
      ptf_ks_sand%max = 0.026_dp
    end if
    if (present(ptf_ks_clay)) then
      ptf_ks_clay%value = ieee_value(ptf_ks_clay%value, ieee_quiet_nan) ! sentinel for derived component value
      ptf_ks_clay%optimize = .false.
      ptf_ks_clay%min = ieee_value(ptf_ks_clay%min, ieee_quiet_nan) ! sentinel for derived component min
      ptf_ks_clay%max = ieee_value(ptf_ks_clay%max, ieee_quiet_nan) ! sentinel for derived component max
      ptf_ks_clay%min = 0.003_dp
      ptf_ks_clay%max = 0.013_dp
    end if
    if (present(root_fraction_forest)) then
      root_fraction_forest%value = ieee_value(root_fraction_forest%value, ieee_quiet_nan) ! sentinel for derived component value
      root_fraction_forest%optimize = .false.
      root_fraction_forest%min = ieee_value(root_fraction_forest%min, ieee_quiet_nan) ! sentinel for derived component min
      root_fraction_forest%max = ieee_value(root_fraction_forest%max, ieee_quiet_nan) ! sentinel for derived component max
      root_fraction_forest%min = 0.9_dp
      root_fraction_forest%max = 0.999_dp
    end if
    if (present(root_fraction_impervious)) then
      root_fraction_impervious%value = ieee_value(root_fraction_impervious%value, ieee_quiet_nan) ! sentinel for derived component value
      root_fraction_impervious%optimize = .false.
      root_fraction_impervious%min = ieee_value(root_fraction_impervious%min, ieee_quiet_nan) ! sentinel for derived component min
      root_fraction_impervious%max = ieee_value(root_fraction_impervious%max, ieee_quiet_nan) ! sentinel for derived component max
      root_fraction_impervious%min = 0.9_dp
      root_fraction_impervious%max = 0.95_dp
    end if
    if (present(root_fraction_pervious)) then
      root_fraction_pervious%value = ieee_value(root_fraction_pervious%value, ieee_quiet_nan) ! sentinel for derived component value
      root_fraction_pervious%optimize = .false.
      root_fraction_pervious%min = ieee_value(root_fraction_pervious%min, ieee_quiet_nan) ! sentinel for derived component min
      root_fraction_pervious%max = ieee_value(root_fraction_pervious%max, ieee_quiet_nan) ! sentinel for derived component max
      root_fraction_pervious%min = 0.001_dp
      root_fraction_pervious%max = 0.09_dp
    end if
    if (present(infiltration_shape_factor)) then
      infiltration_shape_factor%value = ieee_value(infiltration_shape_factor%value, ieee_quiet_nan) ! sentinel for derived component value
      infiltration_shape_factor%optimize = .false.
      infiltration_shape_factor%min = ieee_value(infiltration_shape_factor%min, ieee_quiet_nan) ! sentinel for derived component min
      infiltration_shape_factor%max = ieee_value(infiltration_shape_factor%max, ieee_quiet_nan) ! sentinel for derived component max
      infiltration_shape_factor%min = 1.0_dp
      infiltration_shape_factor%max = 4.0_dp
    end if
  end function nml_soil_moisture_1_init_type


  !> \brief Read soil_moisture_1 namelist from file
  integer function nml_soil_moisture_1_from_file(this, file, errmsg) result(status)
    class(nml_soil_moisture_1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    type(parameter_t) :: organic_matter_forest
    type(parameter_t) :: organic_matter_impervious
    type(parameter_t) :: organic_matter_pervious
    type(parameter_t) :: ptf_lower_66_5_constant
    type(parameter_t) :: ptf_lower_66_5_clay
    type(parameter_t) :: ptf_lower_66_5_bulk_density
    type(parameter_t) :: ptf_upper_66_5_constant
    type(parameter_t) :: ptf_upper_66_5_clay
    type(parameter_t) :: ptf_upper_66_5_bulk_density
    type(parameter_t) :: ptf_ks_constant
    type(parameter_t) :: ptf_ks_sand
    type(parameter_t) :: ptf_ks_clay
    type(parameter_t) :: root_fraction_forest
    type(parameter_t) :: root_fraction_impervious
    type(parameter_t) :: root_fraction_pervious
    type(parameter_t) :: infiltration_shape_factor
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /soil_moisture_1/ &
      organic_matter_forest, &
      organic_matter_impervious, &
      organic_matter_pervious, &
      ptf_lower_66_5_constant, &
      ptf_lower_66_5_clay, &
      ptf_lower_66_5_bulk_density, &
      ptf_upper_66_5_constant, &
      ptf_upper_66_5_clay, &
      ptf_upper_66_5_bulk_density, &
      ptf_ks_constant, &
      ptf_ks_sand, &
      ptf_ks_clay, &
      root_fraction_forest, &
      root_fraction_impervious, &
      root_fraction_pervious, &
      infiltration_shape_factor

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    organic_matter_forest = this%organic_matter_forest
    organic_matter_impervious = this%organic_matter_impervious
    organic_matter_pervious = this%organic_matter_pervious
    ptf_lower_66_5_constant = this%ptf_lower_66_5_constant
    ptf_lower_66_5_clay = this%ptf_lower_66_5_clay
    ptf_lower_66_5_bulk_density = this%ptf_lower_66_5_bulk_density
    ptf_upper_66_5_constant = this%ptf_upper_66_5_constant
    ptf_upper_66_5_clay = this%ptf_upper_66_5_clay
    ptf_upper_66_5_bulk_density = this%ptf_upper_66_5_bulk_density
    ptf_ks_constant = this%ptf_ks_constant
    ptf_ks_sand = this%ptf_ks_sand
    ptf_ks_clay = this%ptf_ks_clay
    root_fraction_forest = this%root_fraction_forest
    root_fraction_impervious = this%root_fraction_impervious
    root_fraction_pervious = this%root_fraction_pervious
    infiltration_shape_factor = this%infiltration_shape_factor

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("soil_moisture_1", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=soil_moisture_1, iostat=iostat, iomsg=iomsg)
    if (iostat /= 0) then
      status = NML_ERR_READ
      if (present(errmsg)) errmsg = trim(iomsg)
      close_status = nml%close()
      return
    end if
    close_status = nml%close(errmsg=errmsg)
    if (close_status /= NML_OK) then
      status = close_status
      return
    end if

    ! assign values
    this%organic_matter_forest = organic_matter_forest
    this%organic_matter_impervious = organic_matter_impervious
    this%organic_matter_pervious = organic_matter_pervious
    this%ptf_lower_66_5_constant = ptf_lower_66_5_constant
    this%ptf_lower_66_5_clay = ptf_lower_66_5_clay
    this%ptf_lower_66_5_bulk_density = ptf_lower_66_5_bulk_density
    this%ptf_upper_66_5_constant = ptf_upper_66_5_constant
    this%ptf_upper_66_5_clay = ptf_upper_66_5_clay
    this%ptf_upper_66_5_bulk_density = ptf_upper_66_5_bulk_density
    this%ptf_ks_constant = ptf_ks_constant
    this%ptf_ks_sand = ptf_ks_sand
    this%ptf_ks_clay = ptf_ks_clay
    this%root_fraction_forest = root_fraction_forest
    this%root_fraction_impervious = root_fraction_impervious
    this%root_fraction_pervious = root_fraction_pervious
    this%infiltration_shape_factor = infiltration_shape_factor

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_soil_moisture_1_from_file

  !> \brief Set soil_moisture_1 values
  integer function nml_soil_moisture_1_set(this, &
    organic_matter_forest, &
    organic_matter_impervious, &
    organic_matter_pervious, &
    ptf_lower_66_5_constant, &
    ptf_lower_66_5_clay, &
    ptf_lower_66_5_bulk_density, &
    ptf_upper_66_5_constant, &
    ptf_upper_66_5_clay, &
    ptf_upper_66_5_bulk_density, &
    ptf_ks_constant, &
    ptf_ks_sand, &
    ptf_ks_clay, &
    root_fraction_forest, &
    root_fraction_impervious, &
    root_fraction_pervious, &
    infiltration_shape_factor, &
    errmsg) result(status)

    class(nml_soil_moisture_1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    type(parameter_t), intent(in) :: organic_matter_forest !< Organic matter content for forest
    type(parameter_t), intent(in) :: organic_matter_impervious !< Organic matter content for impervious areas
    type(parameter_t), intent(in) :: organic_matter_pervious !< Organic matter content for pervious areas
    type(parameter_t), intent(in) :: ptf_lower_66_5_constant !< PTF constant below 66.5 percent sand
    type(parameter_t), intent(in) :: ptf_lower_66_5_clay !< PTF clay multiplier below 66.5 percent sand
    type(parameter_t), intent(in) :: ptf_lower_66_5_bulk_density !< PTF bulk-density multiplier below 66.5 percent sand
    type(parameter_t), intent(in) :: ptf_upper_66_5_constant !< PTF constant above 66.5 percent sand
    type(parameter_t), intent(in) :: ptf_upper_66_5_clay !< PTF clay multiplier above 66.5 percent sand
    type(parameter_t), intent(in) :: ptf_upper_66_5_bulk_density !< PTF bulk-density multiplier above 66.5 percent sand
    type(parameter_t), intent(in) :: ptf_ks_constant !< PTF constant for saturated hydraulic conductivity
    type(parameter_t), intent(in) :: ptf_ks_sand !< PTF sand multiplier for saturated hydraulic conductivity
    type(parameter_t), intent(in) :: ptf_ks_clay !< PTF clay multiplier for saturated hydraulic conductivity
    type(parameter_t), intent(in) :: root_fraction_forest !< Root-fraction coefficient for forest
    type(parameter_t), intent(in) :: root_fraction_impervious !< Root-fraction coefficient for impervious areas
    type(parameter_t), intent(in) :: root_fraction_pervious !< Root-fraction coefficient for pervious areas
    type(parameter_t), intent(in) :: infiltration_shape_factor !< Shape factor partitioning effective precipitation into runoff and infiltration

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    this%organic_matter_forest = organic_matter_forest
    this%organic_matter_impervious = organic_matter_impervious
    this%organic_matter_pervious = organic_matter_pervious
    this%ptf_lower_66_5_constant = ptf_lower_66_5_constant
    this%ptf_lower_66_5_clay = ptf_lower_66_5_clay
    this%ptf_lower_66_5_bulk_density = ptf_lower_66_5_bulk_density
    this%ptf_upper_66_5_constant = ptf_upper_66_5_constant
    this%ptf_upper_66_5_clay = ptf_upper_66_5_clay
    this%ptf_upper_66_5_bulk_density = ptf_upper_66_5_bulk_density
    this%ptf_ks_constant = ptf_ks_constant
    this%ptf_ks_sand = ptf_ks_sand
    this%ptf_ks_clay = ptf_ks_clay
    this%root_fraction_forest = root_fraction_forest
    this%root_fraction_impervious = root_fraction_impervious
    this%root_fraction_pervious = root_fraction_pervious
    this%infiltration_shape_factor = infiltration_shape_factor

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_soil_moisture_1_set

  !> \brief Check whether a namelist value was set
  integer function nml_soil_moisture_1_is_set(this, name, idx, errmsg) result(status)
    class(nml_soil_moisture_1_t), intent(in) :: this !< namelist instance
    character(len=*), intent(in) :: name !< field name
    integer, intent(in), optional :: idx(:) !< optional field index values
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (.not. this%is_configured) then
      status = NML_ERR_NOT_SET
      if (present(errmsg)) errmsg = "namelist not configured; call set or from_file"
      return
    end if
    select case (to_lower(trim(name)))
    case ("organic_matter_forest%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'organic_matter_forest'"
        return
      end if
      if (ieee_is_nan(this%organic_matter_forest%value)) status = NML_ERR_NOT_SET
    case ("organic_matter_forest%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'organic_matter_forest'"
        return
      end if
    case ("organic_matter_forest%min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'organic_matter_forest'"
        return
      end if
    case ("organic_matter_forest%max")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'organic_matter_forest'"
        return
      end if
    case ("organic_matter_forest")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'organic_matter_forest'"
        return
      end if
      if (ieee_is_nan(this%organic_matter_forest%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("organic_matter_impervious%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'organic_matter_impervious'"
        return
      end if
      if (ieee_is_nan(this%organic_matter_impervious%value)) status = NML_ERR_NOT_SET
    case ("organic_matter_impervious%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'organic_matter_impervious'"
        return
      end if
    case ("organic_matter_impervious%min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'organic_matter_impervious'"
        return
      end if
    case ("organic_matter_impervious%max")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'organic_matter_impervious'"
        return
      end if
    case ("organic_matter_impervious")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'organic_matter_impervious'"
        return
      end if
      if (ieee_is_nan(this%organic_matter_impervious%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("organic_matter_pervious%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'organic_matter_pervious'"
        return
      end if
      if (ieee_is_nan(this%organic_matter_pervious%value)) status = NML_ERR_NOT_SET
    case ("organic_matter_pervious%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'organic_matter_pervious'"
        return
      end if
    case ("organic_matter_pervious%min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'organic_matter_pervious'"
        return
      end if
    case ("organic_matter_pervious%max")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'organic_matter_pervious'"
        return
      end if
    case ("organic_matter_pervious")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'organic_matter_pervious'"
        return
      end if
      if (ieee_is_nan(this%organic_matter_pervious%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("ptf_lower_66_5_constant%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_lower_66_5_constant'"
        return
      end if
      if (ieee_is_nan(this%ptf_lower_66_5_constant%value)) status = NML_ERR_NOT_SET
    case ("ptf_lower_66_5_constant%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_lower_66_5_constant'"
        return
      end if
    case ("ptf_lower_66_5_constant%min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_lower_66_5_constant'"
        return
      end if
    case ("ptf_lower_66_5_constant%max")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_lower_66_5_constant'"
        return
      end if
    case ("ptf_lower_66_5_constant")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_lower_66_5_constant'"
        return
      end if
      if (ieee_is_nan(this%ptf_lower_66_5_constant%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("ptf_lower_66_5_clay%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_lower_66_5_clay'"
        return
      end if
      if (ieee_is_nan(this%ptf_lower_66_5_clay%value)) status = NML_ERR_NOT_SET
    case ("ptf_lower_66_5_clay%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_lower_66_5_clay'"
        return
      end if
    case ("ptf_lower_66_5_clay%min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_lower_66_5_clay'"
        return
      end if
    case ("ptf_lower_66_5_clay%max")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_lower_66_5_clay'"
        return
      end if
    case ("ptf_lower_66_5_clay")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_lower_66_5_clay'"
        return
      end if
      if (ieee_is_nan(this%ptf_lower_66_5_clay%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("ptf_lower_66_5_bulk_density%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_lower_66_5_bulk_density'"
        return
      end if
      if (ieee_is_nan(this%ptf_lower_66_5_bulk_density%value)) status = NML_ERR_NOT_SET
    case ("ptf_lower_66_5_bulk_density%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_lower_66_5_bulk_density'"
        return
      end if
    case ("ptf_lower_66_5_bulk_density%min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_lower_66_5_bulk_density'"
        return
      end if
    case ("ptf_lower_66_5_bulk_density%max")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_lower_66_5_bulk_density'"
        return
      end if
    case ("ptf_lower_66_5_bulk_density")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_lower_66_5_bulk_density'"
        return
      end if
      if (ieee_is_nan(this%ptf_lower_66_5_bulk_density%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("ptf_upper_66_5_constant%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_upper_66_5_constant'"
        return
      end if
      if (ieee_is_nan(this%ptf_upper_66_5_constant%value)) status = NML_ERR_NOT_SET
    case ("ptf_upper_66_5_constant%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_upper_66_5_constant'"
        return
      end if
    case ("ptf_upper_66_5_constant%min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_upper_66_5_constant'"
        return
      end if
    case ("ptf_upper_66_5_constant%max")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_upper_66_5_constant'"
        return
      end if
    case ("ptf_upper_66_5_constant")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_upper_66_5_constant'"
        return
      end if
      if (ieee_is_nan(this%ptf_upper_66_5_constant%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("ptf_upper_66_5_clay%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_upper_66_5_clay'"
        return
      end if
      if (ieee_is_nan(this%ptf_upper_66_5_clay%value)) status = NML_ERR_NOT_SET
    case ("ptf_upper_66_5_clay%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_upper_66_5_clay'"
        return
      end if
    case ("ptf_upper_66_5_clay%min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_upper_66_5_clay'"
        return
      end if
    case ("ptf_upper_66_5_clay%max")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_upper_66_5_clay'"
        return
      end if
    case ("ptf_upper_66_5_clay")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_upper_66_5_clay'"
        return
      end if
      if (ieee_is_nan(this%ptf_upper_66_5_clay%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("ptf_upper_66_5_bulk_density%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_upper_66_5_bulk_density'"
        return
      end if
      if (ieee_is_nan(this%ptf_upper_66_5_bulk_density%value)) status = NML_ERR_NOT_SET
    case ("ptf_upper_66_5_bulk_density%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_upper_66_5_bulk_density'"
        return
      end if
    case ("ptf_upper_66_5_bulk_density%min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_upper_66_5_bulk_density'"
        return
      end if
    case ("ptf_upper_66_5_bulk_density%max")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_upper_66_5_bulk_density'"
        return
      end if
    case ("ptf_upper_66_5_bulk_density")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_upper_66_5_bulk_density'"
        return
      end if
      if (ieee_is_nan(this%ptf_upper_66_5_bulk_density%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("ptf_ks_constant%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_ks_constant'"
        return
      end if
      if (ieee_is_nan(this%ptf_ks_constant%value)) status = NML_ERR_NOT_SET
    case ("ptf_ks_constant%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_ks_constant'"
        return
      end if
    case ("ptf_ks_constant%min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_ks_constant'"
        return
      end if
    case ("ptf_ks_constant%max")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_ks_constant'"
        return
      end if
    case ("ptf_ks_constant")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_ks_constant'"
        return
      end if
      if (ieee_is_nan(this%ptf_ks_constant%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("ptf_ks_sand%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_ks_sand'"
        return
      end if
      if (ieee_is_nan(this%ptf_ks_sand%value)) status = NML_ERR_NOT_SET
    case ("ptf_ks_sand%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_ks_sand'"
        return
      end if
    case ("ptf_ks_sand%min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_ks_sand'"
        return
      end if
    case ("ptf_ks_sand%max")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_ks_sand'"
        return
      end if
    case ("ptf_ks_sand")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_ks_sand'"
        return
      end if
      if (ieee_is_nan(this%ptf_ks_sand%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("ptf_ks_clay%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_ks_clay'"
        return
      end if
      if (ieee_is_nan(this%ptf_ks_clay%value)) status = NML_ERR_NOT_SET
    case ("ptf_ks_clay%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_ks_clay'"
        return
      end if
    case ("ptf_ks_clay%min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_ks_clay'"
        return
      end if
    case ("ptf_ks_clay%max")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_ks_clay'"
        return
      end if
    case ("ptf_ks_clay")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_ks_clay'"
        return
      end if
      if (ieee_is_nan(this%ptf_ks_clay%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("root_fraction_forest%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'root_fraction_forest'"
        return
      end if
      if (ieee_is_nan(this%root_fraction_forest%value)) status = NML_ERR_NOT_SET
    case ("root_fraction_forest%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'root_fraction_forest'"
        return
      end if
    case ("root_fraction_forest%min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'root_fraction_forest'"
        return
      end if
    case ("root_fraction_forest%max")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'root_fraction_forest'"
        return
      end if
    case ("root_fraction_forest")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'root_fraction_forest'"
        return
      end if
      if (ieee_is_nan(this%root_fraction_forest%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("root_fraction_impervious%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'root_fraction_impervious'"
        return
      end if
      if (ieee_is_nan(this%root_fraction_impervious%value)) status = NML_ERR_NOT_SET
    case ("root_fraction_impervious%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'root_fraction_impervious'"
        return
      end if
    case ("root_fraction_impervious%min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'root_fraction_impervious'"
        return
      end if
    case ("root_fraction_impervious%max")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'root_fraction_impervious'"
        return
      end if
    case ("root_fraction_impervious")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'root_fraction_impervious'"
        return
      end if
      if (ieee_is_nan(this%root_fraction_impervious%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("root_fraction_pervious%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'root_fraction_pervious'"
        return
      end if
      if (ieee_is_nan(this%root_fraction_pervious%value)) status = NML_ERR_NOT_SET
    case ("root_fraction_pervious%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'root_fraction_pervious'"
        return
      end if
    case ("root_fraction_pervious%min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'root_fraction_pervious'"
        return
      end if
    case ("root_fraction_pervious%max")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'root_fraction_pervious'"
        return
      end if
    case ("root_fraction_pervious")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'root_fraction_pervious'"
        return
      end if
      if (ieee_is_nan(this%root_fraction_pervious%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("infiltration_shape_factor%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'infiltration_shape_factor'"
        return
      end if
      if (ieee_is_nan(this%infiltration_shape_factor%value)) status = NML_ERR_NOT_SET
    case ("infiltration_shape_factor%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'infiltration_shape_factor'"
        return
      end if
    case ("infiltration_shape_factor%min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'infiltration_shape_factor'"
        return
      end if
    case ("infiltration_shape_factor%max")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'infiltration_shape_factor'"
        return
      end if
    case ("infiltration_shape_factor")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'infiltration_shape_factor'"
        return
      end if
      if (ieee_is_nan(this%infiltration_shape_factor%value)) then
        status = NML_ERR_NOT_SET
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_soil_moisture_1_is_set

  !> \brief Validate required values and constraints
  integer function nml_soil_moisture_1_is_valid(this, errmsg) result(status)
    class(nml_soil_moisture_1_t), intent(in) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    integer :: istat

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (.not. this%is_configured) then
      status = NML_ERR_NOT_SET
      if (present(errmsg)) errmsg = "namelist not configured; call set or from_file"
      return
    end if

    ! required parameters
    istat = this%is_set("organic_matter_forest", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: organic_matter_forest"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("organic_matter_impervious", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: organic_matter_impervious"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("organic_matter_pervious", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: organic_matter_pervious"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("ptf_lower_66_5_constant", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: ptf_lower_66_5_constant"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("ptf_lower_66_5_clay", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: ptf_lower_66_5_clay"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("ptf_lower_66_5_bulk_density", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: ptf_lower_66_5_bulk_density"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("ptf_upper_66_5_constant", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: ptf_upper_66_5_constant"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("ptf_upper_66_5_clay", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: ptf_upper_66_5_clay"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("ptf_upper_66_5_bulk_density", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: ptf_upper_66_5_bulk_density"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("ptf_ks_constant", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: ptf_ks_constant"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("ptf_ks_sand", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: ptf_ks_sand"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("ptf_ks_clay", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: ptf_ks_clay"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("root_fraction_forest", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: root_fraction_forest"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("root_fraction_impervious", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: root_fraction_impervious"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("root_fraction_pervious", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: root_fraction_pervious"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("infiltration_shape_factor", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: infiltration_shape_factor"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
  end function nml_soil_moisture_1_is_valid

end module nml_soil_moisture_1
