!> \file nml_parameter_soilmoisture1.f90
!> \copydoc nml_soilmoisture1

!> \brief Soil moisture - Case 1
!> \details Feddes equation for ET reduction, multi-layer infiltration capacity approach, Brooks-Corey like
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_soilmoisture1
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
    to_lower, &
    NML_ERR_PARTLY_SET
  use ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
  ! kind specifiers listed in the nml-tools configuration file
  use mo_kind, only: &
    dp

  implicit none

  ! default values
  real(dp), parameter, public :: ptf_ks_curveslope__default(5) = [60.96_dp, 60.96_dp, 60.96_dp, 0.0_dp, 1.0_dp]

  !> \class nml_soilmoisture1_t
  !> \brief Soil moisture - Case 1
  !> \details Feddes equation for ET reduction, multi-layer infiltration capacity approach, Brooks-Corey like
  type, public :: nml_soilmoisture1_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    real(dp), dimension(5) :: orgmattercontent_forest !< Organic matter content for forest
    real(dp), dimension(5) :: orgmattercontent_impervious !< Organic matter content for impervious
    real(dp), dimension(5) :: orgmattercontent_pervious !< Organic matter content for pervious
    real(dp), dimension(5) :: ptf_lower66_5_constant !< Zacharias PTF parameters below 66.5 % sand content
    real(dp), dimension(5) :: ptf_lower66_5_clay !< Multiplier for clay constant below 66.5 % sand content
    real(dp), dimension(5) :: ptf_lower66_5_db !< Multiplier for mineral bulk density below 66.5 % sand content
    real(dp), dimension(5) :: ptf_higher66_5_constant !< Zacharias PTF parameters above 66.5 % sand content
    real(dp), dimension(5) :: ptf_higher66_5_clay !< Multiplier for clay constant above 66.5 % sand content
    real(dp), dimension(5) :: ptf_higher66_5_db !< Multiplier for mineral bulk density above 66.5 % sand content
    real(dp), dimension(5) :: ptf_ks_constant !< PTF constant for saturated hydraulic conductivity
    real(dp), dimension(5) :: ptf_ks_sand !< Multiplier for sand for saturated hydraulic conductivity
    real(dp), dimension(5) :: ptf_ks_clay !< Multiplier for clay for saturated hydraulic conductivity
    real(dp), dimension(5) :: ptf_ks_curveslope !< Unit conversion factor
    real(dp), dimension(5) :: rootfractioncoefficient_forest !< Root fraction coefficient for forest
    real(dp), dimension(5) :: rootfractioncoefficient_impervious !< Root fraction coefficient for impervious
    real(dp), dimension(5) :: rootfractioncoefficient_pervious !< Root fraction coefficient for pervious
    real(dp), dimension(5) :: infiltrationshapefactor !< Infiltration shape factor
  contains
    procedure :: init => nml_soilmoisture1_init
    procedure :: from_file => nml_soilmoisture1_from_file
    procedure :: set => nml_soilmoisture1_set
    procedure :: is_set => nml_soilmoisture1_is_set
    procedure :: is_valid => nml_soilmoisture1_is_valid
  end type nml_soilmoisture1_t

contains

  !> \brief Initialize defaults and sentinels for soilmoisture1
  integer function nml_soilmoisture1_init(this, errmsg) result(status)
    class(nml_soilmoisture1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! sentinel values for required/optional parameters
    this%orgmattercontent_forest = ieee_value(this%orgmattercontent_forest, ieee_quiet_nan) ! sentinel for required real array
    this%orgmattercontent_impervious = ieee_value(this%orgmattercontent_impervious, ieee_quiet_nan) ! sentinel for required real array
    this%orgmattercontent_pervious = ieee_value(this%orgmattercontent_pervious, ieee_quiet_nan) ! sentinel for required real array
    this%ptf_lower66_5_constant = ieee_value(this%ptf_lower66_5_constant, ieee_quiet_nan) ! sentinel for required real array
    this%ptf_lower66_5_clay = ieee_value(this%ptf_lower66_5_clay, ieee_quiet_nan) ! sentinel for required real array
    this%ptf_lower66_5_db = ieee_value(this%ptf_lower66_5_db, ieee_quiet_nan) ! sentinel for required real array
    this%ptf_higher66_5_constant = ieee_value(this%ptf_higher66_5_constant, ieee_quiet_nan) ! sentinel for required real array
    this%ptf_higher66_5_clay = ieee_value(this%ptf_higher66_5_clay, ieee_quiet_nan) ! sentinel for required real array
    this%ptf_higher66_5_db = ieee_value(this%ptf_higher66_5_db, ieee_quiet_nan) ! sentinel for required real array
    this%ptf_ks_constant = ieee_value(this%ptf_ks_constant, ieee_quiet_nan) ! sentinel for required real array
    this%ptf_ks_sand = ieee_value(this%ptf_ks_sand, ieee_quiet_nan) ! sentinel for required real array
    this%ptf_ks_clay = ieee_value(this%ptf_ks_clay, ieee_quiet_nan) ! sentinel for required real array
    this%rootfractioncoefficient_forest = ieee_value(this%rootfractioncoefficient_forest, ieee_quiet_nan) ! sentinel for required real array
    this%rootfractioncoefficient_impervious = ieee_value(this%rootfractioncoefficient_impervious, ieee_quiet_nan) ! sentinel for required real array
    this%rootfractioncoefficient_pervious = ieee_value(this%rootfractioncoefficient_pervious, ieee_quiet_nan) ! sentinel for required real array
    this%infiltrationshapefactor = ieee_value(this%infiltrationshapefactor, ieee_quiet_nan) ! sentinel for required real array
    ! default values
    this%ptf_ks_curveslope = ptf_ks_curveslope__default
  end function nml_soilmoisture1_init


  !> \brief Read soilmoisture1 namelist from file
  integer function nml_soilmoisture1_from_file(this, file, errmsg) result(status)
    class(nml_soilmoisture1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    real(dp), dimension(5) :: orgmattercontent_forest
    real(dp), dimension(5) :: orgmattercontent_impervious
    real(dp), dimension(5) :: orgmattercontent_pervious
    real(dp), dimension(5) :: ptf_lower66_5_constant
    real(dp), dimension(5) :: ptf_lower66_5_clay
    real(dp), dimension(5) :: ptf_lower66_5_db
    real(dp), dimension(5) :: ptf_higher66_5_constant
    real(dp), dimension(5) :: ptf_higher66_5_clay
    real(dp), dimension(5) :: ptf_higher66_5_db
    real(dp), dimension(5) :: ptf_ks_constant
    real(dp), dimension(5) :: ptf_ks_sand
    real(dp), dimension(5) :: ptf_ks_clay
    real(dp), dimension(5) :: ptf_ks_curveslope
    real(dp), dimension(5) :: rootfractioncoefficient_forest
    real(dp), dimension(5) :: rootfractioncoefficient_impervious
    real(dp), dimension(5) :: rootfractioncoefficient_pervious
    real(dp), dimension(5) :: infiltrationshapefactor
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /soilmoisture1/ &
      orgmattercontent_forest, &
      orgmattercontent_impervious, &
      orgmattercontent_pervious, &
      ptf_lower66_5_constant, &
      ptf_lower66_5_clay, &
      ptf_lower66_5_db, &
      ptf_higher66_5_constant, &
      ptf_higher66_5_clay, &
      ptf_higher66_5_db, &
      ptf_ks_constant, &
      ptf_ks_sand, &
      ptf_ks_clay, &
      ptf_ks_curveslope, &
      rootfractioncoefficient_forest, &
      rootfractioncoefficient_impervious, &
      rootfractioncoefficient_pervious, &
      infiltrationshapefactor

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    orgmattercontent_forest = this%orgmattercontent_forest
    orgmattercontent_impervious = this%orgmattercontent_impervious
    orgmattercontent_pervious = this%orgmattercontent_pervious
    ptf_lower66_5_constant = this%ptf_lower66_5_constant
    ptf_lower66_5_clay = this%ptf_lower66_5_clay
    ptf_lower66_5_db = this%ptf_lower66_5_db
    ptf_higher66_5_constant = this%ptf_higher66_5_constant
    ptf_higher66_5_clay = this%ptf_higher66_5_clay
    ptf_higher66_5_db = this%ptf_higher66_5_db
    ptf_ks_constant = this%ptf_ks_constant
    ptf_ks_sand = this%ptf_ks_sand
    ptf_ks_clay = this%ptf_ks_clay
    ptf_ks_curveslope = this%ptf_ks_curveslope
    rootfractioncoefficient_forest = this%rootfractioncoefficient_forest
    rootfractioncoefficient_impervious = this%rootfractioncoefficient_impervious
    rootfractioncoefficient_pervious = this%rootfractioncoefficient_pervious
    infiltrationshapefactor = this%infiltrationshapefactor

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("soilmoisture1", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=soilmoisture1, iostat=iostat, iomsg=iomsg)
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
    this%orgmattercontent_forest = orgmattercontent_forest
    this%orgmattercontent_impervious = orgmattercontent_impervious
    this%orgmattercontent_pervious = orgmattercontent_pervious
    this%ptf_lower66_5_constant = ptf_lower66_5_constant
    this%ptf_lower66_5_clay = ptf_lower66_5_clay
    this%ptf_lower66_5_db = ptf_lower66_5_db
    this%ptf_higher66_5_constant = ptf_higher66_5_constant
    this%ptf_higher66_5_clay = ptf_higher66_5_clay
    this%ptf_higher66_5_db = ptf_higher66_5_db
    this%ptf_ks_constant = ptf_ks_constant
    this%ptf_ks_sand = ptf_ks_sand
    this%ptf_ks_clay = ptf_ks_clay
    this%ptf_ks_curveslope = ptf_ks_curveslope
    this%rootfractioncoefficient_forest = rootfractioncoefficient_forest
    this%rootfractioncoefficient_impervious = rootfractioncoefficient_impervious
    this%rootfractioncoefficient_pervious = rootfractioncoefficient_pervious
    this%infiltrationshapefactor = infiltrationshapefactor

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_soilmoisture1_from_file

  !> \brief Set soilmoisture1 values
  integer function nml_soilmoisture1_set(this, &
    orgmattercontent_forest, &
    orgmattercontent_impervious, &
    orgmattercontent_pervious, &
    ptf_lower66_5_constant, &
    ptf_lower66_5_clay, &
    ptf_lower66_5_db, &
    ptf_higher66_5_constant, &
    ptf_higher66_5_clay, &
    ptf_higher66_5_db, &
    ptf_ks_constant, &
    ptf_ks_sand, &
    ptf_ks_clay, &
    rootfractioncoefficient_forest, &
    rootfractioncoefficient_impervious, &
    rootfractioncoefficient_pervious, &
    infiltrationshapefactor, &
    ptf_ks_curveslope, &
    errmsg) result(status)

    class(nml_soilmoisture1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    real(dp), dimension(:), intent(in) :: orgmattercontent_forest !< Organic matter content for forest
    real(dp), dimension(:), intent(in) :: orgmattercontent_impervious !< Organic matter content for impervious
    real(dp), dimension(:), intent(in) :: orgmattercontent_pervious !< Organic matter content for pervious
    real(dp), dimension(:), intent(in) :: ptf_lower66_5_constant !< Zacharias PTF parameters below 66.5 % sand content
    real(dp), dimension(:), intent(in) :: ptf_lower66_5_clay !< Multiplier for clay constant below 66.5 % sand content
    real(dp), dimension(:), intent(in) :: ptf_lower66_5_db !< Multiplier for mineral bulk density below 66.5 % sand content
    real(dp), dimension(:), intent(in) :: ptf_higher66_5_constant !< Zacharias PTF parameters above 66.5 % sand content
    real(dp), dimension(:), intent(in) :: ptf_higher66_5_clay !< Multiplier for clay constant above 66.5 % sand content
    real(dp), dimension(:), intent(in) :: ptf_higher66_5_db !< Multiplier for mineral bulk density above 66.5 % sand content
    real(dp), dimension(:), intent(in) :: ptf_ks_constant !< PTF constant for saturated hydraulic conductivity
    real(dp), dimension(:), intent(in) :: ptf_ks_sand !< Multiplier for sand for saturated hydraulic conductivity
    real(dp), dimension(:), intent(in) :: ptf_ks_clay !< Multiplier for clay for saturated hydraulic conductivity
    real(dp), dimension(:), intent(in) :: rootfractioncoefficient_forest !< Root fraction coefficient for forest
    real(dp), dimension(:), intent(in) :: rootfractioncoefficient_impervious !< Root fraction coefficient for impervious
    real(dp), dimension(:), intent(in) :: rootfractioncoefficient_pervious !< Root fraction coefficient for pervious
    real(dp), dimension(:), intent(in) :: infiltrationshapefactor !< Infiltration shape factor
    real(dp), dimension(:), intent(in), optional :: ptf_ks_curveslope !< Unit conversion factor
    integer :: &
      lb__1, &
      ub__1

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    if (size(orgmattercontent_forest, 1) > size(this%orgmattercontent_forest, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'orgmattercontent_forest'"
      return
    end if
    lb__1 = lbound(this%orgmattercontent_forest, 1)
    ub__1 = lb__1 + size(orgmattercontent_forest, 1) - 1
    this%orgmattercontent_forest(lb__1:ub__1) = orgmattercontent_forest
    if (size(orgmattercontent_impervious, 1) > size(this%orgmattercontent_impervious, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'orgmattercontent_impervious'"
      return
    end if
    lb__1 = lbound(this%orgmattercontent_impervious, 1)
    ub__1 = lb__1 + size(orgmattercontent_impervious, 1) - 1
    this%orgmattercontent_impervious(lb__1:ub__1) = orgmattercontent_impervious
    if (size(orgmattercontent_pervious, 1) > size(this%orgmattercontent_pervious, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'orgmattercontent_pervious'"
      return
    end if
    lb__1 = lbound(this%orgmattercontent_pervious, 1)
    ub__1 = lb__1 + size(orgmattercontent_pervious, 1) - 1
    this%orgmattercontent_pervious(lb__1:ub__1) = orgmattercontent_pervious
    if (size(ptf_lower66_5_constant, 1) > size(this%ptf_lower66_5_constant, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'ptf_lower66_5_constant'"
      return
    end if
    lb__1 = lbound(this%ptf_lower66_5_constant, 1)
    ub__1 = lb__1 + size(ptf_lower66_5_constant, 1) - 1
    this%ptf_lower66_5_constant(lb__1:ub__1) = ptf_lower66_5_constant
    if (size(ptf_lower66_5_clay, 1) > size(this%ptf_lower66_5_clay, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'ptf_lower66_5_clay'"
      return
    end if
    lb__1 = lbound(this%ptf_lower66_5_clay, 1)
    ub__1 = lb__1 + size(ptf_lower66_5_clay, 1) - 1
    this%ptf_lower66_5_clay(lb__1:ub__1) = ptf_lower66_5_clay
    if (size(ptf_lower66_5_db, 1) > size(this%ptf_lower66_5_db, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'ptf_lower66_5_db'"
      return
    end if
    lb__1 = lbound(this%ptf_lower66_5_db, 1)
    ub__1 = lb__1 + size(ptf_lower66_5_db, 1) - 1
    this%ptf_lower66_5_db(lb__1:ub__1) = ptf_lower66_5_db
    if (size(ptf_higher66_5_constant, 1) > size(this%ptf_higher66_5_constant, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'ptf_higher66_5_constant'"
      return
    end if
    lb__1 = lbound(this%ptf_higher66_5_constant, 1)
    ub__1 = lb__1 + size(ptf_higher66_5_constant, 1) - 1
    this%ptf_higher66_5_constant(lb__1:ub__1) = ptf_higher66_5_constant
    if (size(ptf_higher66_5_clay, 1) > size(this%ptf_higher66_5_clay, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'ptf_higher66_5_clay'"
      return
    end if
    lb__1 = lbound(this%ptf_higher66_5_clay, 1)
    ub__1 = lb__1 + size(ptf_higher66_5_clay, 1) - 1
    this%ptf_higher66_5_clay(lb__1:ub__1) = ptf_higher66_5_clay
    if (size(ptf_higher66_5_db, 1) > size(this%ptf_higher66_5_db, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'ptf_higher66_5_db'"
      return
    end if
    lb__1 = lbound(this%ptf_higher66_5_db, 1)
    ub__1 = lb__1 + size(ptf_higher66_5_db, 1) - 1
    this%ptf_higher66_5_db(lb__1:ub__1) = ptf_higher66_5_db
    if (size(ptf_ks_constant, 1) > size(this%ptf_ks_constant, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'ptf_ks_constant'"
      return
    end if
    lb__1 = lbound(this%ptf_ks_constant, 1)
    ub__1 = lb__1 + size(ptf_ks_constant, 1) - 1
    this%ptf_ks_constant(lb__1:ub__1) = ptf_ks_constant
    if (size(ptf_ks_sand, 1) > size(this%ptf_ks_sand, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'ptf_ks_sand'"
      return
    end if
    lb__1 = lbound(this%ptf_ks_sand, 1)
    ub__1 = lb__1 + size(ptf_ks_sand, 1) - 1
    this%ptf_ks_sand(lb__1:ub__1) = ptf_ks_sand
    if (size(ptf_ks_clay, 1) > size(this%ptf_ks_clay, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'ptf_ks_clay'"
      return
    end if
    lb__1 = lbound(this%ptf_ks_clay, 1)
    ub__1 = lb__1 + size(ptf_ks_clay, 1) - 1
    this%ptf_ks_clay(lb__1:ub__1) = ptf_ks_clay
    if (size(rootfractioncoefficient_forest, 1) > size(this%rootfractioncoefficient_forest, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'rootfractioncoefficient_forest'"
      return
    end if
    lb__1 = lbound(this%rootfractioncoefficient_forest, 1)
    ub__1 = lb__1 + size(rootfractioncoefficient_forest, 1) - 1
    this%rootfractioncoefficient_forest(lb__1:ub__1) = rootfractioncoefficient_forest
    if (size(rootfractioncoefficient_impervious, 1) > size(this%rootfractioncoefficient_impervious, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'rootfractioncoefficient_impervious'"
      return
    end if
    lb__1 = lbound(this%rootfractioncoefficient_impervious, 1)
    ub__1 = lb__1 + size(rootfractioncoefficient_impervious, 1) - 1
    this%rootfractioncoefficient_impervious(lb__1:ub__1) = rootfractioncoefficient_impervious
    if (size(rootfractioncoefficient_pervious, 1) > size(this%rootfractioncoefficient_pervious, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'rootfractioncoefficient_pervious'"
      return
    end if
    lb__1 = lbound(this%rootfractioncoefficient_pervious, 1)
    ub__1 = lb__1 + size(rootfractioncoefficient_pervious, 1) - 1
    this%rootfractioncoefficient_pervious(lb__1:ub__1) = rootfractioncoefficient_pervious
    if (size(infiltrationshapefactor, 1) > size(this%infiltrationshapefactor, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'infiltrationshapefactor'"
      return
    end if
    lb__1 = lbound(this%infiltrationshapefactor, 1)
    ub__1 = lb__1 + size(infiltrationshapefactor, 1) - 1
    this%infiltrationshapefactor(lb__1:ub__1) = infiltrationshapefactor
    ! override with provided values
    if (present(ptf_ks_curveslope)) then
      if (size(ptf_ks_curveslope, 1) > size(this%ptf_ks_curveslope, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'ptf_ks_curveslope'"
        return
      end if
      lb__1 = lbound(this%ptf_ks_curveslope, 1)
      ub__1 = lb__1 + size(ptf_ks_curveslope, 1) - 1
      this%ptf_ks_curveslope(lb__1:ub__1) = ptf_ks_curveslope
    end if

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_soilmoisture1_set

  !> \brief Check whether a namelist value was set
  integer function nml_soilmoisture1_is_set(this, name, idx, errmsg) result(status)
    class(nml_soilmoisture1_t), intent(in) :: this !< namelist instance
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
    case ("orgmattercontent_forest")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%orgmattercontent_forest), ubound(this%orgmattercontent_forest), &
          "orgMatterContent_forest", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%orgmattercontent_forest(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%orgmattercontent_forest))) status = NML_ERR_NOT_SET
      end if
    case ("orgmattercontent_impervious")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%orgmattercontent_impervious), ubound(this%orgmattercontent_impervious), &
          "orgMatterContent_impervious", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%orgmattercontent_impervious(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%orgmattercontent_impervious))) status = NML_ERR_NOT_SET
      end if
    case ("orgmattercontent_pervious")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%orgmattercontent_pervious), ubound(this%orgmattercontent_pervious), &
          "orgMatterContent_pervious", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%orgmattercontent_pervious(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%orgmattercontent_pervious))) status = NML_ERR_NOT_SET
      end if
    case ("ptf_lower66_5_constant")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%ptf_lower66_5_constant), ubound(this%ptf_lower66_5_constant), &
          "PTF_lower66_5_constant", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%ptf_lower66_5_constant(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%ptf_lower66_5_constant))) status = NML_ERR_NOT_SET
      end if
    case ("ptf_lower66_5_clay")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%ptf_lower66_5_clay), ubound(this%ptf_lower66_5_clay), &
          "PTF_lower66_5_clay", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%ptf_lower66_5_clay(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%ptf_lower66_5_clay))) status = NML_ERR_NOT_SET
      end if
    case ("ptf_lower66_5_db")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%ptf_lower66_5_db), ubound(this%ptf_lower66_5_db), &
          "PTF_lower66_5_Db", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%ptf_lower66_5_db(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%ptf_lower66_5_db))) status = NML_ERR_NOT_SET
      end if
    case ("ptf_higher66_5_constant")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%ptf_higher66_5_constant), ubound(this%ptf_higher66_5_constant), &
          "PTF_higher66_5_constant", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%ptf_higher66_5_constant(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%ptf_higher66_5_constant))) status = NML_ERR_NOT_SET
      end if
    case ("ptf_higher66_5_clay")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%ptf_higher66_5_clay), ubound(this%ptf_higher66_5_clay), &
          "PTF_higher66_5_clay", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%ptf_higher66_5_clay(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%ptf_higher66_5_clay))) status = NML_ERR_NOT_SET
      end if
    case ("ptf_higher66_5_db")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%ptf_higher66_5_db), ubound(this%ptf_higher66_5_db), &
          "PTF_higher66_5_Db", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%ptf_higher66_5_db(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%ptf_higher66_5_db))) status = NML_ERR_NOT_SET
      end if
    case ("ptf_ks_constant")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%ptf_ks_constant), ubound(this%ptf_ks_constant), &
          "PTF_Ks_constant", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%ptf_ks_constant(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%ptf_ks_constant))) status = NML_ERR_NOT_SET
      end if
    case ("ptf_ks_sand")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%ptf_ks_sand), ubound(this%ptf_ks_sand), &
          "PTF_Ks_sand", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%ptf_ks_sand(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%ptf_ks_sand))) status = NML_ERR_NOT_SET
      end if
    case ("ptf_ks_clay")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%ptf_ks_clay), ubound(this%ptf_ks_clay), &
          "PTF_Ks_clay", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%ptf_ks_clay(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%ptf_ks_clay))) status = NML_ERR_NOT_SET
      end if
    case ("ptf_ks_curveslope")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%ptf_ks_curveslope), ubound(this%ptf_ks_curveslope), &
          "PTF_Ks_curveSlope", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("rootfractioncoefficient_forest")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%rootfractioncoefficient_forest), ubound(this%rootfractioncoefficient_forest), &
          "rootFractionCoefficient_forest", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%rootfractioncoefficient_forest(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%rootfractioncoefficient_forest))) status = NML_ERR_NOT_SET
      end if
    case ("rootfractioncoefficient_impervious")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%rootfractioncoefficient_impervious), ubound(this%rootfractioncoefficient_impervious), &
          "rootFractionCoefficient_impervious", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%rootfractioncoefficient_impervious(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%rootfractioncoefficient_impervious))) status = NML_ERR_NOT_SET
      end if
    case ("rootfractioncoefficient_pervious")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%rootfractioncoefficient_pervious), ubound(this%rootfractioncoefficient_pervious), &
          "rootFractionCoefficient_pervious", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%rootfractioncoefficient_pervious(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%rootfractioncoefficient_pervious))) status = NML_ERR_NOT_SET
      end if
    case ("infiltrationshapefactor")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%infiltrationshapefactor), ubound(this%infiltrationshapefactor), &
          "infiltrationShapeFactor", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%infiltrationshapefactor(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%infiltrationshapefactor))) status = NML_ERR_NOT_SET
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_soilmoisture1_is_set

  !> \brief Validate required values and constraints
  integer function nml_soilmoisture1_is_valid(this, errmsg) result(status)
    class(nml_soilmoisture1_t), intent(in) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    integer :: istat

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (.not. this%is_configured) then
      status = NML_ERR_NOT_SET
      if (present(errmsg)) errmsg = "namelist not configured; call set or from_file"
      return
    end if

    ! required arrays
    if (all(ieee_is_nan(this%orgmattercontent_forest(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: orgMatterContent_forest"
      return
    end if
    if (any(ieee_is_nan(this%orgmattercontent_forest(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: orgMatterContent_forest"
      return
    end if
    if (all(ieee_is_nan(this%orgmattercontent_impervious(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: orgMatterContent_impervious"
      return
    end if
    if (any(ieee_is_nan(this%orgmattercontent_impervious(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: orgMatterContent_impervious"
      return
    end if
    if (all(ieee_is_nan(this%orgmattercontent_pervious(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: orgMatterContent_pervious"
      return
    end if
    if (any(ieee_is_nan(this%orgmattercontent_pervious(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: orgMatterContent_pervious"
      return
    end if
    if (all(ieee_is_nan(this%ptf_lower66_5_constant(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: PTF_lower66_5_constant"
      return
    end if
    if (any(ieee_is_nan(this%ptf_lower66_5_constant(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: PTF_lower66_5_constant"
      return
    end if
    if (all(ieee_is_nan(this%ptf_lower66_5_clay(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: PTF_lower66_5_clay"
      return
    end if
    if (any(ieee_is_nan(this%ptf_lower66_5_clay(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: PTF_lower66_5_clay"
      return
    end if
    if (all(ieee_is_nan(this%ptf_lower66_5_db(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: PTF_lower66_5_Db"
      return
    end if
    if (any(ieee_is_nan(this%ptf_lower66_5_db(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: PTF_lower66_5_Db"
      return
    end if
    if (all(ieee_is_nan(this%ptf_higher66_5_constant(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: PTF_higher66_5_constant"
      return
    end if
    if (any(ieee_is_nan(this%ptf_higher66_5_constant(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: PTF_higher66_5_constant"
      return
    end if
    if (all(ieee_is_nan(this%ptf_higher66_5_clay(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: PTF_higher66_5_clay"
      return
    end if
    if (any(ieee_is_nan(this%ptf_higher66_5_clay(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: PTF_higher66_5_clay"
      return
    end if
    if (all(ieee_is_nan(this%ptf_higher66_5_db(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: PTF_higher66_5_Db"
      return
    end if
    if (any(ieee_is_nan(this%ptf_higher66_5_db(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: PTF_higher66_5_Db"
      return
    end if
    if (all(ieee_is_nan(this%ptf_ks_constant(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: PTF_Ks_constant"
      return
    end if
    if (any(ieee_is_nan(this%ptf_ks_constant(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: PTF_Ks_constant"
      return
    end if
    if (all(ieee_is_nan(this%ptf_ks_sand(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: PTF_Ks_sand"
      return
    end if
    if (any(ieee_is_nan(this%ptf_ks_sand(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: PTF_Ks_sand"
      return
    end if
    if (all(ieee_is_nan(this%ptf_ks_clay(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: PTF_Ks_clay"
      return
    end if
    if (any(ieee_is_nan(this%ptf_ks_clay(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: PTF_Ks_clay"
      return
    end if
    if (all(ieee_is_nan(this%rootfractioncoefficient_forest(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: rootFractionCoefficient_forest"
      return
    end if
    if (any(ieee_is_nan(this%rootfractioncoefficient_forest(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: rootFractionCoefficient_forest"
      return
    end if
    if (all(ieee_is_nan(this%rootfractioncoefficient_impervious(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: rootFractionCoefficient_impervious"
      return
    end if
    if (any(ieee_is_nan(this%rootfractioncoefficient_impervious(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: rootFractionCoefficient_impervious"
      return
    end if
    if (all(ieee_is_nan(this%rootfractioncoefficient_pervious(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: rootFractionCoefficient_pervious"
      return
    end if
    if (any(ieee_is_nan(this%rootfractioncoefficient_pervious(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: rootFractionCoefficient_pervious"
      return
    end if
    if (all(ieee_is_nan(this%infiltrationshapefactor(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: infiltrationShapeFactor"
      return
    end if
    if (any(ieee_is_nan(this%infiltrationshapefactor(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: infiltrationShapeFactor"
      return
    end if
  end function nml_soilmoisture1_is_valid

end module nml_soilmoisture1
