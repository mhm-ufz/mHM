!> \file nml_parameter_neutrons2.f90
!> \copydoc nml_neutrons2

!> \brief Neutrons - Case 2
!> \details Parameters for neutrons (case 2 - COSMIC).
!! Ground albedo neutrons - COSMIC version.
!! THIS IS WORK IN PROGRESS, DO NOT USE FOR RESEARCH
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    Jan 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_neutrons2
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

  !> \class nml_neutrons2_t
  !> \brief Neutrons - Case 2
  !> \details Parameters for neutrons (case 2 - COSMIC).
  !! Ground albedo neutrons - COSMIC version.
  !! THIS IS WORK IN PROGRESS, DO NOT USE FOR RESEARCH
  type, public :: nml_neutrons2_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    real(dp), dimension(5) :: cosmic_n0 !< Cosmic N0 parameter
    real(dp), dimension(5) :: cosmic_n1 !< Cosmic N1 parameter
    real(dp), dimension(5) :: cosmic_n2 !< Cosmic N2 parameter
    real(dp), dimension(5) :: cosmic_alpha0 !< Cosmic alpha0 parameter
    real(dp), dimension(5) :: cosmic_alpha1 !< Cosmic alpha1 parameter
    real(dp), dimension(5) :: cosmic_l30 !< Cosmic L30 parameter
    real(dp), dimension(5) :: cosmic_l31 !< Cosmic L31 parameter
    real(dp), dimension(5) :: cosmic_lw0 !< Cosmic LW0 parameter
    real(dp), dimension(5) :: cosmic_lw1 !< Cosmic LW1 parameter
  contains
    procedure :: init => nml_neutrons2_init
    procedure :: from_file => nml_neutrons2_from_file
    procedure :: set => nml_neutrons2_set
    procedure :: is_set => nml_neutrons2_is_set
    procedure :: is_valid => nml_neutrons2_is_valid
  end type nml_neutrons2_t

contains

  !> \brief Initialize defaults and sentinels for neutrons2
  integer function nml_neutrons2_init(this, errmsg) result(status)
    class(nml_neutrons2_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! sentinel values for required/optional parameters
    this%cosmic_n0 = ieee_value(this%cosmic_n0, ieee_quiet_nan) ! sentinel for required real array
    this%cosmic_n1 = ieee_value(this%cosmic_n1, ieee_quiet_nan) ! sentinel for required real array
    this%cosmic_n2 = ieee_value(this%cosmic_n2, ieee_quiet_nan) ! sentinel for required real array
    this%cosmic_alpha0 = ieee_value(this%cosmic_alpha0, ieee_quiet_nan) ! sentinel for required real array
    this%cosmic_alpha1 = ieee_value(this%cosmic_alpha1, ieee_quiet_nan) ! sentinel for required real array
    this%cosmic_l30 = ieee_value(this%cosmic_l30, ieee_quiet_nan) ! sentinel for required real array
    this%cosmic_l31 = ieee_value(this%cosmic_l31, ieee_quiet_nan) ! sentinel for required real array
    this%cosmic_lw0 = ieee_value(this%cosmic_lw0, ieee_quiet_nan) ! sentinel for required real array
    this%cosmic_lw1 = ieee_value(this%cosmic_lw1, ieee_quiet_nan) ! sentinel for required real array
  end function nml_neutrons2_init


  !> \brief Read neutrons2 namelist from file
  integer function nml_neutrons2_from_file(this, file, errmsg) result(status)
    class(nml_neutrons2_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    real(dp), dimension(5) :: cosmic_n0
    real(dp), dimension(5) :: cosmic_n1
    real(dp), dimension(5) :: cosmic_n2
    real(dp), dimension(5) :: cosmic_alpha0
    real(dp), dimension(5) :: cosmic_alpha1
    real(dp), dimension(5) :: cosmic_l30
    real(dp), dimension(5) :: cosmic_l31
    real(dp), dimension(5) :: cosmic_lw0
    real(dp), dimension(5) :: cosmic_lw1
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /neutrons2/ &
      cosmic_n0, &
      cosmic_n1, &
      cosmic_n2, &
      cosmic_alpha0, &
      cosmic_alpha1, &
      cosmic_l30, &
      cosmic_l31, &
      cosmic_lw0, &
      cosmic_lw1

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    cosmic_n0 = this%cosmic_n0
    cosmic_n1 = this%cosmic_n1
    cosmic_n2 = this%cosmic_n2
    cosmic_alpha0 = this%cosmic_alpha0
    cosmic_alpha1 = this%cosmic_alpha1
    cosmic_l30 = this%cosmic_l30
    cosmic_l31 = this%cosmic_l31
    cosmic_lw0 = this%cosmic_lw0
    cosmic_lw1 = this%cosmic_lw1

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("neutrons2", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=neutrons2, iostat=iostat, iomsg=iomsg)
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
    this%cosmic_n0 = cosmic_n0
    this%cosmic_n1 = cosmic_n1
    this%cosmic_n2 = cosmic_n2
    this%cosmic_alpha0 = cosmic_alpha0
    this%cosmic_alpha1 = cosmic_alpha1
    this%cosmic_l30 = cosmic_l30
    this%cosmic_l31 = cosmic_l31
    this%cosmic_lw0 = cosmic_lw0
    this%cosmic_lw1 = cosmic_lw1

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_neutrons2_from_file

  !> \brief Set neutrons2 values
  integer function nml_neutrons2_set(this, &
    cosmic_n0, &
    cosmic_n1, &
    cosmic_n2, &
    cosmic_alpha0, &
    cosmic_alpha1, &
    cosmic_l30, &
    cosmic_l31, &
    cosmic_lw0, &
    cosmic_lw1, &
    errmsg) result(status)

    class(nml_neutrons2_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    real(dp), dimension(:), intent(in) :: cosmic_n0 !< Cosmic N0 parameter
    real(dp), dimension(:), intent(in) :: cosmic_n1 !< Cosmic N1 parameter
    real(dp), dimension(:), intent(in) :: cosmic_n2 !< Cosmic N2 parameter
    real(dp), dimension(:), intent(in) :: cosmic_alpha0 !< Cosmic alpha0 parameter
    real(dp), dimension(:), intent(in) :: cosmic_alpha1 !< Cosmic alpha1 parameter
    real(dp), dimension(:), intent(in) :: cosmic_l30 !< Cosmic L30 parameter
    real(dp), dimension(:), intent(in) :: cosmic_l31 !< Cosmic L31 parameter
    real(dp), dimension(:), intent(in) :: cosmic_lw0 !< Cosmic LW0 parameter
    real(dp), dimension(:), intent(in) :: cosmic_lw1 !< Cosmic LW1 parameter
    integer :: &
      lb__1, &
      ub__1

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    if (size(cosmic_n0, 1) > size(this%cosmic_n0, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'cosmic_n0'"
      return
    end if
    lb__1 = lbound(this%cosmic_n0, 1)
    ub__1 = lb__1 + size(cosmic_n0, 1) - 1
    this%cosmic_n0(lb__1:ub__1) = cosmic_n0
    if (size(cosmic_n1, 1) > size(this%cosmic_n1, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'cosmic_n1'"
      return
    end if
    lb__1 = lbound(this%cosmic_n1, 1)
    ub__1 = lb__1 + size(cosmic_n1, 1) - 1
    this%cosmic_n1(lb__1:ub__1) = cosmic_n1
    if (size(cosmic_n2, 1) > size(this%cosmic_n2, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'cosmic_n2'"
      return
    end if
    lb__1 = lbound(this%cosmic_n2, 1)
    ub__1 = lb__1 + size(cosmic_n2, 1) - 1
    this%cosmic_n2(lb__1:ub__1) = cosmic_n2
    if (size(cosmic_alpha0, 1) > size(this%cosmic_alpha0, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'cosmic_alpha0'"
      return
    end if
    lb__1 = lbound(this%cosmic_alpha0, 1)
    ub__1 = lb__1 + size(cosmic_alpha0, 1) - 1
    this%cosmic_alpha0(lb__1:ub__1) = cosmic_alpha0
    if (size(cosmic_alpha1, 1) > size(this%cosmic_alpha1, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'cosmic_alpha1'"
      return
    end if
    lb__1 = lbound(this%cosmic_alpha1, 1)
    ub__1 = lb__1 + size(cosmic_alpha1, 1) - 1
    this%cosmic_alpha1(lb__1:ub__1) = cosmic_alpha1
    if (size(cosmic_l30, 1) > size(this%cosmic_l30, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'cosmic_l30'"
      return
    end if
    lb__1 = lbound(this%cosmic_l30, 1)
    ub__1 = lb__1 + size(cosmic_l30, 1) - 1
    this%cosmic_l30(lb__1:ub__1) = cosmic_l30
    if (size(cosmic_l31, 1) > size(this%cosmic_l31, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'cosmic_l31'"
      return
    end if
    lb__1 = lbound(this%cosmic_l31, 1)
    ub__1 = lb__1 + size(cosmic_l31, 1) - 1
    this%cosmic_l31(lb__1:ub__1) = cosmic_l31
    if (size(cosmic_lw0, 1) > size(this%cosmic_lw0, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'cosmic_lw0'"
      return
    end if
    lb__1 = lbound(this%cosmic_lw0, 1)
    ub__1 = lb__1 + size(cosmic_lw0, 1) - 1
    this%cosmic_lw0(lb__1:ub__1) = cosmic_lw0
    if (size(cosmic_lw1, 1) > size(this%cosmic_lw1, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'cosmic_lw1'"
      return
    end if
    lb__1 = lbound(this%cosmic_lw1, 1)
    ub__1 = lb__1 + size(cosmic_lw1, 1) - 1
    this%cosmic_lw1(lb__1:ub__1) = cosmic_lw1

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_neutrons2_set

  !> \brief Check whether a namelist value was set
  integer function nml_neutrons2_is_set(this, name, idx, errmsg) result(status)
    class(nml_neutrons2_t), intent(in) :: this !< namelist instance
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
    case ("cosmic_n0")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%cosmic_n0), ubound(this%cosmic_n0), &
          "COSMIC_N0", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%cosmic_n0(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%cosmic_n0))) status = NML_ERR_NOT_SET
      end if
    case ("cosmic_n1")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%cosmic_n1), ubound(this%cosmic_n1), &
          "COSMIC_N1", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%cosmic_n1(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%cosmic_n1))) status = NML_ERR_NOT_SET
      end if
    case ("cosmic_n2")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%cosmic_n2), ubound(this%cosmic_n2), &
          "COSMIC_N2", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%cosmic_n2(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%cosmic_n2))) status = NML_ERR_NOT_SET
      end if
    case ("cosmic_alpha0")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%cosmic_alpha0), ubound(this%cosmic_alpha0), &
          "COSMIC_alpha0", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%cosmic_alpha0(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%cosmic_alpha0))) status = NML_ERR_NOT_SET
      end if
    case ("cosmic_alpha1")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%cosmic_alpha1), ubound(this%cosmic_alpha1), &
          "COSMIC_alpha1", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%cosmic_alpha1(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%cosmic_alpha1))) status = NML_ERR_NOT_SET
      end if
    case ("cosmic_l30")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%cosmic_l30), ubound(this%cosmic_l30), &
          "COSMIC_L30", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%cosmic_l30(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%cosmic_l30))) status = NML_ERR_NOT_SET
      end if
    case ("cosmic_l31")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%cosmic_l31), ubound(this%cosmic_l31), &
          "COSMIC_L31", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%cosmic_l31(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%cosmic_l31))) status = NML_ERR_NOT_SET
      end if
    case ("cosmic_lw0")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%cosmic_lw0), ubound(this%cosmic_lw0), &
          "COSMIC_LW0", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%cosmic_lw0(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%cosmic_lw0))) status = NML_ERR_NOT_SET
      end if
    case ("cosmic_lw1")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%cosmic_lw1), ubound(this%cosmic_lw1), &
          "COSMIC_LW1", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%cosmic_lw1(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%cosmic_lw1))) status = NML_ERR_NOT_SET
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_neutrons2_is_set

  !> \brief Validate required values and constraints
  integer function nml_neutrons2_is_valid(this, errmsg) result(status)
    class(nml_neutrons2_t), intent(in) :: this !< namelist instance
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
    if (all(ieee_is_nan(this%cosmic_n0(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: COSMIC_N0"
      return
    end if
    if (any(ieee_is_nan(this%cosmic_n0(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: COSMIC_N0"
      return
    end if
    if (all(ieee_is_nan(this%cosmic_n1(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: COSMIC_N1"
      return
    end if
    if (any(ieee_is_nan(this%cosmic_n1(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: COSMIC_N1"
      return
    end if
    if (all(ieee_is_nan(this%cosmic_n2(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: COSMIC_N2"
      return
    end if
    if (any(ieee_is_nan(this%cosmic_n2(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: COSMIC_N2"
      return
    end if
    if (all(ieee_is_nan(this%cosmic_alpha0(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: COSMIC_alpha0"
      return
    end if
    if (any(ieee_is_nan(this%cosmic_alpha0(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: COSMIC_alpha0"
      return
    end if
    if (all(ieee_is_nan(this%cosmic_alpha1(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: COSMIC_alpha1"
      return
    end if
    if (any(ieee_is_nan(this%cosmic_alpha1(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: COSMIC_alpha1"
      return
    end if
    if (all(ieee_is_nan(this%cosmic_l30(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: COSMIC_L30"
      return
    end if
    if (any(ieee_is_nan(this%cosmic_l30(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: COSMIC_L30"
      return
    end if
    if (all(ieee_is_nan(this%cosmic_l31(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: COSMIC_L31"
      return
    end if
    if (any(ieee_is_nan(this%cosmic_l31(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: COSMIC_L31"
      return
    end if
    if (all(ieee_is_nan(this%cosmic_lw0(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: COSMIC_LW0"
      return
    end if
    if (any(ieee_is_nan(this%cosmic_lw0(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: COSMIC_LW0"
      return
    end if
    if (all(ieee_is_nan(this%cosmic_lw1(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: COSMIC_LW1"
      return
    end if
    if (any(ieee_is_nan(this%cosmic_lw1(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: COSMIC_LW1"
      return
    end if
  end function nml_neutrons2_is_valid

end module nml_neutrons2
