!> \file nml_parameter_percolation1.f90
!> \copydoc nml_percolation1

!> \brief Percolation - Case 1
!> \details Parameters for percolation.
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    Jan 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_percolation1
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

  !> \class nml_percolation1_t
  !> \brief Percolation - Case 1
  !> \details Parameters for percolation.
  type, public :: nml_percolation1_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    real(dp), dimension(5) :: rechargecoefficient !< Recharge coefficient
    real(dp), dimension(5) :: rechargefactor_karstic !< Recharge factor karstic
    real(dp), dimension(5) :: gain_loss_gwreservoir_karstic !< Gain/loss GW reservoir karstic
  contains
    procedure :: init => nml_percolation1_init
    procedure :: from_file => nml_percolation1_from_file
    procedure :: set => nml_percolation1_set
    procedure :: is_set => nml_percolation1_is_set
    procedure :: is_valid => nml_percolation1_is_valid
  end type nml_percolation1_t

contains

  !> \brief Initialize defaults and sentinels for percolation1
  integer function nml_percolation1_init(this, errmsg) result(status)
    class(nml_percolation1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! sentinel values for required/optional parameters
    this%rechargecoefficient = ieee_value(this%rechargecoefficient, ieee_quiet_nan) ! sentinel for required real array
    this%rechargefactor_karstic = ieee_value(this%rechargefactor_karstic, ieee_quiet_nan) ! sentinel for required real array
    this%gain_loss_gwreservoir_karstic = ieee_value(this%gain_loss_gwreservoir_karstic, ieee_quiet_nan) ! sentinel for required real array
  end function nml_percolation1_init


  !> \brief Read percolation1 namelist from file
  integer function nml_percolation1_from_file(this, file, errmsg) result(status)
    class(nml_percolation1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    real(dp), dimension(5) :: rechargecoefficient
    real(dp), dimension(5) :: rechargefactor_karstic
    real(dp), dimension(5) :: gain_loss_gwreservoir_karstic
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /percolation1/ &
      rechargecoefficient, &
      rechargefactor_karstic, &
      gain_loss_gwreservoir_karstic

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    rechargecoefficient = this%rechargecoefficient
    rechargefactor_karstic = this%rechargefactor_karstic
    gain_loss_gwreservoir_karstic = this%gain_loss_gwreservoir_karstic

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("percolation1", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=percolation1, iostat=iostat, iomsg=iomsg)
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
    this%rechargecoefficient = rechargecoefficient
    this%rechargefactor_karstic = rechargefactor_karstic
    this%gain_loss_gwreservoir_karstic = gain_loss_gwreservoir_karstic

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_percolation1_from_file

  !> \brief Set percolation1 values
  integer function nml_percolation1_set(this, &
    rechargecoefficient, &
    rechargefactor_karstic, &
    gain_loss_gwreservoir_karstic, &
    errmsg) result(status)

    class(nml_percolation1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    real(dp), dimension(:), intent(in) :: rechargecoefficient !< Recharge coefficient
    real(dp), dimension(:), intent(in) :: rechargefactor_karstic !< Recharge factor karstic
    real(dp), dimension(:), intent(in) :: gain_loss_gwreservoir_karstic !< Gain/loss GW reservoir karstic
    integer :: &
      lb__1, &
      ub__1

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    if (size(rechargecoefficient, 1) > size(this%rechargecoefficient, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'rechargecoefficient'"
      return
    end if
    lb__1 = lbound(this%rechargecoefficient, 1)
    ub__1 = lb__1 + size(rechargecoefficient, 1) - 1
    this%rechargecoefficient(lb__1:ub__1) = rechargecoefficient
    if (size(rechargefactor_karstic, 1) > size(this%rechargefactor_karstic, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'rechargefactor_karstic'"
      return
    end if
    lb__1 = lbound(this%rechargefactor_karstic, 1)
    ub__1 = lb__1 + size(rechargefactor_karstic, 1) - 1
    this%rechargefactor_karstic(lb__1:ub__1) = rechargefactor_karstic
    if (size(gain_loss_gwreservoir_karstic, 1) > size(this%gain_loss_gwreservoir_karstic, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'gain_loss_gwreservoir_karstic'"
      return
    end if
    lb__1 = lbound(this%gain_loss_gwreservoir_karstic, 1)
    ub__1 = lb__1 + size(gain_loss_gwreservoir_karstic, 1) - 1
    this%gain_loss_gwreservoir_karstic(lb__1:ub__1) = gain_loss_gwreservoir_karstic

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_percolation1_set

  !> \brief Check whether a namelist value was set
  integer function nml_percolation1_is_set(this, name, idx, errmsg) result(status)
    class(nml_percolation1_t), intent(in) :: this !< namelist instance
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
    case ("rechargecoefficient")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%rechargecoefficient), ubound(this%rechargecoefficient), &
          "rechargeCoefficient", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%rechargecoefficient(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%rechargecoefficient))) status = NML_ERR_NOT_SET
      end if
    case ("rechargefactor_karstic")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%rechargefactor_karstic), ubound(this%rechargefactor_karstic), &
          "rechargeFactor_karstic", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%rechargefactor_karstic(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%rechargefactor_karstic))) status = NML_ERR_NOT_SET
      end if
    case ("gain_loss_gwreservoir_karstic")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%gain_loss_gwreservoir_karstic), ubound(this%gain_loss_gwreservoir_karstic), &
          "gain_loss_GWreservoir_karstic", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%gain_loss_gwreservoir_karstic(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%gain_loss_gwreservoir_karstic))) status = NML_ERR_NOT_SET
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_percolation1_is_set

  !> \brief Validate required values and constraints
  integer function nml_percolation1_is_valid(this, errmsg) result(status)
    class(nml_percolation1_t), intent(in) :: this !< namelist instance
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
    if (all(ieee_is_nan(this%rechargecoefficient(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: rechargeCoefficient"
      return
    end if
    if (any(ieee_is_nan(this%rechargecoefficient(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: rechargeCoefficient"
      return
    end if
    if (all(ieee_is_nan(this%rechargefactor_karstic(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: rechargeFactor_karstic"
      return
    end if
    if (any(ieee_is_nan(this%rechargefactor_karstic(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: rechargeFactor_karstic"
      return
    end if
    if (all(ieee_is_nan(this%gain_loss_gwreservoir_karstic(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: gain_loss_GWreservoir_karstic"
      return
    end if
    if (any(ieee_is_nan(this%gain_loss_gwreservoir_karstic(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: gain_loss_GWreservoir_karstic"
      return
    end if
  end function nml_percolation1_is_valid

end module nml_percolation1
