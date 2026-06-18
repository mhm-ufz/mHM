!> \file nml_parameter_petm2.f90
!> \copydoc nml_petm2

!> \brief PET - Case -2
!> \details Parameters for PET (case -2 - aspect correction).
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    Jan 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_petm2
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

  !> \class nml_petm2_t
  !> \brief PET - Case -2
  !> \details Parameters for PET (case -2 - aspect correction).
  type, public :: nml_petm2_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    real(dp), dimension(5) :: mincorrectionfactorpet !< minimum factor for PET correction with aspect
    real(dp), dimension(5) :: maxcorrectionfactorpet !< maximum factor for PET correction with aspect
    real(dp), dimension(5) :: aspecttresholdpet !< aspect threshold
  contains
    procedure :: init => nml_petm2_init
    procedure :: from_file => nml_petm2_from_file
    procedure :: set => nml_petm2_set
    procedure :: is_set => nml_petm2_is_set
    procedure :: is_valid => nml_petm2_is_valid
  end type nml_petm2_t

contains

  !> \brief Initialize defaults and sentinels for petm2
  integer function nml_petm2_init(this, errmsg) result(status)
    class(nml_petm2_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! sentinel values for required/optional parameters
    this%mincorrectionfactorpet = ieee_value(this%mincorrectionfactorpet, ieee_quiet_nan) ! sentinel for required real array
    this%maxcorrectionfactorpet = ieee_value(this%maxcorrectionfactorpet, ieee_quiet_nan) ! sentinel for required real array
    this%aspecttresholdpet = ieee_value(this%aspecttresholdpet, ieee_quiet_nan) ! sentinel for required real array
  end function nml_petm2_init


  !> \brief Read petm2 namelist from file
  integer function nml_petm2_from_file(this, file, errmsg) result(status)
    class(nml_petm2_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    real(dp), dimension(5) :: mincorrectionfactorpet
    real(dp), dimension(5) :: maxcorrectionfactorpet
    real(dp), dimension(5) :: aspecttresholdpet
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /petm2/ &
      mincorrectionfactorpet, &
      maxcorrectionfactorpet, &
      aspecttresholdpet

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    mincorrectionfactorpet = this%mincorrectionfactorpet
    maxcorrectionfactorpet = this%maxcorrectionfactorpet
    aspecttresholdpet = this%aspecttresholdpet

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("petm2", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=petm2, iostat=iostat, iomsg=iomsg)
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
    this%mincorrectionfactorpet = mincorrectionfactorpet
    this%maxcorrectionfactorpet = maxcorrectionfactorpet
    this%aspecttresholdpet = aspecttresholdpet

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_petm2_from_file

  !> \brief Set petm2 values
  integer function nml_petm2_set(this, &
    mincorrectionfactorpet, &
    maxcorrectionfactorpet, &
    aspecttresholdpet, &
    errmsg) result(status)

    class(nml_petm2_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    real(dp), dimension(:), intent(in) :: mincorrectionfactorpet !< minimum factor for PET correction with aspect
    real(dp), dimension(:), intent(in) :: maxcorrectionfactorpet !< maximum factor for PET correction with aspect
    real(dp), dimension(:), intent(in) :: aspecttresholdpet !< aspect threshold
    integer :: &
      lb__1, &
      ub__1

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    if (size(mincorrectionfactorpet, 1) > size(this%mincorrectionfactorpet, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'mincorrectionfactorpet'"
      return
    end if
    lb__1 = lbound(this%mincorrectionfactorpet, 1)
    ub__1 = lb__1 + size(mincorrectionfactorpet, 1) - 1
    this%mincorrectionfactorpet(lb__1:ub__1) = mincorrectionfactorpet
    if (size(maxcorrectionfactorpet, 1) > size(this%maxcorrectionfactorpet, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'maxcorrectionfactorpet'"
      return
    end if
    lb__1 = lbound(this%maxcorrectionfactorpet, 1)
    ub__1 = lb__1 + size(maxcorrectionfactorpet, 1) - 1
    this%maxcorrectionfactorpet(lb__1:ub__1) = maxcorrectionfactorpet
    if (size(aspecttresholdpet, 1) > size(this%aspecttresholdpet, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'aspecttresholdpet'"
      return
    end if
    lb__1 = lbound(this%aspecttresholdpet, 1)
    ub__1 = lb__1 + size(aspecttresholdpet, 1) - 1
    this%aspecttresholdpet(lb__1:ub__1) = aspecttresholdpet

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_petm2_set

  !> \brief Check whether a namelist value was set
  integer function nml_petm2_is_set(this, name, idx, errmsg) result(status)
    class(nml_petm2_t), intent(in) :: this !< namelist instance
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
    case ("mincorrectionfactorpet")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%mincorrectionfactorpet), ubound(this%mincorrectionfactorpet), &
          "minCorrectionFactorPET", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%mincorrectionfactorpet(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%mincorrectionfactorpet))) status = NML_ERR_NOT_SET
      end if
    case ("maxcorrectionfactorpet")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%maxcorrectionfactorpet), ubound(this%maxcorrectionfactorpet), &
          "maxCorrectionFactorPET", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%maxcorrectionfactorpet(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%maxcorrectionfactorpet))) status = NML_ERR_NOT_SET
      end if
    case ("aspecttresholdpet")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%aspecttresholdpet), ubound(this%aspecttresholdpet), &
          "aspectTresholdPET", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%aspecttresholdpet(idx(1)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%aspecttresholdpet))) status = NML_ERR_NOT_SET
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_petm2_is_set

  !> \brief Validate required values and constraints
  integer function nml_petm2_is_valid(this, errmsg) result(status)
    class(nml_petm2_t), intent(in) :: this !< namelist instance
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
    if (all(ieee_is_nan(this%mincorrectionfactorpet(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: minCorrectionFactorPET"
      return
    end if
    if (any(ieee_is_nan(this%mincorrectionfactorpet(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: minCorrectionFactorPET"
      return
    end if
    if (all(ieee_is_nan(this%maxcorrectionfactorpet(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: maxCorrectionFactorPET"
      return
    end if
    if (any(ieee_is_nan(this%maxcorrectionfactorpet(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: maxCorrectionFactorPET"
      return
    end if
    if (all(ieee_is_nan(this%aspecttresholdpet(:)))) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) errmsg = "required field not set: aspectTresholdPET"
      return
    end if
    if (any(ieee_is_nan(this%aspecttresholdpet(:)))) then
      status = NML_ERR_PARTLY_SET
      if (present(errmsg)) errmsg = "array partly set: aspectTresholdPET"
      return
    end if
  end function nml_petm2_is_valid

end module nml_petm2
