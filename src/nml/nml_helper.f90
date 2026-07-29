!> \dir nml
!> \copydoc f_namelists

!> \defgroup   f_namelists Namelists - Fortran modules
!> \brief      Modules for namelist handling in mHM.
!> \details    These modules provide functionalities to read, write, and manage namelists used in the mHM hydrological model.

!> \file nml_helper.f90
!> \copydoc nml_helper

!> \brief Helper module for namelist file operations
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_helper
  ! kind specifiers used by locally generated derived types
  use mo_kind, only: &
    dp

  !> \brief Buffer length for reading lines
  integer, public :: nml_line_buffer = 1024
  !> \brief Status code: success
  integer, parameter, public :: NML_OK = 0
  !> \brief Status code: file not found
  integer, parameter, public :: NML_ERR_FILE_NOT_FOUND = 1
  !> \brief Status code: failed to open file
  integer, parameter, public :: NML_ERR_OPEN = 2
  !> \brief Status code: file not open
  integer, parameter, public :: NML_ERR_NOT_OPEN = 3
  !> \brief Status code: namelist not found
  integer, parameter, public :: NML_ERR_NML_NOT_FOUND = 4
  !> \brief Status code: read error
  integer, parameter, public :: NML_ERR_READ = 5
  !> \brief Status code: close error
  integer, parameter, public :: NML_ERR_CLOSE = 6
  !> \brief Status code: required value missing
  integer, parameter, public :: NML_ERR_REQUIRED = 10
  !> \brief Status code: enum validation failed
  integer, parameter, public :: NML_ERR_ENUM = 11
  !> \brief Status code: value not set
  integer, parameter, public :: NML_ERR_NOT_SET = 12
  !> \brief Status code: array partially set
  integer, parameter, public :: NML_ERR_PARTLY_SET = 13
  !> \brief Status code: bounds validation failed
  integer, parameter, public :: NML_ERR_BOUNDS = 14
  !> \brief Status code: invalid field name
  integer, parameter, public :: NML_ERR_INVALID_NAME = 20
  !> \brief Status code: invalid index
  integer, parameter, public :: NML_ERR_INVALID_INDEX = 21
  !> \brief Status code: zero opaque handle
  integer, parameter, public :: NML_ERR_INVALID_HANDLE = 22

  !> \brief Shared constants for generated namelist modules
  integer, parameter, public :: buf = 256 !< String buffer length.
  integer, parameter, public :: max_layers__default = 10 !< Maximum number of soil-layer entries stored per domain in a namelist file.
  integer, parameter, public :: n_domains__default = 100 !< Number of domains.
  integer, parameter, public :: n_geo_units__default = 25 !< Number of geological units.

  !> \class parameter_t
  !> \brief Calibration parameter
  !> \details A model parameter with optional calibration metadata.
  type, public :: parameter_t
    real(dp) :: value !< Parameter value
    logical :: optimize !< Include parameter in calibration
    real(dp) :: lower_bound !< Lower calibration bound
    real(dp) :: upper_bound !< Upper calibration bound
  end type parameter_t

  !> \class nml_file_t
  !> \brief Type for namelist file operations
  type, public :: nml_file_t
    integer :: unit = 0
    logical :: is_open = .false.
  contains
    procedure :: open => nml_open
    procedure :: find => nml_find
    procedure :: close => nml_close
  end type nml_file_t

contains

  !> \brief Open a namelist file
  integer function nml_open(this, file, errmsg) result(status)
    class(nml_file_t), intent(inout) :: this
    character(len=*), intent(in) :: file
    character(len=*), intent(out), optional :: errmsg
    integer :: iostat
    logical :: exists
    character(len=nml_line_buffer) :: iomsg
    if (present(errmsg)) errmsg = ""
    status = this%close()
    inquire(file=file, exist=exists)
    if (.not. exists) then
      this%is_open = .false.
      this%unit = 0
      status = NML_ERR_FILE_NOT_FOUND
      if (present(errmsg)) errmsg = "file not found: " // trim(file)
      return
    end if
    open(newunit=this%unit, file=file, status='old', action='read', &
      iostat=iostat, iomsg=iomsg)
    this%is_open = (iostat == 0)
    if (.not. this%is_open) then
      this%unit = 0
      status = NML_ERR_OPEN
      if (present(errmsg)) errmsg = trim(iomsg)
      return
    end if
    status = NML_OK
  end function nml_open

  !> \brief Find a namelist in the opened file
  integer function nml_find(this, nml, errmsg) result(status)
    class(nml_file_t), intent(inout) :: this
    character(len=*), intent(in) :: nml
    character(len=*), intent(out), optional :: errmsg
    integer :: iostat
    character(len=nml_line_buffer) :: line
    character(len=nml_line_buffer) :: iomsg
    if (present(errmsg)) errmsg = ""
    status = NML_ERR_NML_NOT_FOUND
    if (.not. this%is_open) then
      status = NML_ERR_NOT_OPEN
      if (present(errmsg)) errmsg = "file not open"
      return
    end if
    rewind(unit=this%unit)
    do
      read(this%unit, '(A)', iostat=iostat, iomsg=iomsg) line
      if (iostat < 0) exit
      if (iostat > 0) then
        status = NML_ERR_READ
        if (present(errmsg)) errmsg = trim(iomsg)
        return
      end if
      if (index(to_lower(line), '&' // to_lower(trim(nml))) /= 0) then
        status = NML_OK
        backspace(this%unit)
        return
      end if
    end do
    if (present(errmsg)) errmsg = "namelist not found: " // trim(nml)
  end function nml_find

  !> \brief Close the namelist file
  integer function nml_close(this, errmsg) result(status)
    class(nml_file_t), intent(inout) :: this
    character(len=*), intent(out), optional :: errmsg
    integer :: iostat
    character(len=nml_line_buffer) :: iomsg
    if (present(errmsg)) errmsg = ""
    status = NML_OK
    if (this%is_open) then
      close(unit=this%unit, iostat=iostat, iomsg=iomsg)
      if (iostat /= 0) then
        status = NML_ERR_CLOSE
        if (present(errmsg)) errmsg = trim(iomsg)
      end if
    end if
    this%is_open = .false.
    this%unit = 0
  end function nml_close

  !> \brief Convert string to lower case
  pure function to_lower(string) result(lower_string)
    character(len=*), intent(in) :: string
    character(len=len(string)) :: lower_string
    integer, parameter :: shift=iachar('a')-iachar('A'), upA=iachar('A'), upZ=iachar('Z')
    integer :: k, i
    do i = 1, len(string)
      k = ichar(string(i:i))
      if (k>=upA .and. k<=upZ) k = k + shift
      lower_string(i:i) = char(k)
    end do
  end function to_lower

  !> \brief Validate index bounds for array access
  integer function idx_check(idx, lower, upper, field, errmsg) result(status)
    integer, intent(in) :: idx(:)
    integer, intent(in) :: lower(:)
    integer, intent(in) :: upper(:)
    character(len=*), intent(in) :: field
    character(len=*), intent(out), optional :: errmsg

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (size(idx) /= size(lower)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "index rank mismatch for '" // trim(field) // "'"
    else if (any(idx < lower) .or. any(idx > upper)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "index out of bounds for '" // trim(field) // "'"
    end if
  end function idx_check

end module nml_helper
