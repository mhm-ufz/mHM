!> \file nml_parameter_neutrons_2.f90
!> \copydoc nml_neutrons_2

!> \brief Neutrons - Case 2
!> \details Parameters for the experimental COSMIC neutron formulation.
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_neutrons_2
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
    parameter_t
  use ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
  ! kind specifiers listed in the nml-tools configuration file
  use mo_kind, only: &
    dp

  implicit none

  !> \class nml_neutrons_2_t
  !> \brief Neutrons - Case 2
  !> \details Parameters for the experimental COSMIC neutron formulation.
  type, public :: nml_neutrons_2_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    type(parameter_t) :: cosmic_n0 !< COSMIC N0 parameter
    type(parameter_t) :: cosmic_n1 !< COSMIC N1 parameter
    type(parameter_t) :: cosmic_n2 !< COSMIC N2 parameter
    type(parameter_t) :: cosmic_alpha0 !< COSMIC alpha0 parameter
    type(parameter_t) :: cosmic_alpha1 !< COSMIC alpha1 parameter
    type(parameter_t) :: cosmic_l30 !< COSMIC L30 parameter
    type(parameter_t) :: cosmic_l31 !< COSMIC L31 parameter
    type(parameter_t) :: cosmic_lw0 !< COSMIC lattice-water parameter LW0
    type(parameter_t) :: cosmic_lw1 !< COSMIC lattice-water parameter LW1
  contains
    procedure :: init => nml_neutrons_2_init
    procedure :: init_type => nml_neutrons_2_init_type
    procedure :: from_file => nml_neutrons_2_from_file
    procedure :: set => nml_neutrons_2_set
    procedure :: is_set => nml_neutrons_2_is_set
    procedure :: is_valid => nml_neutrons_2_is_valid
  end type nml_neutrons_2_t

contains

  !> \brief Initialize defaults and sentinels for neutrons_2
  integer function nml_neutrons_2_init(this, errmsg) result(status)
    class(nml_neutrons_2_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! derived values
    status = this%init_type( &
      cosmic_n0=this%cosmic_n0, &
      cosmic_n1=this%cosmic_n1, &
      cosmic_n2=this%cosmic_n2, &
      cosmic_alpha0=this%cosmic_alpha0, &
      cosmic_alpha1=this%cosmic_alpha1, &
      cosmic_l30=this%cosmic_l30, &
      cosmic_l31=this%cosmic_l31, &
      cosmic_lw0=this%cosmic_lw0, &
      cosmic_lw1=this%cosmic_lw1, &
      errmsg=errmsg)
    if (status /= NML_OK) return
  end function nml_neutrons_2_init

  !> \brief Initialize derived values with their field-specific defaults
  integer function nml_neutrons_2_init_type(this, &
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
    class(nml_neutrons_2_t), intent(in) :: this !< parent namelist instance
    type(parameter_t), intent(inout), optional :: cosmic_n0 !< COSMIC N0 parameter
    type(parameter_t), intent(inout), optional :: cosmic_n1 !< COSMIC N1 parameter
    type(parameter_t), intent(inout), optional :: cosmic_n2 !< COSMIC N2 parameter
    type(parameter_t), intent(inout), optional :: cosmic_alpha0 !< COSMIC alpha0 parameter
    type(parameter_t), intent(inout), optional :: cosmic_alpha1 !< COSMIC alpha1 parameter
    type(parameter_t), intent(inout), optional :: cosmic_l30 !< COSMIC L30 parameter
    type(parameter_t), intent(inout), optional :: cosmic_l31 !< COSMIC L31 parameter
    type(parameter_t), intent(inout), optional :: cosmic_lw0 !< COSMIC lattice-water parameter LW0
    type(parameter_t), intent(inout), optional :: cosmic_lw1 !< COSMIC lattice-water parameter LW1
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (present(cosmic_n0)) then
      cosmic_n0%value = ieee_value(cosmic_n0%value, ieee_quiet_nan) ! sentinel for derived component value
      cosmic_n0%optimize = .false.
      cosmic_n0%lower_bound = ieee_value(cosmic_n0%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      cosmic_n0%upper_bound = ieee_value(cosmic_n0%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      cosmic_n0%lower_bound = 300.0_dp
      cosmic_n0%upper_bound = 2000.0_dp
    end if
    if (present(cosmic_n1)) then
      cosmic_n1%value = ieee_value(cosmic_n1%value, ieee_quiet_nan) ! sentinel for derived component value
      cosmic_n1%optimize = .false.
      cosmic_n1%lower_bound = ieee_value(cosmic_n1%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      cosmic_n1%upper_bound = ieee_value(cosmic_n1%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      cosmic_n1%lower_bound = 0.01_dp
      cosmic_n1%upper_bound = 10.0_dp
    end if
    if (present(cosmic_n2)) then
      cosmic_n2%value = ieee_value(cosmic_n2%value, ieee_quiet_nan) ! sentinel for derived component value
      cosmic_n2%optimize = .false.
      cosmic_n2%lower_bound = ieee_value(cosmic_n2%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      cosmic_n2%upper_bound = ieee_value(cosmic_n2%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      cosmic_n2%lower_bound = 0.01_dp
      cosmic_n2%upper_bound = 10.0_dp
    end if
    if (present(cosmic_alpha0)) then
      cosmic_alpha0%value = ieee_value(cosmic_alpha0%value, ieee_quiet_nan) ! sentinel for derived component value
      cosmic_alpha0%optimize = .false.
      cosmic_alpha0%lower_bound = ieee_value(cosmic_alpha0%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      cosmic_alpha0%upper_bound = ieee_value(cosmic_alpha0%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      cosmic_alpha0%lower_bound = 0.01_dp
      cosmic_alpha0%upper_bound = 10.0_dp
    end if
    if (present(cosmic_alpha1)) then
      cosmic_alpha1%value = ieee_value(cosmic_alpha1%value, ieee_quiet_nan) ! sentinel for derived component value
      cosmic_alpha1%optimize = .false.
      cosmic_alpha1%lower_bound = ieee_value(cosmic_alpha1%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      cosmic_alpha1%upper_bound = ieee_value(cosmic_alpha1%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      cosmic_alpha1%lower_bound = 0.01_dp
      cosmic_alpha1%upper_bound = 10.0_dp
    end if
    if (present(cosmic_l30)) then
      cosmic_l30%value = ieee_value(cosmic_l30%value, ieee_quiet_nan) ! sentinel for derived component value
      cosmic_l30%optimize = .false.
      cosmic_l30%lower_bound = ieee_value(cosmic_l30%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      cosmic_l30%upper_bound = ieee_value(cosmic_l30%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      cosmic_l30%lower_bound = 26.56_dp
      cosmic_l30%upper_bound = 424.78_dp
    end if
    if (present(cosmic_l31)) then
      cosmic_l31%value = ieee_value(cosmic_l31%value, ieee_quiet_nan) ! sentinel for derived component value
      cosmic_l31%optimize = .false.
      cosmic_l31%lower_bound = ieee_value(cosmic_l31%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      cosmic_l31%upper_bound = ieee_value(cosmic_l31%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      cosmic_l31%lower_bound = -118.3_dp
      cosmic_l31%upper_bound = 200.28_dp
    end if
    if (present(cosmic_lw0)) then
      cosmic_lw0%value = ieee_value(cosmic_lw0%value, ieee_quiet_nan) ! sentinel for derived component value
      cosmic_lw0%optimize = .false.
      cosmic_lw0%lower_bound = ieee_value(cosmic_lw0%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      cosmic_lw0%upper_bound = ieee_value(cosmic_lw0%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      cosmic_lw0%lower_bound = 0.0_dp
      cosmic_lw0%upper_bound = 0.2_dp
    end if
    if (present(cosmic_lw1)) then
      cosmic_lw1%value = ieee_value(cosmic_lw1%value, ieee_quiet_nan) ! sentinel for derived component value
      cosmic_lw1%optimize = .false.
      cosmic_lw1%lower_bound = ieee_value(cosmic_lw1%lower_bound, ieee_quiet_nan) ! sentinel for derived component lower_bound
      cosmic_lw1%upper_bound = ieee_value(cosmic_lw1%upper_bound, ieee_quiet_nan) ! sentinel for derived component upper_bound
      cosmic_lw1%lower_bound = 0.0_dp
      cosmic_lw1%upper_bound = 0.05_dp
    end if
  end function nml_neutrons_2_init_type


  !> \brief Read neutrons_2 namelist from file
  integer function nml_neutrons_2_from_file(this, file, errmsg) result(status)
    class(nml_neutrons_2_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    type(parameter_t) :: cosmic_n0
    type(parameter_t) :: cosmic_n1
    type(parameter_t) :: cosmic_n2
    type(parameter_t) :: cosmic_alpha0
    type(parameter_t) :: cosmic_alpha1
    type(parameter_t) :: cosmic_l30
    type(parameter_t) :: cosmic_l31
    type(parameter_t) :: cosmic_lw0
    type(parameter_t) :: cosmic_lw1
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /neutrons_2/ &
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

    status = nml%find("neutrons_2", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=neutrons_2, iostat=iostat, iomsg=iomsg)
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
  end function nml_neutrons_2_from_file

  !> \brief Set neutrons_2 values
  integer function nml_neutrons_2_set(this, &
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

    class(nml_neutrons_2_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    type(parameter_t), intent(in) :: cosmic_n0 !< COSMIC N0 parameter
    type(parameter_t), intent(in) :: cosmic_n1 !< COSMIC N1 parameter
    type(parameter_t), intent(in) :: cosmic_n2 !< COSMIC N2 parameter
    type(parameter_t), intent(in) :: cosmic_alpha0 !< COSMIC alpha0 parameter
    type(parameter_t), intent(in) :: cosmic_alpha1 !< COSMIC alpha1 parameter
    type(parameter_t), intent(in) :: cosmic_l30 !< COSMIC L30 parameter
    type(parameter_t), intent(in) :: cosmic_l31 !< COSMIC L31 parameter
    type(parameter_t), intent(in) :: cosmic_lw0 !< COSMIC lattice-water parameter LW0
    type(parameter_t), intent(in) :: cosmic_lw1 !< COSMIC lattice-water parameter LW1

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
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
  end function nml_neutrons_2_set

  !> \brief Check whether a namelist value was set
  integer function nml_neutrons_2_is_set(this, name, idx, errmsg) result(status)
    class(nml_neutrons_2_t), intent(in) :: this !< namelist instance
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
    case ("cosmic_n0%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_n0'"
        return
      end if
      if (ieee_is_nan(this%cosmic_n0%value)) status = NML_ERR_NOT_SET
    case ("cosmic_n0%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_n0'"
        return
      end if
    case ("cosmic_n0%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_n0'"
        return
      end if
    case ("cosmic_n0%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_n0'"
        return
      end if
    case ("cosmic_n0")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_n0'"
        return
      end if
      if (ieee_is_nan(this%cosmic_n0%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("cosmic_n1%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_n1'"
        return
      end if
      if (ieee_is_nan(this%cosmic_n1%value)) status = NML_ERR_NOT_SET
    case ("cosmic_n1%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_n1'"
        return
      end if
    case ("cosmic_n1%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_n1'"
        return
      end if
    case ("cosmic_n1%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_n1'"
        return
      end if
    case ("cosmic_n1")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_n1'"
        return
      end if
      if (ieee_is_nan(this%cosmic_n1%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("cosmic_n2%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_n2'"
        return
      end if
      if (ieee_is_nan(this%cosmic_n2%value)) status = NML_ERR_NOT_SET
    case ("cosmic_n2%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_n2'"
        return
      end if
    case ("cosmic_n2%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_n2'"
        return
      end if
    case ("cosmic_n2%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_n2'"
        return
      end if
    case ("cosmic_n2")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_n2'"
        return
      end if
      if (ieee_is_nan(this%cosmic_n2%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("cosmic_alpha0%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_alpha0'"
        return
      end if
      if (ieee_is_nan(this%cosmic_alpha0%value)) status = NML_ERR_NOT_SET
    case ("cosmic_alpha0%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_alpha0'"
        return
      end if
    case ("cosmic_alpha0%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_alpha0'"
        return
      end if
    case ("cosmic_alpha0%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_alpha0'"
        return
      end if
    case ("cosmic_alpha0")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_alpha0'"
        return
      end if
      if (ieee_is_nan(this%cosmic_alpha0%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("cosmic_alpha1%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_alpha1'"
        return
      end if
      if (ieee_is_nan(this%cosmic_alpha1%value)) status = NML_ERR_NOT_SET
    case ("cosmic_alpha1%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_alpha1'"
        return
      end if
    case ("cosmic_alpha1%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_alpha1'"
        return
      end if
    case ("cosmic_alpha1%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_alpha1'"
        return
      end if
    case ("cosmic_alpha1")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_alpha1'"
        return
      end if
      if (ieee_is_nan(this%cosmic_alpha1%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("cosmic_l30%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_l30'"
        return
      end if
      if (ieee_is_nan(this%cosmic_l30%value)) status = NML_ERR_NOT_SET
    case ("cosmic_l30%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_l30'"
        return
      end if
    case ("cosmic_l30%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_l30'"
        return
      end if
    case ("cosmic_l30%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_l30'"
        return
      end if
    case ("cosmic_l30")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_l30'"
        return
      end if
      if (ieee_is_nan(this%cosmic_l30%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("cosmic_l31%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_l31'"
        return
      end if
      if (ieee_is_nan(this%cosmic_l31%value)) status = NML_ERR_NOT_SET
    case ("cosmic_l31%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_l31'"
        return
      end if
    case ("cosmic_l31%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_l31'"
        return
      end if
    case ("cosmic_l31%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_l31'"
        return
      end if
    case ("cosmic_l31")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_l31'"
        return
      end if
      if (ieee_is_nan(this%cosmic_l31%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("cosmic_lw0%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_lw0'"
        return
      end if
      if (ieee_is_nan(this%cosmic_lw0%value)) status = NML_ERR_NOT_SET
    case ("cosmic_lw0%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_lw0'"
        return
      end if
    case ("cosmic_lw0%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_lw0'"
        return
      end if
    case ("cosmic_lw0%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_lw0'"
        return
      end if
    case ("cosmic_lw0")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_lw0'"
        return
      end if
      if (ieee_is_nan(this%cosmic_lw0%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("cosmic_lw1%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_lw1'"
        return
      end if
      if (ieee_is_nan(this%cosmic_lw1%value)) status = NML_ERR_NOT_SET
    case ("cosmic_lw1%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_lw1'"
        return
      end if
    case ("cosmic_lw1%lower_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_lw1'"
        return
      end if
    case ("cosmic_lw1%upper_bound")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_lw1'"
        return
      end if
    case ("cosmic_lw1")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'cosmic_lw1'"
        return
      end if
      if (ieee_is_nan(this%cosmic_lw1%value)) then
        status = NML_ERR_NOT_SET
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_neutrons_2_is_set

  !> \brief Validate required values and constraints
  integer function nml_neutrons_2_is_valid(this, errmsg) result(status)
    class(nml_neutrons_2_t), intent(in) :: this !< namelist instance
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
    istat = this%is_set("cosmic_n0", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: cosmic_n0"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("cosmic_n1", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: cosmic_n1"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("cosmic_n2", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: cosmic_n2"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("cosmic_alpha0", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: cosmic_alpha0"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("cosmic_alpha1", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: cosmic_alpha1"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("cosmic_l30", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: cosmic_l30"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("cosmic_l31", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: cosmic_l31"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("cosmic_lw0", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: cosmic_lw0"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("cosmic_lw1", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: cosmic_lw1"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
  end function nml_neutrons_2_is_valid

end module nml_neutrons_2
