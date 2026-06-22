!> \file nml_config_optimize.f90
!> \copydoc nml_config_optimize

!> \brief Optimization configuration
!> \details Optimization configuration
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_config_optimize
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
  ! kind specifiers listed in the nml-tools configuration file
  use mo_kind, only: &
    i4, &
    dp

  implicit none

  ! default values
  logical, parameter, public :: optimize__default = .false.
  logical, parameter, public :: optimize_restart__default = .false.
  integer(i4), parameter, public :: opti_method__default = 1_i4
  integer(i4), parameter, public :: opti_function__default = 10_i4
  integer(i4), parameter, public :: iterations__default = 7_i4
  integer(i4), parameter, public :: seed__default = 1235876_i4
  real(dp), parameter, public :: dds_r__default = 0.2_dp
  real(dp), parameter, public :: sa_temp__default = -9.0_dp
  integer(i4), parameter, public :: sce_ngs__default = 2_i4
  integer(i4), parameter, public :: sce_npg__default = -9_i4
  integer(i4), parameter, public :: sce_nps__default = -9_i4
  logical, parameter, public :: mcmc_opti__default = .false.
  real(dp), parameter, public :: mcmc_error_params__default(2) = [0.01_dp, 0.6_dp]
  logical, parameter, public :: bfi_calc__default = .true.

  !> \class nml_config_optimize_t
  !> \brief Optimization configuration
  !> \details Optimization configuration
  type, public :: nml_config_optimize_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    logical :: optimize !< Enable optimization
    logical :: optimize_restart !< Restart optimization
    integer(i4) :: opti_method !< Optimization method
    integer(i4) :: opti_function !< Objective function to be optimized
    integer(i4) :: iterations !< Number of iterations
    integer(i4) :: seed !< Seed
    real(dp) :: dds_r !< DDS perturbation rate
    real(dp) :: sa_temp !< Initial Temperature
    integer(i4) :: sce_ngs !< Number of SCE Complexes
    integer(i4) :: sce_npg !< Points per SCE Complex
    integer(i4) :: sce_nps !< Points per SCE Sub-Complex
    logical :: mcmc_opti !< MCMC for optimisation
    real(dp), dimension(2) :: mcmc_error_params !< Parameters of MCMC error model.
    logical :: bfi_calc !< Calculate BFI from discharge time series with the Eckhardt filter
  contains
    procedure :: init => nml_config_optimize_init
    procedure :: from_file => nml_config_optimize_from_file
    procedure :: set => nml_config_optimize_set
    procedure :: is_set => nml_config_optimize_is_set
    procedure :: is_valid => nml_config_optimize_is_valid
  end type nml_config_optimize_t

contains

  !> \brief Initialize defaults and sentinels for config_optimize
  integer function nml_config_optimize_init(this, errmsg) result(status)
    class(nml_config_optimize_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! default values
    this%optimize = optimize__default ! bool values always need a default
    this%optimize_restart = optimize_restart__default ! bool values always need a default
    this%opti_method = opti_method__default
    this%opti_function = opti_function__default
    this%iterations = iterations__default
    this%seed = seed__default
    this%dds_r = dds_r__default
    this%sa_temp = sa_temp__default
    this%sce_ngs = sce_ngs__default
    this%sce_npg = sce_npg__default
    this%sce_nps = sce_nps__default
    this%mcmc_opti = mcmc_opti__default ! bool values always need a default
    this%mcmc_error_params = mcmc_error_params__default
    this%bfi_calc = bfi_calc__default ! bool values always need a default
  end function nml_config_optimize_init


  !> \brief Read config_optimize namelist from file
  integer function nml_config_optimize_from_file(this, file, errmsg) result(status)
    class(nml_config_optimize_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    logical :: optimize
    logical :: optimize_restart
    integer(i4) :: opti_method
    integer(i4) :: opti_function
    integer(i4) :: iterations
    integer(i4) :: seed
    real(dp) :: dds_r
    real(dp) :: sa_temp
    integer(i4) :: sce_ngs
    integer(i4) :: sce_npg
    integer(i4) :: sce_nps
    logical :: mcmc_opti
    real(dp), dimension(2) :: mcmc_error_params
    logical :: bfi_calc
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /config_optimize/ &
      optimize, &
      optimize_restart, &
      opti_method, &
      opti_function, &
      iterations, &
      seed, &
      dds_r, &
      sa_temp, &
      sce_ngs, &
      sce_npg, &
      sce_nps, &
      mcmc_opti, &
      mcmc_error_params, &
      bfi_calc

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    optimize = this%optimize
    optimize_restart = this%optimize_restart
    opti_method = this%opti_method
    opti_function = this%opti_function
    iterations = this%iterations
    seed = this%seed
    dds_r = this%dds_r
    sa_temp = this%sa_temp
    sce_ngs = this%sce_ngs
    sce_npg = this%sce_npg
    sce_nps = this%sce_nps
    mcmc_opti = this%mcmc_opti
    mcmc_error_params = this%mcmc_error_params
    bfi_calc = this%bfi_calc

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("config_optimize", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=config_optimize, iostat=iostat, iomsg=iomsg)
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
    this%optimize = optimize
    this%optimize_restart = optimize_restart
    this%opti_method = opti_method
    this%opti_function = opti_function
    this%iterations = iterations
    this%seed = seed
    this%dds_r = dds_r
    this%sa_temp = sa_temp
    this%sce_ngs = sce_ngs
    this%sce_npg = sce_npg
    this%sce_nps = sce_nps
    this%mcmc_opti = mcmc_opti
    this%mcmc_error_params = mcmc_error_params
    this%bfi_calc = bfi_calc

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_config_optimize_from_file

  !> \brief Set config_optimize values
  integer function nml_config_optimize_set(this, &
    optimize, &
    optimize_restart, &
    opti_method, &
    opti_function, &
    iterations, &
    seed, &
    dds_r, &
    sa_temp, &
    sce_ngs, &
    sce_npg, &
    sce_nps, &
    mcmc_opti, &
    mcmc_error_params, &
    bfi_calc, &
    errmsg) result(status)

    class(nml_config_optimize_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    logical, intent(in), optional :: optimize !< Enable optimization
    logical, intent(in), optional :: optimize_restart !< Restart optimization
    integer(i4), intent(in), optional :: opti_method !< Optimization method
    integer(i4), intent(in), optional :: opti_function !< Objective function to be optimized
    integer(i4), intent(in), optional :: iterations !< Number of iterations
    integer(i4), intent(in), optional :: seed !< Seed
    real(dp), intent(in), optional :: dds_r !< DDS perturbation rate
    real(dp), intent(in), optional :: sa_temp !< Initial Temperature
    integer(i4), intent(in), optional :: sce_ngs !< Number of SCE Complexes
    integer(i4), intent(in), optional :: sce_npg !< Points per SCE Complex
    integer(i4), intent(in), optional :: sce_nps !< Points per SCE Sub-Complex
    logical, intent(in), optional :: mcmc_opti !< MCMC for optimisation
    real(dp), dimension(:), intent(in), optional :: mcmc_error_params !< Parameters of MCMC error model.
    logical, intent(in), optional :: bfi_calc !< Calculate BFI from discharge time series with the Eckhardt filter
    integer :: &
      lb__1, &
      ub__1

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    ! override with provided values
    if (present(optimize)) this%optimize = optimize
    if (present(optimize_restart)) this%optimize_restart = optimize_restart
    if (present(opti_method)) this%opti_method = opti_method
    if (present(opti_function)) this%opti_function = opti_function
    if (present(iterations)) this%iterations = iterations
    if (present(seed)) this%seed = seed
    if (present(dds_r)) this%dds_r = dds_r
    if (present(sa_temp)) this%sa_temp = sa_temp
    if (present(sce_ngs)) this%sce_ngs = sce_ngs
    if (present(sce_npg)) this%sce_npg = sce_npg
    if (present(sce_nps)) this%sce_nps = sce_nps
    if (present(mcmc_opti)) this%mcmc_opti = mcmc_opti
    if (present(mcmc_error_params)) then
      if (size(mcmc_error_params, 1) > size(this%mcmc_error_params, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'mcmc_error_params'"
        return
      end if
      lb__1 = lbound(this%mcmc_error_params, 1)
      ub__1 = lb__1 + size(mcmc_error_params, 1) - 1
      this%mcmc_error_params(lb__1:ub__1) = mcmc_error_params
    end if
    if (present(bfi_calc)) this%bfi_calc = bfi_calc

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_config_optimize_set

  !> \brief Check whether a namelist value was set
  integer function nml_config_optimize_is_set(this, name, idx, errmsg) result(status)
    class(nml_config_optimize_t), intent(in) :: this !< namelist instance
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
    case ("optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'optimize'"
        return
      end if
    case ("optimize_restart")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'optimize_restart'"
        return
      end if
    case ("opti_method")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'opti_method'"
        return
      end if
    case ("opti_function")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'opti_function'"
        return
      end if
    case ("iterations")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'iterations'"
        return
      end if
    case ("seed")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'seed'"
        return
      end if
    case ("dds_r")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'dds_r'"
        return
      end if
    case ("sa_temp")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'sa_temp'"
        return
      end if
    case ("sce_ngs")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'sce_ngs'"
        return
      end if
    case ("sce_npg")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'sce_npg'"
        return
      end if
    case ("sce_nps")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'sce_nps'"
        return
      end if
    case ("mcmc_opti")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'mcmc_opti'"
        return
      end if
    case ("mcmc_error_params")
      if (present(idx)) then
        status = idx_check(idx, lbound(this%mcmc_error_params), ubound(this%mcmc_error_params), &
          "mcmc_error_params", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("bfi_calc")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'BFI_calc'"
        return
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_config_optimize_is_set

  !> \brief Validate required values and constraints
  integer function nml_config_optimize_is_valid(this, errmsg) result(status)
    class(nml_config_optimize_t), intent(in) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    integer :: istat

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (.not. this%is_configured) then
      status = NML_ERR_NOT_SET
      if (present(errmsg)) errmsg = "namelist not configured; call set or from_file"
      return
    end if

  end function nml_config_optimize_is_valid

end module nml_config_optimize
