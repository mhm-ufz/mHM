!> \file    mo_main_config.F90
!> \brief   \copybrief mo_main_config
!> \details \copydetails mo_main_config

!> \brief   Named process-parameter registry for the exchange runtime.
!> \details Components register their own parameter definitions during
!! configuration. The registry provides one deterministic flattened vector for
!! calibration and named slices for process-local consumption.
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jul 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_exchange
#include "logging.h"
module mo_main_config
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use mo_kind, only: i4, dp
  use mo_logging
  use nml_helper, only: parameter_t, to_lower

  implicit none
  private

  integer, parameter :: parameter_name_length = 256

  !> \class process_parameters_t
  !> \brief Registry metadata for one named process.
  type :: process_parameters_t
    character(:), allocatable :: name !< normalized process name
    character(:), allocatable :: group !< selected parameter namelist group
    character(parameter_name_length), allocatable :: local_names(:) !< names used inside the namelist group
    integer(i4) :: first = 0_i4 !< first entry in the flattened registry
    integer(i4) :: count = 0_i4 !< number of process parameters
  end type process_parameters_t

  !> \class parameters_t
  !> \brief Named, process-independent parameter registry.
  type, public :: parameters_t
    private
    type(process_parameters_t), allocatable :: processes(:)
    type(parameter_t), allocatable :: definitions(:)
    character(parameter_name_length), allocatable :: names(:)
    real(dp), allocatable :: values(:)
    logical :: configuring = .false.
    logical :: sealed = .false.
  contains
    procedure :: begin_configuration => parameters_begin_configuration
    procedure :: add_process => parameters_add_process
    procedure :: seal => parameters_seal
    procedure :: reset => parameters_reset
    procedure :: set => parameters_set
    procedure :: as_array => parameters_as_array
    procedure :: get_process => parameters_get_process
    procedure :: is_default => parameters_is_default
    procedure :: get_optimizer_data => parameters_get_optimizer_data
    procedure :: write_namelist => parameters_write_namelist
    procedure :: print => parameters_print
    procedure, private :: require_sealed => parameters_require_sealed
    procedure, private :: validate_values => parameters_validate_values
  end type parameters_t

contains

  !> \brief Clear the registry and open a new configuration phase.
  subroutine parameters_begin_configuration(self)
    class(parameters_t), intent(inout), target :: self

    if (allocated(self%processes)) deallocate(self%processes)
    if (allocated(self%definitions)) deallocate(self%definitions)
    if (allocated(self%names)) deallocate(self%names)
    if (allocated(self%values)) deallocate(self%values)
    allocate(self%processes(0))
    allocate(self%definitions(0))
    allocate(self%names(0))
    allocate(self%values(0))
    self%configuring = .true.
    self%sealed = .false.
  end subroutine parameters_begin_configuration

  !> \brief Register an ordered parameter list for one process.
  subroutine parameters_add_process(self, process, values, names, group)
    class(parameters_t), intent(inout), target :: self
    character(*), intent(in) :: process !< process name
    type(parameter_t), intent(in) :: values(:) !< process parameter definitions
    character(*), intent(in) :: names(:) !< process-local parameter names
    character(*), intent(in) :: group !< selected parameter namelist group
    type(process_parameters_t), allocatable :: new_processes(:)
    type(parameter_t), allocatable :: new_definitions(:)
    character(parameter_name_length), allocatable :: new_names(:)
    character(:), allocatable :: process_name
    character(:), allocatable :: group_name
    character(:), allocatable :: local_name
    character(:), allocatable :: qualified_name
    integer(i4) :: old_process_count
    integer(i4) :: old_parameter_count
    integer(i4) :: i
    integer(i4) :: j

    if (.not.self%configuring .or. self%sealed) then
      log_fatal(*) "parameters%add_process: registry is not open for configuration."
      error stop 1
    end if
    if (size(values) /= size(names)) then
      log_fatal(*) "parameters%add_process: values and names have different sizes for process '", trim(process), "'."
      error stop 1
    end if

    process_name = normalize_name(process)
    group_name = normalize_name(group)
    if (len(process_name) == 0) then
      log_fatal(*) "parameters%add_process: process name must not be empty."
      error stop 1
    end if
    if (len(group_name) == 0) then
      log_fatal(*) "parameters%add_process: namelist group must not be empty for process '", process_name, "'."
      error stop 1
    end if
    do i = 1_i4, size(self%processes, kind=i4)
      if (self%processes(i)%name == process_name) then
        log_fatal(*) "parameters%add_process: process '", process_name, "' is already registered."
        error stop 1
      end if
      if (self%processes(i)%group == group_name) then
        log_fatal(*) "parameters%add_process: namelist group '", group_name, "' is already registered."
        error stop 1
      end if
    end do

    do i = 1_i4, size(names, kind=i4)
      local_name = normalize_name(names(i))
      if (len(local_name) == 0) then
        log_fatal(*) "parameters%add_process: parameter name must not be empty for process '", process_name, "'."
        error stop 1
      end if
      qualified_name = process_name // "." // local_name
      if (len(qualified_name) > parameter_name_length) then
        log_fatal(*) "parameters%add_process: qualified name is too long: '", qualified_name, "'."
        error stop 1
      end if
      do j = 1_i4, i - 1_i4
        if (normalize_name(names(j)) == local_name) then
          log_fatal(*) "parameters%add_process: duplicate parameter '", local_name, "' in process '", process_name, "'."
          error stop 1
        end if
      end do
    end do

    old_process_count = size(self%processes, kind=i4)
    old_parameter_count = size(self%definitions, kind=i4)

    allocate(new_processes(old_process_count + 1_i4))
    if (old_process_count > 0_i4) new_processes(:old_process_count) = self%processes
    new_processes(old_process_count + 1_i4)%name = process_name
    new_processes(old_process_count + 1_i4)%group = group_name
    new_processes(old_process_count + 1_i4)%first = old_parameter_count + 1_i4
    new_processes(old_process_count + 1_i4)%count = size(values, kind=i4)
    allocate(new_processes(old_process_count + 1_i4)%local_names(size(names)))
    do i = 1_i4, size(names, kind=i4)
      new_processes(old_process_count + 1_i4)%local_names(i) = normalize_name(names(i))
    end do
    call move_alloc(new_processes, self%processes)

    allocate(new_definitions(old_parameter_count + size(values, kind=i4)))
    if (old_parameter_count > 0_i4) new_definitions(:old_parameter_count) = self%definitions
    if (size(values) > 0) new_definitions(old_parameter_count + 1_i4:) = values
    call move_alloc(new_definitions, self%definitions)

    allocate(new_names(old_parameter_count + size(names, kind=i4)))
    if (old_parameter_count > 0_i4) new_names(:old_parameter_count) = self%names
    do i = 1_i4, size(names, kind=i4)
      new_names(old_parameter_count + i) = process_name // "." // normalize_name(names(i))
    end do
    call move_alloc(new_names, self%names)
  end subroutine parameters_add_process

  !> \brief Validate definitions and close the configuration phase.
  subroutine parameters_seal(self)
    class(parameters_t), intent(inout), target :: self
    integer(i4) :: i

    if (.not.self%configuring .or. self%sealed) then
      log_fatal(*) "parameters%seal: registry is not open for configuration."
      error stop 1
    end if

    do i = 1_i4, size(self%definitions, kind=i4)
      call validate_definition(self%definitions(i), trim(self%names(i)))
    end do
    self%values = self%definitions%value
    self%configuring = .false.
    self%sealed = .true.
  end subroutine parameters_seal

  !> \brief Restore all current values to their configured inputs.
  subroutine parameters_reset(self)
    class(parameters_t), intent(inout), target :: self

    call self%require_sealed("reset")
    self%values = self%definitions%value
  end subroutine parameters_reset

  !> \brief Replace current values from a complete flattened parameter vector.
  subroutine parameters_set(self, values)
    class(parameters_t), intent(inout), target :: self
    real(dp), intent(in) :: values(:) !< complete flattened parameter vector

    call self%validate_values(values, "set")
    self%values = values
  end subroutine parameters_set

  !> \brief Return a copy of the complete flattened parameter vector.
  function parameters_as_array(self) result(values)
    class(parameters_t), intent(in), target :: self
    real(dp), allocatable :: values(:)

    call self%require_sealed("as_array")
    values = self%values
  end function parameters_as_array

  !> \brief Return current values for a named process.
  function parameters_get_process(self, process) result(values)
    class(parameters_t), intent(in), target :: self
    character(*), intent(in) :: process !< process name
    real(dp), allocatable :: values(:)
    character(:), allocatable :: process_name
    integer(i4) :: first
    integer(i4) :: count
    integer(i4) :: i

    call self%require_sealed("get_process")
    process_name = normalize_name(process)
    do i = 1_i4, size(self%processes, kind=i4)
      if (self%processes(i)%name == process_name) then
        first = self%processes(i)%first
        count = self%processes(i)%count
        allocate(values(count))
        if (count > 0_i4) values = self%values(first:first + count - 1_i4)
        return
      end if
    end do
    log_fatal(*) "parameters%get_process: process '", process_name, "' is not registered."
    error stop 1
  end function parameters_get_process

  !> \brief Return whether a process still uses all configured input values.
  logical function parameters_is_default(self, process)
    class(parameters_t), intent(in), target :: self
    character(*), intent(in) :: process !< process name
    character(:), allocatable :: process_name
    integer(i4) :: first
    integer(i4) :: count
    integer(i4) :: i

    call self%require_sealed("is_default")
    process_name = normalize_name(process)
    do i = 1_i4, size(self%processes, kind=i4)
      if (self%processes(i)%name == process_name) then
        first = self%processes(i)%first
        count = self%processes(i)%count
        parameters_is_default = all(self%values(first:first + count - 1_i4) == &
          self%definitions(first:first + count - 1_i4)%value)
        return
      end if
    end do
    log_fatal(*) "parameters%is_default: process '", process_name, "' is not registered."
    error stop 1
  end function parameters_is_default

  !> \brief Return complete calibration metadata without invoking an optimizer.
  subroutine parameters_get_optimizer_data(self, initial, bounds, mask, names)
    class(parameters_t), intent(in), target :: self
    real(dp), allocatable, intent(out) :: initial(:)
    real(dp), allocatable, intent(out) :: bounds(:, :)
    logical, allocatable, intent(out) :: mask(:)
    character(parameter_name_length), allocatable, intent(out) :: names(:)
    integer(i4) :: count

    call self%require_sealed("get_optimizer_data")
    count = size(self%definitions, kind=i4)
    allocate(initial(count), bounds(count, 2), mask(count), names(count))
    if (count == 0_i4) return
    initial = self%definitions%value
    bounds(:, 1) = self%definitions%lower_bound
    bounds(:, 2) = self%definitions%upper_bound
    mask = self%definitions%optimize
    names = self%names
  end subroutine parameters_get_optimizer_data

  !> \brief Write registered current or supplied values as a forward-run parameter namelist.
  subroutine parameters_write_namelist(self, file, values)
    class(parameters_t), intent(in), target :: self
    character(*), intent(in) :: file !< destination parameter namelist path
    real(dp), intent(in), optional :: values(:) !< complete flattened parameter vector to write
    real(dp), allocatable :: output_values(:)
    character(1024) :: iomsg
    integer(i4) :: i
    integer(i4) :: j
    integer(i4) :: value_index
    integer :: unit
    integer :: status

    call self%require_sealed("write_namelist")
    if (len_trim(file) == 0) then
      log_fatal(*) "parameters%write_namelist: destination path must not be empty."
      error stop 1
    end if
    if (present(values)) then
      output_values = values
    else
      output_values = self%values
    end if
    call self%validate_values(output_values, "write_namelist")

    open(newunit=unit, file=file, status="replace", action="write", form="formatted", iostat=status, iomsg=iomsg)
    if (status /= 0) then
      log_fatal(*) "parameters%write_namelist: could not open '", trim(file), "': ", trim(iomsg)
      error stop 1
    end if

    call write_text_line(unit, "! Final v6 parameter values for forward simulation.", file)
    do i = 1_i4, size(self%processes, kind=i4)
      call write_text_line(unit, "", file)
      call write_text_line(unit, "&" // self%processes(i)%group, file)
      do j = 1_i4, self%processes(i)%count
        value_index = self%processes(i)%first + j - 1_i4
        write(unit, '("  ",a," = ",es24.16e3)', iostat=status, iomsg=iomsg) &
          trim(self%processes(i)%local_names(j)), output_values(value_index)
        if (status /= 0) then
          log_fatal(*) "parameters%write_namelist: could not write '", trim(file), "': ", trim(iomsg)
          error stop 1
        end if
      end do
      call write_text_line(unit, "/", file)
    end do

    close(unit, iostat=status, iomsg=iomsg)
    if (status /= 0) then
      log_fatal(*) "parameters%write_namelist: could not close '", trim(file), "': ", trim(iomsg)
      error stop 1
    end if
  end subroutine parameters_write_namelist

  !> \brief Print current named parameter values.
  subroutine parameters_print(self)
    class(parameters_t), intent(in), target :: self
    integer(i4) :: i

    call self%require_sealed("print")
    do i = 1_i4, size(self%values, kind=i4)
      log_info(*) trim(self%names(i)), " = ", self%values(i)
    end do
  end subroutine parameters_print

  !> \brief Validate one configured parameter definition.
  subroutine validate_definition(definition, name)
    type(parameter_t), intent(in) :: definition
    character(*), intent(in) :: name

    if (.not.ieee_is_finite(definition%value) .or. &
        .not.ieee_is_finite(definition%lower_bound) .or. &
        .not.ieee_is_finite(definition%upper_bound)) then
      log_fatal(*) "parameters%seal: definition for '", trim(name), "' contains a non-finite value."
      error stop 1
    end if
    if (definition%lower_bound > definition%upper_bound) then
      log_fatal(*) "parameters%seal: invalid bounds for '", trim(name), "'."
      error stop 1
    end if
    if (definition%optimize .and. definition%lower_bound >= definition%upper_bound) then
      log_fatal(*) "parameters%seal: optimized parameter '", trim(name), "' requires strict bounds."
      error stop 1
    end if
    if (definition%value < definition%lower_bound .or. definition%value > definition%upper_bound) then
      log_fatal(*) "parameters%seal: initial value for '", trim(name), "' is outside its bounds."
      error stop 1
    end if
  end subroutine validate_definition

  !> \brief Fail unless the registry has completed configuration.
  subroutine parameters_require_sealed(self, method)
    class(parameters_t), intent(in), target :: self
    character(*), intent(in) :: method

    if (.not.self%sealed) then
      log_fatal(*) "parameters%", trim(method), ": registry is not sealed."
      error stop 1
    end if
  end subroutine parameters_require_sealed

  !> \brief Validate a complete flattened parameter vector.
  subroutine parameters_validate_values(self, values, method)
    class(parameters_t), intent(in), target :: self
    real(dp), intent(in) :: values(:)
    character(*), intent(in) :: method
    integer(i4) :: i

    call self%require_sealed(method)
    if (size(values) /= size(self%values)) then
      log_fatal(*) "parameters%", trim(method), ": expected ", size(self%values), " values, got ", size(values), "."
      error stop 1
    end if
    do i = 1_i4, size(values, kind=i4)
      if (.not.ieee_is_finite(values(i))) then
        log_fatal(*) "parameters%", trim(method), ": value for '", trim(self%names(i)), "' is not finite."
        error stop 1
      end if
      if (values(i) < self%definitions(i)%lower_bound .or. values(i) > self%definitions(i)%upper_bound) then
        log_fatal(*) "parameters%", trim(method), ": value for '", trim(self%names(i)), "' is outside its bounds."
        error stop 1
      end if
    end do
  end subroutine parameters_validate_values

  !> \brief Write one text line and fail with file context on an I/O error.
  subroutine write_text_line(unit, line, file)
    integer, intent(in) :: unit
    character(*), intent(in) :: line
    character(*), intent(in) :: file
    character(1024) :: iomsg
    integer :: status

    write(unit, '(a)', iostat=status, iomsg=iomsg) line
    if (status /= 0) then
      log_fatal(*) "parameters%write_namelist: could not write '", trim(file), "': ", trim(iomsg)
      error stop 1
    end if
  end subroutine write_text_line

  !> \brief Normalize registry identifiers for case-insensitive lookup.
  pure function normalize_name(name) result(normalized)
    character(*), intent(in) :: name
    character(:), allocatable :: normalized

    normalized = trim(adjustl(to_lower(name)))
  end function normalize_name

end module mo_main_config
