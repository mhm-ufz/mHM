program test_exchange_failures
  use, intrinsic :: iso_fortran_env, only: error_unit
  use mo_exchange_type, only: exchange_t, var_dp, l0, single_layer
  use mo_kind, only: dp, i4, i8
  use mo_grid, only: grid_t

  implicit none

  character(64) :: scenario

  if (command_argument_count() /= 1) then
    write(error_unit, '(a)') "usage: test_exchange_failures <scenario>"
    error stop 2
  end if
  call get_command_argument(1, scenario)

  select case (trim(scenario))
  case ("duplicate-provider")
    call duplicate_provider()
  case ("missing-provider")
    call missing_provider()
  case ("storage-shape-change")
    call storage_shape_change()
  case ("missing-data")
    call missing_data()
  case ("storage-copy-shape")
    call storage_copy_shape()
  case ("static-stepping")
    call static_stepping()
  case ("alias-unconnected")
    call alias_unconnected()
  case ("storage-rank")
    call storage_rank()
  case ("data-rank")
    call data_rank()
  case default
    write(error_unit, '(a)') "unknown exchange failure scenario: " // trim(scenario)
    error stop 2
  end select

  write(error_unit, '(a)') "exchange failure scenario unexpectedly succeeded: " // trim(scenario)

contains

  subroutine duplicate_provider()
    type(var_dp) :: variable

    variable = var_dp(name="test")
    call variable%provide("Input")
    call variable%provide("Meteo")
  end subroutine duplicate_provider

  subroutine missing_provider()
    type(var_dp) :: variable

    variable = var_dp(name="test")
    call variable%check_provided("mHM")
  end subroutine missing_provider

  subroutine storage_shape_change()
    type(var_dp) :: variable

    variable = var_dp(name="test")
    call variable%provide("Input")
    call variable%prepare_storage("Input", 1_i4, [1_i8])
    call variable%prepare_storage("Input", 1_i4, [2_i8])
  end subroutine storage_shape_change

  subroutine missing_data()
    type(exchange_t) :: exchange

    call exchange%runoff_total%provide("external")
    call exchange%check_data(exchange%runoff_total, "mRM")
  end subroutine missing_data

  subroutine storage_copy_shape()
    type(var_dp) :: variable

    variable = var_dp(name="test", static=.true.)
    call variable%provide("Input")
    call variable%prepare_storage("Input", 0_i4, [2_i8])
    call variable%set_storage([1.0_dp])
  end subroutine storage_copy_shape

  subroutine static_stepping()
    type(var_dp) :: variable

    variable = var_dp(name="test", static=.true.)
    call variable%provide("Input")
    call variable%prepare_data("Input", 1_i4)
  end subroutine static_stepping

  subroutine alias_unconnected()
    type(var_dp) :: source, alias

    source = var_dp(name="source")
    alias = var_dp(name="alias")
    call source%provide("test")
    call alias%provide("test")
    call alias%associate_alias("test", source)
  end subroutine alias_unconnected


  subroutine storage_rank()
    type(var_dp) :: variable

    variable = var_dp(name="test", static=.true.)
    call variable%provide("Input")
    call variable%prepare_storage("Input", 0_i4, [2_i8, 1_i8])
  end subroutine storage_rank

  subroutine data_rank()
    type(exchange_t), target :: exchange
    type(grid_t), target :: grid
    type(var_dp) :: variable
    real(dp), target :: data(2)

    call grid%init(2_i4, 1_i4, cellsize=1.0_dp)
    exchange%level0 => grid
    data = 1.0_dp
    variable = var_dp(name="test", grid=l0, layers=single_layer)
    call variable%provide("test")
    call variable%prepare_data("test", 1_i4)
    variable%data => data
    call exchange%check_data(variable, "MPR")
  end subroutine data_rank


end program test_exchange_failures
