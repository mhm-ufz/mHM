program test_exchange_failures
  use, intrinsic :: iso_fortran_env, only: error_unit
  use mo_exchange_type, only: exchange_t, var_dp
  use mo_kind, only: dp, i4

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
  case ("duplicate-binding")
    call duplicate_binding()
  case ("missing-data")
    call missing_data()
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

  subroutine duplicate_binding()
    type(var_dp) :: variable
    real(dp), target :: data(1)

    data = 1.0_dp
    variable = var_dp(name="test")
    call variable%provide("Input")
    call variable%publish_local("Input", data, 1_i4)
    call variable%publish_local("Input", data, 1_i4)
  end subroutine duplicate_binding

  subroutine missing_data()
    type(exchange_t) :: exchange

    call exchange%runoff_total%provide("external")
    call exchange%check_data(exchange%runoff_total, "mRM")
  end subroutine missing_data

end program test_exchange_failures
