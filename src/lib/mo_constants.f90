!> \file mo_constants.f90

!> \brief Provides computational, mathematical, physical, and file constants

!> \details Provides computational constants like epsilon, mathematical constants such as Pi,
!> physical constants such as the Stefan-Boltzmann constant, and file units for some standard streams
!> such as standard in.

!> \author Matthias Cuntz
!> \date Nov 2011

module mo_constants

  !  This module contains basic and derived constants
  !
  !  Written  Nov 2011, Matthias Cuntz
  !  Modified Mar 2014, Matthias Cuntz - iso_fortran_env

  ! License
  ! -------
  ! This file is part of the UFZ Fortran library.

  ! The UFZ Fortran library is free software: you can redistribute it and/or modify
  ! it under the terms of the GNU Lesser General Public License as published by
  ! the Free Software Foundation, either version 3 of the License, or
  ! (at your option) any later version.

  ! The UFZ Fortran library is distributed in the hope that it will be useful,
  ! but WITHOUT ANY WARRANTY; without even the implied warranty of
  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  ! GNU Lesser General Public License for more details.

  ! You should have received a copy of the GNU Lesser General Public License
  ! along with the UFZ Fortran library (cf. gpl.txt and lgpl.txt).
  ! If not, see <http://www.gnu.org/licenses/>.

  ! Copyright 2011-2014 Matthias Cuntz

  use mo_kind, only : sp, dp, i4, i8
  use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit

  implicit none

  ! Mathematical
  !> Pi in double precision
  real(dp), parameter :: PI_dp = 3.141592653589793238462643383279502884197_dp    ! Pi
  !> Pi in single precision
  real(sp), parameter :: PI_sp = 3.141592653589793238462643383279502884197_sp
  !> Pi/2 in double precision
  real(dp), parameter :: PIO2_dp = 1.57079632679489661923132169163975144209858_dp  ! Pi/2
  !> Pi/2 in single precision
  real(sp), parameter :: PIO2_sp = 1.57079632679489661923132169163975144209858_sp
  !> 2*Pi in double precision
  real(dp), parameter :: TWOPI_dp = 6.283185307179586476925286766559005768394_dp    ! 2*Pi
  !> 2*Pi in single precision
  real(sp), parameter :: TWOPI_sp = 6.283185307179586476925286766559005768394_sp
  !> Square root of 2 in double precision
  real(dp), parameter :: SQRT2_dp = 1.41421356237309504880168872420969807856967_dp  ! Sqrt(2)
  !> Square root of 2 in single precision
  real(sp), parameter :: SQRT2_sp = 1.41421356237309504880168872420969807856967_sp
  !> 2/3 in double precision
  real(dp), parameter :: TWOTHIRD_dp = 0.6666666666666666666666666666666666667_dp      ! 2/3
  !> 2/3 in single precision
  real(sp), parameter :: TWOTHIRD_sp = 0.6666666666666666666666666666666666667_sp
  !> degree to radian conversion (pi/180) in double precision
  real(dp), parameter :: deg2rad_dp = PI_dp / 180._dp                       ! deg2rad
  !> degree to radian conversion (pi/180) in double precision
  real(sp), parameter :: deg2rad_sp = PI_sp / 180._sp
  !> radian to conversion (180/pi) in double precision
  real(dp), parameter :: rad2deg_dp = 180._dp / PI_dp                       ! rad2deg
  !> radian to degree conversion (180/pi) in single precision
  real(sp), parameter :: rad2deg_sp = 180._sp / PI_sp

  !> Time conversion
  !> Seconds per day [s] in single precision
  real(sp), public, parameter :: secday_sp = 86400.0_sp
  real(dp), public, parameter :: secday_dp = 86400.0_dp  ! secday [s]
  real(dp), public, parameter :: DayHours = 24.0_dp  ! hours per day
  real(dp), public, parameter :: YearMonths = 12.0_dp  ! months per year
  real(dp), public, parameter :: YearDays = 365.0_dp  ! days in a year
  real(dp), public, parameter :: DaySecs = 86400.0_dp  ! sec in a day
  real(dp), public, parameter :: HourSecs = 3600.0_dp  ! seconds per hour

  ! Physical
  !> Psychrometric constant [kPa K^-1] in double precision
  real(dp), parameter :: Psychro_dp = 0.0646_dp                 ! psychrometric constant [kPa C-1]
  !> Psychrometric constant [kPa K^-1] in sibgle precision
  real(sp), parameter :: Psychro_sp = 0.0646_sp
  !> Gravity accelaration [m^2 s^-1] in double precision
  real(dp), parameter :: Gravity_dp = 9.81_dp                      ! Gravity acceleration [m^2/s]
  !> Gravity accelaration [m^2 s^-1] in single precision
  real(sp), parameter :: Gravity_sp = 9.81_sp
  !>  Solar constant in [J m^-2 s^-1] in double precision
  real(dp), parameter :: SolarConst_dp = 1367._dp                   ! Solar constant in [W m-2 = kg s-3]
  !>  Solar constant in [J m^-2 s^-1] in single precision
  real(sp), parameter :: SolarConst_sp = 1367._sp
  !> Specific heat for vaporization of water in [J m-2 mm-1] in double precision
  real(dp), parameter :: SpecHeatET_dp = 2.45e06_dp                ! Specific heat in [W s m-2 mm-1 = kg s-2 mm-1]
  !> Specific heat for vaporization of water in [J m-2 mm-1] in single precision
  real(sp), parameter :: SpecHeatET_sp = 2.45e06_sp
  !> Standard temperature [K] in double precision
  real(dp), parameter :: T0_dp = 273.15_dp                    ! Celcius <-> Kelvin [K]
  !> Standard temperature [K] in single precision
  real(sp), parameter :: T0_sp = 273.15_sp
  !> Stefan-Boltzmann constant [W m^-2 K^-4] in double precision
  real(dp), parameter :: sigma_dp = 5.67e-08_dp                  ! Stefan-Boltzmann constant [W/m^2/K^4]
  !> Stefan-Boltzmann constant [W m^-2 K^-4] in single precision
  real(sp), parameter :: sigma_sp = 5.67e-08_sp
  ! Earth radius [m] in single precision
  real(sp), parameter :: RadiusEarth_sp = 6371228._sp
  ! Earth radius [m] in double precision
  real(dp), parameter :: RadiusEarth_dp = 6371228._dp

  !> standard atmospehere
  !> Standard pressure [Pa] in double precision
  real(dp), parameter :: P0_dp = 101325._dp                   ! Standard pressure [Pa]
  !> Standard pressure [Pa] in single precision
  real(sp), parameter :: P0_sp = 101325._sp
  !> standard density  [kg m^-3] in double precision
  real(dp), parameter :: rho0_dp = 1.225_dp                    ! Standard air density
  !> standard density  [kg m^-3] in single precision
  real(sp), parameter :: rho0_sp = 1.225_sp
  !> specific heat capacity of air [J kg^-1 K^-1] in double precision
  real(dp), parameter :: cp0_dp = 1005.0_dp                   ! Standard specific heat of air
  !> specific heat capacity of air [J kg^-1 K^-1] in single precision
  real(sp), parameter :: cp0_sp = 1005.0_sp
  !> specific heat capacity of water [J kg^-1 K^-1] in double precision
  real(dp), parameter :: cp_w_dp = 4.19_dp
  !> specific heat capacity of water [J kg^-1 K^-1] in single precision
  real(sp), parameter :: cp_w_sp = 4.19_sp

  !> Pi in double precision
  real(dp), parameter :: PI_D = 3.141592653589793238462643383279502884197_dp      ! Pi
  !> Pi in single precision
  real(sp), parameter :: PI = 3.141592653589793238462643383279502884197_sp
  !> Pi/2 in double precision
  real(dp), parameter :: PIO2_D = 1.57079632679489661923132169163975144209858_dp    ! Pi/2
  !> Pi/2 in single precision
  real(sp), parameter :: PIO2 = 1.57079632679489661923132169163975144209858_sp
  !> 2*Pi in double precision
  real(dp), parameter :: TWOPI_D = 6.283185307179586476925286766559005768394_dp      ! 2*Pi
  !> 2*Pi in single precision
  real(sp), parameter :: TWOPI = 6.283185307179586476925286766559005768394_sp
  !> Square root of 2 in double precision
  real(dp), parameter :: SQRT2_D = 1.41421356237309504880168872420969807856967_dp    ! Sqrt(2)
  !> Square root of 2 in single precision
  real(sp), parameter :: SQRT2 = 1.41421356237309504880168872420969807856967_sp
  !> Euler''s constant in double precision
  real(dp), parameter :: EULER_D = 0.5772156649015328606065120900824024310422_dp     ! Euler
  !> Euler''s constant in single precision
  real(sp), parameter :: EULER = 0.5772156649015328606065120900824024310422_sp

  ! Standard file units
  !> Standard input file unit
  ! integer, parameter :: nin  = 5   ! standard input stream
  integer, parameter :: nin = input_unit   ! standard input stream
  !> Standard output file unit
  ! integer, parameter :: nout = 6   ! standard output stream
  integer, parameter :: nout = output_unit   ! standard output stream
  !> Standard error file unit
  ! integer, parameter :: nerr = 0   ! error output stream
  integer, parameter :: nerr = error_unit   ! error output stream
  !> Standard file unit for namelist
  integer, parameter :: nnml = 100 ! namelist unit

  ! computational, these values need to be the same!!!
  real(dp), public, parameter :: nodata_dp = -9999.0_dp ! [-]     global no data value
  integer(i4), public, parameter :: nodata_i4 = int(nodata_dp)  ! [-]     global no data value
  integer(i8), public, parameter :: nodata_i8 = int(nodata_dp, kind=i8)  ! [-]     global no data value
  !> epsilon(1.0) in double precision
  real(dp), public, parameter :: eps_dp = epsilon(1.0_dp)
  real(sp), public, parameter :: eps_sp = epsilon(1.0_sp)

end module mo_constants
