!> \file mo_constants.f90

!> \brief Provides computational, mathematical, physical, and file constants

!> \details Provides computational constants like epsilon, mathematical constants such as Pi,
!> physical constants such as the Stefan-Boltzmann constant, and file units for some standard streams
!> such as standard in.

!> \author Matthias Cuntz
!> \date Nov 2011

MODULE mo_constants

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

  USE mo_kind, ONLY: sp, dp
  use, intrinsic :: iso_fortran_env, only: input_unit, output_unit, error_unit

  IMPLICIT NONE

  ! Computational
  !> epsilon(1.0) in double precision
  REAL(dp), PARAMETER :: eps_dp = epsilon(1.0_dp)
  !> epsilon(1.0) in single precision
  REAL(sp), PARAMETER :: eps_sp = epsilon(1.0_sp)

  ! Mathematical
  !> Pi in double precision
  REAL(dp), PARAMETER :: PI_dp    = 3.141592653589793238462643383279502884197_dp    ! Pi
  !> Pi in single precision
  REAL(sp), PARAMETER :: PI_sp    = 3.141592653589793238462643383279502884197_sp
  !> Pi/2 in double precision
  REAL(dp), PARAMETER :: PIO2_dp  = 1.57079632679489661923132169163975144209858_dp  ! Pi/2
  !> Pi/2 in single precision
  REAL(sp), PARAMETER :: PIO2_sp  = 1.57079632679489661923132169163975144209858_sp
  !> 2*Pi in double precision
  REAL(dp), PARAMETER :: TWOPI_dp = 6.283185307179586476925286766559005768394_dp    ! 2*Pi
  !> 2*Pi in single precision
  REAL(sp), PARAMETER :: TWOPI_sp = 6.283185307179586476925286766559005768394_sp
  !> Square root of 2 in double precision
  REAL(dp), PARAMETER :: SQRT2_dp = 1.41421356237309504880168872420969807856967_dp  ! Sqrt(2)
  !> Square root of 2 in single precision
  REAL(sp), PARAMETER :: SQRT2_sp = 1.41421356237309504880168872420969807856967_sp
  !> 2/3 in double precision
  REAL(dp), PARAMETER :: TWOTHIRD_dp = 0.6666666666666666666666666666666666667_dp      ! 2/3
  !> 2/3 in single precision
  REAL(sp), PARAMETER :: TWOTHIRD_sp = 0.6666666666666666666666666666666666667_sp
  !> degree to radian conversion (pi/180) in double precision
  REAL(dp), PARAMETER :: deg2rad_dp = PI_dp/180._dp                       ! deg2rad
  !> degree to radian conversion (pi/180) in double precision
  REAL(sp), PARAMETER :: deg2rad_sp = PI_sp/180._sp
  !> radian to degree conversion (180/pi) in double precision
  REAL(dp), PARAMETER :: rad2deg_dp = 180._dp/PI_dp                       ! rad2deg
  !> radian to degree conversion (180/pi) in single precision
  REAL(sp), PARAMETER :: rad2deg_sp = 180._sp/PI_sp
  !> arc to radian conversion (180/pi) in single precision
  REAL(dp), PARAMETER :: arc2rad_dp = PI_dp/12._dp                        !arc2rad
  !> arc to radian conversion (180/pi) in single precision
  REAL(sp), PARAMETER :: arc2rad_sp = PI_dp/12._dp
  ! seconds per radian of hour angle
  REAL(sp), PARAMETER :: SECPERRAD = 13750.9871_dp


  ! Physical
  !> Seconds per day [s] in double precision
  REAL(dp), PARAMETER :: secday_dp = 86400._dp                       ! secday [s]
  !> Seconds per day [s] in single precision
  REAL(sp), PARAMETER :: secday_sp = 86400._sp
  !> Psychrometric constant [kPa K^-1] in double precision
  REAL(dp), PARAMETER :: Psychro_dp      = 0.0646_dp                 ! psychrometric constant [kPa C-1]
  !> Psychrometric constant [kPa K^-1] in sibgle precision
  REAL(sp), PARAMETER :: Psychro_sp      = 0.0646_sp
  !> Gravity accelaration [m^2 s^-1] in double precision
  REAL(dp), PARAMETER :: Gravity_dp   = 9.81_dp                      ! Gravity acceleration [m^2/s]
  !> Gravity accelaration [m^2 s^-1] in single precision
  REAL(sp), PARAMETER :: Gravity_sp   = 9.81_sp
  !>  Solar constant in [J m^-2 s^-1] in double precision
  REAL(dp), PARAMETER :: SolarConst_dp  = 1367._dp                   ! Solar constant in [W m-2 = kg s-3]
  !>  Solar constant in [J m^-2 s^-1] in single precision
  REAL(sp), PARAMETER :: SolarConst_sp  = 1367._sp
  !> Specific heat for vaporization of water in [J m-2 mm-1] in double precision
  REAL(dp), PARAMETER :: SpecHeatET_dp   = 2.45e06_dp                ! Specific heat in [W s m-2 mm-1 = kg s-2 mm-1]
  !> Specific heat for vaporization of water in [J m-2 mm-1] in single precision
  REAL(sp), PARAMETER :: SpecHeatET_sp   = 2.45e06_sp
  !> Standard temperature [K] in double precision
  REAL(dp), PARAMETER :: T0_dp        = 273.15_dp                    ! Celcius <-> Kelvin [K]
  !> Standard temperature [K] in single precision
  REAL(sp), PARAMETER :: T0_sp        = 273.15_sp
  !> Stefan-Boltzmann constant [W m^-2 K^-4] in double precision
  REAL(dp), PARAMETER :: sigma_dp     = 5.67e-08_dp                  ! Stefan-Boltzmann constant [W/m^2/K^4]
  !> Stefan-Boltzmann constant [W m^-2 K^-4] in single precision
  REAL(sp), PARAMETER :: sigma_sp     = 5.67e-08_sp
  ! Earth radius [m] in single precision
  REAL(sp), PARAMETER   :: RadiusEarth_sp  = 6371228._sp
  ! Earth radius [m] in double precision
  REAL(dp), PARAMETER   :: RadiusEarth_dp  = 6371228._dp
  ! radians of Earth orbit per julian day
  REAL(dp), PARAMETER   :: RADPERDAY = 0.017214
  ! minimum declination (radians)
  REAL(dp), PARAMETER   :: MINDECL   = -0.4092797_dp
  ! julian day offset of winter solstice
  REAL(dp), PARAMETER   ::  DAYSOFF  = 11.25_dp
  ! (kg mol-1) molecular weight of air
  REAL(dp), PARAMETER   ::  MA       = 28.9644e-3_dp
  ! (kg mol-1) molecular weight of water
  REAL(dp), PARAMETER   ::  MW       = 18.0148e-3_dp
  ! (MW/MA) unitless ratio of molec weights
  REAL(dp), PARAMETER   ::  EPS      = MW/MA
  ! (m3 Pa mol-1 K-1) gas law constant
  REAL(dp), PARAMETER   ::  R        = 8.3143_dp
  ! (-K m-1) standard temperature lapse rate
  REAL(dp), PARAMETER   ::  LR_STD   = 0.0065


  !> standard atmospehere
  !> Standard pressure [Pa] in double precision
  REAL(dp), PARAMETER :: P0_dp        = 101325._dp                   ! Standard pressure [Pa]
  !> Standard pressure [Pa] in single precision
  REAL(sp), PARAMETER :: P0_sp        = 101325._sp
  !> standard density  [kg m^-3] in double precision
  REAL(dp), PARAMETER :: rho0_dp       = 1.225_dp                    ! Standard air density
  !> standard density  [kg m^-3] in single precision
  REAL(sp), PARAMETER :: rho0_sp       = 1.225_sp
  !> specific heat capacity of air [J kg^-1 K^-1] in double precision
  REAL(dp), PARAMETER :: cp0_dp        = 1005.0_dp                   ! Standard specific heat of air
  !> specific heat capacity of air [J kg^-1 K^-1] in single precision
  REAL(sp), PARAMETER :: cp0_sp        = 1005.0_sp
  !> (K) standard temp at 0.0 m elevation
  REAL(sp), PARAMETER :: T_STD = 288.15

  !> Pi in double precision
  REAL(dp), PARAMETER :: PI_D    = 3.141592653589793238462643383279502884197_dp      ! Pi
  !> Pi in single precision
  REAL(sp), PARAMETER :: PI      = 3.141592653589793238462643383279502884197_sp
  !> Pi/2 in double precision
  REAL(dp), PARAMETER :: PIO2_D  = 1.57079632679489661923132169163975144209858_dp    ! Pi/2
  !> Pi/2 in single precision
  REAL(sp), PARAMETER :: PIO2    = 1.57079632679489661923132169163975144209858_sp
  !> 2*Pi in double precision
  REAL(dp), PARAMETER :: TWOPI_D = 6.283185307179586476925286766559005768394_dp      ! 2*Pi
  !> 2*Pi in single precision
  REAL(sp), PARAMETER :: TWOPI   = 6.283185307179586476925286766559005768394_sp
  !> Square root of 2 in double precision
  REAL(dp), PARAMETER :: SQRT2_D = 1.41421356237309504880168872420969807856967_dp    ! Sqrt(2)
  !> Square root of 2 in single precision
  REAL(sp), PARAMETER :: SQRT2   = 1.41421356237309504880168872420969807856967_sp
  !> Euler''s constant in double precision
  REAL(dp), PARAMETER :: EULER_D = 0.5772156649015328606065120900824024310422_dp     ! Euler
  !> Euler''s constant in single precision
  REAL(sp), PARAMETER :: EULER   = 0.5772156649015328606065120900824024310422_sp

  ! Standard file units
  !> Standard input file unit
  ! INTEGER, PARAMETER :: nin  = 5   ! standard input stream
  INTEGER, PARAMETER :: nin  = input_unit   ! standard input stream
  !> Standard output file unit
  ! INTEGER, PARAMETER :: nout = 6   ! standard output stream
  INTEGER, PARAMETER :: nout = output_unit   ! standard output stream
  !> Standard error file unit
  ! INTEGER, PARAMETER :: nerr = 0   ! error output stream
  INTEGER, PARAMETER :: nerr = error_unit   ! error output stream
  !> Standard file unit for namelist
  INTEGER, PARAMETER :: nnml = 100 ! namelist unit

END MODULE mo_constants
