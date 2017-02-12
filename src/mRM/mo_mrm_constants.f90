!> \file mo_mrm_constants.f90
!> \brief Provides mRM specific constants
!> \details Provides mRM specific constants such as flood plain elevation.
!> \author Stephan Thober
!> \date Aug 2015
module mo_mrm_constants
  use mo_kind, only: i4, dp
  implicit none
  ! maximum number of outputs (fluxes states) for mHM
  integer(i4), public, parameter :: nOutFlxState       = 1_i4     ! max. number of outputs to write into a netcdf file
  ! computational
  integer(i4), public, parameter :: nodata_i4 = -9999_i4  ! [-]     global no data value
  real(dp),    public, parameter :: nodata_dp = -9999._dp ! [-]     global no data value
  !
  integer(i4), public, parameter :: nRoutingStates = 2 ! Dimension of the auxiliary vectors
  !                                                    ! which store current and past states of
  !                                                    ! incoming and outgoing of discharge at
  !                                                    ! a given node 
  !                                                    ! (1 - past)
  !                                                    ! (2 - current)
  integer(i4), public, parameter :: nColPars = 5_i4 ! number of properties of the global variables
#ifdef CYGWIN
  integer(i4), public, parameter :: maxNoGauges = 50_i4 ! maximal number of gauges allowed
#else
  integer(i4), public, parameter :: maxNoGauges = 200_i4 ! maximal number of gauges allowed
#endif
  integer(i4), public, parameter :: maxNoBasins = 50_i4 ! maximum number of allowed basins
  integer(i4), public, parameter :: maxNLcovers = 50_i4 ! maximum number of allowed LCover scenes
  ! temporal
  real(dp), public, parameter :: HourSecs = 3600.0_dp ! seconds per hour
  ! default inital values for states and fluxes as well as parameter fields
  real(dp), public, parameter :: P1_InitStateFluxes = 0.00_dp
  !
  ! parameters for routing
  real(dp), public, parameter :: rout_space_weight = 0._dp ! space weighting of routing is set to 0._dp
  !                                                        ! This parameter has no effect on the routing
  !                                                        ! results, see Thober et al. 2017
  !
  ! hydrological modelling
  real(dp), public, parameter :: deltaH = 5.000_dp  ! [m]     flood plain elevation, transept, above riverbed
  !
  ! timesteps in [s] that can be selected by adaptive time step
  ! these are multiples of 1 hour and 24 hours
  real(dp), dimension(19), parameter :: given_TS = &
       (/ 60._dp, 120._dp, 180._dp, 240._dp, 300._dp, 360._dp, &
       600._dp, 720._dp, 900._dp, 1200._dp, 1800._dp, 3600._dp, &
       7200._dp, 10800._dp, 14400._dp, 21600._dp, 28800._dp, 43200._dp, &
       86400._dp/)

end module mo_mrm_constants
