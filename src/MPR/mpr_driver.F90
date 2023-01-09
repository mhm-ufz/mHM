!> \file mpr_driver.f90
!> \brief Distributed precipitation-runoff model mHM
#ifdef MPR_STANDALONE
!> \details \copydetails mpr_driver

!> \brief Distributed precipitation-runoff model mHM
!> \details This is the main driver of mHM, which calls
!! one instance of mHM for a multiple domains and a given period.
!! \image html  mhm5-logo.png "Typical mHM cell"
!! \image latex mhm5-logo.pdf "Typical mHM cell" width=10cm
!> \changelog
!! - Robert Schweppe Jun 2018
!!   - refactored from mhm codebase
!> \authors Luis Samaniego & Rohini Kumar (UFZ)
!> \date Dec 2015
!> \version 0.1
!> \copyright (c)2005-2019, Helmholtz-Zentrum fuer Umweltforschung GmbH - UFZ.
!! All rights reserved.
!!
!! This code is a property of:
!!
!! ----------------------------------------------------------
!!
!! Helmholtz-Zentrum fuer Umweltforschung GmbH - UFZ
!! Registered Office: Leipzig
!! Registration Office: Amtsgericht Leipzig
!! Trade Register: Nr. B 4703
!! Chairman of the Supervisory Board: MinDirig Wilfried Kraus
!! Scientific Director: Prof. Dr. Georg Teutsch
!! Administrative Director: Dr. Heike Grassmann
!!
!! ----------------------------------------------------------
!!
!! NEITHER UFZ NOR THE DEVELOPERS MAKES ANY WARRANTY,
!! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE
!! OF THIS SOFTWARE. If software is modified to produce
!! derivative works, such modified software should be
!! clearly marked, so as not to confuse it with the version
!! available from UFZ.  This code can be used for research
!! purposes ONLY provided that the following sources are
!! acknowledged:
!!
!! Samaniego L., Kumar R., Attinger S. (2010): Multiscale
!! parameter regionalization of a grid-based hydrologic
!! model at the mesoscale.  Water Resour. Res., 46,
!! W05523, doi:10.1029/2008WR007327.
!!
!! Kumar, R., L. Samaniego, and S. Attinger (2013), Implications
!! of distributed hydrologic model parameterization on water
!! fluxes at multiple scales and locations, Water Resour. Res.,
!! 49, doi:10.1029/2012WR012195.
!!
!! For commercial applications you have to consult the
!! authorities of the UFZ.
!!
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mpr
program mpr_driver

  use mo_message, only : message
  use mo_string_utils, only : separator
  use mo_mpr_eval, only : mpr_eval
  use mo_read_wrapper, only : read_data
  use mo_mpr_read_config, only : mpr_read_config
  USE mo_common_read_config, ONLY : common_read_config                    ! Read main configuration files
  use mo_common_variables, only : mhmFileRestartOut, write_restart
  use mo_timer, only : timers_init
  use mo_mpr_startup, only : mpr_initialize
  use mo_mpr_restart, only : write_mpr_restart_files
  use mo_mpr_file, only : &
          file_namelist_mpr_param, unamelist_mpr_param, &      ! filename of namelist: mhm model parameter
          file_namelist_mpr, unamelist_mpr ! file containing main configurations
  use mo_kind, only: dp

  implicit none

  ! --------------------------------------------------------------------------
  ! INITIALIZE
  ! --------------------------------------------------------------------------
  call common_read_config(file_namelist_mpr, unamelist_mpr)
  call mpr_read_config(file_namelist_mpr, unamelist_mpr, file_namelist_mpr_param, unamelist_mpr_param)

  ! Start timings
  call timers_init

  call read_data()

  call mpr_initialize()

  ! -----------------------------------------------------------------------
  ! EXECUTION
  ! -----------------------------------------------------------------------
  call mpr_eval()

  ! --------------------------------------------------------------------------
  ! WRITE OUTPUT
  ! --------------------------------------------------------------------------
  if (write_restart) then
    call write_mpr_restart_files(mhmFileRestartOut)
  end if
  ! --------------------------------------------------------------------------
  ! FINISH UP
  ! --------------------------------------------------------------------------
  call message(separator)
  call message('MPR: Finished!')
  call message(separator)
end program mpr_driver
#else

!> \brief dummy module such that this file is never empty for compilation
module dummy_mpr
  implicit none
end module dummy_mpr
#endif
