!> \file mpr_driver.f90
! --------------------------------------------------------------------------
!> \authors   Luis Samaniego & Rohini Kumar (UFZ)
!  CONTACT    luis.samaniego@ufz.de / rohini.kumar@ufz.de
!> \version   0.1
!> \date      Dec 2015

!  PURPOSE
!>            \brief Distributed precipitation-runoff model mHM

!>            \details This is the main driver of mHM, which calls
!>             one instance of mHM for a multiple basins and a given period.

!>              \image html  mhm5-logo.png "Typical mHM cell"
!>              \image latex mhm5-logo.pdf "Typical mHM cell" width=10cm

!>  \copyright (c)2005-2018, Helmholtz-Zentrum fuer Umweltforschung GmbH - UFZ.
!>             All rights reserved.
!>
!>             This code is a property of:
!>
!>             ----------------------------------------------------------
!>
!>             Helmholtz-Zentrum fuer Umweltforschung GmbH - UFZ\n
!>             Registered Office: Leipzig\n
!>             Registration Office: Amtsgericht Leipzig\n
!>             Trade Register: Nr. B 4703\n
!>             Chairman of the Supervisory Board: MinDirig Wilfried Kraus\n
!>             Scientific Director: Prof. Dr. Georg Teutsch\n
!>             Administrative Director: Dr. Heike Grassmann\n
!>
!>             ----------------------------------------------------------
!>
!>             NEITHER UFZ NOR THE DEVELOPERS MAKES ANY WARRANTY,
!>             EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE
!>             OF THIS SOFTWARE. If software is modified to produce
!>             derivative works, such modified software should be
!>             clearly marked, so as not to confuse it with the version
!>             available from UFZ.  This code can be used for research
!>             purposes ONLY provided that the following sources are
!>             acknowledged:
!>
!>                Samaniego L., Kumar R., Attinger S. (2010): Multiscale
!>                parameter regionalization of a grid-based hydrologic
!>                model at the mesoscale.  Water Resour. Res., 46,
!>                W05523, doi:10.1029/2008WR007327.
!>
!>                Kumar, R., L. Samaniego, and S. Attinger (2013), Implications
!>                of distributed hydrologic model parameterization on water
!>                fluxes at multiple scales and locations, Water Resour. Res.,
!>                49, doi:10.1029/2012WR012195.
!>
!>             For commercial applications you have to consult the
!>             authorities of the UFZ.


! REDISTRIBUTION
!             Redistribution and use in source and binary forms,
!             with or without modification, are permitted provided that
!             the following conditions are met: * Redistributions of
!             source code must retain the above copyright notice, this
!             list of conditions and the following disclaimer.  *
!             Redistributions in binary form must reproduce the above
!             copyright notice, this list of conditions and the
!             following disclaimer in the documentation and/or other
!             materials provided with the distribution.  * Neither the
!             name of Helmholtz-Zentrum fuer Umweltforschung GmbH - UFZ,
!             nor the names of its contributors may be used to endorse
!             or promote products derived from this software without
!             specific prior written permission.

! DISCLAIM
!             THIS SOFTWARE IS PROVIDED BY HELMHOLTZ-ZENTRUM FUER
!             UMWELTFORSCHUNG GMBH - UFZ AND CONTRIBUTORS "AS IS" AND
!             ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!             LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
!             FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
!             EVENT SHALL THE HELMHOLTZ-ZENTRUM FUER UMWELTFORSCHUNG
!             GMBH - UFZ AND CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
!             INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!             CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
!             PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
!             DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!             CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!             CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
!             OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!             SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
!             DAMAGE.

!             This program is distributed in the hope that it will be
!             useful,but WITHOUT ANY WARRANTY; without even the implied
!             warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
!             PURPOSE. UFZ and the DEVELOPERS of this code do not take
!             any liabilities on the aplication of this code.

! --------------------------------------------------------------------------
!  HISTORY
!         Written          Sa 24.11.2005    v1.0 main structure
!         Modified  Robert Schweppe, Jun 2018 - refactored from mhm codebase
!
! --------------------------------------------------------------------------
#ifdef MPR_STANDALONE
program mpr_driver

  use mo_finish, only : finish
  use mo_mpr_eval, only : mpr_eval
  use mo_read_wrapper, only : read_data
  use mo_mpr_read_config, only : mpr_read_config
  USE mo_common_read_config, ONLY : common_read_config                    ! Read main configuration files
  use mo_common_variables, only : dirRestartOut, write_restart
  use mo_timer, only : timers_init
  use mo_mpr_startup, only : mpr_initialize
  use mo_mpr_restart, only : write_mpr_restart_files
  use mo_mpr_file, only : &
          file_namelist_mpr_param, unamelist_mpr_param, &      ! filename of namelist: mhm model parameter
          file_namelist_mpr, unamelist_mpr ! file containing main configurations
  use mo_mpr_global_variables, only : c2TSTu

  implicit none

  ! --------------------------------------------------------------------------
  ! INITIALIZE
  ! --------------------------------------------------------------------------
  call common_read_config(file_namelist_mpr, unamelist_mpr)
  call mpr_read_config(file_namelist_mpr, unamelist_mpr, file_namelist_mpr_param, unamelist_mpr_param)

  ! Start timings
  call timers_init

  call read_data()

  ! TODO: this might become part of the namelist in MPR-STANDALONE?!
  c2TSTu = 1.0_dp / 24.0_dp
  call mpr_initialize()

  ! -----------------------------------------------------------------------
  ! EXECUTION
  ! -----------------------------------------------------------------------
  call mpr_eval()

  ! --------------------------------------------------------------------------
  ! WRITE OUTPUT
  ! --------------------------------------------------------------------------
  if (write_restart) then
    call write_mpr_restart_files(dirRestartOut)
  end if
  ! --------------------------------------------------------------------------
  ! FINISH UP
  ! --------------------------------------------------------------------------
  call finish('MPR', 'Finished!')
end program mpr_driver
#else
! dummy module such that this file is never empty for compilation
module dummy_mpr
  implicit none
end module dummy_mpr
#endif
