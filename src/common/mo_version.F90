!> \file mo_version.F90
!> \brief \copybrief mo_version
!> \details \copydetails mo_version

#ifndef MHMVERSION
#define MHMVERSION "0.0.0-dev0"
#endif
#ifndef MHMDATE
#define MHMDATE __DATE__
#endif
!> \brief Compiled mHM version metadata.
!> \details Extracted from mo_file without its legacy file configuration.
!! Provides the compiled version and release date shared by the driver and
!! Python wrapper. Values are supplied by MHMVERSION and MHMDATE at build time.
!> \changelog
!! - Sebastian Mueller Sep 2020: set the version with the preprocessor from the version file.
!> \authors Matthias Cuntz, Sebastian Mueller
!> \date Jan 2012
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_common
module mo_version
  implicit none
  private
  character(*), parameter, public :: version = MHMVERSION !< Compiled model version.
  character(*), parameter, public :: version_date = MHMDATE !< Model version date.
end module mo_version
