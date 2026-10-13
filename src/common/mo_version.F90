!> \brief Compiled mHM version metadata, extracted from mo_file (Matthias Cuntz, Sebastian Mueller).
#ifndef MHMVERSION
#define MHMVERSION "0.0.0-dev0"
#endif
#ifndef MHMDATE
#define MHMDATE __DATE__
#endif
module mo_version
  implicit none
  private
  character(*), parameter, public :: version = MHMVERSION !< Compiled model version.
  character(*), parameter, public :: version_date = MHMDATE !< Model version date.
end module mo_version
