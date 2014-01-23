!> \file mo_set_netcdf_outputs.f90

!> \brief Defines the structure of the netCDF to write the output in.

!> \details All output variables are initialized for the NetCDF.

!> \authors Matthias Zink
!> \date Apr 2013

MODULE mo_set_netcdf_outputs

  ! This module defines the NetCDF variables.

  ! Written Matthias Zink, Apr 2013

  USE mo_kind, ONLY: i4

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: set_netCDF

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !      NAME
  !          WriteFluxStateInit

  !>        \brief Runs mhm with a specific parameter set and returns required variables, e.g. runoff.

  !>        \details Runs mhm with a specific parameter set and returns required variables, e.g. runoff.

  !     INTENT(IN)
  !         None

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !

  !     EXAMPLE
  !        

  !     LITERATURE

  !     HISTORY
  !>        \author Matthias Zink
  !>        \date Apr 2013
  !         Modified, 
  subroutine set_netCDF(NoNetcdfVars, nrows, ncols)
    !
    use mo_NCWrite    , only: V, Dnc, ndims, nVars
    use netcdf        , only: NF90_UNLIMITED, NF90_CHAR, NF90_INT, NF90_DOUBLE

    implicit none
    !
    integer(i4),                         intent(in) :: NoNetcdfVars
    integer(i4),                         intent(in) :: nrows
    integer(i4),                         intent(in) :: ncols
    !
    integer(i4)                                     :: i
    !
    ! define output file
    !
    ! define parameters
    nDims  = 3                                                ! nr. dim types
    nVars  = 5 + NoNetcdfVars                                 ! total nr. var to print
    !
    ! allocate arrays
    if (allocated(Dnc)   ) deallocate ( Dnc)
    if (allocated(V)     ) deallocate ( V  )
    allocate ( Dnc(nDims)         )
    allocate ( V(nVars)           )
    !  
    ! define dimensions (coordinates) => corresponding variable must be created 
    i              = 1
    Dnc(i)%name      = "easting"
    Dnc(i)%len       = nrows
    !
    i              = 2
    Dnc(i)%name      = "northing"
    Dnc(i)%len       = ncols
    !
    i              = 3
    Dnc(i)%name      = "time"
    Dnc(i)%len       = NF90_UNLIMITED

    !
    !
    ! DIMENSION VARIABLES
    i                =  1
    V(i)%name        =  Dnc(i)%name
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  0
    V(i)%nSubs       =  0
    V(i)%nDims       =  1
    V(i)%dimTypes    =  (/i,0,0,0,0/)
    V(i)%wFlag       =  .true.
    ! attributes (other possibilities: add_offset, valid_min, valid_max)  
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "units"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values  = "m"
    !
    V(i)%att(2)%name   = "long_name"
    V(i)%att(2)%xType  = NF90_CHAR
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = "x-coordinate in cartesian coordinates GK4"   
    !
    i                =  2
    V(i)%name        =  Dnc(i)%name
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  0
    V(i)%nSubs       =  0
    V(i)%nDims       =  1
    V(i)%dimTypes    =  (/i,0,0,0,0/)
    ! printing
    V(i)%wFlag       =  .true.
    ! attributes
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "units"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values  = "m"
    !
    V(i)%att(2)%name   = "long_name"
    V(i)%att(2)%xType  = NF90_CHAR
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = "y-coordinate in cartesian coordinates GK4"
    !
    i                =  3
    V(i)%name        =  Dnc(i)%name
    V(i)%xType       =  NF90_INT
    V(i)%nLvls       =  0
    V(i)%nSubs       =  0
    V(i)%nDims       =  1
    V(i)%dimTypes    =  (/i,0,0,0,0/)
    V(i)%wFlag       =  .true.
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "units"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    !
    V(i)%att(2)%name   = "long_name"
    V(i)%att(2)%xType  = NF90_CHAR
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = "time"
    !
    i                =  4
    V(i)%name        =  "lon"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    V(i)%wFlag       =  .true.
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "units"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "degrees_east"
    !
    V(i)%att(2)%name   = "long_name"
    V(i)%att(2)%xType  = NF90_CHAR
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = "longitude"
    !
    i                =  5
    V(i)%name        =  "lat"
    V(i)%xType       =  NF90_DOUBLE
    V(i)%nLvls       =  1
    V(i)%nSubs       =  1
    V(i)%nDims       =  2
    V(i)%dimTypes    =  (/1,2,0,0,0/)
    V(i)%wFlag       =  .true.
    ! attributes 
    V(i)%nAtt          = 2
    !
    V(i)%att(1)%name   = "units"
    V(i)%att(1)%xType  = NF90_CHAR
    V(i)%att(1)%nValues= 1
    V(i)%att(1)%values = "degrees_north"
    !
    V(i)%att(2)%name   = "long_name"
    V(i)%att(2)%xType  = NF90_CHAR
    V(i)%att(2)%nValues= 1
    V(i)%att(2)%values = "latitude"
    !
    ! ****************************** 
    !  dynamic variables data
    ! ******************************
    do i = 6,  NoNetcdfVars + 5
       V(i)%xType       =  NF90_DOUBLE
       V(i)%nLvls       =  1
       V(i)%nSubs       =  1
       V(i)%nDims       =  3
       V(i)%dimTypes    =  (/1,2,3,0,0/)
       ! printing
       V(i)%wFlag       =  .true.
       ! pointer      
       !                   during running time 
       ! attributes 
       V(i)%nAtt          = 6
       !
       V(i)%att(1)%name   = "units"
       V(i)%att(1)%xType  = NF90_CHAR
       V(i)%att(1)%nValues= 1
       !
       V(i)%att(2)%name   = "long_name"
       V(i)%att(2)%xType  = NF90_CHAR
       V(i)%att(2)%nValues= 1
       !
       V(i)%att(3)%name   = "scale_factor"
       V(i)%att(3)%xType  = NF90_DOUBLE
       V(i)%att(3)%nValues= 1
       !
       V(i)%att(4)%name   = "_FillValue"
       V(i)%att(4)%xType  = NF90_DOUBLE
       V(i)%att(4)%nValues= 1
       !
       V(i)%att(5)%name   = "missing_value"
       V(i)%att(5)%xType  = NF90_DOUBLE
       V(i)%att(5)%nValues= 1
       ! !
       V(i)%att(6)%name   = "coordinates"
       V(i)%att(6)%xType  = NF90_CHAR
       V(i)%att(6)%nValues= 1
    end do
    !
    !
  end subroutine set_netCDF

END MODULE mo_set_netcdf_outputs
