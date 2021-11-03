!> \file mo_netcdf.f90

!> \brief NetCDF Fortran 90 interface wrapper

!> \details A wrapper around the NetCDF Fortran 90 interface.
!> \copied from https://github.com/schaefed/mo_netcdf (commit d903514da5c96aee263c2a47c3d6b9df1dbf6ef7)
!
!> \authors David Schaefer
!> \date Jun 2015

module mo_netcdf

  ! This module provides a thin wrapper around the NetCDF Fortran 90 interface,
  ! following a object-oriented approach.

  ! Written  David Schaefer, Jun 2015
  ! Modified Matthias Cuntz, Jan 2016 - compiled with PGI Fortran rev 15.9 - no automatic allocation of left-hand-side
  ! Modified Ricardo Torres, Feb 2017 - add derived type NcGroup and NcAttributable. NcAttributable is the base derived type,
  !                                     NcGroup and NcVariable are extended from it. NcDataset is extended from NcGroup. No more
  !                                     duplicated routines to set attributes.
  ! Modified Robert Schweppe, Nov 2018 - adapted imports: changed from mo_types to mo_kind

  ! License
  ! -------
  ! GNU Lesser General Public License http://www.gnu.org/licenses/

  use mo_kind, only: i1, i2, i4, i8, sp, dp
  use netcdf,          only: &
       nf90_open, nf90_close, nf90_strerror, nf90_def_dim, nf90_def_var,   &
       nf90_put_var, nf90_get_var, nf90_put_att, nf90_get_att,             &
       nf90_inquire, nf90_inq_dimid, nf90_inquire_dimension,               &
       nf90_inq_varid, nf90_inq_varids, nf90_inquire_variable, nf90_inquire_attribute,     &
       nf90_inq_ncid, nf90_inq_grp_parent, nf90_inq_grpname, nf90_def_grp, &
       nf90_rename_dim, nf90_rename_var, nf90_rename_att, nf90_sync,       &
       NF90_OPEN, NF90_NETCDF4, NF90_CREATE, NF90_WRITE, NF90_NOWRITE,     &
       NF90_BYTE, NF90_SHORT, NF90_INT, NF90_INT64, NF90_FLOAT, NF90_DOUBLE,           &
       NF90_FILL_BYTE, NF90_FILL_SHORT, NF90_FILL_INT, NF90_FILL_FLOAT , NF90_FILL_DOUBLE, &
       NF90_NOERR, NF90_UNLIMITED, NF90_GLOBAL, NF90_SHARE, NF90_HDF5, &
       NF90_64BIT_OFFSET, NF90_CLASSIC_MODEL


  implicit none

  ! --------------------------------------------------------------------------------------

  type, abstract :: NcBase

     integer(i4) :: id

   contains

     procedure(getNameInterface), deferred   :: getName
     procedure(getParentInterface), deferred :: getParent

  end type NcBase

  type, abstract, extends(NcBase) :: NcAttributable

   contains

     procedure, public  :: hasAttribute
     procedure, public  :: renameAttribute

     procedure, private :: getAttributableIds
     procedure, private :: setAttributeChar
     procedure, private :: setAttributeI8
     procedure, private :: setAttributeI16
     procedure, private :: setAttributeI32
     procedure, private :: setAttributeI64
     procedure, private :: setAttributeF32
     procedure, private :: setAttributeF64

     procedure, private :: getAttributeChar
     procedure, private :: getAttributeI8
     procedure, private :: getAttributeI16
     procedure, private :: getAttributeI32
     procedure, private :: getAttributeI64
     procedure, private :: getAttributeF32
     procedure, private :: getAttributeF64

     generic, public :: setAttribute => &
          setAttributeChar,  &
          setAttributeI8,    &
          setAttributeI16,   &
          setAttributeI32,   &
          setAttributeI64,   &
          setAttributeF32,   &
          setAttributeF64

     generic, public :: getAttribute => &
          getAttributeChar,  &
          getAttributeI8,    &
          getAttributeI16,   &
          getAttributeI32,   &
          getAttributeI64,   &
          getAttributeF32,   &
          getAttributeF64
  end type NcAttributable

  ! --------------------------------------------------------------------------------------

  type, extends(NcAttributable) :: NcGroup

  contains

     ! getter
     procedure, private :: getVariableIds
     procedure, public  :: getVariables
     procedure, public  :: getUnlimitedDimension
     procedure, public  :: getNoVariables

     procedure, private :: getDimensionByName
     procedure, private :: getDimensionById

     procedure, public  :: getParent => getGroupParent
     procedure, public  :: getName => getGroupName
     procedure, public  :: getGroup => getGroupByName
     procedure, public  :: getVariable => getVariableByName
     generic,   public  :: getDimension => &
          getDimensionById, &
          getDimensionByName

     ! checker
     procedure, public  :: hasVariable
     procedure, public  :: hasDimension
     procedure, public  :: hasGroup
     procedure, public  :: isUnlimited => isDatasetUnlimited

     ! setter
     procedure, public  :: setGroup
     !procedure, private :: setDimension_1Dbounds
     !procedure, private :: setDimension_2Dbounds
     procedure, private :: setDimension_
     procedure, public  :: setDimension
     procedure, private :: setVariableWithTypes
     procedure, private :: setVariableWithNames
     procedure, private :: setVariableWithIds

     generic,   public  :: setVariable => &
          setVariableWithNames, &
          setVariableWithTypes, &
          setVariableWithIds

  end type NcGroup

  interface NcGroup
     procedure newNcGroup
  end interface NcGroup

  ! --------------------------------------------------------------------------------------

  type, extends(NcGroup) :: NcDataset

     character(256) :: fname !> Filename of the opened dataset
     character(1)   :: mode  !> File open mode

   contains

     procedure, public :: sync
     procedure, public :: close

  end type NcDataset

  interface NcDataset
     procedure newNcDataset
  end interface NcDataset

! --------------------------------------------------------------------------------------

  type, extends(NcBase) :: NcDimension

     type(NcGroup)     :: parent  !> The dimension's parent

   contains

     procedure, public :: renameDimension
     procedure, public :: getParent => getDimensionParent
     procedure, public :: getName => getDimensionName
     procedure, public :: getLength => getDimensionLength

     procedure, public :: isUnlimited => isUnlimitedDimension

  end type NcDimension

  interface NcDimension
     procedure newNcDimension
  end interface NcDimension

! --------------------------------------------------------------------------------------

  type, extends(NcAttributable) :: NcVariable

     type(NcGroup)      :: parent   !> The variables's parent

   contains


     procedure, public  :: renameVariable
     procedure, public  :: getParent => getVariableParent
     procedure, public  :: getName => getVariableName
     procedure, private :: getSlicingShape

     procedure, private :: setDataScalarI8
     procedure, private :: setData1dI8
     procedure, private :: setData2dI8
     procedure, private :: setData3dI8
     procedure, private :: setData4dI8
     procedure, private :: setData5dI8
     procedure, private :: setDataScalarI16
     procedure, private :: setData1dI16
     procedure, private :: setData2dI16
     procedure, private :: setData3dI16
     procedure, private :: setData4dI16
     procedure, private :: setData5dI16
     procedure, private :: setDataScalarI32
     procedure, private :: setData1dI32
     procedure, private :: setData2dI32
     procedure, private :: setData3dI32
     procedure, private :: setData4dI32
     procedure, private :: setData5dI32
     procedure, private :: setDataScalarI64
     procedure, private :: setData1dI64
     procedure, private :: setData2dI64
     procedure, private :: setData3dI64
     procedure, private :: setData4dI64
     procedure, private :: setData5dI64
     procedure, private :: setDataScalarF32
     procedure, private :: setData1dF32
     procedure, private :: setData2dF32
     procedure, private :: setData3dF32
     procedure, private :: setData4dF32
     procedure, private :: setData5dF32
     procedure, private :: setDataScalarF64
     procedure, private :: setData1dF64
     procedure, private :: setData2dF64
     procedure, private :: setData3dF64
     procedure, private :: setData4dF64
     procedure, private :: setData5dF64

     procedure, private :: getDataScalarI8
     procedure, private :: getData1dI8
     procedure, private :: getData2dI8
     procedure, private :: getData3dI8
     procedure, private :: getData4dI8
     procedure, private :: getData5dI8
     procedure, private :: getDataScalarI16
     procedure, private :: getData1dI16
     procedure, private :: getData2dI16
     procedure, private :: getData3dI16
     procedure, private :: getData4dI16
     procedure, private :: getData5dI16
     procedure, private :: getDataScalarI32
     procedure, private :: getData1dI32
     procedure, private :: getData2dI32
     procedure, private :: getData3dI32
     procedure, private :: getData4dI32
     procedure, private :: getData5dI32
     procedure, private :: getDataScalarI64
     procedure, private :: getData1dI64
     procedure, private :: getData2dI64
     procedure, private :: getData3dI64
     procedure, private :: getData4dI64
     procedure, private :: getData5dI64
     procedure, private :: getDataScalarF32
     procedure, private :: getData1dF32
     procedure, private :: getData2dF32
     procedure, private :: getData3dF32
     procedure, private :: getData4dF32
     procedure, private :: getData5dF32
     procedure, private :: getDataScalarF64
     procedure, private :: getData1dF64
     procedure, private :: getData2dF64
     procedure, private :: getData3dF64
     procedure, private :: getData4dF64
     procedure, private :: getData5dF64

     procedure, private :: setVariableFillValueI8
     procedure, private :: setVariableFillValueI16
     procedure, private :: setVariableFillValueI32
     procedure, private :: setVariableFillValueI64
     procedure, private :: setVariableFillValueF32
     procedure, private :: setVariableFillValueF64

     procedure, private :: getVariableFillValueI8
     procedure, private :: getVariableFillValueI16
     procedure, private :: getVariableFillValueI32
     procedure, private :: getVariableFillValueI64
     procedure, private :: getVariableFillValueF32
     procedure, private :: getVariableFillValueF64

     procedure, public  :: getNoDimensions

     procedure, public  :: getDimensions => getVariableDimensions

     procedure, public  :: getRank       => getVariableRank

     procedure, public  :: getShape      => getVariableShape

     procedure, public  :: getDtype      => getVariableDtype

     procedure, public  :: isUnlimited   => isUnlimitedVariable

     generic, public :: setData => &
          setDataScalarI8, &
          setData1dI8, &
          setData2dI8, &
          setData3dI8, &
          setData4dI8, &
          setData5dI8, &
          setDataScalarI16, &
          setData1dI16, &
          setData2dI16, &
          setData3dI16, &
          setData4dI16, &
          setData5dI16, &
          setDataScalarI32, &
          setData1dI32, &
          setData2dI32, &
          setData3dI32, &
          setData4dI32, &
          setData5dI32, &
          setDataScalarI64, &
          setData1dI64, &
          setData2dI64, &
          setData3dI64, &
          setData4dI64, &
          setData5dI64, &
          setDataScalarF32, &
          setData1dF32, &
          setData2dF32, &
          setData3dF32, &
          setData4dF32, &
          setData5dF32, &
          setDataScalarF64, &
          setData1dF64, &
          setData2dF64, &
          setData3dF64, &
          setData4dF64, &
          setData5dF64

     generic, public :: getData => &
          getDataScalarI8, &
          getData1dI8, &
          getData2dI8, &
          getData3dI8, &
          getData4dI8, &
          getData5dI8, &
          getDataScalarI16, &
          getData1dI16, &
          getData2dI16, &
          getData3dI16, &
          getData4dI16, &
          getData5dI16, &
          getDataScalarI32, &
          getData1dI32, &
          getData2dI32, &
          getData3dI32, &
          getData4dI32, &
          getData5dI32, &
          getDataScalarI64, &
          getData1dI64, &
          getData2dI64, &
          getData3dI64, &
          getData4dI64, &
          getData5dI64, &
          getDataScalarF32, &
          getData1dF32, &
          getData2dF32, &
          getData3dF32, &
          getData4dF32, &
          getData5dF32, &
          getDataScalarF64, &
          getData1dF64, &
          getData2dF64, &
          getData3dF64, &
          getData4dF64, &
          getData5dF64

     generic, public :: setFillValue => &
          setVariableFillValueI8,  &
          setVariableFillValueI16, &
          setVariableFillValueI32, &
          setVariableFillValueI64, &
          setVariableFillValueF32, &
          setVariableFillValueF64

     generic, public :: getFillValue => &
          getVariableFillValueI8,  &
          getVariableFillValueI16, &
          getVariableFillValueI32, &
          getVariableFillValueI64, &
          getVariableFillValueF32, &
          getVariableFillValueF64

  end type NcVariable

  interface NcVariable
    procedure newNcVariable
  end interface NcVariable
  ! --------------------------------------------------------------------------------------

  ! abstract interfaces
  interface
     function getNameInterface(self)
       import NcBase
       class(NcBase), intent(in) :: self
       character(len=256)        :: getNameInterface
     end function getNameInterface

     function getParentInterface(self)
       import NcBase, NcGroup
       class(NcBase), intent(in) :: self
       type(NcGroup)             :: getParentInterface
     end function getParentInterface
  end interface

  interface operator (==)
     procedure equalNcBases
  end interface operator (==)

contains

  function newNcDataset(fname, fmode, cmode) result(out)
    character(*), intent(in)              :: fname
    character(1), intent(in)              :: fmode
    character(*), intent(inout), optional :: cmode
    integer(i4)                          :: status
    type(NcDataset)                       :: out

    select case(fmode)
    case("w")
       status = nf90_create(trim(fname), getCreationMode(cmode), out%id)
    case("r")
       status = nf90_open(trim(fname), NF90_NOWRITE, out%id)
    case("a")
       status = nf90_open(trim(fname), NF90_WRITE, out%id)
    case default
       write(*,*) "Mode argument must be in 'w','r','a' ! "
       stop 1
    end select
    call check(status,"Failed to open file: " // fname)

    out%fname = fname
    out%mode  = fmode
  end function newNcDataset

  function newNcVariable(id, parent) result(out)
    integer(i4) , intent(in) :: id
    type(NcGroup), intent(in) :: parent
    type(NcVariable)          :: out

    out%id     = id
    out%parent = parent
  end function newNcVariable

  function newNcDimension(id, parent) result(out)
    integer(i4) , intent(in) :: id
    type(NcGroup), intent(in) :: parent
    type(NcDimension)         :: out

    out%id     = id
    out%parent = parent
  end function newNcDimension

  function newNcGroup(id) result(out)
    integer(i4)    , intent(in) :: id
    type(NcGroup)                :: out

    out%id = id
  end function newNcGroup

  subroutine sync(self)
    class(NcDataset) :: self

    call check(nf90_sync(self%id), "Failed to sync file: "//self%fname)
  end subroutine sync

  subroutine close(self)
    class(NcDataset) :: self

    call check(nf90_close(self%id), "Failed to close file: "//self%fname)
  end subroutine close

  function setGroup(self, name)
    class(NcGroup), intent(inout) :: self
    character(*)  , intent(in)    :: name
    integer(i4)                  :: id
    type(NcGroup)                 :: setGroup

    call check(nf90_def_grp(self%id, name, id), "Failed to create new group: " // name)
    setGroup = NcGroup(id)
  end function setGroup

  function getGroupParent(self)
    class(NcGroup), intent(in) :: self
    integer(i4)               :: id
    type(NcGroup)              :: getGroupParent

    call check(nf90_inq_grp_parent(self%id, id), "Failed to get parent group of: " // self%getName())
    getGroupParent = NcGroup(id)
  end function getGroupParent

  function getGroupName(self)
    class(NcGroup), intent(in) :: self
    character(256)             :: getGroupName

    call check(nf90_inq_grpname(self%id, getGroupName), "Failed to inquire group name")
  end function getGroupName

  function getNoVariables(self)
    class(NcGroup), intent(in) :: self
    integer(i4)               :: getNoVariables

    call check(nf90_inquire(self%id, nvariables=getNoVariables), "Failed inquire number of variables")
  end function getNoVariables

  function getDimensionParent(self)
    class(NcDimension), intent(in) :: self
    type(NcGroup)                  :: getDimensionParent

    getDimensionParent = self%parent
  end function getDimensionParent

  function getVariableParent(self)
    class(NcVariable), intent(in) :: self
    type(NcGroup)                  :: getVariableParent

    getVariableParent = self%parent
  end function getVariableParent

  function getVariableIds(self)
    class(NcGroup), intent(in)              :: self
    integer(i4), dimension(:), allocatable :: getVariableIds
    integer(i4)                            :: tmp

    allocate(getVariableIds(self%getNoVariables()))
    call check(nf90_inq_varids(self%id, tmp, getVariableIds), "Failed to inquire variable ids")
  end function getVariableIds

  function getVariables(self)
    class(NcGroup), intent(in)                  :: self
    type(NcVariable), dimension(:), allocatable :: getVariables
    integer(i4), dimension(:), allocatable     :: varids
    integer(i4)                                :: i, nvars

    nvars = self%getNoVariables()
    allocate(getVariables(nvars), varids(nvars))

    varids = self%getVariableIds()
    do i=1,size(varids)
       getVariables(i) = NcVariable(varids(i), self)
    end do

  end function getVariables

  function getDimensionName(self)
    class(NcDimension), intent(in) :: self
    character(len=256)             :: getDimensionName

    call check(nf90_inquire_dimension(self%parent%id, self%id, name=getDimensionName), &
         "Failed to inquire dimension name")
  end function getDimensionName

  function getDimensionLength(self)
    class(NcDimension), intent(in) :: self
    integer(i4)                    :: getDimensionLength

    call check(nf90_inquire_dimension(self%parent%id,self%id,len=getDimensionLength),&
         "Failed to inquire dimension: "//self%getName())
  end function getDimensionLength

  function isDatasetUnlimited(self)
    class(NcGroup), intent(in)   :: self
    logical                      :: isDatasetUnlimited
    integer(i4)                  :: dimid

    call check(nf90_inquire(self%id,unlimitedDimId=dimid), &
         "Failed to inquire group "//self%getName())
    isDatasetUnlimited = (dimid .ne. -1)
  end function isDatasetUnlimited

  function getUnlimitedDimension(self)
    class(NcGroup), intent(in)   :: self
    type(NcDimension)            :: getUnlimitedDimension
    integer(i4)                  :: dimid

    call check(nf90_inquire(self%id,unlimitedDimId=dimid), &
         "Failed to inquire group "//self%getName())

    if (dimid .eq. -1) then
       write(*,*) "Dataset has no unlimited dimension"
       stop 1
    end if

    getUnlimitedDimension = self%getDimension(dimid)
  end function getUnlimitedDimension

  function equalNcBases(left, right) result(out)
    class(NcBase), intent(in) :: left, right
    logical                   :: out

    out = (left%id .eq. right%id)
  end function equalNcBases

  function isUnlimitedDimension(self)
    class(NcDimension), intent(in) :: self
    logical                        :: isUnlimitedDimension

    isUnlimitedDimension = .false.
    if (self%parent%isUnlimited()) then
       isUnlimitedDimension = (self == self%parent%getUnlimitedDimension())
    end if
  end function isUnlimitedDimension

  function setDimension_(self, name, length)
    class(NcGroup), intent(in) :: self
    character(*)  , intent(in) :: name
    integer(i4)   , intent(in), optional :: length

    type(NcDimension)          :: setDimension_
    integer(i4)                :: id, dimlength

    dimlength = NF90_UNLIMITED
    if (present(length)) then
      if (length .le. 0) then
       dimlength = NF90_UNLIMITED
      else
        dimlength = length
      end if
    end if

    call check(nf90_def_dim(self%id, name, dimlength, id), &
         "Failed to create dimension: " // name)

    setDimension_ = NcDimension(id,self)

  end function setDimension_

  function setDimension(self, name, length, bounds, reference, attribute_names, attribute_values)
    class(NcGroup), intent(in) :: self
    character(*)  , intent(in) :: name
    integer(i4)   , intent(in), optional :: length
    real(dp)      , intent(in), optional, dimension(:) :: bounds
    integer(i4)   , intent(in), optional :: reference
    character(256) , intent(in), optional, dimension(:) :: attribute_names
    character(2048) , intent(in), optional, dimension(:) :: attribute_values

    type(NcDimension)          :: setDimension, bnds_dim
    type(NcVariable)           :: nc_var
    integer(i4)                :: dimlength, reference_default, iAtt
    character(256)              :: dim_bound_name
    real(dp), allocatable, dimension(:, :) :: bound_data

    ! set the new ncDimension (integer values and name)
    setDimension = self%setDimension_(name, length)

    if (present(bounds)) then
      ! init
      dimlength = size(bounds)
      reference_default = 1_i4
      if (present(reference)) then
        reference_default = reference
      end if
      ! here we set the reference to ncDimension for labelled ncDimension which in fact is a variable
      nc_var = self%setVariable(name, "f64", [setDimension])
      ! write the data based on the type of reference
      select case(reference_default)
      case(0_i4)
        ! set the start values
        call nc_var%setData(bounds(1:dimlength - 1))
        call nc_var%setAttribute("cell-reference", "start-bound")
      case(1_i4)
        ! set the center values
        call nc_var%setData((bounds(2:dimlength) + bounds(1:dimlength-1)) / 2.0_dp)
        call nc_var%setAttribute("cell-reference", "center-bound")
      case(2_4)
        ! set the end values
        call nc_var%setData(bounds(2:dimlength))
        call nc_var%setAttribute("cell-reference", "end-bound")
      case default
        write(*,*) "reference id for set_Dimension is unknown"
        stop 1
      end select
      ! set attributes
      ! already set attributes
      if (present(attribute_names) .and. present(attribute_values)) then
        do iAtt = 1, size(attribute_names)
          call nc_var%setAttribute(trim(attribute_names(iAtt)), &
                  trim(attribute_values(iAtt)))
        end do
      end if
      ! --- bounds ---
      ! allocate array for data
      allocate(bound_data(dimlength-1, 2_i4))
      ! create dimension name for bounds
      dim_bound_name  = trim(name) // "_bnds"
      ! set the dimensions used for the bounds array
      if (self%hasDimension("bnds")) then
        ! add it to our bounds of ncDimensions for the current array
        bnds_dim = self%getDimension("bnds")
      else
        bnds_dim = self%setDimension_("bnds", 2)
      end if
      nc_var = self%setVariable(dim_bound_name, "f64", [setDimension, bnds_dim])
      bound_data(:, 1) = bounds(1 : dimlength - 1)
      bound_data(:, 2) = bounds(2 : dimlength)
      call nc_var%setData(bound_data)
      deallocate(bound_data)

    end if

  end function setDimension

  ! function setDimension_2Dbounds(self, name, length, bounds, reference, attribute_names, attribute_values)
  !   class(NcGroup), intent(in) :: self
  !   character(*)  , intent(in) :: name
  !   integer(i4)   , intent(in), optional :: length
  !   real(dp)      , intent(in), optional, dimension(:, :) :: bounds
  !   integer(i4)   , intent(in), optional :: reference
  !   character(256) , intent(in), optional, dimension(:) :: attribute_names
  !   character(2048) , intent(in), optional, dimension(:) :: attribute_values
  !
  !   type(NcDimension)          :: setDimension_2Dbounds, bnds_dim
  !   type(NcVariable)           :: nc_var
  !   integer(i4)                :: dimlength, reference_default, iAtt
  !   character(256)              :: dim_bound_name
  !
  !   ! set the new ncDimension (integer values and name)
  !   setDimension_2Dbounds = setDimension_(self, name, length)
  !
  !   if (present(bounds)) then
  !     ! init
  !     dimlength = size(bounds, 1)
  !     reference_default = 1_i4
  !     if (present(reference)) then
  !       reference_default = reference
  !     end if
  !     ! here we set the reference to ncDimension for labelled ncDimension which in fact is a variable
  !     nc_var = self%setVariable(name, "f64", [setDimension_2Dbounds])
  !     ! write the data based on the type of reference
  !     select case(reference_default)
  !     case(0_i4)
  !       ! set the start values
  !       call nc_var%setData(bounds(:, 1))
  !     case(1_i4)
  !       ! set the center values
  !       call nc_var%setData((bounds(:, 1) + bounds(:, 2)) / 2.0_dp)
  !     case(2_4)
  !       ! set the end values
  !       call nc_var%setData(bounds(:, 2))
  !     case default
  !       write(*,*) "reference id for set_Dimension is unknown"
  !       stop 1
  !     end select
  !     ! set attributes
  !     ! already set attributes
  !     if (present(attribute_names) .and. present(attribute_values)) then
  !       do iAtt = 1, size(attribute_names)
  !         call nc_var%setAttribute(trim(attribute_names(iAtt)), &
  !                 trim(attribute_values(iAtt)))
  !       end do
  !     end if
  !     ! --- bounds ---
  !     ! create dimension name for bounds
  !     dim_bound_name  = trim(name) // "_bnds"
  !     ! set the dimensions used for the bounds array
  !     if (self%hasDimension("bnds")) then
  !       ! add it to our bounds of ncDimensions for the current array
  !       bnds_dim = self%getDimension("bnds")
  !     else
  !       bnds_dim = self%setDimension("bnds", 2)
  !     end if
  !     nc_var = self%setVariable(dim_bound_name, "f64", [setDimension_2Dbounds, bnds_dim])
  !     call nc_var%setData(bounds)
  !
  !   end if
  !
  ! end function setDimension_2Dbounds

  function hasVariable(self, name)
    class(NcGroup), intent(in) :: self
    character(*)    , intent(in) :: name
    logical                      :: hasVariable
    integer(i4)                  :: tmpid

    hasVariable = (nf90_inq_varid(self%id,name,tmpid) .eq. NF90_NOERR)
  end function hasVariable

  function hasDimension(self, name)
    class(NcGroup), intent(in) :: self
    character(*)    , intent(in) :: name
    logical                      :: hasDimension
    integer(i4)                  :: tmpid

    hasDimension = (nf90_inq_dimid(self%id,name,tmpid) .eq. NF90_NOERR)
  end function hasDimension

  function hasGroup(self, name)
    class(NcGroup), intent(in) :: self
    character(*)    , intent(in) :: name
    logical                      :: hasGroup
    integer(i4)                  :: tmpid

    hasGroup = (nf90_inq_ncid(self%id,name,tmpid) .eq. NF90_NOERR)
  end function hasGroup

  function setVariableWithIds(self, name, dtype, dimensions, contiguous, &
       chunksizes, deflate_level, shuffle, fletcher32, endianness, &
       cache_size, cache_nelems, cache_preemption)
    class(NcGroup), intent(in)           :: self
    character(*)    , intent(in)           :: name
    character(*)    , intent(in)           :: dtype
    integer(i4)     , intent(in)           :: dimensions(:)
    logical         , intent(in), optional :: contiguous,shuffle, fletcher32
    integer(i4)     , intent(in), optional :: endianness,deflate_level,cache_size, &
         cache_nelems, cache_preemption, chunksizes(:)
    type(NcVariable)                       :: setVariableWithIds
    integer(i4)                            :: varid, status

    status = nf90_def_var(self%id, name, getDtypeFromString(dtype), dimensions, varid, contiguous, &
         chunksizes, deflate_level, shuffle, fletcher32, endianness, &
         cache_size, cache_nelems, cache_preemption)
    call check(status, "Failed to create variable: " // name)
    setVariableWithIds = NcVariable(varid, self)
  end function setVariableWithIds

  function setVariableWithNames(self, name, dtype, dimensions, contiguous, &
       chunksizes, deflate_level, shuffle, fletcher32, endianness, &
       cache_size, cache_nelems, cache_preemption)

    class(NcGroup), intent(in)              :: self
    character(*)    , intent(in)              :: name
    character(*)    , intent(in)              :: dtype
    character(*)    , intent(in)              :: dimensions(:)
    logical         , intent(in), optional    :: contiguous,shuffle, fletcher32
    integer(i4)     , intent(in), optional    :: endianness,deflate_level,cache_size, &
         cache_nelems, cache_preemption, chunksizes(:)
    type(NcVariable)                          :: setVariableWithNames
    type(NcDimension)                         :: dim
    integer(i4)                               :: i, dimids(size(dimensions))

    do i = 1,size(dimensions)
       dim = self%getDimension(dimensions(i))
       dimids(i) = dim%id
    end do

    setVariableWithNames = setVariableWithIds(self, name, dtype, dimids, contiguous, &
         chunksizes, deflate_level, shuffle, fletcher32, endianness, &
         cache_size, cache_nelems, cache_preemption)
  end function setVariableWithNames

  function setVariableWithTypes(self, name, dtype, dimensions, contiguous, &
       chunksizes, deflate_level, shuffle, fletcher32, endianness, &
       cache_size, cache_nelems, cache_preemption)
    class(NcGroup)   , intent(in)              :: self
    character(*)     , intent(in)              :: name
    character(*)     , intent(in)              :: dtype
    type(NcDimension), intent(in)              :: dimensions(:)
    logical          , intent(in), optional    :: contiguous,shuffle, fletcher32
    integer(i4)      , intent(in), optional    :: endianness,deflate_level,cache_size, &
         cache_nelems, cache_preemption, chunksizes(:)
    type(NcVariable)                           :: setVariableWithTypes
    type(NcDimension)                          :: dim
    integer(i4)                                :: i, dimids(size(dimensions))

    do i = 1,size(dimensions)
       dim = dimensions(i)
       dimids(i) = dim%id
    end do

    setVariableWithTypes = setVariableWithIds(self, name, dtype, dimids, contiguous, &
         chunksizes, deflate_level, shuffle, fletcher32, endianness, &
         cache_size, cache_nelems, cache_preemption)
  end function setVariableWithTypes

  function getDimensionById(self, id)
    class(NcGroup), intent(in) :: self
    integer(i4)                  :: id
    type(NcDimension)            :: getDimensionById
    character(32)                :: msg, name

    write(msg,*) id
    call check(nf90_inquire_dimension(self%id,id,name), &
         "Could not inquire dimension: " // msg)
    getDimensionById = NcDimension(id,self)
  end function getDimensionById

  function getDimensionByName(self, name)
    class(NcGroup), intent(in) :: self
    character(*)                 :: name
    type(NcDimension)            :: getDimensionByName
    integer(i4)                  :: id

    call check(nf90_inq_dimid(self%id,name,id), &
         "Could not inquire dimension: " // name)
    getDimensionByName = self%getDimensionById(id)
  end function getDimensionByName

  function getGroupByName(self, name)
    class(NcGroup), intent(in) :: self
    character(*)  , intent(in) :: name
    type(NcGroup)              :: getGroupByName
    integer(i4)                :: id

    call check(nf90_inq_ncid(self%id, name, id), &
         "Could not inquire variable: " // name)
    getGroupByName = NcGroup(id)
  end function getGroupByName

  function getVariableByName(self, name)
    class(NcGroup), intent(in) :: self
    character(*)    , intent(in) :: name
    type(NcVariable)             :: getVariableByName
    integer(i4)                  :: id

    call check(nf90_inq_varid(self%id, name, id), &
         "Could not inquire variable: " // name)
    getVariableByName = NcVariable(id, self)

  end function getVariableByName

  function getVariableName(self)
    class(NcVariable), intent(in) :: self
    character(len=256)            :: getVariableName

    call check(nf90_inquire_variable(self%parent%id, self%id, name=getVariableName), &
         "Could not inquire variable name")
  end function getVariableName

  function getNoDimensions(self)
    class(NcVariable), intent(in) :: self
    integer(i4)                   :: getNoDimensions

    call check(nf90_inquire_variable(self%parent%id,self%id,ndims=getNoDimensions), &
         "Could not inquire variable: " // self%getName())
  end function getNoDimensions

  function getVariableDimensions(self)
    class(NcVariable), intent(in)  :: self
    type(NcDimension), allocatable :: getVariableDimensions(:)
    integer(i4)      , allocatable :: dimids(:)
    integer(i4)                    :: ii , ndims

    ndims = self%getNoDimensions()
    allocate(dimids(ndims), getVariableDimensions(ndims))
    call check(nf90_inquire_variable(self%parent%id,self%id,dimids=dimids), &
         "Could not inquire variable: " // self%getName())

    do ii=1,ndims
       getVariableDimensions (ii) = self%parent%getDimension(dimids(ii))
    end do
  end function getVariableDimensions

  function getVariableShape(self)
    class(NcVariable), intent(in)  :: self
    integer(i4)     , allocatable :: getVariableShape(:)
    type(NcDimension), allocatable :: dims(:)
    integer(i4)                   :: ii, ndims

    ndims = self%getNoDimensions()
    allocate(getVariableShape(ndims), dims(ndims))

    dims = self%getDimensions()
    do ii = 1,size(dims)
       getVariableShape(ii) = dims(ii)%getLength()
    end do
  end function getVariableShape

  function getVariableRank(self)
    class(NcVariable), intent(in)  :: self
    integer(i4)                   :: getVariableRank

    getVariableRank = size(self%getDimensions())
  end function getVariableRank

  function getVariableDtype(self)
    class(NcVariable), intent(in) :: self
    integer(i4)                   :: dtype
    character(3)                  :: getVariableDtype

    call check(nf90_inquire_variable(self%parent%id,self%id,xtype=dtype),&
         "Could not inquire variable: " // self%getName())
    getVariableDtype = getDtypeFromInteger(dtype)
  end function getVariableDtype

  function isUnlimitedVariable(self)
    class(NcVariable), intent(in)  :: self
    logical                        :: isUnlimitedVariable
    type(NcDimension), allocatable :: dims(:)
    type(NcDimension)              :: dim
    integer(i4)                    :: ii

    allocate(dims(self%getNoDimensions()))

    isUnlimitedVariable = .false.
    dims = self%getDimensions()

    do ii = 1,size(dims)
       dim = dims(ii)
       if (dim%isUnlimited()) then
          isUnlimitedVariable = .true.
       end if
    end do
  end function isUnlimitedVariable

  logical function hasAttribute(self,name)
    class(NcAttributable), intent(inout) :: self
    character(*)     , intent(in) :: name
    integer(i4)                   :: status

    select type (self)
    class is (NcGroup)
        status = nf90_inquire_attribute(self%id, NF90_GLOBAL, name)
    class is (NcVariable)
        status = nf90_inquire_attribute(self%parent%id, self%id, name)
    end select

    hasAttribute = (status .eq. NF90_NOERR)
  end function hasAttribute

  subroutine setAttributeChar(self, name, data)
    class(NcAttributable), intent(in) :: self
    character(*)         , intent(in) :: name
    character(*)         , intent(in) :: data
    integer(i4)                      :: ids(2)

    ids = self%getAttributableIds()
    call check(nf90_put_att(ids(1), ids(2), name, data), &
         "Failed to write attribute: " // name)

  end subroutine setAttributeChar

  subroutine setAttributeI8(self, name, data)
    class(NcAttributable), intent(in) :: self
    character(*)         , intent(in) :: name
    integer(i1)          , intent(in) :: data
    integer(i4)                      :: ids(2)

    ids = self%getAttributableIds()
    call check(nf90_put_att(ids(1), ids(2), name, data), &
         "Failed to write attribute: " // name)

  end subroutine setAttributeI8

  subroutine setAttributeI16(self, name, data)
    class(NcAttributable), intent(in) :: self
    character(*)         , intent(in) :: name
    integer(i2)         , intent(in) :: data
    integer(i4)                      :: ids(2)

    ids = self%getAttributableIds()
    call check(nf90_put_att(ids(1), ids(2), name, data), &
         "Failed to write attribute: " // name)

  end subroutine setAttributeI16

  subroutine setAttributeI32(self, name, data)
    class(NcAttributable), intent(in) :: self
    character(*)         , intent(in) :: name
    integer(i4)         , intent(in) :: data
    integer(i4)                      :: ids(2)

    ids = self%getAttributableIds()
    call check(nf90_put_att(ids(1), ids(2), name, data), &
         "Failed to write attribute: " // name)

  end subroutine setAttributeI32

  subroutine setAttributeI64(self, name, data)
    class(NcAttributable), intent(in) :: self
    character(*)         , intent(in) :: name
    integer(i8)         , intent(in) :: data
    integer(i4)                      :: ids(2)

    ids = self%getAttributableIds()
    call check(nf90_put_att(ids(1), ids(2), name, data), &
         "Failed to write attribute: " // name)

  end subroutine setAttributeI64

  subroutine setAttributeF32(self, name, data)
    class(NcAttributable), intent(in) :: self
    character(*)         , intent(in) :: name
    real(sp)            , intent(in) :: data
    integer(i4)                      :: ids(2)

    ids = self%getAttributableIds()
    call check(nf90_put_att(ids(1), ids(2), name, data), &
         "Failed to write attribute: " // name)

  end subroutine setAttributeF32

  subroutine setAttributeF64(self, name, data)
    class(NcAttributable), intent(in) :: self
    character(*)         , intent(in) :: name
    real(dp)            , intent(in) :: data
    integer(i4)                      :: ids(2)

    ids = self%getAttributableIds()
    call check(nf90_put_att(ids(1), ids(2), name, data), &
         "Failed to write attribute: " // name)

  end subroutine setAttributeF64

  subroutine getAttributeChar(self, name, avalue)
    class(NcAttributable), intent(in)  :: self
    character(*)         , intent(in)  :: name
    character(*)         , intent(out) :: avalue
    integer(i4)                       :: length, ids(2)

    ids = self%getAttributableIds()
    call check(nf90_inquire_attribute(ids(1), ids(2), name, len=length),&
         "Could not inquire attribute " // name)
    call check(nf90_get_att(ids(1), ids(2), name, avalue), &
         "Could not read attribute "//name)
  end subroutine getAttributeChar

  subroutine getAttributeI8(self, name, avalue)
    class(NcAttributable) , intent(in)  :: self
    character(*)          , intent(in)  :: name
    integer(i1)           , intent(out) :: avalue
    integer(i4)                        :: length, ids(2)

    ids = self%getAttributableIds()
    call check(nf90_inquire_attribute(ids(1), ids(2), name, len=length),&
         "Could not inquire attribute " // name)
    call check(nf90_get_att(ids(1), ids(2), name, avalue), &
         "Could not read attribute "//name)

  end subroutine getAttributeI8

  subroutine getAttributeI16(self, name, avalue)
    class(NcAttributable) , intent(in)  :: self
    character(*)          , intent(in)  :: name
    integer(i2)          , intent(out) :: avalue
    integer(i4)                        :: length, ids(2)

    ids = self%getAttributableIds()
    call check(nf90_inquire_attribute(ids(1), ids(2), name, len=length),&
         "Could not inquire attribute " // name)
    call check(nf90_get_att(ids(1), ids(2), name, avalue), &
         "Could not read attribute "//name)

  end subroutine getAttributeI16

  subroutine getAttributeI32(self, name, avalue)
    class(NcAttributable) , intent(in)  :: self
    character(*)          , intent(in)  :: name
    integer(i4)          , intent(out) :: avalue
    integer(i4)                        :: length, ids(2)

    ids = self%getAttributableIds()
    call check(nf90_inquire_attribute(ids(1), ids(2), name, len=length),&
         "Could not inquire attribute " // name)
    call check(nf90_get_att(ids(1), ids(2), name, avalue), &
         "Could not read attribute "//name)

  end subroutine getAttributeI32

  subroutine getAttributeI64(self, name, avalue)
    class(NcAttributable) , intent(in)  :: self
    character(*)          , intent(in)  :: name
    integer(i8)          , intent(out) :: avalue
    integer(i4)                        :: length, ids(2)

    ids = self%getAttributableIds()
    call check(nf90_inquire_attribute(ids(1), ids(2), name, len=length),&
         "Could not inquire attribute " // name)
    call check(nf90_get_att(ids(1), ids(2), name, avalue), &
         "Could not read attribute "//name)

  end subroutine getAttributeI64

  subroutine getAttributeF32(self, name, avalue)
    class(NcAttributable) , intent(in)  :: self
    character(*)          , intent(in)  :: name
    real(sp)             , intent(out) :: avalue
    integer(i4)                        :: length, ids(2)

    ids = self%getAttributableIds()
    call check(nf90_inquire_attribute(ids(1), ids(2), name, len=length),&
         "Could not inquire attribute " // name)
    call check(nf90_get_att(ids(1), ids(2), name, avalue), &
         "Could not read attribute "//name)

  end subroutine getAttributeF32

  subroutine getAttributeF64(self, name, avalue)
    class(NcAttributable) , intent(in)  :: self
    character(*)          , intent(in)  :: name
    real(dp)             , intent(out) :: avalue
    integer(i4)                        :: length, ids(2)

    ids = self%getAttributableIds()
    call check(nf90_inquire_attribute(ids(1), ids(2), name, len=length),&
         "Could not inquire attribute " // name)
    call check(nf90_get_att(ids(1), ids(2), name, avalue), &
         "Could not read attribute "//name)

 end subroutine getAttributeF64

  function getAttributableIds(self)
    class(NcAttributable), intent(in) :: self
    integer(i4)                      :: getAttributableIds(2)
    select type(self)
    class is (NcGroup)
       getAttributableIds(1) = self%id
       getAttributableIds(2) = NF90_GLOBAL
    class is (NcVariable)
       getAttributableIds(1) = self%parent%id
       getAttributableIds(2) = self%id
    end select
  end function getAttributableIds

  subroutine renameAttribute(self, oldname, newname)
    class(NcAttributable), intent(inout) :: self
    character(len=*), intent(in)         :: oldname, newname
    integer(i4)                         :: ids(2)
    ids = self%getAttributableIds()
    call check(nf90_rename_att(ids(1), ids(2), oldname, newname), "Failed to rename attribute: " // oldname)
  end subroutine renameAttribute

  subroutine renameVariable(self, name)
    class(NcVariable), intent(inout) :: self
    character(len=*),  intent(in)    :: name
    call check(nf90_rename_var(self%parent%id, self%id, name), "Failed to rename variable: " // self%getName())
  end subroutine renameVariable

  subroutine renameDimension(self, name)
    class(NcDimension), intent(inout) :: self
    character(len=*),  intent(in)     :: name
    call check(nf90_rename_dim(self%parent%id, self%id, name), "Failed to rename dimension: " // self%getName())
  end subroutine renameDimension

  subroutine setVariableFillValueI8(self, fvalue)
    class(NcVariable), intent(inout)  :: self
    integer(i1)      , intent(in)  :: fvalue

    if (.not. self%hasAttribute("_FillValue")) then
       call self%setAttribute("_FillValue",fvalue)
    end if
  end subroutine setVariableFillValueI8

  subroutine setVariableFillValueI16(self, fvalue)
    class(NcVariable), intent(inout)  :: self
    integer(i2)      , intent(in)  :: fvalue

    if (.not. self%hasAttribute("_FillValue")) then
       call self%setAttribute("_FillValue",fvalue)
    end if
  end subroutine setVariableFillValueI16

  subroutine setVariableFillValueI32(self, fvalue)
    class(NcVariable), intent(inout)  :: self
    integer(i4)      , intent(in)  :: fvalue

    if (.not. self%hasAttribute("_FillValue")) then
       call self%setAttribute("_FillValue",fvalue)
    end if
  end subroutine setVariableFillValueI32

  subroutine setVariableFillValueI64(self, fvalue)
    class(NcVariable), intent(inout)  :: self
    integer(i8)      , intent(in)  :: fvalue

    if (.not. self%hasAttribute("_FillValue")) then
       call self%setAttribute("_FillValue",fvalue)
    end if
  end subroutine setVariableFillValueI64

  subroutine setVariableFillValueF32(self, fvalue)
    class(NcVariable), intent(inout)  :: self
    real(sp)         , intent(in)  :: fvalue

    if (.not. self%hasAttribute("_FillValue")) then
       call self%setAttribute("_FillValue",fvalue)
    end if
  end subroutine setVariableFillValueF32

  subroutine setVariableFillValueF64(self, fvalue)
    class(NcVariable), intent(inout)  :: self
    real(dp)         , intent(in)  :: fvalue

    if (.not. self%hasAttribute("_FillValue")) then
       call self%setAttribute("_FillValue",fvalue)
    end if
  end subroutine setVariableFillValueF64

  subroutine getVariableFillValueI8(self, fvalue)
    class(NcVariable), intent(inout)  :: self
    integer(i1)      , intent(out) :: fvalue

    if (self%hasAttribute("_FillValue")) then
       call self%getAttribute("_FillValue", fvalue)
    else
       fvalue = NF90_FILL_BYTE
    end if
  end subroutine getVariableFillValueI8

  subroutine getVariableFillValueI16(self, fvalue)
    class(NcVariable), intent(inout)  :: self
    integer(i2)      , intent(out) :: fvalue

    if (self%hasAttribute("_FillValue")) then
       call self%getAttribute("_FillValue", fvalue)
    else
       fvalue = NF90_FILL_SHORT
    end if
  end subroutine getVariableFillValueI16

  subroutine getVariableFillValueI32(self, fvalue)
    class(NcVariable), intent(inout)  :: self
    integer(i4)      , intent(out) :: fvalue

    if (self%hasAttribute("_FillValue")) then
       call self%getAttribute("_FillValue", fvalue)
    else
       fvalue = NF90_FILL_INT
    end if
  end subroutine getVariableFillValueI32

  subroutine getVariableFillValueI64(self, fvalue)
    class(NcVariable), intent(inout)  :: self
    integer(i8)      , intent(out) :: fvalue

    if (self%hasAttribute("_FillValue")) then
       call self%getAttribute("_FillValue", fvalue)
    else
       fvalue = NF90_FILL_INT
    end if
  end subroutine getVariableFillValueI64

  subroutine getVariableFillValueF32(self, fvalue)
    class(NcVariable), intent(inout)  :: self
    real(sp)         , intent(out) :: fvalue

    if (self%hasAttribute("_FillValue")) then
       call self%getAttribute("_FillValue", fvalue)
    else
       fvalue = NF90_FILL_FLOAT
    end if
  end subroutine getVariableFillValueF32

  subroutine getVariableFillValueF64(self, fvalue)
    class(NcVariable), intent(inout)  :: self
    real(dp)         , intent(out) :: fvalue

    if (self%hasAttribute("_FillValue")) then
       call self%getAttribute("_FillValue", fvalue)
    else
       fvalue = NF90_FILL_DOUBLE
    end if
  end subroutine getVariableFillValueF64

  subroutine setDataScalarI8(self, values, start)
    class(NcVariable), intent(in)           :: self
    integer(i1)      , intent(in)           :: values
    integer(i4)      , intent(in), optional :: start(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setDataScalarI8

  subroutine setData1dI8(self, values, start, cnt, stride, map)
    class(NcVariable), intent(in)           :: self
    integer(i1)      , intent(in)           :: values(:)
    integer(i4)      , intent(in), optional :: start(:), cnt(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setData1dI8

  subroutine setData2dI8(self, values, start, cnt, stride, map)
    class(NcVariable), intent(in)           :: self
    integer(i1)      , intent(in)           :: values(:,:)
    integer(i4)      , intent(in), optional :: start(:), cnt(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setData2dI8

  subroutine setData3dI8(self, values, start, cnt, stride, map)
    class(NcVariable), intent(in)           :: self
    integer(i1)      , intent(in)           :: values(:,:,:)
    integer(i4)      , intent(in), optional :: start(:), cnt(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setData3dI8

  subroutine setData4dI8(self, values, start, cnt, stride, map)
    class(NcVariable), intent(in)           :: self
    integer(i1)      , intent(in)           :: values(:,:,:,:)
    integer(i4)      , intent(in), optional :: start(:), cnt(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setData4dI8

  subroutine setData5dI8(self, values, start, cnt, stride, map)
    class(NcVariable), intent(in)           :: self
    integer(i1)      , intent(in)           :: values(:,:,:,:,:)
    integer(i4)      , intent(in), optional :: start(:), cnt(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setData5dI8

  subroutine setDataScalarI16(self, values, start)
    class(NcVariable), intent(in)           :: self
    integer(i2)      , intent(in)           :: values
    integer(i4)      , intent(in), optional :: start(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setDataScalarI16

  subroutine setData1dI16(self, values, start, cnt, stride, map)
    class(NcVariable), intent(in)           :: self
    integer(i2)      , intent(in)           :: values(:)
    integer(i4)      , intent(in), optional :: start(:), cnt(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setData1dI16

  subroutine setData2dI16(self, values, start, cnt, stride, map)
    class(NcVariable), intent(in)           :: self
    integer(i2)      , intent(in)           :: values(:,:)
    integer(i4)      , intent(in), optional :: start(:), cnt(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setData2dI16

  subroutine setData3dI16(self, values, start, cnt, stride, map)
    class(NcVariable), intent(in)           :: self
    integer(i2)      , intent(in)           :: values(:,:,:)
    integer(i4)      , intent(in), optional :: start(:), cnt(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setData3dI16

  subroutine setData4dI16(self, values, start, cnt, stride, map)
    class(NcVariable), intent(in)           :: self
    integer(i2)      , intent(in)           :: values(:,:,:,:)
    integer(i4)      , intent(in), optional :: start(:), cnt(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setData4dI16

  subroutine setData5dI16(self, values, start, cnt, stride, map)
    class(NcVariable), intent(in)           :: self
    integer(i2)      , intent(in)           :: values(:,:,:,:,:)
    integer(i4)      , intent(in), optional :: start(:), cnt(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setData5dI16

  subroutine setDataScalarI32(self, values, start)
    class(NcVariable), intent(in)           :: self
    integer(i4)      , intent(in)           :: values
    integer(i4)      , intent(in), optional :: start(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setDataScalarI32

  subroutine setData1dI32(self, values, start, cnt, stride, map)
    class(NcVariable), intent(in)           :: self
    integer(i4)      , intent(in)           :: values(:)
    integer(i4)      , intent(in), optional :: start(:), cnt(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setData1dI32

  subroutine setData2dI32(self, values, start, cnt, stride, map)
    class(NcVariable), intent(in)           :: self
    integer(i4)      , intent(in)           :: values(:,:)
    integer(i4)      , intent(in), optional :: start(:), cnt(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setData2dI32

  subroutine setData3dI32(self, values, start, cnt, stride, map)
    class(NcVariable), intent(in)           :: self
    integer(i4)      , intent(in)           :: values(:,:,:)
    integer(i4)      , intent(in), optional :: start(:), cnt(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setData3dI32

  subroutine setData4dI32(self, values, start, cnt, stride, map)
    class(NcVariable), intent(in)           :: self
    integer(i4)      , intent(in)           :: values(:,:,:,:)
    integer(i4)      , intent(in), optional :: start(:), cnt(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setData4dI32

  subroutine setData5dI32(self, values, start, cnt, stride, map)
    class(NcVariable), intent(in)           :: self
    integer(i4)      , intent(in)           :: values(:,:,:,:,:)
    integer(i4)      , intent(in), optional :: start(:), cnt(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setData5dI32

  subroutine setDataScalarI64(self, values, start)
    class(NcVariable), intent(in)           :: self
    integer(i8)      , intent(in)           :: values
    integer(i4)      , intent(in), optional :: start(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setDataScalarI64

  subroutine setData1dI64(self, values, start, cnt, stride, map)
    class(NcVariable), intent(in)           :: self
    integer(i8)      , intent(in)           :: values(:)
    integer(i4)      , intent(in), optional :: start(:), cnt(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setData1dI64

  subroutine setData2dI64(self, values, start, cnt, stride, map)
    class(NcVariable), intent(in)           :: self
    integer(i8)      , intent(in)           :: values(:,:)
    integer(i4)      , intent(in), optional :: start(:), cnt(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setData2dI64

  subroutine setData3dI64(self, values, start, cnt, stride, map)
    class(NcVariable), intent(in)           :: self
    integer(i8)      , intent(in)           :: values(:,:,:)
    integer(i4)      , intent(in), optional :: start(:), cnt(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setData3dI64

  subroutine setData4dI64(self, values, start, cnt, stride, map)
    class(NcVariable), intent(in)           :: self
    integer(i8)      , intent(in)           :: values(:,:,:,:)
    integer(i4)      , intent(in), optional :: start(:), cnt(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setData4dI64

  subroutine setData5dI64(self, values, start, cnt, stride, map)
    class(NcVariable), intent(in)           :: self
    integer(i8)      , intent(in)           :: values(:,:,:,:,:)
    integer(i4)      , intent(in), optional :: start(:), cnt(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setData5dI64

  subroutine setDataScalarF32(self, values, start)
    class(NcVariable), intent(in)           :: self
    real(sp)         , intent(in)           :: values
    integer(i4)      , intent(in), optional :: start(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setDataScalarF32

  subroutine setData1dF32(self, values, start, cnt, stride, map)
    class(NcVariable), intent(in)           :: self
    real(sp)         , intent(in)           :: values(:)
    integer(i4)      , intent(in), optional :: start(:), cnt(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setData1dF32

  subroutine setData2dF32(self, values, start, cnt, stride, map)
    class(NcVariable), intent(in)           :: self
    real(sp)         , intent(in)           :: values(:,:)
    integer(i4)      , intent(in), optional :: start(:), cnt(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setData2dF32

  subroutine setData3dF32(self, values, start, cnt, stride, map)
    class(NcVariable), intent(in)           :: self
    real(sp)         , intent(in)           :: values(:,:,:)
    integer(i4)      , intent(in), optional :: start(:), cnt(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setData3dF32

  subroutine setData4dF32(self, values, start, cnt, stride, map)
    class(NcVariable), intent(in)           :: self
    real(sp)         , intent(in)           :: values(:,:,:,:)
    integer(i4)      , intent(in), optional :: start(:), cnt(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setData4dF32

  subroutine setData5dF32(self, values, start, cnt, stride, map)
    class(NcVariable), intent(in)           :: self
    real(sp)         , intent(in)           :: values(:,:,:,:,:)
    integer(i4)      , intent(in), optional :: start(:), cnt(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setData5dF32

  subroutine setDataScalarF64(self, values, start)
    class(NcVariable), intent(in)           :: self
    real(dp)         , intent(in)           :: values
    integer(i4)      , intent(in), optional :: start(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setDataScalarF64

  subroutine setData1dF64(self, values, start, cnt, stride, map)
    class(NcVariable), intent(in)           :: self
    real(dp)         , intent(in)           :: values(:)
    integer(i4)      , intent(in), optional :: start(:), cnt(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setData1dF64

  subroutine setData2dF64(self, values, start, cnt, stride, map)
    class(NcVariable), intent(in)           :: self
    real(dp)         , intent(in)           :: values(:,:)
    integer(i4)      , intent(in), optional :: start(:), cnt(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setData2dF64

  subroutine setData3dF64(self, values, start, cnt, stride, map)
    class(NcVariable), intent(in)           :: self
    real(dp)         , intent(in)           :: values(:,:,:)
    integer(i4)      , intent(in), optional :: start(:), cnt(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setData3dF64

  subroutine setData4dF64(self, values, start, cnt, stride, map)
    class(NcVariable), intent(in)           :: self
    real(dp)         , intent(in)           :: values(:,:,:,:)
    integer(i4)      , intent(in), optional :: start(:), cnt(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setData4dF64

  subroutine setData5dF64(self, values, start, cnt, stride, map)
    class(NcVariable), intent(in)           :: self
    real(dp)         , intent(in)           :: values(:,:,:,:,:)
    integer(i4)      , intent(in), optional :: start(:), cnt(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
         "Failed to write data into variable: " // trim(self%getName()))
  end subroutine setData5dF64

  subroutine getDataScalarI8(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)            :: self
    integer(i4)     , intent(in) , optional :: start(:), cnt(:), stride(:), map(:)
    integer(i1)      , intent(out)           :: data
    integer(i1)                              :: tmp(1)

    call check (nf90_get_var(self%parent%id, self%id, tmp, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
    data = tmp(1)
  end subroutine getDataScalarI8

  subroutine getData1dI8(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)     , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i1)      , intent(out), allocatable :: data(:)
    integer(i4)                  , allocatable :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, size(shape(data)))

    allocate(data(datashape(1)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData1dI8

  subroutine getData2dI8(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)     , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i1)      , intent(out), allocatable :: data(:,:)
    integer(i4)                  , allocatable :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, size(shape(data)))

    allocate(data(datashape(1), datashape(2)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData2dI8

  subroutine getData3dI8(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)     , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i1)      , intent(out), allocatable :: data(:,:,:)
    integer(i4)                  , allocatable :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, size(shape(data)))

    allocate(data(datashape(1), datashape(2), datashape(3)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData3dI8

  subroutine getData4dI8(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)     , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i1)      , intent(out), allocatable :: data(:,:,:,:)
    integer(i4)                  , allocatable :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, size(shape(data)))

    allocate(data(datashape(1), datashape(2), datashape(3), datashape(4)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData4dI8

  subroutine getData5dI8(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)     , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i1)      , intent(out), allocatable :: data(:,:,:,:,:)
    integer(i4)                  , allocatable :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, size(shape(data)))

    allocate(data(datashape(1), datashape(2), datashape(3), datashape(4), datashape(5)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData5dI8

  subroutine getDataScalarI16(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i2)      , intent(out)              :: data
    integer(i2)                                 :: tmp(1)

    call check (nf90_get_var(self%parent%id, self%id, tmp, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
    data = tmp(1)
  end subroutine getDataScalarI16

  subroutine getData1dI16(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)     , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i2)     , intent(out), allocatable :: data(:)
    integer(i4)                  , allocatable :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, size(shape(data)))

    allocate(data(datashape(1)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData1dI16

  subroutine getData2dI16(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)     , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i2)     , intent(out), allocatable :: data(:,:)
    integer(i4)                  , allocatable :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, size(shape(data)))

    allocate(data(datashape(1), datashape(2)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData2dI16

  subroutine getData3dI16(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)     , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i2)     , intent(out), allocatable :: data(:,:,:)
    integer(i4)                  , allocatable :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, size(shape(data)))

    allocate(data(datashape(1), datashape(2), datashape(3)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData3dI16

  subroutine getData4dI16(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)     , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i2)     , intent(out), allocatable :: data(:,:,:,:)
    integer(i4)                  , allocatable :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape( start, cnt, stride)
    datashape = getReadShape(slcshape, size(shape(data)))

    allocate(data(datashape(1), datashape(2), datashape(3), datashape(4)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData4dI16

  subroutine getData5dI16(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)     , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i2)     , intent(out), allocatable :: data(:,:,:,:,:)
    integer(i4)                  , allocatable :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, size(shape(data)))

    allocate(data(datashape(1), datashape(2), datashape(3), datashape(4), datashape(5)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData5dI16

  subroutine getDataScalarI32(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)     , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i4)     , intent(out)              :: data
    integer(i4)                                :: tmp(1)

    call check (nf90_get_var(self%parent%id, self%id, tmp, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
    data = tmp(1)
  end subroutine getDataScalarI32

  subroutine getData1dI32(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)     , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i4)     , intent(out), allocatable :: data(:)
    integer(i4)                  , allocatable :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, size(shape(data)))

    allocate(data(datashape(1)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData1dI32

  subroutine getData2dI32(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)     , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i4)     , intent(out), allocatable :: data(:,:)
    integer(i4)                  , allocatable :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, size(shape(data)))

    allocate(data(datashape(1), datashape(2)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData2dI32

  subroutine getData3dI32(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)     , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i4)     , intent(out), allocatable :: data(:,:,:)
    integer(i4)                  , allocatable :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, size(shape(data)))

    allocate(data(datashape(1), datashape(2), datashape(3)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData3dI32

  subroutine getData4dI32(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)     , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i4)     , intent(out), allocatable :: data(:,:,:,:)
    integer(i4)                  , allocatable :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, size(shape(data)))

    allocate(data(datashape(1), datashape(2), datashape(3), datashape(4)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData4dI32

  subroutine getData5dI32(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)     , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i4)     , intent(out), allocatable :: data(:,:,:,:,:)
    integer(i4)                  , allocatable :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, size(shape(data)))

    allocate(data(datashape(1), datashape(2), datashape(3), datashape(4), datashape(5)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData5dI32

  subroutine getDataScalarI64(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)     , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i8)     , intent(out)              :: data
    integer(i4)                                :: tmp(1)

    call check (nf90_get_var(self%parent%id, self%id, tmp, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
    data = tmp(1)
  end subroutine getDataScalarI64

  subroutine getData1dI64(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)     , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i8)     , intent(out), allocatable :: data(:)
    integer(i4)                  , allocatable :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, size(shape(data)))

    allocate(data(datashape(1)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData1dI64

  subroutine getData2dI64(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)     , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i8)     , intent(out), allocatable :: data(:,:)
    integer(i4)                  , allocatable :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, size(shape(data)))

    allocate(data(datashape(1), datashape(2)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData2dI64

  subroutine getData3dI64(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)     , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i8)     , intent(out), allocatable :: data(:,:,:)
    integer(i4)                  , allocatable :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, size(shape(data)))

    allocate(data(datashape(1), datashape(2), datashape(3)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData3dI64

  subroutine getData4dI64(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)     , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i8)     , intent(out), allocatable :: data(:,:,:,:)
    integer(i4)                  , allocatable :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, size(shape(data)))

    allocate(data(datashape(1), datashape(2), datashape(3), datashape(4)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData4dI64

  subroutine getData5dI64(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)     , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i8)     , intent(out), allocatable :: data(:,:,:,:,:)
    integer(i4)                  , allocatable :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, size(shape(data)))

    allocate(data(datashape(1), datashape(2), datashape(3), datashape(4), datashape(5)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData5dI64

  subroutine getDataScalarF32(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)             :: self
    integer(i4)     , intent(in) , optional  :: start(:), cnt(:), stride(:), map(:)
    real(sp)        , intent(out)            :: data
    real(sp)                                 :: tmp(1)

    call check (nf90_get_var(self%parent%id, self%id, tmp, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
    data = tmp(1)
  end subroutine getDataScalarF32

  subroutine getData1dF32(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)     , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    real(sp)        , intent(out), allocatable :: data(:)
    integer(i4)                  , allocatable :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, size(shape(data)))

    allocate(data(datashape(1)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData1dF32

  subroutine getData2dF32(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)     , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    real(sp)        , intent(out), allocatable :: data(:,:)
    integer(i4)                  , allocatable :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, size(shape(data)))

    allocate(data(datashape(1), datashape(2)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData2dF32

  subroutine getData3dF32(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)     , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    real(sp)        , intent(out), allocatable :: data(:,:,:)
    integer(i4)                  , allocatable :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, size(shape(data)))

    allocate(data(datashape(1), datashape(2), datashape(3)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData3dF32

  subroutine getData4dF32(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)     , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    real(sp)        , intent(out), allocatable :: data(:,:,:,:)
    integer(i4)                  , allocatable :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, size(shape(data)))

    allocate(data(datashape(1), datashape(2), datashape(3), datashape(4)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData4dF32

  subroutine getData5dF32(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)     , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    real(sp)        , intent(out), allocatable :: data(:,:,:,:,:)
    integer(i4)                  , allocatable :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, size(shape(data)))

    allocate(data(datashape(1), datashape(2), datashape(3), datashape(4), datashape(5)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData5dF32

  subroutine getDataScalarF64(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)             :: self
    integer(i4)     , intent(in) , optional  :: start(:), cnt(:), stride(:), map(:)
    real(dp)        , intent(out)            :: data
    real(dp)                                 :: tmp(1)

    call check (nf90_get_var(self%parent%id, self%id, tmp, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
    data = tmp(1)
  end subroutine getDataScalarF64

  subroutine getData1dF64(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)     , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    real(dp)        , intent(out), allocatable :: data(:)
    integer(i4)                  , allocatable :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, size(shape(data)))

    allocate(data(datashape(1)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData1dF64

  subroutine getData2dF64(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)     , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    real(dp)        , intent(out), allocatable :: data(:,:)
    integer(i4)                  , allocatable :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, size(shape(data)))

    allocate(data(datashape(1), datashape(2)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData2dF64

  subroutine getData3dF64(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)     , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    real(dp)        , intent(out), allocatable :: data(:,:,:)
    integer(i4)                  , allocatable :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, size(shape(data)))

    allocate(data(datashape(1), datashape(2), datashape(3)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData3dF64

  subroutine getData4dF64(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)     , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    real(dp)        , intent(out), allocatable :: data(:,:,:,:)
    integer(i4)                  , allocatable :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, size(shape(data)))

    allocate(data(datashape(1), datashape(2), datashape(3), datashape(4)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData4dF64

  subroutine getData5dF64(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)     , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    real(dp)        , intent(out), allocatable :: data(:,:,:,:,:)
    integer(i4)                  , allocatable :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, size(shape(data)))

    allocate(data(datashape(1), datashape(2), datashape(3), datashape(4), datashape(5)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData5dF64

  function getSlicingShape(self, instart, incnt, instride) result(out)
    class(NcVariable), intent(in)           :: self
    integer(i4)     , intent(in), optional :: instart(:), incnt(:), instride(:)
    integer(i4)     , allocatable          :: out(:)

    out = self%getShape()

    if (present(incnt)) then
       out(:size(incnt)) = incnt
       ! out = incnt
    else
       if (present(instart)) then
          out(:size(instart)) = out(:size(instart)) - (instart - 1)
       end if
       if (present(instride)) then
          out(:size(instride)) = out(:size(instride)) / instride
       end if
    end if

  end function getSlicingShape

  function getReadShape(slcshape, outrank) result(out)
    integer(i4), intent(in)   :: slcshape(:)
    integer(i4), intent(in)   :: outrank
    integer(i4)               :: naxis
    integer(i4), allocatable  :: out(:)

    naxis = count(slcshape .gt. 1)

    if (all(slcshape .eq. 1)) then
       ! return 1-element array
      allocate(out(size(slcshape)))
       out(:) = 1
    else if (size(slcshape) .eq. outrank) then
       ! sizes fit
       out = slcshape
    else if (naxis .eq. outrank) then
       out = pack(slcshape, slcshape .gt. 1)
    ! else if (naxis .lt. outrank) then
       ! would be nice...
    else
       write(*,*) "Given indices do not match output variable rank!"
       stop 1
    end if
  end function getReadShape

  function getDtypeFromString(dtype)
    integer(i4)          :: getDtypeFromString
    character(*)         :: dtype

    select case(dtype)
    case("f32")
       getDtypeFromString = NF90_FLOAT
    case("f64")
       getDtypeFromString = NF90_DOUBLE
    case("i8")
       getDtypeFromString = NF90_BYTE
    case("i16")
       getDtypeFromString = NF90_SHORT
    case("i32")
       getDtypeFromString = NF90_INT
    case("i64")
       getDtypeFromString = NF90_INT64
    case default
       write(*,*) "Datatype not understood: ", dtype
       stop 1
    end select
  end function getDtypeFromString

  function getDtypeFromInteger(dtype)
    character(3) :: getDtypeFromInteger
    integer(i4)  :: dtype

    select case(dtype)
    case(NF90_FLOAT)
       getDtypeFromInteger = "f32"
    case(NF90_DOUBLE)
       getDtypeFromInteger = "f64"
    case(NF90_BYTE)
       getDtypeFromInteger = "i8"
    case(NF90_SHORT)
       getDtypeFromInteger = "i16"
    case(NF90_INT)
       getDtypeFromInteger = "i32"
    case(NF90_INT64)
       getDtypeFromInteger = "i64"
    case default
       write(*,*) "Datatype not understood: ", dtype
       stop 1
    end select
  end function getDtypeFromInteger

  function getCreationMode(cmode)
    character(*), intent(in), optional :: cmode
    integer(i4)                       :: getCreationMode
    character(256)                     :: mode

    if (.not. (present(cmode))) then
       mode = "NETCDF4"
    else
       mode = cmode
    end if

    select case(trim(mode))
    case ("NETCDF4")
       getCreationMode = NF90_NETCDF4
    case ("SHARE")
       getCreationMode = NF90_SHARE
    case ("CLASSIC")
       getCreationMode = NF90_CLASSIC_MODEL
    case ("HDF5")
       getCreationMode = NF90_HDF5
    case ("64BIT_OFFSET")
       getCreationMode = NF90_64BIT_OFFSET
    case default
       print*, "Creation mode not understood: " // trim(mode)
       stop 1
    end select

  end function getCreationMode

  subroutine check(status, msg)
    integer(i4) , intent(in) :: status
    character(*), intent(in) :: msg

    if (status .ne. NF90_NOERR) then
       write(*,*) msg
       write(*,*) nf90_strerror(status)
       stop 1
    end if
  end subroutine check

end module mo_netcdf

