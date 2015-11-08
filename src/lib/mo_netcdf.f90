!> \file mo_netcdf.f90

!> \brief NetCDF Fortran 90 interface wrapper

!> \details A thin wrapper around the NetCDF Fortran 90 interface.
!>          Provided are currently 3 user facing derived Types:
!>             1. NcDataset
!>             2. NcDimension
!>             3. NcVariable
!
!> \authors David Schaefer
!> \date Jun 2015


module mo_netcdf

  ! This module provides a thin wrapper around the NetCDF Fortran 90 interface,
  ! following a somehow object-oriented approach.

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

  ! Copyright 2011-2015 David Schaefer

  use mo_kind,         only: i2, i4, sp, dp
  use netcdf,          only: &
       nf90_open, nf90_close, nf90_strerror, nf90_def_dim, nf90_def_var,   &
       nf90_put_var, nf90_get_var, nf90_put_att, nf90_get_att,             &
       nf90_inquire, nf90_inq_dimid, nf90_inquire_dimension,               &
       nf90_inq_varid, nf90_inquire_variable, nf90_inquire_attribute,      &
       NF90_OPEN, NF90_NETCDF4, NF90_CREATE, NF90_WRITE, NF90_NOWRITE,     &
       NF90_SHORT, NF90_INT, NF90_FLOAT, NF90_DOUBLE,                      &
       NF90_FILL_SHORT, NF90_FILL_INT, NF90_FILL_FLOAT , NF90_FILL_DOUBLE, &
       NF90_NOERR, NF90_UNLIMITED, NF90_GLOBAL


  implicit none

  ! --------------------------------------------------------------------------------------
  !
  !     NAME
  !         NcDataset
  !
  !     PURPOSE
  !>        \brief Provides basic file modification functionality
  !>        \details Bound to this derived type is the basic file level create/retrieve
  !>                 functionality, i.e. functions/subroutines to create/retrieve
  !>                 dimensions, variables and global attributes.
  !>                 All files created by this derived type and its procedures are
  !>                 are NF90_NETCDF4 only.
  !>                 The supported modes are:
  !>                     r: read
  !>                     w: write/create
  !<                     a: alter
  !
  !     INTENT(IN)
  !>        \param[in] "character(*) :: fname"
  !>        \param[in] "character(1) :: mode"
  !
  !     INTENT(INOUT)
  !         None
  !
  !     INTENT(OUT)
  !         None
  !
  !     RETURN
  !>        \return "type(NcDataset)"
  !
  !     EXAMPLE
  !         See test file
  !
  !     HISTORY
  !         \author David Schaefer
  !         \date June 2015
  !
  ! --------------------------------------------------------------------------------------
  type NcDataset

     character(256)                       :: fname !> Filename of the opened dataset
     character(1)                         :: mode  !> File open mode
     integer(i4)                          :: id    !> NetCDF id

   contains

     procedure, public :: initNcDataset

     procedure, private :: setGlobalAttributeChar
     procedure, private :: setGlobalAttributeI2
     procedure, private :: setGlobalAttributeI4
     procedure, private :: setGlobalAttributeSp
     procedure, private :: setGlobalAttributeDp

     procedure, private :: getGlobalAttributeChar
     procedure, private :: getGlobalAttributeI2
     procedure, private :: getGlobalAttributeI4
     procedure, private :: getGlobalAttributeSp
     procedure, private :: getGlobalAttributeDp

     procedure, private :: getDimensionByName
     procedure, private :: getDimensionById

     procedure, private :: setVariableWithTypes
     procedure, private :: setVariableWithNames
     procedure, private :: setVariableWithIds

     procedure, private :: getVariableByName

     ! -----------------------------------------------------------------------------------
     !
     !     NAME
     !         hasVariable
     !
     !     PURPOSE
     !>        \brief Check if variable exists
     !>        \details Returns true if a variable with the given name exists, false
     !>                 otherwise.
     !
     !     INTENT(IN)
     !>        \param[in] "character(*) :: name"
     !
     !     INTENT(INOUT)
     !         None
     !
     !     INTENT(OUT)
     !         None
     !
     !     RETURN
     !>        \return "logical"
     !
     !     EXAMPLE
     !         See test file
     !
     !     HISTORY
     !>        \author David Schaefer
     !>        \date June 2015
     !
     ! -----------------------------------------------------------------------------------
     procedure, public  :: hasVariable

     ! -----------------------------------------------------------------------------------
     !
     !     NAME
     !         hasDimension
     !
     !     PURPOSE
     !>        \brief Check if dimension exists
     !>        \details Returns true if a dimension with the given name exists, false
     !>                 otherwise.
     !
     !     INTENT(IN)
     !>        \param[in] "character(*) :: name"
     !
     !     INTENT(INOUT)
     !         None
     !
     !     INTENT(OUT)
     !         None
     !
     !     RETURN
     !>        \return "logical"
     !
     !     EXAMPLE
     !         See test file
     !
     !     HISTORY
     !>        \author David Schaefer
     !>        \date June 2015
     !
     ! -----------------------------------------------------------------------------------
     procedure, public  :: hasDimension

     ! -----------------------------------------------------------------------------------
     !
     !     NAME
     !         getUnlimitedDimension
     !
     !     PURPOSE
     !>        \brief Return the unlimited dimension of the dataset
     !>        \details Returns the NcDimension derived type of the unlimited dimension.
     !>                 The program will teminate abruptly if no such dimension exists.
     !
     !     INTENT(IN)
     !>        \param[in] "character(*) :: name"
     !
     !     INTENT(INOUT)
     !         None
     !
     !     INTENT(OUT)
     !         None
     !
     !     RETURN
     !>        \return "logical"
     !
     !     EXAMPLE
     !         See test file
     !
     !     HISTORY
     !>        \author David Schaefer
     !>        \date June 2015
     !
     ! -----------------------------------------------------------------------------------
     procedure, public  :: getUnlimitedDimension

     ! -----------------------------------------------------------------------------------
     !
     !     NAME
     !         isUnlimited
     !
     !     PURPOSE
     !>        \brief Check if the dataset is unlimited
     !>        \details Returns true if the dataset contains an unlimited dimension,
     !>                 false otherwise.
     !
     !     INTENT(IN)
     !         None
     !
     !     INTENT(INOUT)
     !         None
     !
     !     INTENT(OUT)
     !         None
     !
     !     RETURN
     !>        \return "logical"
     !
     !     EXAMPLE
     !         See test file
     !
     !     HISTORY
     !>        \author David Schaefer
     !>        \date June 2015
     !
     ! -----------------------------------------------------------------------------------
     procedure, public  :: isUnlimited => isDatasetUnlimited

     ! -----------------------------------------------------------------------------------
     !
     !     NAME
     !         close
     !
     !     PURPOSE
     !>        \brief Close the datset
     !>        \details Close the NetCDF datset. The program will teminate abruptly if
     !>                 the file cannot be closed correctly.
     !
     !     INTENT(IN)
     !         None
     !
     !     INTENT(INOUT)
     !         None
     !
     !     INTENT(OUT)
     !         None
     !
     !     RETURN
     !         None
     !
     !     EXAMPLE
     !         See test file
     !
     !     HISTORY
     !>        \author David Schaefer
     !>        \date June 2015
     !
     ! -----------------------------------------------------------------------------------
     procedure, public  :: close

     ! -----------------------------------------------------------------------------------
     !
     !     NAME
     !         setDimension
     !
     !     PURPOSE
     !>        \brief Create a new dimension
     !>        \details Create a new dimension of given length. A length < 0 indicates an
     !>                 unlimited dimension. The program will teminate abruptly if the
     !>                 dimension cannot be created.
     !
     !     INTENT(IN)
     !>        \param[in] "character(*) :: name"
     !>        \param[in] "integer(i4)  :: length"
     !
     !     INTENT(INOUT)
     !         None
     !
     !     INTENT(OUT)
     !         None
     !
     !     RETURN
     !>        \return NcDimension
     !
     !     EXAMPLE
     !         See test file
     !
     !     HISTORY
     !>        \author David Schaefer
     !>        \date June 2015
     !
     ! -----------------------------------------------------------------------------------
     procedure, public  :: setDimension

     ! -----------------------------------------------------------------------------------
     !
     !     NAME
     !         setAttribute
     !
     !     PURPOSE
     !>        \brief Create a new global Attribute
     !>        \details Create a new global attribute from given name and value.
     !>                 The program will teminate abruptly if the attribute cannot be
     !>                 created.
     !
     !     INTENT(IN)
     !>        \param[in] "character(*)                               :: name"
     !>        \param[in] "character(*)/integer(i4)/real(sp)/real(dp) :: value"
     !
     !     INTENT(INOUT)
     !         None
     !
     !     INTENT(OUT)
     !         None
     !
     !     RETURN
     !         None
     !
     !     EXAMPLE
     !         See test file
     !
     !     HISTORY
     !>        \author David Schaefer
     !>        \date June 2015
     !
     ! -----------------------------------------------------------------------------------
     generic,   public  :: setAttribute => &
          setGlobalAttributeChar, &
          setGlobalAttributeI2,   &
          setGlobalAttributeI4,   &
          setGlobalAttributeSp,   &
          setGlobalAttributeDp

     ! -----------------------------------------------------------------------------------
     !
     !     NAME
     !         getAttribute
     !
     !     PURPOSE
     !>        \brief Retrieve global attribute value
     !>        \details Retrieve the value for a global attribute specified by its name.
     !>                 The program will teminate abruptly if the attribute does not
     !>                 exist.
     !
     !     INTENT(IN)
     !>        \param[in] "character(*) :: name"
     !
     !     INTENT(INOUT)
     !         None
     !
     !     INTENT(OUT)
     !>        \param[out] "character(*)/integer(i4)/real(sp)/real(dp) :: value"
     !
     !     RETURN
     !         None
     !
     !     EXAMPLE
     !         See test file
     !
     !     HISTORY
     !>        \author David Schaefer
     !>        \date June 2015
     !
     ! -----------------------------------------------------------------------------------
     generic,   public  :: getAttribute => &
          getGlobalAttributeChar, &
          getGlobalAttributeI2,   &
          getGlobalAttributeI4,   &
          getGlobalAttributeSp,   &
          getGlobalAttributeDp

     ! -----------------------------------------------------------------------------------
     !
     !     NAME
     !         getDimension
     !
     !     PURPOSE
     !>        \brief Retrieve NcDimension
     !>        \details Retrieve the NcDimension derived type for the dimension specified by
     !>                 its name or id. The program will teminate abruptly if no such
     !>                 dimension exists.
     !
     !     INTENT(IN)
     !>        \param[in] "character(*)/integer(i4) :: name/id"
     !
     !     INTENT(INOUT)
     !         None
     !
     !     INTENT(OUT)
     !         None
     !
     !     RETURN
     !>        \return NcDimension
     !
     !     EXAMPLE
     !         See test file
     !
     !     HISTORY
     !>        \author David Schaefer
     !>        \date June 2015
     !
     ! -----------------------------------------------------------------------------------
     generic,   public  :: getDimension => &
          getDimensionById, &
          getDimensionByName

     ! -----------------------------------------------------------------------------------
     !
     !     NAME
     !         setVariable
     !
     !     PURPOSE
     !>        \brief Create a NetCDF variable
     !>        \details Create a NetCDF Variable with given name, data type and dimensions.
     !>                 All optional arguments to the nf90_def_var function are supported.
     !>                 The program will teminate abruptly if the variable cannot be
     !>                 created.
     !>                 Supported data types and their string encodings:
     !>                     NF90_SHORT  -> "i16"
     !>                     NF90_INT    -> "i32"
     !>                     NF90_FLOAT  -> "f32"
     !>                     NF90_DOUBLE -> "f64"
     !
     !     INTENT(IN)
     !>        \param[in] "character(*) :: name"
     !>        \param[in] "character(3) :: dtype"
     !>        \param[in] "integer(i4)/character(*)/type(NcDataset) :: dimensions(:)"
     !
     !     INTENT(IN), OPTIONAL
     !>       \param[in] "logical     :: contiguous"
     !>       \param[in] "integer(i4) :: chunksizes(:)"
     !>       \param[in] "integer(i4) :: deflate_level"
     !>       \param[in] "logical     :: shuffle"
     !>       \param[in] "logical     :: fletcher32"
     !>       \param[in] "integer(i4) :: endianess"
     !>       \param[in] "integer(i4) :: cache_size"
     !>       \param[in] "integer(i4) :: cache_nelems"
     !>       \param[in] "integer(i4) :: cache_preemption"
     !
     !     INTENT(INOUT)
     !         None
     !
     !     INTENT(OUT)
     !         None
     !
     !     RETURN
     !>        \return NcVariable
     !
     !     EXAMPLE
     !         See test file
     !
     !     HISTORY
     !>        \author David Schaefer
     !>        \date June 2015
     !
     ! -----------------------------------------------------------------------------------
     generic,   public  :: setVariable => &
          setVariableWithNames, &
          setVariableWithTypes, &
          setVariableWithIds

     ! -----------------------------------------------------------------------------------
     !
     !     NAME
     !         getVariable
     !
     !     PURPOSE
     !>        \brief Retrieve NcVariable
     !>        \details Retrieve the NcVariable derived type for the variable specified by its
     !>                 name. The program will teminate abruptly if no such dimension
     !>                 exists.
     !
     !     INTENT(IN)
     !>        \param[in] "character(*) :: name"
     !
     !     INTENT(INOUT)
     !         None
     !
     !     INTENT(OUT)
     !         None
     !
     !     RETURN
     !>        \return NcVariable
     !
     !     EXAMPLE
     !         See test file
     !
     !     HISTORY
     !>        \author David Schaefer
     !>        \date June 2015
     !
     ! -----------------------------------------------------------------------------------
     generic,  public  :: getVariable => &
          getVariableByName

  end type NcDataset

  interface NcDataset
     procedure newNcDataset
  end interface NcDataset

  ! --------------------------------------------------------------------------------------
  !
  !     NAME
  !         NcDimension
  !
  !     PURPOSE
  !>        \brief Provides the dimension access functionality
  !>        \details Bound to this derived type is some necessary inquire functionality.
  !>                 This type is not to be instantiated directly! Use the
  !>                 getDimension/setDimension functions of a NcDataset instance as a
  !>                 "constructor".
  !
  !     INTENT(IN)
  !         None
  !
  !     INTENT(INOUT)
  !         None
  !
  !     INTENT(OUT)
  !         None
  !
  !     RETURN
  !         None
  !
  !     EXAMPLE
  !         See test file
  !
  !     HISTORY
  !         \author David Schaefer
  !         \date June 2015
  !
  ! --------------------------------------------------------------------------------------
  type NcDimension

     integer(i4)     :: id      !> The NetCDF dimension id
     type(NcDataset) :: parent  !> The dimension's parent
     character(16)   :: name    !> The dimension's name

   contains

     procedure, public :: initNcDimension

     ! -----------------------------------------------------------------------------------
     !
     !     NAME
     !         getLength
     !
     !     PURPOSE
     !>        \brief Retrieve the dimension length
     !>        \details Return the length of the dimension
     !
     !     INTENT(IN)
     !         None
     !
     !     INTENT(INOUT)
     !         None
     !
     !     INTENT(OUT)
     !         None
     !
     !     RETURN
     !         \return "integer(i4)"
     !
     !     EXAMPLE
     !         See test file
     !
     !     HISTORY
     !         \author David Schaefer
     !         \date June 2015
     !
     ! -----------------------------------------------------------------------------------
     procedure, public :: getLength   => getDimensionLength

     ! -----------------------------------------------------------------------------------
     !
     !     NAME
     !         isUnlimited
     !
     !     PURPOSE
     !>        \brief Check if the dimension is unlimited
     !>        \details Returns true if the dimension is unlimited, false otherwise.
     !
     !     INTENT(IN)
     !         None
     !
     !     INTENT(INOUT)
     !         None
     !
     !     INTENT(OUT)
     !         None
     !
     !     RETURN
     !>        \return "logical"
     !
     !     EXAMPLE
     !         See test file
     !
     !     HISTORY
     !>        \author David Schaefer
     !>        \date June 2015
     !
     ! -----------------------------------------------------------------------------------
     procedure, public :: isUnlimited => isUnlimitedDimension

  end type NcDimension

  interface NcDimension
     procedure newNcDimension
  end interface NcDimension

  interface operator (==)
     procedure equalNcDimensions
  end interface operator (==)

  ! --------------------------------------------------------------------------------------
  !
  !     NAME
  !         NcVariable
  !
  !     PURPOSE
  !>        \brief Provides the variable releated read/write functionality
  !>        \details Bound to this derived type is the main variable related NetCDF
  !>                 functionality, i.e. reading/writing data and attributes, as well
  !>                 as some inquire procedures.
  !>                 This type is not to be instatiated directly, instead the
  !>                 getVariable/setVariable functions of a NcDataset instance should
  !>                 be used as a "constructor".
  !
  !     INTENT(IN)
  !         None
  !
  !     INTENT(INOUT)
  !         None
  !
  !     INTENT(OUT)
  !         None
  !
  !     RETURN
  !         None
  !
  !     EXAMPLE
  !         See test file
  !
  !     HISTORY
  !         \author David Schaefer
  !         \date June 2015
  !
  ! --------------------------------------------------------------------------------------
  type NcVariable

     character(32)   :: name     !> Variable name
     integer(i4)     :: id       !> Variable id
     type(NcDataset) :: parent   !> The variables's parent

   contains

     procedure, public :: initNcVariable

     procedure, private :: setVariableAttributeChar
     procedure, private :: setVariableAttributeI2
     procedure, private :: setVariableAttributeI4
     procedure, private :: setVariableAttributeSp
     procedure, private :: setVariableAttributeDp

     procedure, private :: getVariableAttributeChar
     procedure, private :: getVariableAttributeI2
     procedure, private :: getVariableAttributeI4
     procedure, private :: getVariableAttributeSp
     procedure, private :: getVariableAttributeDp

     procedure, private :: setDataScalarI2
     procedure, private :: setData1dI2
     procedure, private :: setData2dI2
     procedure, private :: setData3dI2
     procedure, private :: setData4dI2
     procedure, private :: setData5dI2
     procedure, private :: setDataScalarI4
     procedure, private :: setData1dI4
     procedure, private :: setData2dI4
     procedure, private :: setData3dI4
     procedure, private :: setData4dI4
     procedure, private :: setData5dI4
     procedure, private :: setDataScalarSp
     procedure, private :: setData1dSp
     procedure, private :: setData2dSp
     procedure, private :: setData3dSp
     procedure, private :: setData4dSp
     procedure, private :: setData5dSp
     procedure, private :: setDataScalarDp
     procedure, private :: setData1dDp
     procedure, private :: setData2dDp
     procedure, private :: setData3dDp
     procedure, private :: setData4dDp
     procedure, private :: setData5dDp

     procedure, private :: getDataScalarI2
     procedure, private :: getData1dI2
     procedure, private :: getData2dI2
     procedure, private :: getData3dI2
     procedure, private :: getData4dI2
     procedure, private :: getData5dI2
     procedure, private :: getDataScalarI4
     procedure, private :: getData1dI4
     procedure, private :: getData2dI4
     procedure, private :: getData3dI4
     procedure, private :: getData4dI4
     procedure, private :: getData5dI4
     procedure, private :: getDataScalarSp
     procedure, private :: getData1dSp
     procedure, private :: getData2dSp
     procedure, private :: getData3dSp
     procedure, private :: getData4dSp
     procedure, private :: getData5dSp
     procedure, private :: getDataScalarDp
     procedure, private :: getData1dDp
     procedure, private :: getData2dDp
     procedure, private :: getData3dDp
     procedure, private :: getData4dDp
     procedure, private :: getData5dDp

     procedure, private :: setVariableFillValueI2
     procedure, private :: setVariableFillValueI4
     procedure, private :: setVariableFillValueSp
     procedure, private :: setVariableFillValueDp

     procedure, private :: getVariableFillValueI2
     procedure, private :: getVariableFillValueI4
     procedure, private :: getVariableFillValueSp
     procedure, private :: getVariableFillValueDp

     ! -----------------------------------------------------------------------------------
     !
     !     NAME
     !         getNoDimensions
     !
     !     PURPOSE
     !>        \brief Retrieve the number of dimensions
     !>        \details Return the number of dimensions associated with variable
     !
     !     INTENT(IN)
     !         None
     !
     !     INTENT(INOUT)
     !         None
     !
     !     INTENT(OUT)
     !         None
     !
     !     RETURN
     !         \return "integer(i4)"
     !
     !     EXAMPLE
     !         See test file
     !
     !     HISTORY
     !         \author David Schaefer
     !         \date June 2015
     !
     ! -----------------------------------------------------------------------------------
     procedure, public  :: getNoDimensions

     ! -----------------------------------------------------------------------------------
     !
     !     NAME
     !         getDimensions
     !
     !     PURPOSE
     !>        \brief Retrieve the variable dimensions
     !>        \details Return the ids of the dimensions associated with variable.
     !
     !     INTENT(IN)
     !         None
     !
     !     INTENT(INOUT)
     !         None
     !
     !     INTENT(OUT)
     !         None
     !
     !     RETURN
     !         \return "integer(i4), alloctable, dimension(:)"
     !
     !     EXAMPLE
     !         See test file
     !
     !     HISTORY
     !         \author David Schaefer
     !         \date June 2015
     !
     ! -----------------------------------------------------------------------------------
     procedure, public  :: getDimensions => getVariableDimensions

     ! -----------------------------------------------------------------------------------
     !
     !     NAME
     !         getShape
     !
     !     PURPOSE
     !>        \brief Retrieve the shape of the variable
     !>        \details Return the shape of the variable.
     !
     !     INTENT(IN)
     !         None
     !
     !     INTENT(INOUT)
     !         None
     !
     !     INTENT(OUT)
     !         None
     !
     !     RETURN
     !         \return "integer(i4), alloctable, dimension(:)"
     !
     !     EXAMPLE
     !         See test file
     !
     !     HISTORY
     !         \author David Schaefer
     !         \date June 2015
     !
     ! -----------------------------------------------------------------------------------
     procedure, public  :: getShape      => getVariableShape

     ! -----------------------------------------------------------------------------------
     !
     !     NAME
     !         getDtype
     !
     !     PURPOSE
     !>        \brief Retrieve the variable data type
     !>        \details Return the encoded data type of the variable.
     !>                 Data type encodeings
     !>                     "f32" -> NF90_FLOAT
     !>                     "f64" -> NF90_DOUBLE
     !>                     "i16" -> NF90_SHORT
     !>                     "i32" -> NF90_INT
     !>                     "i64" -> NF90_INT64
     !
     !     INTENT(IN)
     !         None
     !
     !     INTENT(INOUT)
     !         None
     !
     !     INTENT(OUT)
     !         None
     !
     !     RETURN
     !         \return "character(3)"
     !
     !     EXAMPLE
     !         See test file
     !
     !     HISTORY
     !         \author David Schaefer
     !         \date June 2015
     !
     ! -----------------------------------------------------------------------------------
     procedure, public  :: getDtype      => getVariableDtype

     ! -----------------------------------------------------------------------------------
     !
     !     NAME
     !         hasAttribute
     !
     !     PURPOSE
     !>        \brief Check if attribute exists
     !>        \details Returns true if an attribute with the given name exists, false
     !>                 otherwise.
     !
     !     INTENT(IN)
     !         \param[in] "character(*) :: name"
     !
     !     INTENT(INOUT)
     !         None
     !
     !     INTENT(OUT)
     !         None
     !
     !     RETURN
     !         \return "logical"
     !
     !     EXAMPLE
     !         See test file
     !
     !     HISTORY
     !         \author David Schaefer
     !         \date June 2015
     !
     ! -----------------------------------------------------------------------------------
     procedure, public  :: hasAttribute

     ! -----------------------------------------------------------------------------------
     !
     !     NAME
     !         isUnlimited
     !
     !     PURPOSE
     !>        \brief Check if the variable is unlimited
     !>        \details Returns true if the variable has an unlimited dimension,
     !>                 false otherwise.
     !
     !     INTENT(IN)
     !         None
     !
     !     INTENT(INOUT)
     !         None
     !
     !     INTENT(OUT)
     !         None
     !
     !     RETURN
     !>        \return "logical"
     !
     !     EXAMPLE
     !         See test file
     !
     !     HISTORY
     !>        \author David Schaefer
     !>        \date June 2015
     !
     ! -----------------------------------------------------------------------------------
     procedure, public  :: isUnlimited   => isUnlimitedVariable

     ! -----------------------------------------------------------------------------------
     !
     !     NAME
     !         setData
     !
     !     PURPOSE
     !>        \brief Write data to variable
     !>        \details Write the given data into the variable at an optionally given
     !>                 position.
     !>                 All optional arguments to the nf90_put_var function are supported.
     !>                 A write error will result in abrupt program termination.
     !
     !     INTENT(IN)
     !>        \param[in] "integer(i4)/real(sp)/real(dp) :: &
     !>                                   values/(:)/(:,:)/(:,:,:)/(:,:,:,:)/(:,:,:,:,:)"
     !
     !     INTENT(IN), OPTIONAL
     !>        \param[in]  "integer(i4) :: start(:), count(:), stride(:), map(:)"
     !
     !     INTENT(INOUT)
     !         None
     !
     !     INTENT(OUT)
     !         None
     !
     !     RETURN
     !         None
     !
     !     EXAMPLE
     !         See test file
     !
     !     HISTORY
     !>        \author David Schaefer
     !>        \date June 2015
     !
     ! -----------------------------------------------------------------------------------
     generic, public :: setData => &
          setDataScalarI2, &
          setData1dI2, &
          setData2dI2, &
          setData3dI2, &
          setData4dI2, &
          setData5dI2, &
          setDataScalarI4, &
          setData1dI4, &
          setData2dI4, &
          setData3dI4, &
          setData4dI4, &
          setData5dI4, &
          setDataScalarSp, &
          setData1dSp, &
          setData2dSp, &
          setData3dSp, &
          setData4dSp, &
          setData5dSp, &
          setDataScalarDp, &
          setData1dDp, &
          setData2dDp, &
          setData3dDp, &
          setData4dDp, &
          setData5dDp

     ! -----------------------------------------------------------------------------------
     !
     !     NAME
     !         getData
     !
     !     PURPOSE
     !>        \brief Retrieve data
     !>        \details Read the data from an optionally given position.
     !>                 All optional arguments to the nf90_get_var function are supported.
     !>                 A read error will result in abrupt program termination.
     !
     !     INTENT(IN)
     !         None
     !
     !     INTENT(IN), OPTIONAL
     !>        \param[in] "integer(i4) :: start(:), count(:), stride(:), map(:)"
     !
     !     INTENT(INOUT)
     !         None
     !
     !     INTENT(OUT)
     !>      \param[out] "integer(i4)/real(sp)/real(dp), allocatable :: &
     !>                                   values/(:)/(:,:)/(:,:,:)/(:,:,:,:)/(:,:,:,:,:)"
     !
     !
     !     RETURN
     !         None
     !
     !     EXAMPLE
     !         See test file
     !
     !     HISTORY
     !>        \author David Schaefer
     !>        \date June 2015
     !
     ! -----------------------------------------------------------------------------------
     generic, public :: getData => &
          getDataScalarI2, &
          getData1dI2, &
          getData2dI2, &
          getData3dI2, &
          getData4dI2, &
          getData5dI2, &
          getDataScalarI4, &
          getData1dI4, &
          getData2dI4, &
          getData3dI4, &
          getData4dI4, &
          getData5dI4, &
          getDataScalarSp, &
          getData1dSp, &
          getData2dSp, &
          getData3dSp, &
          getData4dSp, &
          getData5dSp, &
          getDataScalarDp, &
          getData1dDp, &
          getData2dDp, &
          getData3dDp, &
          getData4dDp, &
          getData5dDp


     ! -----------------------------------------------------------------------------------
     !
     !     NAME
     !         setFillValue
     !
     !     PURPOSE
     !>        \brief Set the variable fill value
     !>        \details Define the variable fill value.
     !>                 A write error results in abrupt program temination.
     !>        \note This procedure must be called AFTER the variable was created but
     !>              BEFORE data is first written.
     !
     !     INTENT(IN)
     !>        \param[in] "integer(i4)/real(sp)/real(dp) :: fvalue"
     !
     !     INTENT(INOUT)
     !         None
     !
     !     INTENT(OUT)
     !         None
     !
     !     RETURN
     !         None
     !
     !     EXAMPLE
     !         See test file
     !
     !     HISTORY
     !>        \author David Schaefer
     !>        \date June 2015
     !
     ! -----------------------------------------------------------------------------------
          generic, public :: setFillValue => &
          setVariableFillValueI2, &
          setVariableFillValueI4, &
          setVariableFillValueSp, &
          setVariableFillValueDp

     ! -----------------------------------------------------------------------------------
     !
     !     NAME
     !         getFillValue
     !
     !     PURPOSE
     !>        \brief  Retrieve the variable fill value
     !>        \details Retrieve the variable fill value or a default value if fill value
     !>                 was not explicitly set.
     !>                 A read error results in abrupt program temination.
     !
     !     INTENT(IN)
     !         None
     !
     !     INTENT(INOUT)
     !         None
     !
     !     INTENT(OUT)
     !>        \param[out] "integer(i4)/real(sp)/real(dp) :: fvalue"
     !
     !     RETURN
     !         None
     !
     !     EXAMPLE
     !         See test file
     !
     !     HISTORY
     !>        \author David Schaefer
     !>        \date June 2015
     !
     ! -----------------------------------------------------------------------------------
     generic, public :: getFillValue => &
          getVariableFillValueI2, &
          getVariableFillValueI4, &
          getVariableFillValueSp, &
          getVariableFillValueDp

     ! -----------------------------------------------------------------------------------
     !
     !     NAME
     !         setAttribute
     !
     !     PURPOSE
     !>        \brief Create a new variable attribute
     !>        \details Create a new variable attribute from given name and value.
     !>                 A write error results in abrupt program temination.
     !
     !     INTENT(IN)
     !>        \param[in] "character(*)                               :: name"
     !>        \param[in] "character(*)/integer(i4)/real(sp)/real(dp) :: value"
     !
     !     INTENT(INOUT)
     !         None
     !
     !     INTENT(OUT)
     !         None
     !
     !     RETURN
     !         None
     !
     !     EXAMPLE
     !         See test file
     !
     !     HISTORY
     !>        \author David Schaefer
     !>        \date June 2015
     !
     ! -----------------------------------------------------------------------------------
     generic, public :: setAttribute => &
          setVariableAttributeChar, &
          setVariableAttributeI2,   &
          setVariableAttributeI4,   &
          setVariableAttributeSp,   &
          setVariableAttributeDp

     ! -----------------------------------------------------------------------------------
     !
     !     NAME
     !         getAttribute
     !
     !     PURPOSE
     !>        \brief Retrieve variable attribute value
     !>        \details Retrieve the value for a variable attribute specified by its name.
     !>                 The program will teminate abruptly if the attribute does not
     !>                 exist.
     !
     !     INTENT(IN)
     !>        \param[in] "character(*) :: name"
     !
     !     INTENT(INOUT)
     !         None
     !
     !     INTENT(OUT)
     !>        \param[out] "character(*)/integer(i4)/real(sp)/real(dp) :: value"
     !
     !     RETURN
     !         None
     !
     !     EXAMPLE
     !         See test file
     !
     !     HISTORY
     !>        \author David Schaefer
     !>        \date June 2015
     !
     ! -----------------------------------------------------------------------------------
     generic, public :: getAttribute => &
          getVariableAttributeChar, &
          getVariableAttributeI2,   &
          getVariableAttributeI4,   &
          getVariableAttributeSp,   &
          getVariableAttributeDp

  end type NcVariable

  interface NcVariable
     procedure newNcVariable
  end interface NcVariable

contains

  subroutine initNcVariable(self, id, parent, name)
    class(NcVariable), intent(inout) :: self
    integer(i4)      , intent(in)    :: id
    type(NcDataset)  , intent(in)    :: parent
    character(*)     , intent(in)    :: name

    self%id     = id
    self%parent = parent
    self%name   = name
  end subroutine initNcVariable

  subroutine initNcDimension(self, id, parent, name)
    class(NcDimension), intent(inout) :: self
    integer(i4)       , intent(in)    :: id
    type(NcDataset)   , intent(in)    :: parent
    character(*)      , intent(in)    :: name

    self%id     = id
    self%parent = parent
    self%name   = name
  end subroutine initNcDimension

  subroutine initNcDataset(self, fname, mode)
    class(NcDataset), intent(inout) :: self
    character(*)    , intent(in)    :: fname
    character(1)    , intent(in)    :: mode
    integer(i4)                     :: status

    select case(mode)
    case("w")
       status = nf90_create(trim(fname), NF90_NETCDF4, self%id)
    case("r")
       status = nf90_open(trim(fname), NF90_NOWRITE, self%id)
    case("a")
       status = nf90_open(trim(fname), NF90_WRITE, self%id)
    case default
       write(*,*) "Mode argument must be in 'w','r','a' ! "
       stop 1
    end select
    call check(status,"Failed to open file: " // fname)

    self%fname = fname
    self%mode  = mode
  end subroutine initNcDataset

  type(NcVariable) function newNcVariable(id, parent, name)
    integer(i4)    , intent(in) :: id
    type(NcDataset), intent(in) :: parent
    character(*)   , intent(in) :: name

    call newNcVariable%initNcVariable(id, parent, name)
  end function newNcVariable

  type(NcDimension) function newNcDimension(id, parent, name)
    integer(i4)    , intent(in) :: id
    type(NcDataset), intent(in) :: parent
    character(*)   , intent(in) :: name

    call newNcDimension%initNcDimension(id, parent, name)
  end function newNcDimension

  type(NcDataset) function newNcDataset(fname, mode)
    character(*), intent(in) :: fname
    character(1), intent(in) :: mode

    call newNcDataset%initNcDataset(fname,mode)
  end function newNcDataset

  subroutine close(self)
    class(NcDataset) :: self

    call check(nf90_close(self%id), "Failed to close file: "//self%fname)
  end subroutine close


  function getDimensionLength(self)
    class(NcDimension), intent(in) :: self
    integer(i4)                    :: getDimensionLength

    call check(nf90_inquire_dimension(self%parent%id,self%id,len=getDimensionLength),&
         "Failed to inquire dimension: "//self%name)
  end function getDimensionLength

  function isDatasetUnlimited(self)
    class(NcDataset), intent(in) :: self
    logical                      :: isDatasetUnlimited
    integer(i4)                  :: dimid

    call check(nf90_inquire(self%id,unlimitedDimId=dimid), &
         "Failed to inquire file "//self%fname)
    isDatasetUnlimited = (dimid .ne. -1)
  end function isDatasetUnlimited

  function getUnlimitedDimension(self)
    class(NcDataset), intent(in) :: self
    type(NcDimension)            :: getUnlimitedDimension
    integer(i4)                  :: dimid

    call check(nf90_inquire(self%id,unlimitedDimId=dimid), &
         "Failed to inquire file "//self%fname)

    if (dimid .eq. -1) then
       write(*,*) "Dataset has no unlimited dimension"
       stop 1
    end if

    getUnlimitedDimension = self%getDimension(dimid)
  end function getUnlimitedDimension

  logical function equalNcDimensions(dim1, dim2)
    type(NcDimension), intent(in) :: dim1, dim2

    equalNcDimensions = (dim1%id .eq. dim2%id)
  end function equalNcDimensions

  function isUnlimitedDimension(self)
    class(NcDimension), intent(in) :: self
    logical                        :: isUnlimitedDimension

    isUnlimitedDimension = .false.
    if (self%parent%isUnlimited()) then
       isUnlimitedDimension = (self == self%parent%getUnlimitedDimension())
    end if
  end function isUnlimitedDimension

  function setDimension(self, name, length)
    class(NcDataset), intent(in) :: self
    character(*)    , intent(in) :: name
    integer(i4)     , intent(in) :: length
    type(NcDimension)            :: setDimension
    integer(i4)                  :: id, dimlength

    if (length .le. 0) then
       dimlength = NF90_UNLIMITED
    else
       dimlength = length
    end if

    call check(nf90_def_dim(self%id, name, dimlength, id), &
         "Failed to create dimension: " // name)

    setDimension = NcDimension(id,self,trim(name))
  end function setDimension

  function hasVariable(self, name)
    class(NcDataset), intent(in) :: self
    character(*)    , intent(in) :: name
    logical                      :: hasVariable
    integer(i4)                  :: tmpid

    hasVariable = (nf90_inq_varid(self%id,name,tmpid) .eq. NF90_NOERR)
  end function hasVariable

  function hasDimension(self, name)
    class(NcDataset), intent(in) :: self
    character(*)    , intent(in) :: name
    logical                      :: hasDimension
    integer(i4)                  :: tmpid

    hasDimension = (nf90_inq_dimid(self%id,name,tmpid) .eq. NF90_NOERR)
  end function hasDimension


  function setVariableWithIds(self, name, dtype, dimensions, contiguous, &
       chunksizes, deflate_level, shuffle, fletcher32, endianness, &
       cache_size, cache_nelems, cache_preemption)
    class(NcDataset), intent(in)           :: self
    character(*)    , intent(in)           :: name
    character(3)    , intent(in)           :: dtype
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
    setVariableWithIds = NcVariable(varid,self,name)
  end function setVariableWithIds

  function setVariableWithNames(self, name, dtype, dimensions, contiguous, &
       chunksizes, deflate_level, shuffle, fletcher32, endianness, &
       cache_size, cache_nelems, cache_preemption)

    class(NcDataset), intent(in)              :: self
    character(*)    , intent(in)              :: name
    character(3)    , intent(in)              :: dtype
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
    class(NcDataset) , intent(in)              :: self
    character(*)     , intent(in)              :: name
    character(3)     , intent(in)              :: dtype
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
    class(NcDataset), intent(in) :: self
    integer(i4)                  :: id
    type(NcDimension)            :: getDimensionById
    character(32)                :: msg, name

    write(msg,*) id
    call check(nf90_inquire_dimension(self%id,id,name), &
         "Could not inquire dimension: " // msg)
    getDimensionById = NcDimension(id,self,name)
  end function getDimensionById

  function getDimensionByName(self, name)
    class(NcDataset), intent(in) :: self
    character(*)                 :: name
    type(NcDimension)            :: getDimensionByName
    integer(i4)                  :: id

    call check(nf90_inq_dimid(self%id,name,id), &
         "Could not inquire dimension: " // name)
    getDimensionByName = self%getDimensionById(id)
  end function getDimensionByName

  function getVariableByName(self, name)
    class(NcDataset), intent(in) :: self
    character(*)    , intent(in) :: name
    type(NcVariable)             :: getVariableByName
    integer(i4)                  :: id

    call check(nf90_inq_varid(self%id,name,id), &
         "Could not inquire variable: " // name)
    getVariableByName = NcVariable(id,self,name)
  end function getVariableByName

  function getNoDimensions(self)
    class(NcVariable), intent(in) :: self
    integer(i4)                   :: getNoDimensions

    call check(nf90_inquire_variable(self%parent%id,self%id,ndims=getNoDimensions), &
         "Could not inquire variable: " // self%name)
  end function getNoDimensions

  function getVariableDimensions (self)
    class(NcVariable), intent(in)  :: self
    type(NcDimension), allocatable :: getVariableDimensions(:)
    integer(i4)      , allocatable :: dimids(:)
    integer(i4)                    :: ii , ndims

    ndims = self%getNoDimensions()
    allocate(dimids(ndims),getVariableDimensions (ndims))
    call check(nf90_inquire_variable(self%parent%id,self%id,dimids=dimids), &
         "Could not inquire variable: " // self%name)

    do ii=1,ndims
       getVariableDimensions (ii) = self%parent%getDimension(dimids(ii))
    end do
  end function getVariableDimensions

  function getVariableShape(self)
    class(NcVariable), intent(in)  :: self
    integer(i4)      , allocatable :: getVariableShape(:)
    type(NcDimension), allocatable :: dims(:)
    integer(i4)                    :: ii

    dims = self%getDimensions()
    allocate(getVariableShape(size(dims)))

    do ii = 1,size(dims)
       getVariableShape(ii) = dims(ii)%getLength()
    end do
  end function getVariableShape

  function getVariableDtype(self)
    class(NcVariable), intent(in) :: self
    integer(i4)                   :: dtype
    character(3)                  :: getVariableDtype

    call check(nf90_inquire_variable(self%parent%id,self%id,xtype=dtype),&
         "Could not inquire variable: " // self%name)
    getVariableDtype = getDtypeFromInteger(dtype)
  end function getVariableDtype

  function isUnlimitedVariable(self)
    class(NcVariable), intent(in)  :: self
    logical                        :: isUnlimitedVariable
    type(NcDimension), allocatable :: dims(:)
    type(NcDimension)              :: dim
    integer(i4)                    :: ii

    isUnlimitedVariable = .false.
    dims = self%getDimensions()

    do ii = 1,size(dims)
       dim = dims(ii)
       if (dim%isUnlimited()) then
          isUnlimitedVariable = .true.
       end if
    end do
  end function isUnlimitedVariable

  function hasAttribute(self,name)
    class(NcVariable), intent(in) :: self
    character(*)     , intent(in) :: name
    logical                       :: hasAttribute
    integer(i4)                   :: status

    status = nf90_inquire_attribute(self%parent%id, self%id, name)
    hasAttribute = (status .eq. NF90_NOERR)
  end function hasAttribute

  subroutine setGlobalAttributeChar(self, name, data)
    class(NcDataset), intent(in) :: self
    character(*)    , intent(in) :: name
    character(*)    , intent(in) :: data

    call check(nf90_put_att(self%id,NF90_GLOBAL,name,data), &
         "Failed to write attribute: " // name )
  end subroutine setGlobalAttributeChar

  subroutine setGlobalAttributeI2(self, name, data)
    class(NcDataset), intent(in) :: self
    character(*)    , intent(in) :: name
    integer(i2)     , intent(in) :: data

    call check(nf90_put_att(self%id,NF90_GLOBAL,name,data), &
         "Failed to write attribute: " // name )
  end subroutine setGlobalAttributeI2

  subroutine setGlobalAttributeI4(self, name, data)
    class(NcDataset), intent(in) :: self
    character(*)    , intent(in) :: name
    integer(i4)     , intent(in) :: data

    call check(nf90_put_att(self%id,NF90_GLOBAL,name,data), &
         "Failed to write attribute: " // name )
  end subroutine setGlobalAttributeI4

  subroutine setGlobalAttributeSp(self, name, data)
    class(NcDataset), intent(in) :: self
    character(*)    , intent(in) :: name
    real(sp)        , intent(in) :: data

    call check(nf90_put_att(self%id,NF90_GLOBAL,name,data), &
         "Failed to write attribute: " // name )
  end subroutine setGlobalAttributeSp

  subroutine setGlobalAttributeDp(self, name, data)
    class(NcDataset), intent(in) :: self
    character(*)    , intent(in) :: name
    real(dp)        , intent(in) :: data

    call check(nf90_put_att(self%id,NF90_GLOBAL,name,data), &
         "Failed to write attribute: " // name )
  end subroutine setGlobalAttributeDp

  subroutine getGlobalAttributeChar(self, name, avalue)
    class(NcDataset), intent(in)  :: self
    character(*)    , intent(in)  :: name
    character(*)    , intent(out) :: avalue
    integer(i4)                   :: length

    call check(nf90_inquire_attribute(self%id,NF90_GLOBAL,name,len=length),&
         "Could not inquire attribute "//name)
    call check(nf90_get_att(self%id,NF90_GLOBAL,name,avalue), &
         "Could not read attribute "//name)
  end subroutine getGlobalAttributeChar

  subroutine getGlobalAttributeI2(self, name, avalue)
    class(NcDataset), intent(in)  :: self
    character(*)    , intent(in)  :: name
    integer(i2)     , intent(out) :: avalue
    integer(i4)                   :: length

    call check(nf90_inquire_attribute(self%id,NF90_GLOBAL,name,len=length),&
         "Could not inquire attribute "//name)
    call check(nf90_get_att(self%id,NF90_GLOBAL,name,avalue), &
         "Could not read attribute "//name)
  end subroutine getGlobalAttributeI2

  subroutine getGlobalAttributeI4(self, name, avalue)
    class(NcDataset), intent(in)  :: self
    character(*)    , intent(in)  :: name
    integer(i4)     , intent(out) :: avalue
    integer(i4)                   :: length

    call check(nf90_inquire_attribute(self%id,NF90_GLOBAL,name,len=length),&
         "Could not inquire attribute "//name)
    call check(nf90_get_att(self%id,NF90_GLOBAL,name,avalue), &
         "Could not read attribute "//name)
  end subroutine getGlobalAttributeI4

  subroutine getGlobalAttributeSp(self, name, avalue)
    class(NcDataset), intent(in)  :: self
    character(*)    , intent(in)  :: name
    real(sp)        , intent(out) :: avalue
    integer(i4)                   :: length

    call check(nf90_inquire_attribute(self%id,NF90_GLOBAL,name,len=length),&
         "Could not inquire attribute "//name)
    call check(nf90_get_att(self%id,NF90_GLOBAL,name,avalue), &
         "Could not read attribute "//name)
  end subroutine getGlobalAttributeSp

  subroutine getGlobalAttributeDp(self, name, avalue)
    class(NcDataset), intent(in)  :: self
    character(*)    , intent(in)  :: name
    real(dp)        , intent(out) :: avalue
    integer(i4)                   :: length

    call check(nf90_inquire_attribute(self%id,NF90_GLOBAL,name,len=length),&
         "Could not inquire attribute "//name)
    call check(nf90_get_att(self%id,NF90_GLOBAL,name,avalue), &
         "Could not read attribute "//name)
  end subroutine getGlobalAttributeDp

  subroutine setVariableAttributeChar(self, name, data)
    class(NcVariable), intent(in) :: self
    character(*)     , intent(in) :: name
    character(*)     , intent(in) :: data

    call check(nf90_put_att(self%parent%id,self%id,name,data), &
         "Failed to write attribute: " // name )
  end subroutine setVariableAttributeChar

  subroutine setVariableAttributeI2(self, name, data)
    class(NcVariable), intent(in) :: self
    character(*)     , intent(in) :: name
    integer(i2)      , intent(in) :: data

    call check(nf90_put_att(self%parent%id,self%id,name,data), &
         "Failed to write attribute: " // name )
  end subroutine setVariableAttributeI2

  subroutine setVariableAttributeI4(self, name, data)
    class(NcVariable), intent(in) :: self
    character(*)     , intent(in) :: name
    integer(i4)      , intent(in) :: data

    call check(nf90_put_att(self%parent%id,self%id,name,data), &
         "Failed to write attribute: " // name )
  end subroutine setVariableAttributeI4

  subroutine setVariableAttributeSp(self, name, data)
    class(NcVariable), intent(in) :: self
    character(*)     , intent(in) :: name
    real(sp)         , intent(in) :: data

    call check(nf90_put_att(self%parent%id,self%id,name,data), &
         "Failed to write attribute: " // name )
  end subroutine setVariableAttributeSp

  subroutine setVariableAttributeDp(self, name, data)
    class(NcVariable), intent(in) :: self
    character(*)     , intent(in) :: name
    real(dp)         , intent(in) :: data

    call check(nf90_put_att(self%parent%id,self%id,name,data), &
         "Failed to write attribute: " // name )
  end subroutine setVariableAttributeDp

  subroutine getVariableAttributeChar(self, name, avalue)
    class(NcVariable), intent(in)  :: self
    character(*)     , intent(in)  :: name
    character(*)     , intent(out) :: avalue
    integer(i4)                    :: length

    call check(nf90_inquire_attribute(self%parent%id,self%id,name,len=length),&
         "Could not inquire attribute "//name)
    call check(nf90_get_att(self%parent%id,self%id,name,avalue), &
         "Could not read attribute "//name)
  end subroutine getVariableAttributeChar

  subroutine getVariableAttributeI2(self, name, avalue)
    class(NcVariable), intent(in)  :: self
    character(*)     , intent(in)  :: name
    integer(i2)      , intent(out) :: avalue
    integer(i4)                    :: length

    call check(nf90_inquire_attribute(self%parent%id,self%id,name,len=length),&
         "Could not inquire attribute "//name)
    call check(nf90_get_att(self%parent%id,self%id,name,avalue), &
         "Could not read attribute "//name)
  end subroutine getVariableAttributeI2

  subroutine getVariableAttributeI4(self, name, avalue)
    class(NcVariable), intent(in)  :: self
    character(*)     , intent(in)  :: name
    integer(i4)      , intent(out) :: avalue
    integer(i4)                    :: length

    call check(nf90_inquire_attribute(self%parent%id,self%id,name,len=length),&
         "Could not inquire attribute "//name)
    call check(nf90_get_att(self%parent%id,self%id,name,avalue), &
         "Could not read attribute "//name)
  end subroutine getVariableAttributeI4

  subroutine getVariableAttributeSp(self, name, avalue)
    class(NcVariable), intent(in)  :: self
    character(*)     , intent(in)  :: name
    real(sp)         , intent(out) :: avalue
    integer(i4)                    :: length

    call check(nf90_inquire_attribute(self%parent%id,self%id,name,len=length),&
         "Could not inquire attribute "//name)
    call check(nf90_get_att(self%parent%id,self%id,name,avalue), &
         "Could not read attribute "//name)
  end subroutine getVariableAttributeSp

  subroutine getVariableAttributeDp(self, name, avalue)
    class(NcVariable), intent(in)  :: self
    character(*)     , intent(in)  :: name
    real(dp)         , intent(out) :: avalue
    integer(i4)                    :: length

    call check(nf90_inquire_attribute(self%parent%id,self%id,name,len=length),&
         "Could not inquire attribute "//name)
    call check(nf90_get_att(self%parent%id,self%id,name,avalue), &
         "Could not read attribute "//name)
  end subroutine getVariableAttributeDp

  subroutine setVariableFillValueI2(self, fvalue)
    class(NcVariable), intent(in)  :: self
    integer(i2)      , intent(in)  :: fvalue

    if (.not. self%hasAttribute("_FillValue")) then
       call self%setAttribute("_FillValue",fvalue)
    end if
  end subroutine setVariableFillValueI2

  subroutine setVariableFillValueI4(self, fvalue)
    class(NcVariable), intent(in)  :: self
    integer(i4)      , intent(in)  :: fvalue

    if (.not. self%hasAttribute("_FillValue")) then
       call self%setAttribute("_FillValue",fvalue)
    end if
  end subroutine setVariableFillValueI4

  subroutine setVariableFillValueSp(self, fvalue)
    class(NcVariable), intent(in)  :: self
    real(sp)         , intent(in)  :: fvalue

    if (.not. self%hasAttribute("_FillValue")) then
       call self%setAttribute("_FillValue",fvalue)
    end if
  end subroutine setVariableFillValueSp

  subroutine setVariableFillValueDp(self, fvalue)
    class(NcVariable), intent(in)  :: self
    real(dp)         , intent(in)  :: fvalue

    if (.not. self%hasAttribute("_FillValue")) then
       call self%setAttribute("_FillValue",fvalue)
    end if
  end subroutine setVariableFillValueDp

  subroutine getVariableFillValueI2(self, fvalue)
    class(NcVariable), intent(in)  :: self
    integer(i2)      , intent(out) :: fvalue

    if (self%hasAttribute("_FillValue")) then
       call self%getAttribute("_FillValue", fvalue)
    else
       fvalue = NF90_FILL_SHORT
    end if
  end subroutine getVariableFillValueI2

  subroutine getVariableFillValueI4(self, fvalue)
    class(NcVariable), intent(in)  :: self
    integer(i4)      , intent(out) :: fvalue

    if (self%hasAttribute("_FillValue")) then
       call self%getAttribute("_FillValue", fvalue)
    else
       fvalue = NF90_FILL_INT
    end if
  end subroutine getVariableFillValueI4

  subroutine getVariableFillValueSp(self, fvalue)
    class(NcVariable), intent(in)  :: self
    real(sp)         , intent(out) :: fvalue

    if (self%hasAttribute("_FillValue")) then
       call self%getAttribute("_FillValue", fvalue)
    else
       fvalue = NF90_FILL_FLOAT
    end if
  end subroutine getVariableFillValueSp

  subroutine getVariableFillValueDp(self, fvalue)
    class(NcVariable), intent(in)  :: self
    real(dp)         , intent(out) :: fvalue

    if (self%hasAttribute("_FillValue")) then
       call self%getAttribute("_FillValue", fvalue)
    else
       fvalue = NF90_FILL_DOUBLE
    end if
  end subroutine getVariableFillValueDp

  subroutine setDataScalarI2(self, values, start)
    class(NcVariable), intent(in)           :: self
    integer(i2)      , intent(in)           :: values
    integer(i4)      , intent(in), optional :: start(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start), &
         "Failed to write data into variable: " // trim(self%name))
  end subroutine setDataScalarI2

  subroutine setData1dI2(self, values, start, count, stride, map)
    class(NcVariable), intent(in)           :: self
    integer(i2)      , intent(in)           :: values(:)
    integer(i4)      , intent(in), optional :: start(:), count(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, count, stride, map), &
         "Failed to write data into variable: " // trim(self%name))
  end subroutine setData1dI2

  subroutine setData2dI2(self, values, start, count, stride, map)
    class(NcVariable), intent(in)           :: self
    integer(i2)      , intent(in)           :: values(:,:)
    integer(i4)      , intent(in), optional :: start(:), count(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, count, stride, map), &
         "Failed to write data into variable: " // trim(self%name))
  end subroutine setData2dI2

  subroutine setData3dI2(self, values, start, count, stride, map)
    class(NcVariable), intent(in)           :: self
    integer(i2)      , intent(in)           :: values(:,:,:)
    integer(i4)      , intent(in), optional :: start(:), count(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, count, stride, map), &
         "Failed to write data into variable: " // trim(self%name))
  end subroutine setData3dI2

  subroutine setData4dI2(self, values, start, count, stride, map)
    class(NcVariable), intent(in)           :: self
    integer(i2)      , intent(in)           :: values(:,:,:,:)
    integer(i4)      , intent(in), optional :: start(:), count(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, count, stride, map), &
         "Failed to write data into variable: " // trim(self%name))
  end subroutine setData4dI2

  subroutine setData5dI2(self, values, start, count, stride, map)
    class(NcVariable), intent(in)           :: self
    integer(i2)      , intent(in)           :: values(:,:,:,:,:)
    integer(i4)      , intent(in), optional :: start(:), count(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, count, stride, map), &
         "Failed to write data into variable: " // trim(self%name))
  end subroutine setData5dI2

  subroutine setDataScalarI4(self, values, start)
    class(NcVariable), intent(in)           :: self
    integer(i4)      , intent(in)           :: values
    integer(i4)      , intent(in), optional :: start(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start), &
         "Failed to write data into variable: " // trim(self%name))
  end subroutine setDataScalarI4

  subroutine setData1dI4(self, values, start, count, stride, map)
    class(NcVariable), intent(in)           :: self
    integer(i4)      , intent(in)           :: values(:)
    integer(i4)      , intent(in), optional :: start(:), count(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, count, stride, map), &
         "Failed to write data into variable: " // trim(self%name))
  end subroutine setData1dI4

  subroutine setData2dI4(self, values, start, count, stride, map)
    class(NcVariable), intent(in)           :: self
    integer(i4)      , intent(in)           :: values(:,:)
    integer(i4)      , intent(in), optional :: start(:), count(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, count, stride, map), &
         "Failed to write data into variable: " // trim(self%name))
  end subroutine setData2dI4

  subroutine setData3dI4(self, values, start, count, stride, map)
    class(NcVariable), intent(in)           :: self
    integer(i4)      , intent(in)           :: values(:,:,:)
    integer(i4)      , intent(in), optional :: start(:), count(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, count, stride, map), &
         "Failed to write data into variable: " // trim(self%name))
  end subroutine setData3dI4

  subroutine setData4dI4(self, values, start, count, stride, map)
    class(NcVariable), intent(in)           :: self
    integer(i4)      , intent(in)           :: values(:,:,:,:)
    integer(i4)      , intent(in), optional :: start(:), count(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, count, stride, map), &
         "Failed to write data into variable: " // trim(self%name))
  end subroutine setData4dI4

  subroutine setData5dI4(self, values, start, count, stride, map)
    class(NcVariable), intent(in)           :: self
    integer(i4)      , intent(in)           :: values(:,:,:,:,:)
    integer(i4)      , intent(in), optional :: start(:), count(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, count, stride, map), &
         "Failed to write data into variable: " // trim(self%name))
  end subroutine setData5dI4

  subroutine setDataScalarSp(self, values, start)
    class(NcVariable), intent(in)           :: self
    real(sp)         , intent(in)           :: values
    integer(i4)      , intent(in), optional :: start(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start), &
         "Failed to write data into variable: " // trim(self%name))
  end subroutine setDataScalarSp

  subroutine setData1dSp(self, values, start, count, stride, map)
    class(NcVariable), intent(in)           :: self
    real(sp)         , intent(in)           :: values(:)
    integer(i4)      , intent(in), optional :: start(:), count(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, count, stride, map), &
         "Failed to write data into variable: " // trim(self%name))
  end subroutine setData1dSp

  subroutine setData2dSp(self, values, start, count, stride, map)
    class(NcVariable), intent(in)           :: self
    real(sp)         , intent(in)           :: values(:,:)
    integer(i4)      , intent(in), optional :: start(:), count(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, count, stride, map), &
         "Failed to write data into variable: " // trim(self%name))
  end subroutine setData2dSp

  subroutine setData3dSp(self, values, start, count, stride, map)
    class(NcVariable), intent(in)           :: self
    real(sp)         , intent(in)           :: values(:,:,:)
    integer(i4)      , intent(in), optional :: start(:), count(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, count, stride, map), &
         "Failed to write data into variable: " // trim(self%name))
  end subroutine setData3dSp

  subroutine setData4dSp(self, values, start, count, stride, map)
    class(NcVariable), intent(in)           :: self
    real(sp)         , intent(in)           :: values(:,:,:,:)
    integer(i4)      , intent(in), optional :: start(:), count(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, count, stride, map), &
         "Failed to write data into variable: " // trim(self%name))
  end subroutine setData4dSp

  subroutine setData5dSp(self, values, start, count, stride, map)
    class(NcVariable), intent(in)           :: self
    real(sp)         , intent(in)           :: values(:,:,:,:,:)
    integer(i4)      , intent(in), optional :: start(:), count(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, count, stride, map), &
         "Failed to write data into variable: " // trim(self%name))
  end subroutine setData5dSp

  subroutine setDataScalarDp(self, values, start)
    class(NcVariable), intent(in)           :: self
    real(dp)         , intent(in)           :: values
    integer(i4)      , intent(in), optional :: start(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start), &
         "Failed to write data into variable: " // trim(self%name))
  end subroutine setDataScalarDp

  subroutine setData1dDp(self, values, start, count, stride, map)
    class(NcVariable), intent(in)           :: self
    real(dp)         , intent(in)           :: values(:)
    integer(i4)      , intent(in), optional :: start(:), count(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, count, stride, map), &
         "Failed to write data into variable: " // trim(self%name))
  end subroutine setData1dDp

  subroutine setData2dDp(self, values, start, count, stride, map)
    class(NcVariable), intent(in)           :: self
    real(dp)         , intent(in)           :: values(:,:)
    integer(i4)      , intent(in), optional :: start(:), count(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, count, stride, map), &
         "Failed to write data into variable: " // trim(self%name))
  end subroutine setData2dDp

  subroutine setData3dDp(self, values, start, count, stride, map)
    class(NcVariable), intent(in)           :: self
    real(dp)         , intent(in)           :: values(:,:,:)
    integer(i4)      , intent(in), optional :: start(:), count(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, count, stride, map), &
         "Failed to write data into variable: " // trim(self%name))
  end subroutine setData3dDp

  subroutine setData4dDp(self, values, start, count, stride, map)
    class(NcVariable), intent(in)           :: self
    real(dp)         , intent(in)           :: values(:,:,:,:)
    integer(i4)      , intent(in), optional :: start(:), count(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, count, stride, map), &
         "Failed to write data into variable: " // trim(self%name))
  end subroutine setData4dDp

  subroutine setData5dDp(self, values, start, count, stride, map)
    class(NcVariable), intent(in)           :: self
    real(dp)         , intent(in)           :: values(:,:,:,:,:)
    integer(i4)      , intent(in), optional :: start(:), count(:), stride(:), map(:)

    call check( nf90_put_var(self%parent%id, self%id, values, start, count, stride, map), &
         "Failed to write data into variable: " // trim(self%name))
  end subroutine setData5dDp

  subroutine getDataScalarI2(self, data, start, count, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), count(:), stride(:), map(:)
    integer(i2)      , intent(out)              :: data
    integer(i2)                                 :: tmp(1)

    call check (nf90_get_var(self%parent%id, self%id, tmp, start, count, stride, map), &
         "Could not read data from variable: "//trim(self%name))
    data = tmp(1)
  end subroutine getDataScalarI2

  subroutine getData1dI2(self, data, start, count, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), count(:), stride(:), map(:)
    integer(i2)      , intent(out), allocatable :: data(:)
    integer(i4)                   , allocatable :: datashape(:)

    datashape = getReadDataShape(self,1,start,count,stride)
    allocate(data(datashape(1)) )

    call check (nf90_get_var(self%parent%id, self%id, data, start, count, stride, map), &
         "Could not read data from variable: "//trim(self%name))
  end subroutine getData1dI2

  subroutine getData2dI2(self, data, start, count, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), count(:), stride(:), map(:)
    integer(i2)      , intent(out), allocatable :: data(:,:)
    integer(i4)                   , allocatable :: datashape(:)

    datashape = getReadDataShape(self,2,start,count,stride)
    allocate(data(datashape(1),datashape(2)))

    call check (nf90_get_var(self%parent%id, self%id, data, start, count, stride, map), &
         "Could not read data from variable: "//trim(self%name))
  end subroutine getData2dI2

  subroutine getData3dI2(self, data, start, count, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), count(:), stride(:), map(:)
    integer(i2)      , intent(out), allocatable :: data(:,:,:)
    integer(i4)                   , allocatable :: datashape(:)

    datashape = getReadDataShape(self,3,start,count,stride)
    allocate(data(datashape(1),datashape(2),datashape(3)))

    call check (nf90_get_var(self%parent%id, self%id, data, start, count, stride, map), &
         "Could not read data from variable: "//trim(self%name))
  end subroutine getData3dI2

  subroutine getData4dI2(self, data, start, count, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), count(:), stride(:), map(:)
    integer(i2)      , intent(out), allocatable :: data(:,:,:,:)
    integer(i4)                   , allocatable :: datashape(:)

    datashape = getReadDataShape(self,3,start,count,stride)
    allocate(data(datashape(1),datashape(2),datashape(3),datashape(4)))

    call check (nf90_get_var(self%parent%id, self%id, data, start, count, stride, map), &
         "Could not read data from variable: "//trim(self%name))
  end subroutine getData4dI2

  subroutine getData5dI2(self, data, start, count, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), count(:), stride(:), map(:)
    integer(i2)      , intent(out), allocatable :: data(:,:,:,:,:)
    integer(i4)                   , allocatable :: datashape(:)

    datashape = getReadDataShape(self,3,start,count,stride)
    allocate(data(datashape(1),datashape(2),datashape(3),datashape(4),datashape(5)))

    call check (nf90_get_var(self%parent%id, self%id, data, start, count, stride, map), &
         "Could not read data from variable: "//trim(self%name))
  end subroutine getData5dI2

  subroutine getDataScalarI4(self, data, start, count, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), count(:), stride(:), map(:)
    integer(i4)      , intent(out)              :: data
    integer(i4)                                 :: tmp(1)

    call check (nf90_get_var(self%parent%id, self%id, tmp, start, count, stride, map), &
         "Could not read data from variable: "//trim(self%name))
    data = tmp(1)
  end subroutine getDataScalarI4

  subroutine getData1dI4(self, data, start, count, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), count(:), stride(:), map(:)
    integer(i4)      , intent(out), allocatable :: data(:)
    integer(i4)                   , allocatable :: datashape(:)

    datashape = getReadDataShape(self,1,start,count,stride)
    allocate(data(datashape(1)) )

    call check (nf90_get_var(self%parent%id, self%id, data, start, count, stride, map), &
         "Could not read data from variable: "//trim(self%name))
  end subroutine getData1dI4

  subroutine getData2dI4(self, data, start, count, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), count(:), stride(:), map(:)
    integer(i4)      , intent(out), allocatable :: data(:,:)
    integer(i4)                   , allocatable :: datashape(:)

    datashape = getReadDataShape(self,2,start,count,stride)
    allocate(data(datashape(1),datashape(2)))

    call check (nf90_get_var(self%parent%id, self%id, data, start, count, stride, map), &
         "Could not read data from variable: "//trim(self%name))
  end subroutine getData2dI4

  subroutine getData3dI4(self, data, start, count, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), count(:), stride(:), map(:)
    integer(i4)      , intent(out), allocatable :: data(:,:,:)
    integer(i4)                   , allocatable :: datashape(:)

    datashape = getReadDataShape(self,3,start,count,stride)
    allocate(data(datashape(1),datashape(2),datashape(3)))

    call check (nf90_get_var(self%parent%id, self%id, data, start, count, stride, map), &
         "Could not read data from variable: "//trim(self%name))
  end subroutine getData3dI4

  subroutine getData4dI4(self, data, start, count, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), count(:), stride(:), map(:)
    integer(i4)      , intent(out), allocatable :: data(:,:,:,:)
    integer(i4)                   , allocatable :: datashape(:)

    datashape = getReadDataShape(self,3,start,count,stride)
    allocate(data(datashape(1),datashape(2),datashape(3),datashape(4)))

    call check (nf90_get_var(self%parent%id, self%id, data, start, count, stride, map), &
         "Could not read data from variable: "//trim(self%name))
  end subroutine getData4dI4

  subroutine getData5dI4(self, data, start, count, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional   :: start(:), count(:), stride(:), map(:)
    integer(i4)      , intent(out), allocatable :: data(:,:,:,:,:)
    integer(i4)                   , allocatable :: datashape(:)

    datashape = getReadDataShape(self,3,start,count,stride)
    allocate(data(datashape(1),datashape(2),datashape(3),datashape(4),datashape(5)))

    call check (nf90_get_var(self%parent%id, self%id, data, start, count, stride, map), &
         "Could not read data from variable: "//trim(self%name))
  end subroutine getData5dI4

  subroutine getDataScalarSp(self, data, start, count, stride, map)
    class(NcVariable), intent(in)             :: self
    integer(i4)      , intent(in) , optional  :: start(:), count(:), stride(:), map(:)
    real(sp)         , intent(out)            :: data
    real(sp)                                  :: tmp(1)

    call check (nf90_get_var(self%parent%id, self%id, tmp, start, count, stride, map), &
         "Could not read data from variable: "//trim(self%name))
    data = tmp(1)
  end subroutine getDataScalarSp

  subroutine getData1dSp(self, data, start, count, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), count(:), stride(:), map(:)
    real(sp)         , intent(out), allocatable :: data(:)
    integer(i4)                   , allocatable :: datashape(:)

    datashape = getReadDataShape(self,1,start,count,stride)
    allocate(data(datashape(1)))

    call check (nf90_get_var(self%parent%id, self%id, data, start, count, stride, map), &
         "Could not read data from variable: "//trim(self%name))
  end subroutine getData1dSp

  subroutine getData2dSp(self, data, start, count, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), count(:), stride(:), map(:)
    real(sp)         , intent(out), allocatable :: data(:,:)
    integer(i4)                   , allocatable :: datashape(:)

    datashape = getReadDataShape(self,2,start,count,stride)
    allocate(data(datashape(1),datashape(2)))

    call check (nf90_get_var(self%parent%id, self%id, data, start, count, stride, map), &
         "Could not read data from variable: "//trim(self%name))
  end subroutine getData2dSp

  subroutine getData3dSp(self, data, start, count, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), count(:), stride(:), map(:)
    real(sp)         , intent(out), allocatable :: data(:,:,:)
    integer(i4)                   , allocatable :: datashape(:)

    datashape = getReadDataShape(self,3,start,count,stride)
    allocate(data(datashape(1),datashape(2),datashape(3)))

    call check (nf90_get_var(self%parent%id, self%id, data, start, count, stride, map), &
         "Could not read data from variable: "//trim(self%name))
  end subroutine getData3dSp

  subroutine getData4dSp(self, data, start, count, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), count(:), stride(:), map(:)
    real(sp)         , intent(out), allocatable :: data(:,:,:,:)
    integer(i4)                   , allocatable :: datashape(:)

    datashape = getReadDataShape(self,3,start,count,stride)
    allocate(data(datashape(1),datashape(2),datashape(3),datashape(4)))

    call check (nf90_get_var(self%parent%id, self%id, data, start, count, stride, map), &
         "Could not read data from variable: "//trim(self%name))
  end subroutine getData4dSp

  subroutine getData5dSp(self, data, start, count, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), count(:), stride(:), map(:)
    real(sp)         , intent(out), allocatable :: data(:,:,:,:,:)
    integer(i4)                   , allocatable :: datashape(:)

    datashape = getReadDataShape(self,3,start,count,stride)
    allocate(data(datashape(1),datashape(2),datashape(3),datashape(4),datashape(5)))

    call check (nf90_get_var(self%parent%id, self%id, data, start, count, stride, map), &
         "Could not read data from variable: "//trim(self%name))
  end subroutine getData5dSp

  subroutine getDataScalarDp(self, data, start, count, stride, map)
    class(NcVariable), intent(in)             :: self
    integer(i4)      , intent(in) , optional  :: start(:), count(:), stride(:), map(:)
    real(dp)         , intent(out)            :: data
    real(dp)                                  :: tmp(1)

    call check (nf90_get_var(self%parent%id, self%id, tmp, start, count, stride, map), &
         "Could not read data from variable: "//trim(self%name))
    data = tmp(1)
  end subroutine getDataScalarDp

  subroutine getData1dDp(self, data, start, count, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), count(:), stride(:), map(:)
    real(dp)         , intent(out), allocatable :: data(:)
    integer(i4)                   , allocatable :: datashape(:)

    datashape = getReadDataShape(self,1,start,count,stride)
    allocate(data(datashape(1)))

    call check (nf90_get_var(self%parent%id, self%id, data, start, count, stride, map), &
         "Could not read data from variable: "//trim(self%name))
  end subroutine getData1dDp

  subroutine getData2dDp(self, data, start, count, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), count(:), stride(:), map(:)
    real(dp)         , intent(out), allocatable :: data(:,:)
    integer(i4)                   , allocatable :: datashape(:)

    datashape = getReadDataShape(self,2,start,count,stride)
    allocate(data(datashape(1),datashape(2)))

    call check (nf90_get_var(self%parent%id, self%id, data, start, count, stride, map), &
         "Could not read data from variable: "//trim(self%name))
  end subroutine getData2dDp

  subroutine getData3dDp(self, data, start, count, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), count(:), stride(:), map(:)
    real(dp)         , intent(out), allocatable :: data(:,:,:)
    integer(i4)                   , allocatable :: datashape(:)

    datashape = getReadDataShape(self,3,start,count,stride)
    allocate(data(datashape(1),datashape(2),datashape(3)))

    call check (nf90_get_var(self%parent%id, self%id, data, start, count, stride, map), &
         "Could not read data from variable: "//trim(self%name))
  end subroutine getData3dDp

  subroutine getData4dDp(self, data, start, count, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), count(:), stride(:), map(:)
    real(dp)         , intent(out), allocatable :: data(:,:,:,:)
    integer(i4)                   , allocatable :: datashape(:)

    datashape = getReadDataShape(self,3,start,count,stride)
    allocate(data(datashape(1),datashape(2),datashape(3),datashape(4)))

    call check (nf90_get_var(self%parent%id, self%id, data, start, count, stride, map), &
         "Could not read data from variable: "//trim(self%name))
  end subroutine getData4dDp

  subroutine getData5dDp(self, data, start, count, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), count(:), stride(:), map(:)
    real(dp)         , intent(out), allocatable :: data(:,:,:,:,:)
    integer(i4)                   , allocatable :: datashape(:)

    datashape = getReadDataShape(self,3,start,count,stride)
    allocate(data(datashape(1),datashape(2),datashape(3),datashape(4),datashape(5)))

    call check (nf90_get_var(self%parent%id, self%id, data, start, count, stride, map), &
         "Could not read data from variable: "//trim(self%name))
  end subroutine getData5dDp

  function getReadDataShape(var, datarank, instart, incount, instride)
    type(NcVariable), intent(in)           :: var
    integer(i4)     , intent(in)           :: datarank
    integer(i4)     , intent(in), optional :: instart(:), incount(:), instride(:)
    integer(i4)                            :: getReadDataShape(datarank)
    integer(i4)     , allocatable          :: datashape(:)

    datashape = var%getShape()
    if (present(incount)) then
       datashape = incount
    else
       if (present(instart)) then
          datashape(:size(instart)) = datashape(:size(instart)) - (instart - 1)
       end if
       if (present(instride)) then
          datashape(:size(instride)) = datashape(:size(instride)) / instride
       end if
    end if

    if (count(datashape .gt. 1) .ne. datarank) then
       write(*,*) "Given read parameters do not match output variable rank!"
       stop 1
    end if

    getReadDataShape = pack(datashape, datashape .gt. 1)
  end function getReadDataShape

  function getDtypeFromString(dtype)
    integer(i4)          :: getDtypeFromString
    character(*)         :: dtype

    select case(dtype)
    case("f32")
       getDtypeFromString = NF90_FLOAT
    case("f64")
       getDtypeFromString = NF90_DOUBLE
    case("i16")
       getDtypeFromString = NF90_SHORT
    case("i32")
       getDtypeFromString = NF90_INT
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
    case(NF90_SHORT)
       getDtypeFromInteger = "i16"
    case(NF90_INT)
       getDtypeFromInteger = "i32"
    case default
       write(*,*) "Datatype not understood: ", dtype
       stop 1
    end select
  end function getDtypeFromInteger

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

