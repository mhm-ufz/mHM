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

  ! Written  David Schaefer, Jun 2015
  ! Modified Matthias Cuntz, Jan 2016 - compiled with PGI Fortran rev 15.9 - no automatic allocation of left-hand-side
  
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

  use mo_kind,         only: i1, i2, i4, sp, dp
  use netcdf,          only: &
       nf90_open, nf90_close, nf90_strerror, nf90_def_dim, nf90_def_var,   &
       nf90_put_var, nf90_get_var, nf90_put_att, nf90_get_att,             &
       nf90_inquire, nf90_inq_dimid, nf90_inquire_dimension,               &
       nf90_inq_varid, nf90_inq_varids, nf90_inquire_variable, nf90_inquire_attribute,      &
       NF90_OPEN, NF90_NETCDF4, NF90_CREATE, NF90_WRITE, NF90_NOWRITE,     &
       NF90_BYTE, NF90_SHORT, NF90_INT, NF90_FLOAT, NF90_DOUBLE,                      &
       NF90_FILL_BYTE, NF90_FILL_SHORT, NF90_FILL_INT, NF90_FILL_FLOAT , NF90_FILL_DOUBLE, &
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

     procedure, public  :: initNcDataset

     procedure, public  :: getNoVariables
     procedure, public  :: getVariableIds
     procedure, public  :: getVariables
     
     procedure, private :: setGlobalAttributeChar
     procedure, private :: setGlobalAttributeI8
     procedure, private :: setGlobalAttributeI16
     procedure, private :: setGlobalAttributeI32
     procedure, private :: setGlobalAttributeF32
     procedure, private :: setGlobalAttributeF64

     procedure, private :: getGlobalAttributeChar
     procedure, private :: getGlobalAttributeI8
     procedure, private :: getGlobalAttributeI16
     procedure, private :: getGlobalAttributeI32
     procedure, private :: getGlobalAttributeF32
     procedure, private :: getGlobalAttributeF64

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
          setGlobalAttributeI8,   &
          setGlobalAttributeI16,   &
          setGlobalAttributeI32,   &
          setGlobalAttributeF32,   &
          setGlobalAttributeF64

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
          getGlobalAttributeI8,   &
          getGlobalAttributeI16,   &
          getGlobalAttributeI32,   &
          getGlobalAttributeF32,   &
          getGlobalAttributeF64

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
     !>                     NF90_BYTE   -> "i8"
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

   contains

     procedure, public :: initNcDimension
     
     procedure, public :: getName => getDimensionName
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

     integer(i4)     :: id       !> Variable id
     type(NcDataset) :: parent   !> The variables's parent
     ! character(256)  :: name     !> Variable name

   contains

     procedure, public  :: initNcVariable

     procedure, public  :: getName => getVariableName
     
     procedure, private :: setVariableAttributeChar
     procedure, private :: setVariableAttributeI8
     procedure, private :: setVariableAttributeI16
     procedure, private :: setVariableAttributeI32
     procedure, private :: setVariableAttributeF32
     procedure, private :: setVariableAttributeF64

     procedure, private :: getVariableAttributeChar
     procedure, private :: getVariableAttributeI8
     procedure, private :: getVariableAttributeI16
     procedure, private :: getVariableAttributeI32
     procedure, private :: getVariableAttributeF32
     procedure, private :: getVariableAttributeF64

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
     procedure, private :: setVariableFillValueF32
     procedure, private :: setVariableFillValueF64

     procedure, private :: getVariableFillValueI8
     procedure, private :: getVariableFillValueI16
     procedure, private :: getVariableFillValueI32
     procedure, private :: getVariableFillValueF32
     procedure, private :: getVariableFillValueF64

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
     !>                     "i8"  -> NF90_BYTE
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
     !>        \param[in]  "integer(i4) :: start(:), cnt(:), stride(:), map(:)"
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
     !>        \param[in] "integer(i4) :: start(:), cnt(:), stride(:), map(:)"
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
          setVariableFillValueI8, &
          setVariableFillValueI16, &
          setVariableFillValueI32, &
          setVariableFillValueF32, &
          setVariableFillValueF64

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
          getVariableFillValueI8, &
          getVariableFillValueI16, &
          getVariableFillValueI32, &
          getVariableFillValueF32, &
          getVariableFillValueF64

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
          setVariableAttributeI8,   &
          setVariableAttributeI16,   &
          setVariableAttributeI32,   &
          setVariableAttributeF32,   &
          setVariableAttributeF64

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
          getVariableAttributeI8,   &
          getVariableAttributeI16,   &
          getVariableAttributeI32,   &
          getVariableAttributeF32,   &
          getVariableAttributeF64

  end type NcVariable

  interface NcVariable
     procedure newNcVariable
  end interface NcVariable

contains

  subroutine initNcVariable(self, id, parent)
    class(NcVariable), intent(inout) :: self
    integer(i4)      , intent(in)    :: id
    type(NcDataset)  , intent(in)    :: parent

    self%id     = id
    self%parent = parent
  end subroutine initNcVariable

  subroutine initNcDimension(self, id, parent)
    class(NcDimension), intent(inout) :: self
    integer(i4)       , intent(in)    :: id
    type(NcDataset)   , intent(in)    :: parent

    self%id     = id
    self%parent = parent
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

  type(NcVariable) function newNcVariable(id, parent)
    integer(i4)    , intent(in) :: id
    type(NcDataset), intent(in) :: parent

    call newNcVariable%initNcVariable(id, parent)
  end function newNcVariable

  type(NcDimension) function newNcDimension(id, parent)
    integer(i4)    , intent(in) :: id
    type(NcDataset), intent(in) :: parent

    call newNcDimension%initNcDimension(id, parent)
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

  function getNoVariables(self)
    class(NcDataset), intent(in) :: self
    integer(i4)                  :: getNoVariables

    call check(nf90_inquire(self%id, nvariables=getNoVariables), "Failed inquire number of variables")
  end function getNoVariables

  function getVariableIds(self)
    class(NcDataset), intent(in)           :: self
    integer(i4), dimension(:), allocatable :: getVariableIds
    integer(i4)                            :: tmp
    
    allocate(getVariableIds(self%getNoVariables()))
    call check(nf90_inq_varids(self%id, tmp, getVariableIds), "Failed to inquire variable ids")
  end function getVariableIds
  
  function getVariables(self)
    class(NcDataset), intent(in)                :: self
    type(NcVariable), dimension(:), allocatable :: getVariables
    integer(i4), dimension(:), allocatable      :: varids
    integer(i4)                                 :: i, nvars

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

    setDimension = NcDimension(id,self)
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
    setVariableWithIds = NcVariable(varid, self)
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
    getDimensionById = NcDimension(id,self)
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
    integer(i4)      , allocatable :: getVariableShape(:)
    type(NcDimension), allocatable :: dims(:)
    integer(i4)                    :: ii, ndims

    ndims = self%getNoDimensions()
    allocate(getVariableShape(ndims), dims(ndims))

    dims = self%getDimensions()
    do ii = 1,size(dims)
       getVariableShape(ii) = dims(ii)%getLength()
    end do
  end function getVariableShape

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

  subroutine setGlobalAttributeI8(self, name, data)
    class(NcDataset), intent(in) :: self
    character(*)    , intent(in) :: name
    integer(i1)     , intent(in) :: data

    call check(nf90_put_att(self%id,NF90_GLOBAL,name,data), &
         "Failed to write attribute: " // name )
  end subroutine setGlobalAttributeI8

  subroutine setGlobalAttributeI16(self, name, data)
    class(NcDataset), intent(in) :: self
    character(*)    , intent(in) :: name
    integer(i2)     , intent(in) :: data

    call check(nf90_put_att(self%id,NF90_GLOBAL,name,data), &
         "Failed to write attribute: " // name )
  end subroutine setGlobalAttributeI16

  subroutine setGlobalAttributeI32(self, name, data)
    class(NcDataset), intent(in) :: self
    character(*)    , intent(in) :: name
    integer(i4)     , intent(in) :: data

    call check(nf90_put_att(self%id,NF90_GLOBAL,name,data), &
         "Failed to write attribute: " // name )
  end subroutine setGlobalAttributeI32

  subroutine setGlobalAttributeF32(self, name, data)
    class(NcDataset), intent(in) :: self
    character(*)    , intent(in) :: name
    real(sp)        , intent(in) :: data

    call check(nf90_put_att(self%id,NF90_GLOBAL,name,data), &
         "Failed to write attribute: " // name )
  end subroutine setGlobalAttributeF32

  subroutine setGlobalAttributeF64(self, name, data)
    class(NcDataset), intent(in) :: self
    character(*)    , intent(in) :: name
    real(dp)        , intent(in) :: data

    call check(nf90_put_att(self%id,NF90_GLOBAL,name,data), &
         "Failed to write attribute: " // name )
  end subroutine setGlobalAttributeF64

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

  subroutine getGlobalAttributeI8(self, name, avalue)
    class(NcDataset), intent(in)  :: self
    character(*)    , intent(in)  :: name
    integer(i1)     , intent(out) :: avalue
    integer(i4)                   :: length

    call check(nf90_inquire_attribute(self%id,NF90_GLOBAL,name,len=length),&
         "Could not inquire attribute "//name)
    call check(nf90_get_att(self%id,NF90_GLOBAL,name,avalue), &
         "Could not read attribute "//name)
  end subroutine getGlobalAttributeI8

  subroutine getGlobalAttributeI16(self, name, avalue)
    class(NcDataset), intent(in)  :: self
    character(*)    , intent(in)  :: name
    integer(i2)     , intent(out) :: avalue
    integer(i4)                   :: length

    call check(nf90_inquire_attribute(self%id,NF90_GLOBAL,name,len=length),&
         "Could not inquire attribute "//name)
    call check(nf90_get_att(self%id,NF90_GLOBAL,name,avalue), &
         "Could not read attribute "//name)
  end subroutine getGlobalAttributeI16

  subroutine getGlobalAttributeI32(self, name, avalue)
    class(NcDataset), intent(in)  :: self
    character(*)    , intent(in)  :: name
    integer(i4)     , intent(out) :: avalue
    integer(i4)                   :: length

    call check(nf90_inquire_attribute(self%id,NF90_GLOBAL,name,len=length),&
         "Could not inquire attribute "//name)
    call check(nf90_get_att(self%id,NF90_GLOBAL,name,avalue), &
         "Could not read attribute "//name)
  end subroutine getGlobalAttributeI32

  subroutine getGlobalAttributeF32(self, name, avalue)
    class(NcDataset), intent(in)  :: self
    character(*)    , intent(in)  :: name
    real(sp)        , intent(out) :: avalue
    integer(i4)                   :: length

    call check(nf90_inquire_attribute(self%id,NF90_GLOBAL,name,len=length),&
         "Could not inquire attribute "//name)
    call check(nf90_get_att(self%id,NF90_GLOBAL,name,avalue), &
         "Could not read attribute "//name)
  end subroutine getGlobalAttributeF32

  subroutine getGlobalAttributeF64(self, name, avalue)
    class(NcDataset), intent(in)  :: self
    character(*)    , intent(in)  :: name
    real(dp)        , intent(out) :: avalue
    integer(i4)                   :: length

    call check(nf90_inquire_attribute(self%id,NF90_GLOBAL,name,len=length),&
         "Could not inquire attribute "//name)
    call check(nf90_get_att(self%id,NF90_GLOBAL,name,avalue), &
         "Could not read attribute "//name)
  end subroutine getGlobalAttributeF64

  subroutine setVariableAttributeChar(self, name, data)
    class(NcVariable), intent(in) :: self
    character(*)     , intent(in) :: name
    character(*)     , intent(in) :: data

    call check(nf90_put_att(self%parent%id,self%id,name,data), &
         "Failed to write attribute: " // name )
  end subroutine setVariableAttributeChar

  subroutine setVariableAttributeI8(self, name, data)
    class(NcVariable), intent(in) :: self
    character(*)     , intent(in) :: name
    integer(i1)      , intent(in) :: data

    call check(nf90_put_att(self%parent%id,self%id,name,data), &
         "Failed to write attribute: " // name )
  end subroutine setVariableAttributeI8
  
  subroutine setVariableAttributeI16(self, name, data)
    class(NcVariable), intent(in) :: self
    character(*)     , intent(in) :: name
    integer(i2)      , intent(in) :: data

    call check(nf90_put_att(self%parent%id,self%id,name,data), &
         "Failed to write attribute: " // name )
  end subroutine setVariableAttributeI16

  subroutine setVariableAttributeI32(self, name, data)
    class(NcVariable), intent(in) :: self
    character(*)     , intent(in) :: name
    integer(i4)      , intent(in) :: data

    call check(nf90_put_att(self%parent%id,self%id,name,data), &
         "Failed to write attribute: " // name )
  end subroutine setVariableAttributeI32

  subroutine setVariableAttributeF32(self, name, data)
    class(NcVariable), intent(in) :: self
    character(*)     , intent(in) :: name
    real(sp)         , intent(in) :: data

    call check(nf90_put_att(self%parent%id,self%id,name,data), &
         "Failed to write attribute: " // name )
  end subroutine setVariableAttributeF32

  subroutine setVariableAttributeF64(self, name, data)
    class(NcVariable), intent(in) :: self
    character(*)     , intent(in) :: name
    real(dp)         , intent(in) :: data

    call check(nf90_put_att(self%parent%id,self%id,name,data), &
         "Failed to write attribute: " // name )
  end subroutine setVariableAttributeF64

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

  subroutine getVariableAttributeI8(self, name, avalue)
    class(NcVariable), intent(in)  :: self
    character(*)     , intent(in)  :: name
    integer(i1)      , intent(out) :: avalue
    integer(i4)                    :: length

    call check(nf90_inquire_attribute(self%parent%id,self%id,name,len=length),&
         "Could not inquire attribute "//name)
    call check(nf90_get_att(self%parent%id,self%id,name,avalue), &
         "Could not read attribute "//name)
  end subroutine getVariableAttributeI8

  subroutine getVariableAttributeI16(self, name, avalue)
    class(NcVariable), intent(in)  :: self
    character(*)     , intent(in)  :: name
    integer(i2)      , intent(out) :: avalue
    integer(i4)                    :: length

    call check(nf90_inquire_attribute(self%parent%id,self%id,name,len=length),&
         "Could not inquire attribute "//name)
    call check(nf90_get_att(self%parent%id,self%id,name,avalue), &
         "Could not read attribute "//name)
  end subroutine getVariableAttributeI16

  subroutine getVariableAttributeI32(self, name, avalue)
    class(NcVariable), intent(in)  :: self
    character(*)     , intent(in)  :: name
    integer(i4)      , intent(out) :: avalue
    integer(i4)                    :: length

    call check(nf90_inquire_attribute(self%parent%id,self%id,name,len=length),&
         "Could not inquire attribute "//name)
    call check(nf90_get_att(self%parent%id,self%id,name,avalue), &
         "Could not read attribute "//name)
  end subroutine getVariableAttributeI32

  subroutine getVariableAttributeF32(self, name, avalue)
    class(NcVariable), intent(in)  :: self
    character(*)     , intent(in)  :: name
    real(sp)         , intent(out) :: avalue
    integer(i4)                    :: length

    call check(nf90_inquire_attribute(self%parent%id,self%id,name,len=length),&
         "Could not inquire attribute "//name)
    call check(nf90_get_att(self%parent%id,self%id,name,avalue), &
         "Could not read attribute "//name)
  end subroutine getVariableAttributeF32

  subroutine getVariableAttributeF64(self, name, avalue)
    class(NcVariable), intent(in)  :: self
    character(*)     , intent(in)  :: name
    real(dp)         , intent(out) :: avalue
    integer(i4)                    :: length

    call check(nf90_inquire_attribute(self%parent%id,self%id,name,len=length),&
         "Could not inquire attribute "//name)
    call check(nf90_get_att(self%parent%id,self%id,name,avalue), &
         "Could not read attribute "//name)
  end subroutine getVariableAttributeF64

  subroutine setVariableFillValueI8(self, fvalue)
    class(NcVariable), intent(in)  :: self
    integer(i1)      , intent(in)  :: fvalue

    if (.not. self%hasAttribute("_FillValue")) then
       call self%setAttribute("_FillValue",fvalue)
    end if
  end subroutine setVariableFillValueI8
  
  subroutine setVariableFillValueI16(self, fvalue)
    class(NcVariable), intent(in)  :: self
    integer(i2)      , intent(in)  :: fvalue

    if (.not. self%hasAttribute("_FillValue")) then
       call self%setAttribute("_FillValue",fvalue)
    end if
  end subroutine setVariableFillValueI16

  subroutine setVariableFillValueI32(self, fvalue)
    class(NcVariable), intent(in)  :: self
    integer(i4)      , intent(in)  :: fvalue

    if (.not. self%hasAttribute("_FillValue")) then
       call self%setAttribute("_FillValue",fvalue)
    end if
  end subroutine setVariableFillValueI32

  subroutine setVariableFillValueF32(self, fvalue)
    class(NcVariable), intent(in)  :: self
    real(sp)         , intent(in)  :: fvalue

    if (.not. self%hasAttribute("_FillValue")) then
       call self%setAttribute("_FillValue",fvalue)
    end if
  end subroutine setVariableFillValueF32

  subroutine setVariableFillValueF64(self, fvalue)
    class(NcVariable), intent(in)  :: self
    real(dp)         , intent(in)  :: fvalue

    if (.not. self%hasAttribute("_FillValue")) then
       call self%setAttribute("_FillValue",fvalue)
    end if
  end subroutine setVariableFillValueF64

  subroutine getVariableFillValueI8(self, fvalue)
    class(NcVariable), intent(in)  :: self
    integer(i1)      , intent(out) :: fvalue

    if (self%hasAttribute("_FillValue")) then
       call self%getAttribute("_FillValue", fvalue)
    else
       fvalue = NF90_FILL_BYTE
    end if
  end subroutine getVariableFillValueI8

  subroutine getVariableFillValueI16(self, fvalue)
    class(NcVariable), intent(in)  :: self
    integer(i2)      , intent(out) :: fvalue

    if (self%hasAttribute("_FillValue")) then
       call self%getAttribute("_FillValue", fvalue)
    else
       fvalue = NF90_FILL_SHORT
    end if
  end subroutine getVariableFillValueI16

  subroutine getVariableFillValueI32(self, fvalue)
    class(NcVariable), intent(in)  :: self
    integer(i4)      , intent(out) :: fvalue

    if (self%hasAttribute("_FillValue")) then
       call self%getAttribute("_FillValue", fvalue)
    else
       fvalue = NF90_FILL_INT
    end if
  end subroutine getVariableFillValueI32

  subroutine getVariableFillValueF32(self, fvalue)
    class(NcVariable), intent(in)  :: self
    real(sp)         , intent(out) :: fvalue

    if (self%hasAttribute("_FillValue")) then
       call self%getAttribute("_FillValue", fvalue)
    else
       fvalue = NF90_FILL_FLOAT
    end if
  end subroutine getVariableFillValueF32

  subroutine getVariableFillValueF64(self, fvalue)
    class(NcVariable), intent(in)  :: self
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
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i1)      , intent(out)              :: data
    integer(i1)                                 :: tmp(1)

    call check (nf90_get_var(self%parent%id, self%id, tmp, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
    data = tmp(1)
  end subroutine getDataScalarI8

  subroutine getData1dI8(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i1)      , intent(out), allocatable :: data(:)
    integer(i4)                   , allocatable :: datashape(:)

    allocate(datashape(1))
    datashape = getReadDataShape(self, 1, start, cnt, stride)

    allocate(data(datashape(1)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData1dI8

  subroutine getData2dI8(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i1)      , intent(out), allocatable :: data(:,:)
    integer(i4)                   , allocatable :: datashape(:)

    allocate(datashape(2))
    datashape = getReadDataShape(self, 2, start, cnt, stride)

    allocate(data(datashape(1), datashape(2)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData2dI8

  subroutine getData3dI8(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i1)      , intent(out), allocatable :: data(:,:,:)
    integer(i4)                   , allocatable :: datashape(:)

    allocate(datashape(3))
    datashape = getReadDataShape(self, 3, start, cnt, stride)

    allocate(data(datashape(1), datashape(2), datashape(3)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData3dI8

  subroutine getData4dI8(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i1)      , intent(out), allocatable :: data(:,:,:,:)
    integer(i4)                   , allocatable :: datashape(:)

    allocate(datashape(4))
    datashape = getReadDataShape(self, 4, start, cnt, stride)

    allocate(data(datashape(1), datashape(2), datashape(3), datashape(4)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData4dI8

  subroutine getData5dI8(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i1)      , intent(out), allocatable :: data(:,:,:,:,:)
    integer(i4)                   , allocatable :: datashape(:)

    allocate(datashape(5))
    datashape = getReadDataShape(self, 5, start, cnt, stride)

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
    integer(i4)      , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i2)      , intent(out), allocatable :: data(:)
    integer(i4)                   , allocatable :: datashape(:)

    allocate(datashape(1))
    datashape = getReadDataShape(self, 1, start, cnt, stride)

    allocate(data(datashape(1)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData1dI16

  subroutine getData2dI16(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i2)      , intent(out), allocatable :: data(:,:)
    integer(i4)                   , allocatable :: datashape(:)

    allocate(datashape(2))
    datashape = getReadDataShape(self, 2, start, cnt, stride)

    allocate(data(datashape(1), datashape(2)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData2dI16

  subroutine getData3dI16(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i2)      , intent(out), allocatable :: data(:,:,:)
    integer(i4)                   , allocatable :: datashape(:)

    allocate(datashape(3))
    datashape = getReadDataShape(self, 3, start, cnt, stride)

    allocate(data(datashape(1), datashape(2), datashape(3)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData3dI16

  subroutine getData4dI16(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i2)      , intent(out), allocatable :: data(:,:,:,:)
    integer(i4)                   , allocatable :: datashape(:)

    allocate(datashape(4))
    datashape = getReadDataShape(self, 4, start, cnt, stride)

    allocate(data(datashape(1), datashape(2), datashape(3), datashape(4)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData4dI16

  subroutine getData5dI16(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i2)      , intent(out), allocatable :: data(:,:,:,:,:)
    integer(i4)                   , allocatable :: datashape(:)

    allocate(datashape(5))
    datashape = getReadDataShape(self, 5, start, cnt, stride)

    allocate(data(datashape(1), datashape(2), datashape(3), datashape(4), datashape(5)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData5dI16

  subroutine getDataScalarI32(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i4)      , intent(out)              :: data
    integer(i4)                                 :: tmp(1)

    call check (nf90_get_var(self%parent%id, self%id, tmp, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
    data = tmp(1)
  end subroutine getDataScalarI32

  subroutine getData1dI32(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i4)      , intent(out), allocatable :: data(:)
    integer(i4)                   , allocatable :: datashape(:)
    
    allocate(datashape(1))
    datashape = getReadDataShape(self, 1, start, cnt, stride)
    
    allocate(data(datashape(1)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData1dI32

  subroutine getData2dI32(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i4)      , intent(out), allocatable :: data(:,:)
    integer(i4)                   , allocatable :: datashape(:)

    allocate(datashape(2))
    datashape = getReadDataShape(self, 2, start, cnt, stride)

    allocate(data(datashape(1), datashape(2)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData2dI32

  subroutine getData3dI32(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i4)      , intent(out), allocatable :: data(:,:,:)
    integer(i4)                   , allocatable :: datashape(:)

    allocate(datashape(3))
    datashape = getReadDataShape(self, 3, start, cnt, stride)

    allocate(data(datashape(1), datashape(2), datashape(3)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData3dI32

  subroutine getData4dI32(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    integer(i4)      , intent(out), allocatable :: data(:,:,:,:)
    integer(i4)                   , allocatable :: datashape(:)

    allocate(datashape(4))
    datashape = getReadDataShape(self, 4, start, cnt, stride)

    allocate(data(datashape(1), datashape(2), datashape(3), datashape(4)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData4dI32

  subroutine getData5dI32(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional   :: start(:), cnt(:), stride(:), map(:)
    integer(i4)      , intent(out), allocatable :: data(:,:,:,:,:)
    integer(i4)                   , allocatable :: datashape(:)

    allocate(datashape(5))
    datashape = getReadDataShape(self, 5, start, cnt, stride)

    allocate(data(datashape(1), datashape(2), datashape(3), datashape(4), datashape(5)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData5dI32

  subroutine getDataScalarF32(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)             :: self
    integer(i4)      , intent(in) , optional  :: start(:), cnt(:), stride(:), map(:)
    real(sp)         , intent(out)            :: data
    real(sp)                                  :: tmp(1)

    call check (nf90_get_var(self%parent%id, self%id, tmp, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
    data = tmp(1)
  end subroutine getDataScalarF32

  subroutine getData1dF32(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    real(sp)         , intent(out), allocatable :: data(:)
    integer(i4)                   , allocatable :: datashape(:)

    allocate(datashape(1))
    datashape = getReadDataShape(self, 1, start, cnt, stride)

    allocate(data(datashape(1)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData1dF32

  subroutine getData2dF32(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    real(sp)         , intent(out), allocatable :: data(:,:)
    integer(i4)                   , allocatable :: datashape(:)

    allocate(datashape(2))
    datashape = getReadDataShape(self, 2, start, cnt, stride)

    allocate(data(datashape(1), datashape(2)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData2dF32

  subroutine getData3dF32(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    real(sp)         , intent(out), allocatable :: data(:,:,:)
    integer(i4)                   , allocatable :: datashape(:)

    allocate(datashape(3))
    datashape = getReadDataShape(self, 3, start, cnt, stride)

    allocate(data(datashape(1), datashape(2), datashape(3)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData3dF32

  subroutine getData4dF32(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    real(sp)         , intent(out), allocatable :: data(:,:,:,:)
    integer(i4)                   , allocatable :: datashape(:)

    allocate(datashape(4))
    datashape = getReadDataShape(self, 4, start, cnt, stride)

    allocate(data(datashape(1), datashape(2), datashape(3), datashape(4)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData4dF32

  subroutine getData5dF32(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    real(sp)         , intent(out), allocatable :: data(:,:,:,:,:)
    integer(i4)                   , allocatable :: datashape(:)

    allocate(datashape(5))
    datashape = getReadDataShape(self, 5, start, cnt, stride)

    allocate(data(datashape(1), datashape(2), datashape(3), datashape(4), datashape(5)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData5dF32

  subroutine getDataScalarF64(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)             :: self
    integer(i4)      , intent(in) , optional  :: start(:), cnt(:), stride(:), map(:)
    real(dp)         , intent(out)            :: data
    real(dp)                                  :: tmp(1)

    call check (nf90_get_var(self%parent%id, self%id, tmp, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
    data = tmp(1)
  end subroutine getDataScalarF64

  subroutine getData1dF64(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    real(dp)         , intent(out), allocatable :: data(:)
    integer(i4)                   , allocatable :: datashape(:)

    allocate(datashape(1))
    datashape = getReadDataShape(self, 1, start, cnt, stride)
    
    allocate(data(datashape(1)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData1dF64

  subroutine getData2dF64(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    real(dp)         , intent(out), allocatable :: data(:,:)
    integer(i4)                   , allocatable :: datashape(:)

    allocate(datashape(2))
    datashape = getReadDataShape(self, 2, start, cnt, stride)

    allocate(data(datashape(1), datashape(2)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData2dF64

  subroutine getData3dF64(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    real(dp)         , intent(out), allocatable :: data(:,:,:)
    integer(i4)                   , allocatable :: datashape(:)

    allocate(datashape(3))
    datashape = getReadDataShape(self, 3, start, cnt, stride)

    allocate(data(datashape(1), datashape(2), datashape(3)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData3dF64

  subroutine getData4dF64(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    real(dp)         , intent(out), allocatable :: data(:,:,:,:)
    integer(i4)                   , allocatable :: datashape(:)

    allocate(datashape(4))
    datashape = getReadDataShape(self, 4, start, cnt, stride)

    allocate(data(datashape(1), datashape(2), datashape(3), datashape(4)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData4dF64

  subroutine getData5dF64(self, data, start, cnt, stride, map)
    class(NcVariable), intent(in)               :: self
    integer(i4)      , intent(in) , optional    :: start(:), cnt(:), stride(:), map(:)
    real(dp)         , intent(out), allocatable :: data(:,:,:,:,:)
    integer(i4)                   , allocatable :: datashape(:)

    allocate(datashape(5))
    datashape = getReadDataShape(self, 5, start, cnt, stride)

    allocate(data(datashape(1), datashape(2), datashape(3), datashape(4), datashape(5)))
    call check (nf90_get_var(self%parent%id, self%id, data, start, cnt, stride, map), &
         "Could not read data from variable: "//trim(self%getName()))
  end subroutine getData5dF64

  function getReadDataShape(var, datarank, instart, incnt, instride)
    type(NcVariable), intent(in)           :: var
    integer(i4)     , intent(in)           :: datarank
    integer(i4)     , intent(in), optional :: instart(:), incnt(:), instride(:)
    integer(i4)     , allocatable          :: readshape(:)
    integer(i4)                            :: naxis
    integer(i4)                            :: getReadDataShape(datarank)

    readshape = var%getShape()
    
    if (present(incnt)) then
       readshape = incnt
    else
       if (present(instart)) then
          readshape(:size(instart)) = readshape(:size(instart)) - (instart - 1)
       end if
       if (present(instride)) then
          readshape(:size(instride)) = readshape(:size(instride)) / instride
       end if
    end if
    
    naxis = count(readshape .gt. 1)
    
    if (all(readshape .eq. 1)) then
       ! return 1-element array
       getReadDataShape(:) = 1 !readshape(1:datarank+1)
    else if (size(readshape) .eq. datarank) then
       ! sizes fit
       getReadDataShape = readshape
    else if (naxis .eq. datarank) then
       getReadDataShape = pack(readshape, readshape .gt. 1)
    ! else if (naxis .lt. datarank) then
       ! would be nice...
    else
       write(*,*) "Given data reading parameters do not match output variable rank!"
       stop 1
    end if
       
  end function getReadDataShape

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

