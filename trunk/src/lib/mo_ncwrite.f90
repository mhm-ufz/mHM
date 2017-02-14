module mo_ncwrite

  ! This module provides a structure and subroutines for writing netcdf files.

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

  ! Copyright 2011-2012 Luis Samaniego, Stephan Thober, Matthias Cuntz

  use mo_kind,         only: i1, i4, sp, dp
  use mo_string_utils, only: nonull
  use mo_utils,        only: ne

  ! functions and constants of netcdf4 library
  use netcdf,  only: nf90_create, nf90_def_dim, NF90_UNLIMITED, nf90_def_var, &
       NF90_CHAR, nf90_put_att, NF90_INT, NF90_INT, NF90_GLOBAL, &
       nf90_enddef, nf90_put_var, NF90_FLOAT, NF90_DOUBLE, NF90_BYTE, &
       NF90_close, nf90_noerr, nf90_strerror, NF90_CLOBBER, &
       NF90_MAX_NAME, NF90_WRITE, nf90_inq_varid, nf90_inquire_variable, &
       nf90_inquire_dimension, nf90_open, NF90_64BIT_OFFSET, NF90_NETCDF4, &
       nf90_inq_varid, nf90_inq_dimid, nf90_inquire, nf90_get_var, nf90_fill_float, &
       nf90_fill_double, nf90_fill_int

  ! public routines -------------------------------------------------------------------
  public :: close_netcdf         ! save and close the netcdf file
  public :: create_netcdf        ! create the nc file with variables and their attributes, after they were set
  public :: dump_netcdf          ! simple dump of variable into a netcdf file
  public :: var2nc               ! simple dump of multiple variables with attributes into one netcdf file
  public :: write_dynamic_netcdf ! write dynamically (one record after the other) in the file
  public :: write_static_netcdf  ! write static data in the file
  ! public types -------------------------------------------------------------------
  public :: dims
  public :: variable
  public :: attribute

  ! public parameters
  integer(i4), public, parameter :: nMaxDim = 5         ! nr. max dimensions
  integer(i4), public, parameter :: nMaxAtt = 20        ! nr. max attributes
  integer(i4), public, parameter :: maxLen  = 256       ! nr. string length
  integer(i4), public, parameter :: nGAtt   = 20        ! nr. global attributes
  integer(i4), public, parameter :: nAttDim = 2         ! dim array of attribute values

  ! public types -----------------------------------------------------------------
  type dims
     character (len=maxLen)                 :: name                ! dim. name
     integer(i4)                            :: len                 ! dim. lenght, undefined time => NF90_UNLIMITED
     integer(i4)                            :: dimId               ! dim. Id
  end type dims

  type attribute
     character (len=maxLen)                  :: name                ! attribute name
     integer(i4)                             :: xType               ! attribute of the values
     integer(i4)                             :: nValues             ! number of attributes
     character (len=maxLen)                  :: values              ! numbers or "characters" separed by spaces
  end type attribute

  type variable
     character (len=maxLen)                 :: name                ! short name
     integer(i4)                            :: xType               ! NF90 var. type
     integer(i4)                            :: nLvls               ! number of levels
     integer(i4)                            :: nSubs               ! number of subparts
     logical                                :: unlimited           ! time limited
     integer(i4)                            :: variD               ! Id
     integer(i4)                            :: nDims               ! field dimension
     integer(i4), dimension(nMaxDim)        :: dimIds              ! passing var. dimensions
     integer(i4), dimension(nMaxDim)        :: dimTypes            ! type of dimensions
     integer(i4)                            :: nAtt                ! nr. attributes
     type(attribute), dimension(nMaxAtt)    :: att                 ! var. attributes
     integer(i4), dimension(nMaxDim)        :: start               ! starting indices for netcdf
     integer(i4), dimension(nMaxDim)        :: count               ! counter          for netcdf
     logical                                :: wFlag               ! write flag
     integer(i1),                     pointer :: G0_b              ! array pointing model variables
     integer(i1), dimension(:      ), pointer :: G1_b              ! array pointing model variables
     integer(i1), dimension(:,:    ), pointer :: G2_b              ! array pointing model variables
     integer(i1), dimension(:,:,:  ), pointer :: G3_b              ! array pointing model variables
     integer(i1), dimension(:,:,:,:), pointer :: G4_b              ! array pointing model variables
     integer(i4),                     pointer :: G0_i              ! array pointing model variables
     integer(i4), dimension(:      ), pointer :: G1_i              ! array pointing model variables
     integer(i4), dimension(:,:    ), pointer :: G2_i              ! array pointing model variables
     integer(i4), dimension(:,:,:  ), pointer :: G3_i              ! array pointing model variables
     integer(i4), dimension(:,:,:,:), pointer :: G4_i              ! array pointing model variables
     real(sp),                        pointer :: G0_f              ! array pointing model variables
     real(sp),    dimension(:    ),   pointer :: G1_f              ! array pointing model variables
     real(sp),    dimension(:,:    ), pointer :: G2_f              ! array pointing model variables
     real(sp),    dimension(:,:,:  ), pointer :: G3_f              ! array pointing model variables
     real(sp),    dimension(:,:,:,:), pointer :: G4_f              ! array pointing model variables
     real(dp),                        pointer :: G0_d              ! array pointing model variables
     real(dp),    dimension(:      ), pointer :: G1_d              ! array pointing model variables
     real(dp),    dimension(:,:    ), pointer :: G2_d              ! array pointing model variables
     real(dp),    dimension(:,:,:  ), pointer :: G3_d              ! array pointing model variables
     real(dp),    dimension(:,:,:,:), pointer :: G4_d              ! array pointing model variables
  end type variable

  ! public variables -----------------------------------------------------------------
  integer(i4),     public                            :: nVars   ! nr. variables
  integer(i4),     public                            :: nDims   ! nr. dimensions
  type (dims),     public, dimension(:), allocatable :: Dnc     ! dimensions list
  type(variable),  public, dimension(:), allocatable :: V       ! variable list, THIS STRUCTURE WILL BE WRITTEN IN THE FILE
  type(attribute), public, dimension(nGAtt)          :: gatt    ! global attributes for netcdf


  ! ------------------------------------------------------------------

  !     NAME
  !         dump_netcdf

  !     PURPOSE
  !         Simple write of a variable in a netcdf file.

  !         The variabel can be 1 to 5 dimensional and single or double precision.

  !         1D and 2D are dumped as static variables. From 3 to 5 dimension, the last
  !         dimension will be defined as time.
  !         The Variable will be called var.

  !     CALLING SEQUENCE
  !         call dump_netcdf(filename, arr, append, lfs, netcdf4, deflate_level)

  !     INTENT(IN)
  !         character(len=*) :: filename                  name of netcdf output file
  !         real(sp/dp)      :: arr(:[,:[,:[,:[,:]]]])    1D to 5D-array with input numbers

  !     INTENT(IN), OPTIONAL
  !         logical          :: lfs                       True: enable netcdf3 large file support, i.e. 64-bit offset
  !         logical          :: logical,                  True: use netcdf4 format
  !         integer(i4)      :: deflate_level             Compression level in netcdf4 (default: 1)

  !     RESTRICTIONS
  !         If dimension

  !     EXAMPLE
  !         call dump_netcdf('test.nc', myarray)
  !         call dump_netcdf('test.nc', myarray, netcdf4=.true.)
  !         -> see also example in test directory

  !     LITERATURE
  !         None

  !     HISTORY
  !         Written,  Matthias Cuntz, Nov 2012
  !         Modified, Stephan Thober, Nov 2012 - added functions for i4 variables
  !                   Matthias Cuntz and Juliane Mai
  !                                   Nov 2012 - append
  !                                            - fake time dimension for 1D and 2D
  !                                            - make i4 behave exactly as sp and dp
  !                                   Mar 2013 - lfs, netcdf4, deflate_level
  interface dump_netcdf
     module procedure dump_netcdf_1d_sp, dump_netcdf_2d_sp, dump_netcdf_3d_sp, &
          dump_netcdf_4d_sp, dump_netcdf_5d_sp, &
          dump_netcdf_1d_dp, dump_netcdf_2d_dp, dump_netcdf_3d_dp, &
          dump_netcdf_4d_dp, dump_netcdf_5d_dp, &
          dump_netcdf_1d_i4, dump_netcdf_2d_i4, dump_netcdf_3d_i4, &
          dump_netcdf_4d_i4, dump_netcdf_5d_i4
  end interface dump_netcdf


  ! ------------------------------------------------------------------
  !
  !     NAME
  !         var2nc
  !
  !     PURPOSE
  !         Provide an intermediate way to write netcdf files
  !         between dump_netcdf and create_netcdf
  !
  !>        \brief Extended dump_netcdf for multiple variables
  !
  !>        \details Write different variables including attributes to netcdf
  !>        file. The attributes are restricted to long_name, units,
  !>        and missing_value. It is also possible to append variables
  !>        when an unlimited dimension is specified.
  !
  !     CALLING SEQUENCE
  !>        call var2nc(f_name, arr, dnames, v_name, dim_unlimited,
  !>                    long_name, units, missing_value, create)
  !
  !     INTENT(IN)
  !>        \param[in] "character(*) :: f_name"                             filename
  !>        \param[in] "integer(i4)/real(sp,dp) :: arr(:[,:[,:[,:[,:]]]])"  array to write
  !>        \param[in] "character(*) :: dnames(:)"                          dimension names
  !>        \param[in] "character(*) :: v_name"                             variable name
  !
  !     INTENT(IN), OPTIONAL
  !>        \param[in] "integer(i4),             optional :: dim_unlimited"     index of unlimited dimension
  !>        \param[in] "character(*),            optional :: long_name"         attribute
  !>        \param[in] "character(*),            optional :: units"             attribute
  !>        \param[in] "integer(i4)/real(sp,dp), optional :: missing_value"     attribute
  !>        \param[in] "character(256), dimension(:,:), optional :: attributes" two dimensional array of attributes
  !>           size of first dimension equals number of attributes
  !>           first entry of second dimension equals attribute name (e.g. long_name)
  !>           second entry of second dimension equals attribute value (e.g. precipitation)
  !>           every attribute is written as string with the exception of missing_value
  !>        \param[in] "logical,                 optional :: create"            flag - specify whether a
  !>                                                                            output file should be
  !>                                                                            created, default
  !>        \param[inout] "integer(i4)/real(sp,dp), optional :: ncid"   if not given filename will be opened and closed
  !>                                                                    if given and <0 then file will be opened
  !>                                                                    and ncid will return the file unit.
  !>                                                                    if given and >0 then file is assumed open and
  !>                                                                    ncid is used as file unit.
  !>        \param[in] "integer(i4),                optional :: nrec"   if given: start point on unlimited dimension.
  !
  !     RESTRICTIONS
  !>        \note It is not allowed to write the folloing numbers for the indicated type\n
  !>                        number | kind\n
  !>                -2.1474836E+09 | integer(i4)\n
  !>                 9.9692100E+36 | real(sp)\n
  !>        9.9692099683868690E+36 | real(dp)\n
  !>        These numbers are netcdf fortran 90 constants! They are used to determine the 
  !>        chunksize of the already written variable. Hence, this routine cannot append
  !>        correctly to variables when these numbers are used. Only five dimensional 
  !>        variables can be written, only one unlimited dimension can be defined.
  !
  !     EXAMPLE
  !         Let <field> be some three dimensional array
  !           dnames(1) = 'x'
  !           dnames(2) = 'y'
  !           dnames(3) = 'time'
  !
  !         The simplest call to write <field> to a file is
  !           call var2nc('test.nc', field, dnames, 'h')
  !
  !         With attributes it looks like
  !           call var2nc('test.nc', field, dnames, 'h', &
  !                       long_name = 'height', units = '[m]', missing_value = -9999)
  !         or alternatively
  !         character(256), dimension(3,2) :: attributes
  !         attributes(1,1) = 'long_name'
  !         attributes(1,2) = 'precipitation'
  !         attributes(2,1) = 'units'
  !         attributes(2,2) = '[mm/d]'
  !         attributes(3,1) = 'missing_value'
  !         attributes(3,2) = '-9999.'
  !         call var2nc('test.nc', field, dnames, 'h', attributes = attributes, create = .true. )

  !         To be able to dynamically write <field>, an unlimited dimension
  !         needs to be specified (in this example field could also be only two
  !         dimensional)
  !         call var2nc('test.nc', field(:,:,1), dnames, 'h', dim_unlimited=3)
  !         Now one can append an arbitrary number of time steps, e.g., the next 9
  !         and the time has to be added again before
  !           call var2nc('test.nc', (/20,...,100/), dnames(3:3), 'time',
  !                       dim_unlimited = 1 )
  !           call var2nc('test.nc', field(:,:,2:10, dnames, 'h', dim_unlimited=3)
  !
  !         You can also write another variable sharing the same dimensions
  !           call var2nc('test.nc', field_2, dnames(1:2), 'h_2')
  !
  !         The netcdf file can stay open after the first call and subsequent calls can use the file unit
  !           ncid = -1_i4
  !           call var2nc('test.nc', field_1, dnames(1:1), 'h_1', ncid=ncid) ! opens file
  !           call var2nc('test.nc', field_2, dnames(1:2), 'h_2', ncid=ncid) ! uses ncid from last call
  !           call close_netcdf(ncid)
  !
  !         One can also give the start record number (on the unlimited dimension)
  !           ncid = -1_i4
  !           call var2nc('test.nc', time1, dnames(3:3), 'time', dim_unlimited=1_i4, ncid=ncid, create=.true.)
  !           do i=1, n
  !               call var2nc('test.nc', field_2(:,:,i), dnames(1:3), 'h_2', dim_unlimited=3_i4, ncid=ncid, nrec=i)
  !           end do
  !           call close_netcdf(ncid)
  !
  !         That's it, enjoy!
  !           -> see also example in test program
  !
  !    LITERATURE
  !        The manual of the used netcdf fortran library can be found in
  !        Robert Pincus & Ross Rew, The netcdf Fortran 90 Interface Guide
  !
  !    HISTORY
  !>       \author Stephan Thober & Matthias Cuntz
  !>       \date May 2014
  !        Modified Stephan Thober - Jun 2014 added deflate, shuffle, and chunksizes
  !                 Stephan Thober - Jun 2014 automatically append variable at the end,
  !                                           renamed _FillValue to missing_value
  !                 Stephan Thober - Jul 2014 add attributes array, introduced unlimited dimension
  !                                           that is added to the dimensions of the given array
  !                 Stephan Thober - Jan 2015 changed chunk_size convention to one chunk per unit in
  !                                           unlimited dimension (typically time)
  !                 Matthias Cuntz - Feb 2015 d_unlimit was not set in 5d cases
  !                                           use ne from mo_utils for fill value comparisons
  !                                           dummy(1) was sp instead of i4 in var2nc_1d_i4
  !                 Matthias Cuntz - May 2015 ncid for opening the file only once
  !                                           nrec for writing a specific record

  interface var2nc
     module procedure var2nc_1d_i4, var2nc_1d_sp, var2nc_1d_dp, &
          var2nc_2d_i4, var2nc_2d_sp, var2nc_2d_dp, &
          var2nc_3d_i4, var2nc_3d_sp, var2nc_3d_dp, &
          var2nc_4d_i4, var2nc_4d_sp, var2nc_4d_dp, &
          var2nc_5d_i4, var2nc_5d_sp, var2nc_5d_dp
  end interface var2nc

  ! ----------------------------------------------------------------------------

  private

  ! ----------------------------------------------------------------------------

contains

  ! ----------------------------------------------------------------------------

  ! NAME
  !     close_netcdf

  ! PURPOSE
  !     closes a stream of an open netcdf file and saves the file.

  ! CALLING SEQUENCE
  !     call close_netcdf(nc)

  ! INTENT(IN)
  !     integer(i4) :: ncid - stream id of an open netcdf file which shall be closed

  ! RESTRICTIONS
  !     Closes only an already open stream

  ! EXAMPLE
  !     See test_mo_ncwrite

  ! LITERATURE
  !     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90.html

  ! HISTORY
  !     Written  Luis Samaniego  Feb 2011
  !     Modified Stephan Thober  Dec 2011 - added comments and generalized
  !              Matthias Cuntz  Jan 2012 - Info
  !              Matthias Cuntz  Mar 2013 - removed Info

  subroutine close_netcdf(ncId)

    implicit none

    integer(i4), intent(in)           :: ncId

    ! close: save new netcdf dataset
    call check(nf90_close(ncId))

  end subroutine close_netcdf

  ! ------------------------------------------------------------------------------

  ! NAME
  !     create_netcdf

  ! PURPOSE
  !     This subroutine will open a new netcdf file and write the variable
  !     attributes stored in the structure V in the file. Therefore V
  !     has to be already set. See the file set_netcdf in test_mo_ncwrite
  !     for an example.

  ! CALLING SEQUENCE
  !     call create_netcdf(File, ncid, lfs, netcdf4, deflate_level)

  ! INTENT(IN)
  !     character(len=maxLen) :: File             Filename of file to be written

  ! INTENT(IN), OPTIONAL
  !     logical               :: lfs              True: enable netcdf3 large file support, i.e. 64-bit offset
  !     logical               :: logical          True: use netcdf4 format
  !     integer(i4)           :: deflate_level    compression level in netcdf4 (default: 1)

  ! INTENT(OUT)
  !     integer(i4)           :: ncid             integer value of the stream for the opened file

  ! RESTRICTIONS
  !     This routine only writes attributes and variables which have been stored in V
  !     nothing else.

  ! EXAMPLE
  !     see test_mo_ncwrite

  ! LITERATURE
  !     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90.html

  ! HISTORY
  !     Written  Luis Samaniego  Feb 2011
  !     Modified Stephan Thober  Dec 2011 - added comments and generalized
  !              Matthias Cuntz  Jan 2012 - Info
  !              Stephan Thober  Feb 2013 - added flag for large file support
  !              Matthias Cuntz  Mar 2013 - netcdf4, deflate_level
  !              Stephan Thober  Mar 2013 - buffersize
  !              Matthias Cuntz  Mar 2013 - removed Info

  subroutine create_netcdf(Filename, ncid, lfs, netcdf4, deflate_level)

    implicit none

    ! netcdf related variables
    character(len=*), intent(in)           :: Filename
    integer(i4),      intent(out)          :: ncid
    logical,          intent(in), optional :: lfs           ! netcdf3 Large File Support
    logical,          intent(in), optional :: netcdf4       ! netcdf4
    integer(i4),      intent(in), optional :: deflate_level ! compression level in netcdf4

    integer(i4)                               :: i, j, k
    integer(i4), dimension(nAttDim)           :: att_INT
    real(sp),    dimension(nAttDim)           :: att_FLOAT
    real(dp),    dimension(nAttDim)           :: att_DOUBLE
    character(len=maxLen), dimension(nAttDim) :: att_CHAR
    logical                                   :: LargeFile
    logical                                   :: inetcdf4
    integer(i4)                               :: deflate
    integer(i4)                               :: buffersize
    integer(i4), dimension(:), allocatable    :: chunksizes ! Size of chunks in netcdf4 writing

    LargeFile = .false.
    if (present(lfs)) LargeFile = lfs
    inetcdf4 = .true.
    if (present(netcdf4)) inetcdf4 = netcdf4
    deflate = 1
    if (present(deflate_level)) deflate = deflate_level
    ! 1  Create netcdf dataset: enter define mode     ->  get ncId
    if (inetcdf4) then
       call check(nf90_create(trim(Filename), NF90_NETCDF4, ncId))
    else
       if (LargeFile) then
          call check(nf90_create(trim(Filename), NF90_64BIT_OFFSET, ncId, chunksize=buffersize))
       else
          ! let the netcdf library chose a buffersize, that results in lesser write system calls
          call check(nf90_create(trim(Filename), NF90_CLOBBER, ncId, chunksize=buffersize))
       end if
    endif

    ! 2  Define dimensions                                 -> get dimId
    do i=1, nDims
       call check(nf90_def_dim(ncId, Dnc(i)%name , Dnc(i)%len, Dnc(i)%dimId))
    end do

    ! 3 Define dimids array, which is used to pass the dimids of the dimensions of
    ! the netcdf variables
    do i = 1, nVars
       V(i)%unlimited = .false.
       V(i)%dimids    = 0
       V(i)%start     = 1
       V(i)%count     = 1
       do k = 1, V(i)%nDims
          if (Dnc(V(i)%dimTypes(k))%len  == NF90_UNLIMITED) V(i)%unlimited = .true.
          V(i)%dimids(k) = Dnc(V(i)%dimTypes(k))%dimId
       end do
       if (V(i)%unlimited) then
          ! set counts for unlimited files (time is always the last dimension)
          if (V(i)%nDims == 1) cycle
          do k = 1, V(i)%nDims - 1
             V(i)%count(k)  = Dnc(V(i)%dimTypes(k))%len
          end do
       end if
    end do

    ! 4 Define the netcdf variables and attributes                            -> get varId
    allocate(chunksizes(maxval(V(1:nVars)%nDims)))
    do i=1, nVars
       if (.not. V(i)%wFlag) cycle
       if (inetcdf4) then
          chunksizes(1:V(i)%nDims) = Dnc(V(i)%dimTypes(1:V(i)%nDims))%len
          chunksizes(V(i)%nDims)   = 1
          call check(nf90_def_var(ncId, V(i)%name, V(i)%xtype, V(i)%dimids(1:V(i)%nDims), V(i)%varId, &
               chunksizes=chunksizes(1:V(i)%nDims), shuffle=.true., deflate_level=deflate))
       else
          call check(nf90_def_var(ncId, V(i)%name, V(i)%xtype, V(i)%dimids(1:V(i)%nDims), V(i)%varId))
       endif
       do k = 1, V(i)%nAtt
          select case (V(i)%att(k)%xType)
          case (NF90_CHAR)
             ! read(V(i)%att(k)%values, *) (att_CHAR(j), j =1, V(i)%att(k)%nValues)
             read(V(i)%att(k)%values, '(a)')  att_CHAR(1)
             call check(nf90_put_att (ncId, V(i)%varId, V(i)%att(k)%name, att_CHAR(1)))
          case (NF90_INT)
             read(V(i)%att(k)%values, *) (att_INT(j), j =1, V(i)%att(k)%nValues)
             call check(nf90_put_att (ncId, V(i)%varId, V(i)%att(k)%name, att_INT(1:V(i)%att(k)%nValues)))
          case (NF90_FLOAT)
             read(V(i)%att(k)%values, *) (att_FLOAT(j), j =1, V(i)%att(k)%nValues)
             call check(nf90_put_att (ncId, V(i)%varId, V(i)%att(k)%name, att_FLOAT(1:V(i)%att(k)%nValues)))
          case (NF90_DOUBLE)
             read(V(i)%att(k)%values, *) (att_DOUBLE(j), j =1, V(i)%att(k)%nValues)
             call check(nf90_put_att (ncId, V(i)%varId, V(i)%att(k)%name, att_DOUBLE(1:V(i)%att(k)%nValues)))
          end select
       end do
    end do

    ! 5 Global attributes
    do k = 1, nGAtt
       if (nonull(Gatt(k)%name)) then
          call check(nf90_put_att(ncId, NF90_GLOBAL, Gatt(k)%name, Gatt(k)%values))
       endif
    end do

    ! 6 end definitions: leave define mode
    call check(nf90_enddef(ncId))

    deallocate(chunksizes)

  end subroutine create_netcdf


  ! ------------------------------------------------------------------

  subroutine dump_netcdf_1d_sp(filename, arr, append, lfs, netcdf4, deflate_level)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    real(sp),         dimension(:),     intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file
    logical,          optional,         intent(in) :: lfs           ! netcdf3 Large File Support
    logical,          optional,         intent(in) :: netcdf4       ! netcdf4
    integer(i4),      optional,         intent(in) :: deflate_level ! compression level in netcdf4

    integer(i4),      parameter         :: ndim = 1 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim+1) :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim+1) :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+2) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim+1) :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim+1) :: counter  ! length array for write of each time step
    integer(i4),      dimension(ndim+1) :: chunksizes ! Size of chunks in netcdf4 writing
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension
    logical                  :: LargeFile
    logical                  :: inetcdf4
    integer(i4)              :: deflate
    integer(i4)              :: buffersize

    ! append or not
    if (present(append)) then
       if (append) then
          iappend = .true.
       else
          iappend = .false.
       endif
    else
       iappend = .false.
    endif
    LargeFile = .false.
    if (present(lfs)) LargeFile = lfs
    inetcdf4 = .false.
    if (present(netcdf4)) inetcdf4 = netcdf4
    deflate = 1
    if (present(deflate_level)) deflate = deflate_level

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)

    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))

       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim+1)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+2)))
       call check(nf90_inquire_variable(ncid, varid(ndim+2), ndims=idim, dimids=dimid))
       if (idim /= ndim+1) stop "dump_netcdf_1d_sp: number of variable dimensions /= number of file variable dimensions."

       ! inquire dimensions
       do i=1, ndim+1
          call check(nf90_inquire_dimension(ncid, dimid(i), name, dims(i)))
          if (i < ndim+1) then
             if (trim(name) /= dnames(i)) stop "dump_netcdf_1d_sp: dimension name problem."
             if (dims(i) /= size(arr,i)) stop "dump_netcdf_1d_sp: variable dimension /= file variable dimension."
          else
             if (trim(name) /= 'time') stop "dump_netcdf_1d_sp: time name problem."
          endif
       enddo

       ! append
       start(:)        = 1
       counter(:)      = dims
       counter(ndim+1) = 1
       do i=1, 1
          start(ndim+1) = dims(ndim+1) + i
          call check(nf90_put_var(ncid, varid(ndim+1), (/dims(ndim+1)+i/), (/dims(ndim+1)+i/)))
          call check(nf90_put_var(ncid, varid(ndim+2), arr, start, counter))
       end do
    else
       ! open file
       if (inetcdf4) then
          call check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
       else
          if (LargeFile) then
             call check(nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid, chunksize=buffersize))
          else
             call check(nf90_create(trim(filename), NF90_CLOBBER, ncid, chunksize=buffersize))
          end if
       endif

       ! define dims
       dims(1:ndim) = shape(arr)
       do i=1, ndim
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       dims(ndim+1) = 1
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim+1)))

       ! define dim variables
       do i=1, ndim
          if (inetcdf4) then
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          else
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          endif
       end do
       ! define time variable
       if (inetcdf4) then
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim+1), varid(ndim+1)))
       else
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim+1), varid(ndim+1)))
       endif

       ! define variable
       if (inetcdf4) then
          chunksizes(1:ndim) = dims(1:ndim)
          chunksizes(ndim+1) = 1
          call check(nf90_def_var(ncid, 'var', NF90_FLOAT, dimid, varid(ndim+2), &
               chunksizes=chunksizes, shuffle=.true., deflate_level=deflate))
       else
          call check(nf90_def_var(ncid, 'var', NF90_FLOAT, dimid, varid(ndim+2)))
       endif

       ! end define mode
       call check(nf90_enddef(ncid))

       ! write dimensions
       do i=1, ndim
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do

       ! write time and variable
       start(:)        = 1
       counter(:)      = dims
       counter(ndim+1) = 1
       do i=1, 1
          start(ndim+1) = i
          call check(nf90_put_var(ncid, varid(ndim+1), (/i/), (/i/)))
          call check(nf90_put_var(ncid, varid(ndim+2), arr, start, counter))
       end do
    endif

    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_1d_sp


  subroutine dump_netcdf_2d_sp(filename, arr, append, lfs, netcdf4, deflate_level)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    real(sp),         dimension(:,:),   intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file
    logical,          optional,         intent(in) :: lfs           ! netcdf3 Large File Support
    logical,          optional,         intent(in) :: netcdf4       ! netcdf4
    integer(i4),      optional,         intent(in) :: deflate_level ! compression level in netcdf4

    integer(i4),      parameter         :: ndim = 2 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim+1) :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim+1) :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+2) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim+1) :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim+1) :: counter  ! length array for write of each time step
    integer(i4),      dimension(ndim+1) :: chunksizes ! Size of chunks in netcdf4 writing
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension
    logical                  :: LargeFile
    logical                  :: inetcdf4
    integer(i4)              :: deflate
    integer(i4)              :: buffersize

    ! append or not
    if (present(append)) then
       if (append) then
          iappend = .true.
       else
          iappend = .false.
       endif
    else
       iappend = .false.
    endif
    LargeFile = .false.
    if (present(lfs)) LargeFile = lfs
    inetcdf4 = .false.
    if (present(netcdf4)) inetcdf4 = netcdf4
    deflate = 1
    if (present(deflate_level)) deflate = deflate_level

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)

    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))

       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim+1)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+2)))
       call check(nf90_inquire_variable(ncid, varid(ndim+2), ndims=idim, dimids=dimid))
       if (idim /= ndim+1) stop "dump_netcdf_2d_sp: number of variable dimensions /= number of file variable dimensions."

       ! inquire dimensions
       do i=1, ndim+1
          call check(nf90_inquire_dimension(ncid, dimid(i), name, dims(i)))
          if (i < ndim+1) then
             if (trim(name) /= dnames(i)) stop "dump_netcdf_2d_sp: dimension name problem."
             if (dims(i) /= size(arr,i)) stop "dump_netcdf_2d_sp: variable dimension /= file variable dimension."
          else
             if (trim(name) /= 'time') stop "dump_netcdf_2d_sp: time name problem."
          endif
       enddo

       ! append
       start(:)        = 1
       counter(:)      = dims
       counter(ndim+1) = 1
       do i=1, 1
          start(ndim+1) = dims(ndim+1) + i
          call check(nf90_put_var(ncid, varid(ndim+1), (/dims(ndim+1)+i/), (/dims(ndim+1)+i/)))
          call check(nf90_put_var(ncid, varid(ndim+2), arr, start, counter))
       end do
    else
       ! open file
       if (inetcdf4) then
          call check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
       else
          if (LargeFile) then
             call check(nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid, chunksize=buffersize))
          else
             call check(nf90_create(trim(filename), NF90_CLOBBER, ncid, chunksize=buffersize))
          end if
       endif

       ! define dims
       dims(1:ndim) = shape(arr)
       do i=1, ndim
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       dims(ndim+1) = 1
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim+1)))

       ! define dim variables
       do i=1, ndim
          if (inetcdf4) then
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          else
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          endif
       end do
       ! define time variable
       if (inetcdf4) then
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim+1), varid(ndim+1)))
       else
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim+1), varid(ndim+1)))
       endif

       ! define variable
       if (inetcdf4) then
          chunksizes(1:ndim) = dims(1:ndim)
          chunksizes(ndim+1) = 1
          call check(nf90_def_var(ncid, 'var', NF90_FLOAT, dimid, varid(ndim+2), &
               chunksizes=chunksizes, shuffle=.true., deflate_level=deflate))
       else
          call check(nf90_def_var(ncid, 'var', NF90_FLOAT, dimid, varid(ndim+2)))
       endif

       ! end define mode
       call check(nf90_enddef(ncid))

       ! write dimensions
       do i=1, ndim
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do

       ! write time and variable
       start(:)        = 1
       counter(:)      = dims
       counter(ndim+1) = 1
       do i=1, 1
          start(ndim+1) = i
          call check(nf90_put_var(ncid, varid(ndim+1), (/i/), (/i/)))
          call check(nf90_put_var(ncid, varid(ndim+2), arr, start, counter))
       end do
    endif

    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_2d_sp


  subroutine dump_netcdf_3d_sp(filename, arr, append, lfs, netcdf4, deflate_level)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    real(sp),         dimension(:,:,:), intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file
    logical,          optional,         intent(in) :: lfs           ! netcdf3 Large File Support
    logical,          optional,         intent(in) :: netcdf4       ! netcdf4
    integer(i4),      optional,         intent(in) :: deflate_level ! compression level in netcdf4

    integer(i4),      parameter         :: ndim = 3 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim)   :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim)   :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+1) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim)   :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim)   :: counter  ! length array for write of each time step
    integer(i4),      dimension(ndim)   :: chunksizes ! Size of chunks in netcdf4 writing
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension
    logical                  :: LargeFile
    logical                  :: inetcdf4
    integer(i4)              :: deflate
    integer(i4)              :: buffersize

    ! append or not
    if (present(append)) then
       if (append) then
          iappend = .true.
       else
          iappend = .false.
       endif
    else
       iappend = .false.
    endif
    LargeFile = .false.
    if (present(lfs)) LargeFile = lfs
    inetcdf4 = .false.
    if (present(netcdf4)) inetcdf4 = netcdf4
    deflate = 1
    if (present(deflate_level)) deflate = deflate_level

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)

    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))

       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+1)))
       call check(nf90_inquire_variable(ncid, varid(ndim+1), ndims=idim, dimids=dimid))
       if (idim /= ndim) stop "dump_netcdf_3d_sp: number of variable dimensions /= number of file variable dimensions."

       ! inquire dimensions
       do i=1, ndim
          call check(nf90_inquire_dimension(ncid, dimid(i), name, dims(i)))
          if (i < ndim) then
             if (trim(name) /= dnames(i)) stop "dump_netcdf_3d_sp: dimension name problem."
             if (dims(i) /= size(arr,i)) stop "dump_netcdf_3d_sp: variable dimension /= file variable dimension."
          else
             if (trim(name) /= 'time') stop "dump_netcdf_3d_sp: time name problem."
          endif
       enddo

       ! append
       start(:)      = 1
       counter(:)    = dims
       counter(ndim) = 1
       do i=1, size(arr,ndim)
          start(ndim) = dims(ndim) + i
          call check(nf90_put_var(ncid, varid(ndim), (/dims(ndim)+i/), (/dims(ndim)+i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,i), start, counter))
       end do
    else
       ! open file
       if (inetcdf4) then
          call check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
       else
          if (LargeFile) then
             call check(nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid, chunksize=buffersize))
          else
             call check(nf90_create(trim(filename), NF90_CLOBBER, ncid, chunksize=buffersize))
          end if
       endif

       ! define dims
       dims = shape(arr)
       do i=1, ndim-1
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim)))

       ! define dim variables
       do i=1, ndim-1
          if (inetcdf4) then
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          else
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          endif
       end do
       ! define time variable
       if (inetcdf4) then
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       else
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       endif

       ! define variable
       if (inetcdf4) then
          chunksizes(1:ndim-1) = dims(1:ndim-1)
          chunksizes(ndim)     = 1
          call check(nf90_def_var(ncid, 'var', NF90_FLOAT, dimid, varid(ndim+1), &
               chunksizes=chunksizes, shuffle=.true., deflate_level=deflate))
       else
          call check(nf90_def_var(ncid, 'var', NF90_FLOAT, dimid, varid(ndim+1)))
       endif

       ! end define mode
       call check(nf90_enddef(ncid))

       ! write dimensions
       do i=1, ndim-1
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do

       ! write time and variable
       start(:) = 1
       counter(:) = dims
       counter(ndim) = 1
       do i=1, dims(ndim)
          start(ndim) = i
          call check(nf90_put_var(ncid, varid(ndim), (/i/), (/i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,i), start, counter))
       end do
    endif

    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_3d_sp


  subroutine dump_netcdf_4d_sp(filename, arr, append, lfs, netcdf4, deflate_level)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    real(sp),         dimension(:,:,:,:), intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file
    logical,          optional,         intent(in) :: lfs           ! netcdf3 Large File Support
    logical,          optional,         intent(in) :: netcdf4       ! netcdf4
    integer(i4),      optional,         intent(in) :: deflate_level ! compression level in netcdf4

    integer(i4),      parameter         :: ndim = 4 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim)   :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim)   :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+1) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim)   :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim)   :: counter  ! length array for write of each time step
    integer(i4),      dimension(ndim)   :: chunksizes ! Size of chunks in netcdf4 writing
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension
    logical                  :: LargeFile
    logical                  :: inetcdf4
    integer(i4)              :: deflate
    integer(i4)              :: buffersize

    ! append or not
    if (present(append)) then
       if (append) then
          iappend = .true.
       else
          iappend = .false.
       endif
    else
       iappend = .false.
    endif
    LargeFile = .false.
    if (present(lfs)) LargeFile = lfs
    inetcdf4 = .false.
    if (present(netcdf4)) inetcdf4 = netcdf4
    deflate = 1
    if (present(deflate_level)) deflate = deflate_level

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)

    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))

       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+1)))
       call check(nf90_inquire_variable(ncid, varid(ndim+1), ndims=idim, dimids=dimid))
       if (idim /= ndim) stop "dump_netcdf_4d_sp: number of variable dimensions /= number of file variable dimensions."

       ! inquire dimensions
       do i=1, ndim
          call check(nf90_inquire_dimension(ncid, dimid(i), name, dims(i)))
          if (i < ndim) then
             if (trim(name) /= dnames(i)) stop "dump_netcdf_4d_sp: dimension name problem."
             if (dims(i) /= size(arr,i)) stop "dump_netcdf_4d_sp: variable dimension /= file variable dimension."
          else
             if (trim(name) /= 'time') stop "dump_netcdf_4d_sp: time name problem."
          endif
       enddo

       ! append
       start(:)      = 1
       counter(:)    = dims
       counter(ndim) = 1
       do i=1, size(arr,ndim)
          start(ndim) = dims(ndim) + i
          call check(nf90_put_var(ncid, varid(ndim), (/dims(ndim)+i/), (/dims(ndim)+i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,:,i), start, counter))
       end do
    else
       ! open file
       if (inetcdf4) then
          call check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
       else
          if (LargeFile) then
             call check(nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid, chunksize=buffersize))
          else
             call check(nf90_create(trim(filename), NF90_CLOBBER, ncid, chunksize=buffersize))
          end if
       endif

       ! define dims
       dims = shape(arr)
       do i=1, ndim-1
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim)))

       ! define dim variables
       do i=1, ndim-1
          if (inetcdf4) then
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          else
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          endif
       end do
       ! define time variable
       if (inetcdf4) then
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       else
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       endif

       ! define variable
       if (inetcdf4) then
          chunksizes(1:ndim-1) = dims(1:ndim-1)
          chunksizes(ndim)     = 1
          call check(nf90_def_var(ncid, 'var', NF90_FLOAT, dimid, varid(ndim+1), &
               chunksizes=chunksizes, shuffle=.true., deflate_level=deflate))
       else
          call check(nf90_def_var(ncid, 'var', NF90_FLOAT, dimid, varid(ndim+1)))
       endif

       ! end define mode
       call check(nf90_enddef(ncid))

       ! write dimensions
       do i=1, ndim-1
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do

       ! write time and variable
       start(:) = 1
       counter(:) = dims
       counter(ndim) = 1
       do i=1, dims(ndim)
          start(ndim) = i
          call check(nf90_put_var(ncid, varid(ndim), (/i/), (/i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,:,i), start, counter))
       end do
    endif

    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_4d_sp


  subroutine dump_netcdf_5d_sp(filename, arr, append, lfs, netcdf4, deflate_level)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    real(sp),         dimension(:,:,:,:,:), intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file
    logical,          optional,         intent(in) :: lfs           ! netcdf3 Large File Support
    logical,          optional,         intent(in) :: netcdf4       ! netcdf4
    integer(i4),      optional,         intent(in) :: deflate_level ! compression level in netcdf4

    integer(i4),      parameter         :: ndim = 5 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim)   :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim)   :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+1) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim)   :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim)   :: counter  ! length array for write of each time step
    integer(i4),      dimension(ndim)   :: chunksizes ! Size of chunks in netcdf4 writing
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension
    logical                  :: LargeFile
    logical                  :: inetcdf4
    integer(i4)              :: deflate
    integer(i4)              :: buffersize

    ! append or not
    if (present(append)) then
       if (append) then
          iappend = .true.
       else
          iappend = .false.
       endif
    else
       iappend = .false.
    endif
    LargeFile = .false.
    if (present(lfs)) LargeFile = lfs
    inetcdf4 = .false.
    if (present(netcdf4)) inetcdf4 = netcdf4
    deflate = 1
    if (present(deflate_level)) deflate = deflate_level

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)

    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))

       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+1)))
       call check(nf90_inquire_variable(ncid, varid(ndim+1), ndims=idim, dimids=dimid))
       if (idim /= ndim) stop "dump_netcdf_5d_sp: number of variable dimensions /= number of file variable dimensions."

       ! inquire dimensions
       do i=1, ndim
          call check(nf90_inquire_dimension(ncid, dimid(i), name, dims(i)))
          if (i < ndim) then
             if (trim(name) /= dnames(i)) stop "dump_netcdf_5d_sp: dimension name problem."
             if (dims(i) /= size(arr,i)) stop "dump_netcdf_5d_sp: variable dimension /= file variable dimension."
          else
             if (trim(name) /= 'time') stop "dump_netcdf_5d_sp: time name problem."
          endif
       enddo

       ! append
       start(:)      = 1
       counter(:)    = dims
       counter(ndim) = 1
       do i=1, size(arr,ndim)
          start(ndim) = dims(ndim) + i
          call check(nf90_put_var(ncid, varid(ndim), (/dims(ndim)+i/), (/dims(ndim)+i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,:,:,i), start, counter))
       end do
    else
       ! open file
       if (inetcdf4) then
          call check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
       else
          if (LargeFile) then
             call check(nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid, chunksize=buffersize))
          else
             call check(nf90_create(trim(filename), NF90_CLOBBER, ncid, chunksize=buffersize))
          end if
       endif

       ! define dims
       dims = shape(arr)
       do i=1, ndim-1
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim)))

       ! define dim variables
       do i=1, ndim-1
          if (inetcdf4) then
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          else
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          endif
       end do
       ! define time variable
       if (inetcdf4) then
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       else
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       endif

       ! define variable
       if (inetcdf4) then
          chunksizes(1:ndim-1) = dims(1:ndim-1)
          chunksizes(ndim)     = 1
          call check(nf90_def_var(ncid, 'var', NF90_FLOAT, dimid, varid(ndim+1), &
               chunksizes=chunksizes, shuffle=.true., deflate_level=deflate))
       else
          call check(nf90_def_var(ncid, 'var', NF90_FLOAT, dimid, varid(ndim+1)))
       endif

       ! end define mode
       call check(nf90_enddef(ncid))

       ! write dimensions
       do i=1, ndim-1
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do

       ! write time and variable
       start(:) = 1
       counter(:) = dims
       counter(ndim) = 1
       do i=1, dims(ndim)
          start(ndim) = i
          call check(nf90_put_var(ncid, varid(ndim), (/i/), (/i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,:,:,i), start, counter))
       end do
    endif

    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_5d_sp


  subroutine dump_netcdf_1d_dp(filename, arr, append, lfs, netcdf4, deflate_level)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    real(dp),         dimension(:),     intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file
    logical,          optional,         intent(in) :: lfs           ! netcdf3 Large File Support
    logical,          optional,         intent(in) :: netcdf4       ! netcdf4
    integer(i4),      optional,         intent(in) :: deflate_level ! compression level in netcdf4

    integer(i4),      parameter         :: ndim = 1 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim+1) :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim+1) :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+2) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim+1) :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim+1) :: counter  ! length array for write of each time step
    integer(i4),      dimension(ndim+1) :: chunksizes ! Size of chunks in netcdf4 writing
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension
    logical                  :: LargeFile
    logical                  :: inetcdf4
    integer(i4)              :: deflate
    integer(i4)              :: buffersize

    ! append or not
    if (present(append)) then
       if (append) then
          iappend = .true.
       else
          iappend = .false.
       endif
    else
       iappend = .false.
    endif
    LargeFile = .false.
    if (present(lfs)) LargeFile = lfs
    inetcdf4 = .false.
    if (present(netcdf4)) inetcdf4 = netcdf4
    deflate = 1
    if (present(deflate_level)) deflate = deflate_level

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)

    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))

       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim+1)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+2)))
       call check(nf90_inquire_variable(ncid, varid(ndim+2), ndims=idim, dimids=dimid))
       if (idim /= ndim+1) stop "dump_netcdf_1d_dp: number of variable dimensions /= number of file variable dimensions."

       ! inquire dimensions
       do i=1, ndim+1
          call check(nf90_inquire_dimension(ncid, dimid(i), name, dims(i)))
          if (i < ndim+1) then
             if (trim(name) /= dnames(i)) stop "dump_netcdf_1d_dp: dimension name problem."
             if (dims(i) /= size(arr,i)) stop "dump_netcdf_1d_dp: variable dimension /= file variable dimension."
          else
             if (trim(name) /= 'time') stop "dump_netcdf_1d_dp: time name problem."
          endif
       enddo

       ! append
       start(:)        = 1
       counter(:)      = dims
       counter(ndim+1) = 1
       do i=1, 1
          start(ndim+1) = dims(ndim+1) + i
          call check(nf90_put_var(ncid, varid(ndim+1), (/dims(ndim+1)+i/), (/dims(ndim+1)+i/)))
          call check(nf90_put_var(ncid, varid(ndim+2), arr, start, counter))
       end do
    else
       ! open file
       if (inetcdf4) then
          call check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
       else
          if (LargeFile) then
             call check(nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid, chunksize=buffersize))
          else
             call check(nf90_create(trim(filename), NF90_CLOBBER, ncid, chunksize=buffersize))
          end if
       endif

       ! define dims
       dims(1:ndim) = shape(arr)
       do i=1, ndim
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       dims(ndim+1) = 1
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim+1)))

       ! define dim variables
       do i=1, ndim
          if (inetcdf4) then
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          else
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          endif
       end do
       ! define time variable
       if (inetcdf4) then
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim+1), varid(ndim+1)))
       else
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim+1), varid(ndim+1)))
       endif

       ! define variable
       if (inetcdf4) then
          chunksizes(1:ndim) = dims(1:ndim)
          chunksizes(ndim+1) = 1
          call check(nf90_def_var(ncid, 'var', NF90_DOUBLE, dimid, varid(ndim+2), &
               chunksizes=chunksizes, shuffle=.true., deflate_level=deflate))
       else
          call check(nf90_def_var(ncid, 'var', NF90_DOUBLE, dimid, varid(ndim+2)))
       endif

       ! end define mode
       call check(nf90_enddef(ncid))

       ! write dimensions
       do i=1, ndim
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do

       ! write time and variable
       start(:)        = 1
       counter(:)      = dims
       counter(ndim+1) = 1
       do i=1, 1
          start(ndim+1) = i
          call check(nf90_put_var(ncid, varid(ndim+1), (/i/), (/i/)))
          call check(nf90_put_var(ncid, varid(ndim+2), arr, start, counter))
       end do
    endif

    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_1d_dp


  subroutine dump_netcdf_2d_dp(filename, arr, append, lfs, netcdf4, deflate_level)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    real(dp),         dimension(:,:),   intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file
    logical,          optional,         intent(in) :: lfs           ! netcdf3 Large File Support
    logical,          optional,         intent(in) :: netcdf4       ! netcdf4
    integer(i4),      optional,         intent(in) :: deflate_level ! compression level in netcdf4

    integer(i4),      parameter         :: ndim = 2 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim+1) :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim+1) :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+2) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim+1) :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim+1) :: counter  ! length array for write of each time step
    integer(i4),      dimension(ndim+1) :: chunksizes ! Size of chunks in netcdf4 writing
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension
    logical                  :: LargeFile
    logical                  :: inetcdf4
    integer(i4)              :: deflate
    integer(i4)              :: buffersize

    ! append or not
    if (present(append)) then
       if (append) then
          iappend = .true.
       else
          iappend = .false.
       endif
    else
       iappend = .false.
    endif
    LargeFile = .false.
    if (present(lfs)) LargeFile = lfs
    inetcdf4 = .false.
    if (present(netcdf4)) inetcdf4 = netcdf4
    deflate = 1
    if (present(deflate_level)) deflate = deflate_level

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)

    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))

       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim+1)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+2)))
       call check(nf90_inquire_variable(ncid, varid(ndim+2), ndims=idim, dimids=dimid))
       if (idim /= ndim+1) stop "dump_netcdf_2d_dp: number of variable dimensions /= number of file variable dimensions."

       ! inquire dimensions
       do i=1, ndim+1
          call check(nf90_inquire_dimension(ncid, dimid(i), name, dims(i)))
          if (i < ndim+1) then
             if (trim(name) /= dnames(i)) stop "dump_netcdf_2d_dp: dimension name problem."
             if (dims(i) /= size(arr,i)) stop "dump_netcdf_2d_dp: variable dimension /= file variable dimension."
          else
             if (trim(name) /= 'time') stop "dump_netcdf_2d_dp: time name problem."
          endif
       enddo

       ! append
       start(:)        = 1
       counter(:)      = dims
       counter(ndim+1) = 1
       do i=1, 1
          start(ndim+1) = dims(ndim+1) + i
          call check(nf90_put_var(ncid, varid(ndim+1), (/dims(ndim+1)+i/), (/dims(ndim+1)+i/)))
          call check(nf90_put_var(ncid, varid(ndim+2), arr, start, counter))
       end do
    else
       ! open file
       if (inetcdf4) then
          call check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
       else
          if (LargeFile) then
             call check(nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid, chunksize=buffersize))
          else
             call check(nf90_create(trim(filename), NF90_CLOBBER, ncid, chunksize=buffersize))
          end if
       endif

       ! define dims
       dims(1:ndim) = shape(arr)
       do i=1, ndim
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       dims(ndim+1) = 1
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim+1)))

       ! define dim variables
       do i=1, ndim
          if (inetcdf4) then
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          else
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          endif
       end do
       ! define time variable
       if (inetcdf4) then
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim+1), varid(ndim+1)))
       else
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim+1), varid(ndim+1)))
       endif

       ! define variable
       if (inetcdf4) then
          chunksizes(1:ndim) = dims(1:ndim)
          chunksizes(ndim+1) = 1
          call check(nf90_def_var(ncid, 'var', NF90_DOUBLE, dimid, varid(ndim+2), &
               chunksizes=chunksizes, shuffle=.true., deflate_level=deflate))
       else
          call check(nf90_def_var(ncid, 'var', NF90_DOUBLE, dimid, varid(ndim+2)))
       endif

       ! end define mode
       call check(nf90_enddef(ncid))

       ! write dimensions
       do i=1, ndim
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do

       ! write time and variable
       start(:)        = 1
       counter(:)      = dims
       counter(ndim+1) = 1
       do i=1, 1
          start(ndim+1) = i
          call check(nf90_put_var(ncid, varid(ndim+1), (/i/), (/i/)))
          call check(nf90_put_var(ncid, varid(ndim+2), arr, start, counter))
       end do
    endif

    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_2d_dp


  subroutine dump_netcdf_3d_dp(filename, arr, append, lfs, netcdf4, deflate_level)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    real(dp),         dimension(:,:,:), intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file
    logical,          optional,         intent(in) :: lfs           ! netcdf3 Large File Support
    logical,          optional,         intent(in) :: netcdf4       ! netcdf4
    integer(i4),      optional,         intent(in) :: deflate_level ! compression level in netcdf4

    integer(i4),      parameter         :: ndim = 3 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim)   :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim)   :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+1) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim)   :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim)   :: counter  ! length array for write of each time step
    integer(i4),      dimension(ndim)   :: chunksizes ! Size of chunks in netcdf4 writing
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension
    logical                  :: LargeFile
    logical                  :: inetcdf4
    integer(i4)              :: deflate
    integer(i4)              :: buffersize

    ! append or not
    if (present(append)) then
       if (append) then
          iappend = .true.
       else
          iappend = .false.
       endif
    else
       iappend = .false.
    endif
    LargeFile = .false.
    if (present(lfs)) LargeFile = lfs
    inetcdf4 = .false.
    if (present(netcdf4)) inetcdf4 = netcdf4
    deflate = 1
    if (present(deflate_level)) deflate = deflate_level

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)

    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))

       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+1)))
       call check(nf90_inquire_variable(ncid, varid(ndim+1), ndims=idim, dimids=dimid))
       if (idim /= ndim) stop "dump_netcdf_3d_dp: number of variable dimensions /= number of file variable dimensions."

       ! inquire dimensions
       do i=1, ndim
          call check(nf90_inquire_dimension(ncid, dimid(i), name, dims(i)))
          if (i < ndim) then
             if (trim(name) /= dnames(i)) stop "dump_netcdf_3d_dp: dimension name problem."
             if (dims(i) /= size(arr,i)) stop "dump_netcdf_3d_dp: variable dimension /= file variable dimension."
          else
             if (trim(name) /= 'time') stop "dump_netcdf_3d_dp: time name problem."
          endif
       enddo

       ! append
       start(:)      = 1
       counter(:)    = dims
       counter(ndim) = 1
       do i=1, size(arr,ndim)
          start(ndim) = dims(ndim) + i
          call check(nf90_put_var(ncid, varid(ndim), (/dims(ndim)+i/), (/dims(ndim)+i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,i), start, counter))
       end do
    else
       ! open file
       if (inetcdf4) then
          call check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
          ! call check(nf90_set_fill(ncid, NF90_NOFILL, old_fill_mode))
       else
          if (LargeFile) then
             call check(nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid, chunksize=buffersize))
          else
             call check(nf90_create(trim(Filename), NF90_CLOBBER, ncid, chunksize=buffersize))
          end if
       endif

       ! define dims
       dims = shape(arr)
       do i=1, ndim-1
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim)))

       ! define dim variables
       do i=1, ndim-1
          if (inetcdf4) then
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          else
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          endif
       end do
       ! define time variable
       if (inetcdf4) then
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       else
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       endif

       ! define variable
       if (inetcdf4) then
          chunksizes(1:ndim-1) = dims(1:ndim-1)
          chunksizes(ndim)     = 1
          call check(nf90_def_var(ncid, 'var', NF90_DOUBLE, dimid, varid(ndim+1), &
               chunksizes=chunksizes, shuffle=.true., deflate_level=deflate))
       else
          call check(nf90_def_var(ncid, 'var', NF90_DOUBLE, dimid, varid(ndim+1)))
       endif

       ! end define mode
       call check(nf90_enddef(ncid))

       ! write dimensions
       do i=1, ndim-1
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do

       ! write time and variable
       start(:) = 1
       counter(:) = dims
       counter(ndim) = 1
       do i=1, dims(ndim)
          start(ndim) = i
          call check(nf90_put_var(ncid, varid(ndim), (/i/), (/i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,i), start, counter))
       end do
    endif

    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_3d_dp


  subroutine dump_netcdf_4d_dp(filename, arr, append, lfs, netcdf4, deflate_level)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    real(dp),         dimension(:,:,:,:), intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file
    logical,          optional,         intent(in) :: lfs           ! netcdf3 Large File Support
    logical,          optional,         intent(in) :: netcdf4       ! netcdf4
    integer(i4),      optional,         intent(in) :: deflate_level ! compression level in netcdf4

    integer(i4),      parameter         :: ndim = 4 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim)   :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim)   :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+1) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim)   :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim)   :: counter  ! length array for write of each time step
    integer(i4),      dimension(ndim)   :: chunksizes ! Size of chunks in netcdf4 writing
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension
    logical                  :: LargeFile
    logical                  :: inetcdf4
    integer(i4)              :: deflate
    integer(i4)              :: buffersize

    ! append or not
    if (present(append)) then
       if (append) then
          iappend = .true.
       else
          iappend = .false.
       endif
    else
       iappend = .false.
    endif
    LargeFile = .false.
    if (present(lfs)) LargeFile = lfs
    inetcdf4 = .false.
    if (present(netcdf4)) inetcdf4 = netcdf4
    deflate = 1
    if (present(deflate_level)) deflate = deflate_level

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)

    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))

       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+1)))
       call check(nf90_inquire_variable(ncid, varid(ndim+1), ndims=idim, dimids=dimid))
       if (idim /= ndim) stop "dump_netcdf_4d_dp: number of variable dimensions /= number of file variable dimensions."

       ! inquire dimensions
       do i=1, ndim
          call check(nf90_inquire_dimension(ncid, dimid(i), name, dims(i)))
          if (i < ndim) then
             if (trim(name) /= dnames(i)) stop "dump_netcdf_4d_dp: dimension name problem."
             if (dims(i) /= size(arr,i)) stop "dump_netcdf_4d_dp: variable dimension /= file variable dimension."
          else
             if (trim(name) /= 'time') stop "dump_netcdf_4d_dp: time name problem."
          endif
       enddo

       ! append
       start(:)      = 1
       counter(:)    = dims
       counter(ndim) = 1
       do i=1, size(arr,ndim)
          start(ndim) = dims(ndim) + i
          call check(nf90_put_var(ncid, varid(ndim), (/dims(ndim)+i/), (/dims(ndim)+i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,:,i), start, counter))
       end do
    else
       ! open file
       if (inetcdf4) then
          call check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
       else
          if (LargeFile) then
             call check(nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid, chunksize=buffersize))
          else
             call check(nf90_create(trim(filename), NF90_CLOBBER, ncid, chunksize=buffersize))
          end if
       endif

       ! define dims
       dims = shape(arr)
       do i=1, ndim-1
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim)))

       ! define dim variables
       do i=1, ndim-1
          if (inetcdf4) then
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          else
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          endif
       end do
       ! define time variable
       if (inetcdf4) then
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       else
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       endif

       ! define variable
       if (inetcdf4) then
          chunksizes(1:ndim-1) = dims(1:ndim-1)
          chunksizes(ndim)     = 1
          call check(nf90_def_var(ncid, 'var', NF90_DOUBLE, dimid, varid(ndim+1), &
               chunksizes=chunksizes, shuffle=.true., deflate_level=deflate))
       else
          call check(nf90_def_var(ncid, 'var', NF90_DOUBLE, dimid, varid(ndim+1)))
       endif

       ! end define mode
       call check(nf90_enddef(ncid))

       ! write dimensions
       do i=1, ndim-1
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do

       ! write time and variable
       start(:) = 1
       counter(:) = dims
       counter(ndim) = 1
       do i=1, dims(ndim)
          start(ndim) = i
          call check(nf90_put_var(ncid, varid(ndim), (/i/), (/i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,:,i), start, counter))
       end do
    endif

    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_4d_dp


  subroutine dump_netcdf_5d_dp(filename, arr, append, lfs, netcdf4, deflate_level)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    real(dp),         dimension(:,:,:,:,:), intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file
    logical,          optional,         intent(in) :: lfs           ! netcdf3 Large File Support
    logical,          optional,         intent(in) :: netcdf4       ! netcdf4
    integer(i4),      optional,         intent(in) :: deflate_level ! compression level in netcdf4

    integer(i4),      parameter         :: ndim = 5 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim)   :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim)   :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+1) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim)   :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim)   :: counter  ! length array for write of each time step
    integer(i4),      dimension(ndim)   :: chunksizes ! Size of chunks in netcdf4 writing
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension
    logical                  :: LargeFile
    logical                  :: inetcdf4
    integer(i4)              :: deflate
    integer(i4)              :: buffersize

    ! append or not
    if (present(append)) then
       if (append) then
          iappend = .true.
       else
          iappend = .false.
       endif
    else
       iappend = .false.
    endif
    LargeFile = .false.
    if (present(lfs)) LargeFile = lfs
    inetcdf4 = .false.
    if (present(netcdf4)) inetcdf4 = netcdf4
    deflate = 1
    if (present(deflate_level)) deflate = deflate_level

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)

    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))

       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+1)))
       call check(nf90_inquire_variable(ncid, varid(ndim+1), ndims=idim, dimids=dimid))
       if (idim /= ndim) stop "dump_netcdf_5d_dp: number of variable dimensions /= number of file variable dimensions."

       ! inquire dimensions
       do i=1, ndim
          call check(nf90_inquire_dimension(ncid, dimid(i), name, dims(i)))
          if (i < ndim) then
             if (trim(name) /= dnames(i)) stop "dump_netcdf_5d_dp: dimension name problem."
             if (dims(i) /= size(arr,i)) stop "dump_netcdf_5d_dp: variable dimension /= file variable dimension."
          else
             if (trim(name) /= 'time') stop "dump_netcdf_5d_dp: time name problem."
          endif
       enddo

       ! append
       start(:)      = 1
       counter(:)    = dims
       counter(ndim) = 1
       do i=1, size(arr,ndim)
          start(ndim) = dims(ndim) + i
          call check(nf90_put_var(ncid, varid(ndim), (/dims(ndim)+i/), (/dims(ndim)+i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,:,:,i), start, counter))
       end do
    else
       ! open file
       if (inetcdf4) then
          call check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
       else
          if (LargeFile) then
             call check(nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid, chunksize=buffersize))
          else
             call check(nf90_create(trim(filename), NF90_CLOBBER, ncid, chunksize=buffersize))
          end if
       endif

       ! define dims
       dims = shape(arr)
       do i=1, ndim-1
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim)))

       ! define dim variables
       do i=1, ndim-1
          if (inetcdf4) then
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          else
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          endif
       end do
       ! define time variable
       if (inetcdf4) then
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       else
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       endif

       ! define variable
       if (inetcdf4) then
          chunksizes(1:ndim-1) = dims(1:ndim-1)
          chunksizes(ndim)     = 1
          call check(nf90_def_var(ncid, 'var', NF90_DOUBLE, dimid, varid(ndim+1), &
               chunksizes=chunksizes, shuffle=.true., deflate_level=deflate))
       else
          call check(nf90_def_var(ncid, 'var', NF90_DOUBLE, dimid, varid(ndim+1)))
       endif

       ! end define mode
       call check(nf90_enddef(ncid))

       ! write dimensions
       do i=1, ndim-1
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do

       ! write time and variable
       start(:) = 1
       counter(:) = dims
       counter(ndim) = 1
       do i=1, dims(ndim)
          start(ndim) = i
          call check(nf90_put_var(ncid, varid(ndim), (/i/), (/i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,:,:,i), start, counter))
       end do
    endif

    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_5d_dp


  subroutine dump_netcdf_1d_i4(filename, arr, append, lfs, netcdf4, deflate_level)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    integer(i4),         dimension(:),     intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file
    logical,          optional,         intent(in) :: lfs           ! netcdf3 Large File Support
    logical,          optional,         intent(in) :: netcdf4       ! netcdf4
    integer(i4),      optional,         intent(in) :: deflate_level ! compression level in netcdf4

    integer(i4),      parameter         :: ndim = 1 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim+1) :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim+1) :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+2) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim+1) :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim+1) :: counter  ! length array for write of each time step
    integer(i4),      dimension(ndim+1) :: chunksizes ! Size of chunks in netcdf4 writing
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension
    logical                  :: LargeFile
    logical                  :: inetcdf4
    integer(i4)              :: deflate
    integer(i4)              :: buffersize

    ! append or not
    if (present(append)) then
       if (append) then
          iappend = .true.
       else
          iappend = .false.
       endif
    else
       iappend = .false.
    endif
    LargeFile = .false.
    if (present(lfs)) LargeFile = lfs
    inetcdf4 = .false.
    if (present(netcdf4)) inetcdf4 = netcdf4
    deflate = 1
    if (present(deflate_level)) deflate = deflate_level

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)

    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))

       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim+1)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+2)))
       call check(nf90_inquire_variable(ncid, varid(ndim+2), ndims=idim, dimids=dimid))
       if (idim /= ndim+1) stop "dump_netcdf_1d_i4: number of variable dimensions /= number of file variable dimensions."

       ! inquire dimensions
       do i=1, ndim+1
          call check(nf90_inquire_dimension(ncid, dimid(i), name, dims(i)))
          if (i < ndim+1) then
             if (trim(name) /= dnames(i)) stop "dump_netcdf_1d_i4: dimension name problem."
             if (dims(i) /= size(arr,i)) stop "dump_netcdf_1d_i4: variable dimension /= file variable dimension."
          else
             if (trim(name) /= 'time') stop "dump_netcdf_1d_i4: time name problem."
          endif
       enddo

       ! append
       start(:)        = 1
       counter(:)      = dims
       counter(ndim+1) = 1
       do i=1, 1
          start(ndim+1) = dims(ndim+1) + i
          call check(nf90_put_var(ncid, varid(ndim+1), (/dims(ndim+1)+i/), (/dims(ndim+1)+i/)))
          call check(nf90_put_var(ncid, varid(ndim+2), arr, start, counter))
       end do
    else
       ! open file
       if (inetcdf4) then
          call check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
       else
          if (LargeFile) then
             call check(nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid, chunksize=buffersize))
          else
             call check(nf90_create(trim(filename), NF90_CLOBBER, ncid, chunksize=buffersize))
          end if
       endif

       ! define dims
       dims(1:ndim) = shape(arr)
       do i=1, ndim
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       dims(ndim+1) = 1
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim+1)))

       ! define dim variables
       do i=1, ndim
          if (inetcdf4) then
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          else
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          endif
       end do
       ! define time variable
       if (inetcdf4) then
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim+1), varid(ndim+1)))
       else
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim+1), varid(ndim+1)))
       endif

       ! define variable
       if (inetcdf4) then
          chunksizes(1:ndim) = dims(1:ndim)
          chunksizes(ndim+1) = 1
          call check(nf90_def_var(ncid, 'var', NF90_INT, dimid, varid(ndim+2), &
               chunksizes=chunksizes, shuffle=.true., deflate_level=deflate))
       else
          call check(nf90_def_var(ncid, 'var', NF90_INT, dimid, varid(ndim+2)))
       endif

       ! end define mode
       call check(nf90_enddef(ncid))

       ! write dimensions
       do i=1, ndim
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do

       ! write time and variable
       start(:)        = 1
       counter(:)      = dims
       counter(ndim+1) = 1
       do i=1, 1
          start(ndim+1) = i
          call check(nf90_put_var(ncid, varid(ndim+1), (/i/), (/i/)))
          call check(nf90_put_var(ncid, varid(ndim+2), arr, start, counter))
       end do
    endif

    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_1d_i4


  subroutine dump_netcdf_2d_i4(filename, arr, append, lfs, netcdf4, deflate_level)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    integer(i4),         dimension(:,:),   intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file
    logical,          optional,         intent(in) :: lfs           ! netcdf3 Large File Support
    logical,          optional,         intent(in) :: netcdf4       ! netcdf4
    integer(i4),      optional,         intent(in) :: deflate_level ! compression level in netcdf4

    integer(i4),      parameter         :: ndim = 2 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim+1) :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim+1) :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+2) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim+1) :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim+1) :: counter  ! length array for write of each time step
    integer(i4),      dimension(ndim+1) :: chunksizes ! Size of chunks in netcdf4 writing
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension
    logical                  :: LargeFile
    logical                  :: inetcdf4
    integer(i4)              :: deflate
    integer(i4)              :: buffersize

    ! append or not
    if (present(append)) then
       if (append) then
          iappend = .true.
       else
          iappend = .false.
       endif
    else
       iappend = .false.
    endif
    LargeFile = .false.
    if (present(lfs)) LargeFile = lfs
    inetcdf4 = .false.
    if (present(netcdf4)) inetcdf4 = netcdf4
    deflate = 1
    if (present(deflate_level)) deflate = deflate_level

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)

    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))

       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim+1)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+2)))
       call check(nf90_inquire_variable(ncid, varid(ndim+2), ndims=idim, dimids=dimid))
       if (idim /= ndim+1) stop "dump_netcdf_2d_i4: number of variable dimensions /= number of file variable dimensions."

       ! inquire dimensions
       do i=1, ndim+1
          call check(nf90_inquire_dimension(ncid, dimid(i), name, dims(i)))
          if (i < ndim+1) then
             if (trim(name) /= dnames(i)) stop "dump_netcdf_2d_i4: dimension name problem."
             if (dims(i) /= size(arr,i)) stop "dump_netcdf_2d_i4: variable dimension /= file variable dimension."
          else
             if (trim(name) /= 'time') stop "dump_netcdf_2d_i4: time name problem."
          endif
       enddo

       ! append
       start(:)        = 1
       counter(:)      = dims
       counter(ndim+1) = 1
       do i=1, 1
          start(ndim+1) = dims(ndim+1) + i
          call check(nf90_put_var(ncid, varid(ndim+1), (/dims(ndim+1)+i/), (/dims(ndim+1)+i/)))
          call check(nf90_put_var(ncid, varid(ndim+2), arr, start, counter))
       end do
    else
       ! open file
       if (inetcdf4) then
          call check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
       else
          if (LargeFile) then
             call check(nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid, chunksize=buffersize))
          else
             call check(nf90_create(trim(filename), NF90_CLOBBER, ncid, chunksize=buffersize))
          end if
       endif

       ! define dims
       dims(1:ndim) = shape(arr)
       do i=1, ndim
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       dims(ndim+1) = 1
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim+1)))

       ! define dim variables
       do i=1, ndim
          if (inetcdf4) then
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          else
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          endif
       end do
       ! define time variable
       if (inetcdf4) then
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim+1), varid(ndim+1)))
       else
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim+1), varid(ndim+1)))
       endif

       ! define variable
       if (inetcdf4) then
          chunksizes(1:ndim) = dims(1:ndim)
          chunksizes(ndim+1) = 1
          call check(nf90_def_var(ncid, 'var', NF90_INT, dimid, varid(ndim+2), &
               chunksizes=chunksizes, shuffle=.true., deflate_level=deflate))
       else
          call check(nf90_def_var(ncid, 'var', NF90_INT, dimid, varid(ndim+2)))
       endif

       ! end define mode
       call check(nf90_enddef(ncid))

       ! write dimensions
       do i=1, ndim
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do

       ! write time and variable
       start(:)        = 1
       counter(:)      = dims
       counter(ndim+1) = 1
       do i=1, 1
          start(ndim+1) = i
          call check(nf90_put_var(ncid, varid(ndim+1), (/i/), (/i/)))
          call check(nf90_put_var(ncid, varid(ndim+2), arr, start, counter))
       end do
    endif

    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_2d_i4


  subroutine dump_netcdf_3d_i4(filename, arr, append, lfs, netcdf4, deflate_level)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    integer(i4),      dimension(:,:,:), intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file
    logical,          optional,         intent(in) :: lfs           ! netcdf3 Large File Support
    logical,          optional,         intent(in) :: netcdf4       ! netcdf4
    integer(i4),      optional,         intent(in) :: deflate_level ! compression level in netcdf4

    integer(i4),      parameter         :: ndim = 3 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim)   :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim)   :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+1) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim)   :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim)   :: counter  ! length array for write of each time step
    integer(i4),      dimension(ndim)   :: chunksizes ! Size of chunks in netcdf4 writing
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension
    logical                  :: LargeFile
    logical                  :: inetcdf4
    integer(i4)              :: deflate
    integer(i4)              :: buffersize

    ! append or not
    if (present(append)) then
       if (append) then
          iappend = .true.
       else
          iappend = .false.
       endif
    else
       iappend = .false.
    endif
    LargeFile = .false.
    if (present(lfs)) LargeFile = lfs
    inetcdf4 = .false.
    if (present(netcdf4)) inetcdf4 = netcdf4
    deflate = 1
    if (present(deflate_level)) deflate = deflate_level

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)

    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))

       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+1)))
       call check(nf90_inquire_variable(ncid, varid(ndim+1), ndims=idim, dimids=dimid))
       if (idim /= ndim) stop "dump_netcdf_3d_i4: number of variable dimensions /= number of file variable dimensions."

       ! inquire dimensions
       do i=1, ndim
          call check(nf90_inquire_dimension(ncid, dimid(i), name, dims(i)))
          if (i < ndim) then
             if (trim(name) /= dnames(i)) stop "dump_netcdf_3d_i4: dimension name problem."
             if (dims(i) /= size(arr,i)) stop "dump_netcdf_3d_i4: variable dimension /= file variable dimension."
          else
             if (trim(name) /= 'time') stop "dump_netcdf_3d_i4: time name problem."
          endif
       enddo

       ! append
       start(:)      = 1
       counter(:)    = dims
       counter(ndim) = 1
       do i=1, size(arr,ndim)
          start(ndim) = dims(ndim) + i
          call check(nf90_put_var(ncid, varid(ndim), (/dims(ndim)+i/), (/dims(ndim)+i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,i), start, counter))
       end do
    else
       ! open file
       if (inetcdf4) then
          call check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
       else
          if (LargeFile) then
             call check(nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid, chunksize=buffersize))
          else
             call check(nf90_create(trim(filename), NF90_CLOBBER, ncid, chunksize=buffersize))
          end if
       endif

       ! define dims
       dims = shape(arr)
       do i=1, ndim-1
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim)))

       ! define dim variables
       do i=1, ndim-1
          if (inetcdf4) then
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          else
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          endif
       end do
       ! define time variable
       if (inetcdf4) then
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       else
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       endif

       ! define variable
       if (inetcdf4) then
          chunksizes(1:ndim-1) = dims(1:ndim-1)
          chunksizes(ndim)     = 1
          call check(nf90_def_var(ncid, 'var', NF90_INT, dimid, varid(ndim+1), &
               chunksizes=chunksizes, shuffle=.true., deflate_level=deflate))
       else
          call check(nf90_def_var(ncid, 'var', NF90_INT, dimid, varid(ndim+1)))
       endif

       ! end define mode
       call check(nf90_enddef(ncid))

       ! write dimensions
       do i=1, ndim-1
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do

       ! write time and variable
       start(:) = 1
       counter(:) = dims
       counter(ndim) = 1
       do i=1, dims(ndim)
          start(ndim) = i
          call check(nf90_put_var(ncid, varid(ndim), (/i/), (/i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,i), start, counter))
       end do
    endif

    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_3d_i4


  subroutine dump_netcdf_4d_i4(filename, arr, append, lfs, netcdf4, deflate_level)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    integer(i4),         dimension(:,:,:,:), intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file
    logical,          optional,         intent(in) :: lfs           ! netcdf3 Large File Support
    logical,          optional,         intent(in) :: netcdf4       ! netcdf4
    integer(i4),      optional,         intent(in) :: deflate_level ! compression level in netcdf4

    integer(i4),      parameter         :: ndim = 4 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim)   :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim)   :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+1) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim)   :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim)   :: counter  ! length array for write of each time step
    integer(i4),      dimension(ndim)   :: chunksizes ! Size of chunks in netcdf4 writing
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension
    logical                  :: LargeFile
    logical                  :: inetcdf4
    integer(i4)              :: deflate
    integer(i4)              :: buffersize

    ! append or not
    if (present(append)) then
       if (append) then
          iappend = .true.
       else
          iappend = .false.
       endif
    else
       iappend = .false.
    endif
    LargeFile = .false.
    if (present(lfs)) LargeFile = lfs
    inetcdf4 = .false.
    if (present(netcdf4)) inetcdf4 = netcdf4
    deflate = 1
    if (present(deflate_level)) deflate = deflate_level

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)

    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))

       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+1)))
       call check(nf90_inquire_variable(ncid, varid(ndim+1), ndims=idim, dimids=dimid))
       if (idim /= ndim) stop "dump_netcdf_4d_i4: number of variable dimensions /= number of file variable dimensions."

       ! inquire dimensions
       do i=1, ndim
          call check(nf90_inquire_dimension(ncid, dimid(i), name, dims(i)))
          if (i < ndim) then
             if (trim(name) /= dnames(i)) stop "dump_netcdf_4d_i4: dimension name problem."
             if (dims(i) /= size(arr,i)) stop "dump_netcdf_4d_i4: variable dimension /= file variable dimension."
          else
             if (trim(name) /= 'time') stop "dump_netcdf_4d_i4: time name problem."
          endif
       enddo

       ! append
       start(:)      = 1
       counter(:)    = dims
       counter(ndim) = 1
       do i=1, size(arr,ndim)
          start(ndim) = dims(ndim) + i
          call check(nf90_put_var(ncid, varid(ndim), (/dims(ndim)+i/), (/dims(ndim)+i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,:,i), start, counter))
       end do
    else
       ! open file
       if (inetcdf4) then
          call check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
       else
          if (LargeFile) then
             call check(nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid, chunksize=buffersize))
          else
             call check(nf90_create(trim(filename), NF90_CLOBBER, ncid, chunksize=buffersize))
          end if
       endif

       ! define dims
       dims = shape(arr)
       do i=1, ndim-1
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim)))

       ! define dim variables
       do i=1, ndim-1
          if (inetcdf4) then
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          else
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          endif
       end do
       ! define time variable
       if (inetcdf4) then
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       else
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       endif

       ! define variable
       if (inetcdf4) then
          chunksizes(1:ndim-1) = dims(1:ndim-1)
          chunksizes(ndim)     = 1
          call check(nf90_def_var(ncid, 'var', NF90_INT, dimid, varid(ndim+1), &
               chunksizes=chunksizes, shuffle=.true., deflate_level=deflate))
       else
          call check(nf90_def_var(ncid, 'var', NF90_INT, dimid, varid(ndim+1)))
       endif

       ! end define mode
       call check(nf90_enddef(ncid))

       ! write dimensions
       do i=1, ndim-1
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do

       ! write time and variable
       start(:) = 1
       counter(:) = dims
       counter(ndim) = 1
       do i=1, dims(ndim)
          start(ndim) = i
          call check(nf90_put_var(ncid, varid(ndim), (/i/), (/i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,:,i), start, counter))
       end do
    endif

    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_4d_i4


  subroutine dump_netcdf_5d_i4(filename, arr, append, lfs, netcdf4, deflate_level)

    implicit none

    character(len=*),                   intent(in) :: filename ! netcdf file name
    integer(i4),      dimension(:,:,:,:,:), intent(in) :: arr      ! input array
    logical,          optional,         intent(in) :: append   ! append to existing file
    logical,          optional,         intent(in) :: lfs           ! netcdf3 Large File Support
    logical,          optional,         intent(in) :: netcdf4       ! netcdf4
    integer(i4),      optional,         intent(in) :: deflate_level ! compression level in netcdf4

    integer(i4),      parameter         :: ndim = 5 ! Routine for ndim dimensional array
    character(len=1), dimension(4)      :: dnames   ! Common dimension names
    integer(i4),      dimension(ndim)   :: dims     ! Size of each dimension
    integer(i4),      dimension(ndim)   :: dimid    ! netcdf IDs of each dimension
    integer(i4),      dimension(ndim+1) :: varid    ! dimension variables and var id
    integer(i4),      dimension(ndim)   :: start    ! start array for write of each time step
    integer(i4),      dimension(ndim)   :: counter  ! length array for write of each time step
    integer(i4),      dimension(ndim)   :: chunksizes ! Size of chunks in netcdf4 writing
    integer(i4) :: ncid                             ! netcdf file id
    integer(i4) :: i, j
    logical                  :: iappend
    integer(i4)              :: idim    ! read dimension on append
    character(NF90_MAX_NAME) :: name    ! name of dimensions from nf90_inquire_dimension
    logical                  :: LargeFile
    logical                  :: inetcdf4
    integer(i4)              :: deflate
    integer(i4)              :: buffersize

    ! append or not
    if (present(append)) then
       if (append) then
          iappend = .true.
       else
          iappend = .false.
       endif
    else
       iappend = .false.
    endif
    LargeFile = .false.
    if (present(lfs)) LargeFile = lfs
    inetcdf4 = .false.
    if (present(netcdf4)) inetcdf4 = netcdf4
    deflate = 1
    if (present(deflate_level)) deflate = deflate_level

    ! dimension names
    dnames(1:4) = (/ 'x', 'y', 'z', 'l' /)

    if (iappend) then
       ! open file
       call check(nf90_open(trim(filename), NF90_WRITE, ncid))

       ! inquire variables time and var
       call check(nf90_inq_varid(ncid, 'time', varid(ndim)))
       call check(nf90_inq_varid(ncid, 'var',  varid(ndim+1)))
       call check(nf90_inquire_variable(ncid, varid(ndim+1), ndims=idim, dimids=dimid))
       if (idim /= ndim) stop "dump_netcdf_5d_i4: number of variable dimensions /= number of file variable dimensions."

       ! inquire dimensions
       do i=1, ndim
          call check(nf90_inquire_dimension(ncid, dimid(i), name, dims(i)))
          if (i < ndim) then
             if (trim(name) /= dnames(i)) stop "dump_netcdf_5d_i4: dimension name problem."
             if (dims(i) /= size(arr,i)) stop "dump_netcdf_5d_i4: variable dimension /= file variable dimension."
          else
             if (trim(name) /= 'time') stop "dump_netcdf_5d_i4: time name problem."
          endif
       enddo

       ! append
       start(:)      = 1
       counter(:)    = dims
       counter(ndim) = 1
       do i=1, size(arr,ndim)
          start(ndim) = dims(ndim) + i
          call check(nf90_put_var(ncid, varid(ndim), (/dims(ndim)+i/), (/dims(ndim)+i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,:,:,i), start, counter))
       end do
    else
       ! open file
       if (inetcdf4) then
          call check(nf90_create(trim(filename), NF90_NETCDF4, ncid))
       else
          if (LargeFile) then
             call check(nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid, chunksize=buffersize))
          else
             call check(nf90_create(trim(filename), NF90_CLOBBER, ncid, chunksize=buffersize))
          end if
       endif

       ! define dims
       dims = shape(arr)
       do i=1, ndim-1
          call check(nf90_def_dim(ncid, dnames(i), dims(i), dimid(i)))
       end do
       ! define dim time
       call check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid(ndim)))

       ! define dim variables
       do i=1, ndim-1
          if (inetcdf4) then
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          else
             call check(nf90_def_var(ncid, dnames(i), NF90_INT, dimid(i), varid(i)))
          endif
       end do
       ! define time variable
       if (inetcdf4) then
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       else
          call check(nf90_def_var(ncid, 'time', NF90_INT, dimid(ndim), varid(ndim)))
       endif

       ! define variable
       if (inetcdf4) then
          chunksizes(1:ndim-1) = dims(1:ndim-1)
          chunksizes(ndim)     = 1
          call check(nf90_def_var(ncid, 'var', NF90_INT, dimid, varid(ndim+1), &
               chunksizes=chunksizes, shuffle=.true., deflate_level=deflate))
       else
          call check(nf90_def_var(ncid, 'var', NF90_INT, dimid, varid(ndim+1)))
       endif

       ! end define mode
       call check(nf90_enddef(ncid))

       ! write dimensions
       do i=1, ndim-1
          call check(nf90_put_var(ncid, varid(i), (/ (j, j=1,dims(i)) /)))
       end do

       ! write time and variable
       start(:) = 1
       counter(:) = dims
       counter(ndim) = 1
       do i=1, dims(ndim)
          start(ndim) = i
          call check(nf90_put_var(ncid, varid(ndim), (/i/), (/i/)))
          call check(nf90_put_var(ncid, varid(ndim+1), arr(:,:,:,:,i), start, counter))
       end do
    endif

    ! close netcdf file
    call check(nf90_close(ncid))

  end subroutine dump_netcdf_5d_i4


  ! ------------------------------------------------------------------


  subroutine var2nc_1d_i4( f_name, arr, dnames, v_name, dim_unlimited, &
       long_name, units, missing_value, attributes, create, ncid, nrec )
    !
    implicit none
    !
    integer(i4), parameter :: ndim_const = 1
    integer(i4)            :: ndim
    ! input variables
    character(len=*),                           intent(in) :: f_name
    integer(i4),      dimension(:),             intent(in) :: arr
    character(len=*),                           intent(in) :: v_name
    character(len=*), dimension(:),             intent(in) :: dnames
    ! attributes
    integer(i4),                      optional, intent(in) :: dim_unlimited
    character(len=*),                 optional, intent(in) :: long_name
    character(len=*),                 optional, intent(in) :: units
    character(256),   dimension(:,:), optional, intent(in) :: attributes
    integer(i4),                      optional, intent(in) :: missing_value
    logical,                          optional, intent(in) :: create
    integer(i4),                      optional, intent(inout) :: ncid
    integer(i4),                      optional, intent(in) :: nrec
    ! local variables
    logical                                :: create_loc
    character(256)                         :: dummy_name
    integer(i4)                            :: deflate
    integer(i4), dimension(:), allocatable :: chunksizes
    integer(i4), dimension(:), allocatable :: start     ! start array for write
    integer(i4), dimension(:), allocatable :: counter   ! length array for write
    integer(i4)                            :: idim      ! read dimension on append
    integer(i4)                            :: d_unlimit ! index of unlimited dimension
    integer(i4)                            :: u_dimid
    integer(i4)                            :: u_len     ! length of unlimited dimension
    integer(i4)                            :: f_handle
    integer(i4), dimension(:), allocatable :: dims
    integer(i4), dimension(:), allocatable :: dimid     ! netcdf IDs of each dimension
    integer(i4), dimension(:), allocatable :: varid     ! dimension variables and var id
    integer(i4)                            :: i         ! loop indices
    integer(i4), dimension(:), allocatable :: dummy_count
    integer(i4), dimension(1)              :: dummy     ! dummy read
    logical                                :: openfile  ! tmp logical
    !
    ndim = size(dnames, 1)
    ! consistency checks
    d_unlimit    = 0_i4
    if ( present( dim_unlimited ) ) d_unlimit = dim_unlimited
    if ( ( size( dnames, 1) .eq. ndim_const + 1 ) .and. ( d_unlimit .ne. ndim_const+1 ) ) then
       print *, '***ERROR one more dimension name specified than dimension of array, but the last one is not unlimited'
       stop '***ERROR see StdOut'
    end if
    if ( size( dnames, 1) .gt. ndim_const + 1 ) then
       print *, '***ERROR too many dimension name specified, should be atmost ndim_const + 1'
       stop '***ERROR see StdOut'
    end if
    allocate( chunksizes(  ndim ) )
    allocate( start(       ndim ) )
    allocate( counter(     ndim ) )
    allocate( dims(        ndim ) )
    allocate( dimid(       ndim ) )
    allocate( varid(  1 +  ndim ) )
    allocate( dummy_count( ndim ) )
    ! initialize
    deflate      = 1
    if ( ndim .gt. ndim_const ) then
       chunksizes     = (/ size( arr, 1), 1 /)
       dims(1:ndim-1) = shape( arr )
       dims(ndim)     = 1
    else
       chunksizes   = (/ size( arr, 1) /)
       dims(1:ndim_const) = shape( arr )
    end if
    start(:)     = 1
    counter(:)   = dims
    dummy        = nf90_fill_int
    dummy_count  = 1
    ! open the netcdf file
    if (present(ncid)) then
       if (ncid < 0_i4) then
          openfile = .true.
       else
          openfile = .false.
          f_handle = ncid
       endif
    else
       openfile = .true.
    endif
    if (openfile) then
       create_loc = .false.
       if (present(create)) create_loc = create
       f_handle = open_netcdf(f_name, create=create_loc)
    endif
    ! check whether variable exists
    if ( nf90_noerr .eq. nf90_inq_varid( f_handle, v_name, varid(ndim+1)) ) then
       ! append
       call check(nf90_inquire_variable(f_handle, varid(ndim+1), ndims=idim, dimids=dimid))
       if (idim .ne. ndim) stop "var2nc_1d_i4: number of variable dimensions /= number of file variable dimensions."
       ! check unlimited dimension
       call check(nf90_inquire( f_handle, unlimitedDimId = u_dimid ))
       if ( u_dimid .eq. -1 ) stop 'var2nc_1d_i4: cannot append, no unlimited dimension defined'
       ! check for unlimited dimension
       if ( dimid( d_unlimit ) .ne. u_dimid ) stop 'var2nc_1d_i4: unlimited dimension not specified correctly'
       if (present(nrec)) then
          start(d_unlimit) = nrec
       else
          ! get length of unlimited dimension
          call check(nf90_inquire_dimension( f_handle, u_dimid, len = u_len ) )
          ! adapt start, that is find last written chunk
          do i = u_len, 1, -1
             if ( dummy(1) /= nf90_fill_int ) exit
             start(d_unlimit) = i
             call check( nf90_get_var( f_handle, varid(ndim+1), dummy, start, dummy_count ) )
          end do
          start(d_unlimit) = start(d_unlimit) + 1
       endif
    else
       ! define dimension
       do i = 1, ndim
          ! check whether dimension exists
          if ( nf90_noerr .ne. nf90_inq_dimid( f_handle, dnames(i), dimid(i)) ) then
             ! create dimension
             if ( i .eq. d_unlimit ) then
                ! define unlimited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), NF90_UNLIMITED, dimid(i) ))
             else
                ! define limited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), dims(i), dimid(i) ))
             end if
          end if
       end do
       ! define variable
       call check(nf90_def_var(f_handle, v_name, NF90_INT, dimid, varid(ndim+1),&
            chunksizes=chunksizes(:), shuffle=.true., deflate_level=deflate))
       !
       if ( present( long_name   ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'long_name', long_name ) )
       if ( present( units      ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'units', units ) )
       if ( present( missing_value ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'missing_value', missing_value ) )
       if ( present( attributes ) ) then
          do i = 1, size( attributes, dim = 1 )
             if ( trim( attributes(i,1)) .eq. 'missing_value' ) then
                ! write number
                read(attributes(i,2),'(I6)') dummy(1)
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), dummy(1) ) )
             else
                ! write string
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), trim(attributes(i,2)) ) )
             end if
          end do
       end if
       ! end definition
       call check(nf90_enddef(f_handle))
    end if
    ! inquire dimensions
    do i=1, ndim_const
       call check(nf90_inquire_dimension(f_handle, dimid(i), dummy_name, dims(i)))
       if (trim(dummy_name) .ne. dnames(i)) &
            stop "var2nc_1d_i4: dimension name problem."
       if ( (dims(i) .ne. size(arr,i)) .and. ( d_unlimit .ne. i ) ) &
            stop "var2nc_1d_i4: variable dimension /= file variable dimension."
    enddo
    ! write time and variable
    call check(nf90_put_var(f_handle, varid(ndim+1), arr, start, counter ) )
    ! close netcdf_dataset
    if (present(ncid)) then
       if (ncid < 0_i4) ncid = f_handle
    else
       call close_netcdf( f_handle )
    endif
    !
  end subroutine var2nc_1d_i4

  subroutine var2nc_1d_sp( f_name, arr, dnames, v_name, dim_unlimited, &
       long_name, units, missing_value, attributes, create, ncid, nrec )
    !
    implicit none
    !
    integer(i4), parameter :: ndim_const = 1
    integer(i4)            :: ndim
    ! input variables
    character(len=*),                           intent(in) :: f_name
    real(sp),         dimension(:),             intent(in) :: arr
    character(len=*), dimension(:),             intent(in) :: dnames
    character(len=*),                           intent(in) :: v_name
    ! attributes
    integer(i4),                      optional, intent(in) :: dim_unlimited
    character(len=*),                 optional, intent(in) :: long_name
    character(len=*),                 optional, intent(in) :: units
    character(256),   dimension(:,:), optional, intent(in) :: attributes
    real(sp),                         optional, intent(in) :: missing_value
    logical,                          optional, intent(in) :: create
    integer(i4),                      optional, intent(inout) :: ncid
    integer(i4),                      optional, intent(in) :: nrec
    ! local variables
    logical                                :: create_loc
    character(256)                         :: dummy_name
    integer(i4)                            :: deflate
    integer(i4), dimension(:), allocatable :: chunksizes
    integer(i4), dimension(:), allocatable :: start     ! start array for write
    integer(i4), dimension(:), allocatable :: counter   ! length array for write
    integer(i4)                            :: idim      ! read dimension on append
    integer(i4)                            :: d_unlimit ! index of unlimited dimension
    integer(i4)                            :: u_dimid
    integer(i4)                            :: u_len     ! length of unlimited dimension
    integer(i4)                            :: f_handle
    integer(i4), dimension(:), allocatable :: dims
    integer(i4), dimension(:), allocatable :: dimid     ! netcdf IDs of each dimension
    integer(i4), dimension(:), allocatable :: varid     ! dimension variables and var id
    integer(i4)                            :: i         ! loop indices
    integer(i4), dimension(:), allocatable :: dummy_count
    real(sp),    dimension(1)              :: dummy     ! dummy read
    logical                                :: openfile  ! tmp logical
    !
    ndim = size(dnames, 1)
    ! consistency checks
    d_unlimit    = 0_i4
    if ( present( dim_unlimited ) ) d_unlimit = dim_unlimited
    if ( ( size( dnames, 1) .eq. ndim_const + 1 ) .and. ( d_unlimit .ne. ndim_const+1 ) ) then
       print *, '***ERROR one more dimension name specified than dimension of array, but the last one is not unlimited'
       stop '***ERROR see StdOut'
    end if
    if ( size( dnames, 1) .gt. ndim_const + 1 ) then
       print *, '***ERROR too many dimension name specified, should be atmost ndim_const + 1'
       stop '***ERROR see StdOut'
    end if
    allocate( chunksizes(  ndim ) )
    allocate( start(       ndim ) )
    allocate( counter(     ndim ) )
    allocate( dims(        ndim ) )
    allocate( dimid(       ndim ) )
    allocate( varid(  1 +  ndim ) )
    allocate( dummy_count( ndim ) )
    ! initialize
    deflate      = 1
    if ( ndim .gt. ndim_const ) then
       chunksizes     = (/ size( arr, 1), 1 /)
       dims(1:ndim-1) = shape( arr )
       dims(ndim)     = 1
    else
       chunksizes   = (/ size( arr, 1) /)
       dims(1:ndim_const) = shape( arr )
    end if
    start(:)     = 1
    counter(:)   = dims
    dummy_count  = 1
    dummy        = nf90_fill_float
    ! open the netcdf file
    if (present(ncid)) then
       if (ncid < 0_i4) then
          openfile = .true.
       else
          openfile = .false.
          f_handle = ncid
       endif
    else
       openfile = .true.
    endif
    if (openfile) then
       create_loc = .false.
       if (present(create)) create_loc = create
       f_handle = open_netcdf(f_name, create=create_loc)
    endif
    ! check whether variable exists
    if ( nf90_noerr .eq. nf90_inq_varid( f_handle, v_name, varid(ndim+1)) ) then
       ! append
       call check(nf90_inquire_variable(f_handle, varid(ndim+1), ndims=idim, dimids=dimid))
       if (idim .ne. ndim) stop "var2nc_1d_sp: number of variable dimensions /= number of file variable dimensions."
       ! check unlimited dimension
       call check(nf90_inquire( f_handle, unlimitedDimId = u_dimid ))
       if ( u_dimid .eq. -1 ) stop 'var2nc_1d_sp: cannot append, no unlimited dimension defined'
       ! check for unlimited dimension
       if ( dimid( d_unlimit ) .ne. u_dimid ) stop 'var2nc_1d_sp: unlimited dimension not specified correctly'
       if (present(nrec)) then
          start(d_unlimit) = nrec
       else
          ! get length of unlimited dimension
          call check(nf90_inquire_dimension( f_handle, u_dimid, len = u_len ) )
          ! adapt start, that is find last written chunk
          do i = u_len, 1, -1
             if ( ne(dummy(1),nf90_fill_float) ) exit
             start(d_unlimit) = i
             call check( nf90_get_var( f_handle, varid(ndim+1), dummy, start, dummy_count ) )
          end do
          start(d_unlimit) = start(d_unlimit) + 1
       endif
    else
       ! define dimension
       do i = 1, ndim
          ! check whether dimension exists
          if ( nf90_noerr .ne. nf90_inq_dimid( f_handle, dnames(i), dimid(i)) ) then
             ! create dimension
             if ( i .eq. d_unlimit ) then
                ! define unlimited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), NF90_UNLIMITED, dimid(i) ))
             else
                ! define limited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), dims(i), dimid(i) ))
             end if
          end if
       end do
       ! define variable
       call check(nf90_def_var(f_handle, v_name, NF90_FLOAT, dimid, varid(ndim+1),&
            chunksizes=chunksizes(:), shuffle=.true., deflate_level=deflate))
       if ( present( long_name   ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'long_name', long_name ) )
       if ( present( units      ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'units', units ) )
       if ( present( missing_value ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'missing_value', missing_value ) )
       if ( present( attributes ) ) then
          do i = 1, size( attributes, dim = 1 )
             if ( trim( attributes(i,1)) .eq. 'missing_value' ) then
                ! write number
                read(attributes(i,2),'(F10.2)') dummy(1)
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), dummy(1) ) )
             else
                ! write string
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), trim(attributes(i,2)) ) )
             end if
          end do
       end if
       ! end definition
       call check(nf90_enddef(f_handle))
    end if
    ! inquire dimensions
    do i=1, ndim_const
       call check(nf90_inquire_dimension(f_handle, dimid(i), dummy_name, dims(i)))
       if (trim(dummy_name) .ne. dnames(i)) &
            stop "var2nc_1d_sp: dimension name problem."
       if ( (dims(i) .ne. size(arr,i)) .and. ( d_unlimit .ne. i ) ) &
            stop "var2nc_1d_sp: variable dimension /= file variable dimension."
    enddo
    ! write time and variable
    call check(nf90_put_var(f_handle, varid(ndim+1), arr, start, counter ) )
    ! close netcdf_dataset
    if (present(ncid)) then
       if (ncid < 0_i4) ncid = f_handle
    else
       call close_netcdf( f_handle )
    endif
    !
  end subroutine var2nc_1d_sp

  subroutine var2nc_1d_dp( f_name, arr, dnames, v_name, dim_unlimited, &
       long_name, units, missing_value, attributes, create, ncid, nrec )
    !
    implicit none
    !
    integer(i4), parameter :: ndim_const = 1
    integer(i4)            :: ndim
    ! input variables
    character(len=*),                           intent(in) :: f_name
    real(dp),         dimension(:),             intent(in) :: arr
    character(len=*),                           intent(in) :: v_name
    character(len=*), dimension(:),             intent(in) :: dnames
    ! attributes
    integer(i4),                      optional, intent(in) :: dim_unlimited
    character(len=*),                 optional, intent(in) :: long_name
    character(len=*),                 optional, intent(in) :: units
    character(256),   dimension(:,:), optional, intent(in) :: attributes
    real(dp),                         optional, intent(in) :: missing_value
    logical,                          optional, intent(in) :: create
    integer(i4),                      optional, intent(inout) :: ncid
    integer(i4),                      optional, intent(in) :: nrec
    ! local variables
    logical                                :: create_loc
    character(256)                         :: dummy_name
    integer(i4)                            :: deflate
    integer(i4), dimension(:), allocatable :: chunksizes
    integer(i4), dimension(:), allocatable :: start     ! start array for write
    integer(i4), dimension(:), allocatable :: counter   ! length array for write
    integer(i4)                            :: idim      ! read dimension on append
    integer(i4)                            :: d_unlimit ! index of unlimited dimension
    integer(i4)                            :: u_dimid
    integer(i4)                            :: u_len     ! length of unlimited dimension
    integer(i4)                            :: f_handle
    integer(i4), dimension(:), allocatable :: dims
    integer(i4), dimension(:), allocatable :: dimid     ! netcdf IDs of each dimension
    integer(i4), dimension(:), allocatable :: varid     ! dimension variables and var id
    integer(i4)                            :: i         ! loop indices
    integer(i4), dimension(:), allocatable :: dummy_count
    real(dp),    dimension(1)              :: dummy     ! dummy read
    logical                                :: openfile  ! tmp logical
    !
    ndim = size(dnames, 1)
    ! consistency checks
    d_unlimit    = 0_i4
    if ( present( dim_unlimited ) ) d_unlimit = dim_unlimited
    if ( ( size( dnames, 1) .eq. ndim_const + 1 ) .and. ( d_unlimit .ne. ndim_const+1 ) ) then
       print *, '***ERROR one more dimension name specified than dimension of array, but the last one is not unlimited'
       stop '***ERROR see StdOut'
    end if
    if ( size( dnames, 1) .gt. ndim_const + 1 ) then
       print *, '***ERROR too many dimension name specified, should be atmost ndim_const + 1'
       stop '***ERROR see StdOut'
    end if
    allocate( chunksizes(  ndim ) )
    allocate( start(       ndim ) )
    allocate( counter(     ndim ) )
    allocate( dims(        ndim ) )
    allocate( dimid(       ndim ) )
    allocate( varid(  1 +  ndim ) )
    allocate( dummy_count( ndim ) )
    ! initialize
    deflate      = 1
    if ( ndim .gt. ndim_const ) then
       chunksizes     = (/ size( arr, 1), 1 /)
       dims(1:ndim-1) = shape( arr )
       dims(ndim)     = 1
    else
       chunksizes   = (/ size( arr, 1) /)
       dims(1:ndim_const) = shape( arr )
    end if
    start(:)     = 1
    counter(:)   = dims
    dummy_count  = 1
    dummy        = nf90_fill_double
    ! open the netcdf file
    if (present(ncid)) then
       if (ncid < 0_i4) then
          openfile = .true.
       else
          openfile = .false.
          f_handle = ncid
       endif
    else
       openfile = .true.
    endif
    if (openfile) then
       create_loc = .false.
       if (present(create)) create_loc = create
       f_handle = open_netcdf(f_name, create=create_loc)
    endif
    ! check whether variable exists
    if ( nf90_noerr .eq. nf90_inq_varid( f_handle, v_name, varid(ndim+1)) ) then
       ! append
       call check(nf90_inquire_variable(f_handle, varid(ndim+1), ndims=idim, dimids=dimid))
       if (idim .ne. ndim) stop "var2nc_1d_dp: number of variable dimensions /= number of file variable dimensions."
       ! check unlimited dimension
       call check(nf90_inquire( f_handle, unlimitedDimId = u_dimid ))
       if ( u_dimid .eq. -1 ) stop 'var2nc_1d_dp: cannot append, no unlimited dimension defined'
       ! check for unlimited dimension
       if ( dimid( d_unlimit ) .ne. u_dimid ) stop 'var2nc_1d_dp: unlimited dimension not specified correctly'
       if (present(nrec)) then
          start(d_unlimit) = nrec
       else
          ! get length of unlimited dimension
          call check(nf90_inquire_dimension( f_handle, u_dimid, len = u_len ) )
          ! adapt start, that is find last written chunk
          do i = u_len, 1, -1
             if ( ne(dummy(1),nf90_fill_double) ) exit
             start(d_unlimit) = i
             call check( nf90_get_var( f_handle, varid(ndim+1), dummy, start, dummy_count ) )
          end do
          start(d_unlimit) = start(d_unlimit) + 1
       endif
    else
       ! define dimension
       do i = 1, ndim
          ! check whether dimension exists
          if ( nf90_noerr .ne. nf90_inq_dimid( f_handle, dnames(i), dimid(i)) ) then
             ! create dimension
             if ( i .eq. d_unlimit ) then
                ! define unlimited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), NF90_UNLIMITED, dimid(i) ))
             else
                ! define limited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), dims(i), dimid(i) ))
             end if
          end if
       end do
       ! define variable
       call check(nf90_def_var(f_handle, v_name, NF90_DOUBLE, dimid, varid(ndim+1),&
            chunksizes=chunksizes(:), shuffle=.true., deflate_level=deflate))
       if ( present( long_name   ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'long_name', long_name ) )
       if ( present( units      ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'units', units ) )
       if ( present( missing_value ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'missing_value', missing_value ) )
       if ( present( attributes ) ) then
          do i = 1, size( attributes, dim = 1 )
             if ( trim( attributes(i,1)) .eq. 'missing_value' ) then
                ! write number
                read(attributes(i,2),'(F10.2)') dummy(1)
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), dummy(1) ) )
             else
                ! write string
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), trim(attributes(i,2)) ) )
             end if
          end do
       end if
       ! end definition
       call check(nf90_enddef(f_handle))
    end if
    ! inquire dimensions
    do i=1, ndim_const
       call check(nf90_inquire_dimension(f_handle, dimid(i), dummy_name, dims(i)))
       if (trim(dummy_name) .ne. dnames(i)) &
            stop "var2nc_1d_dp: dimension name problem."
       if ( (dims(i) .ne. size(arr,i)) .and. ( d_unlimit .ne. i ) ) &
            stop "var2nc_1d_dp: variable dimension /= file variable dimension."
    enddo
    ! write time and variable
    call check(nf90_put_var(f_handle, varid(ndim+1), arr, start, counter ) )
    ! close netcdf_dataset
    if (present(ncid)) then
       if (ncid < 0_i4) ncid = f_handle
    else
       call close_netcdf( f_handle )
    endif
    !
  end subroutine var2nc_1d_dp

  subroutine var2nc_2d_i4( f_name, arr, dnames, v_name, dim_unlimited, &
       long_name, units, missing_value, attributes, create, ncid, nrec )
    !
    implicit none
    !
    integer(i4), parameter :: ndim_const = 2
    integer(i4)            :: ndim
    ! input variables
    character(len=*),                           intent(in) :: f_name
    integer(i4),      dimension(:,:),           intent(in) :: arr
    character(len=*),                           intent(in) :: v_name
    character(len=*), dimension(:),             intent(in) :: dnames
    ! optional
    integer(i4),                      optional, intent(in) :: dim_unlimited
    character(len=*),                 optional, intent(in) :: long_name
    character(len=*),                 optional, intent(in) :: units
    integer(i4),                      optional, intent(in) :: missing_value
    character(256),   dimension(:,:), optional, intent(in) :: attributes
    logical,                          optional, intent(in) :: create
    integer(i4),                      optional, intent(inout) :: ncid
    integer(i4),                      optional, intent(in) :: nrec
    ! local variables
    logical                                :: create_loc
    character(256)                         :: dummy_name
    integer(i4)                            :: deflate
    integer(i4), dimension(:), allocatable :: chunksizes
    integer(i4), dimension(:), allocatable :: start     ! start array for write
    integer(i4), dimension(:), allocatable :: counter   ! length array for write
    integer(i4)                            :: idim      ! read dimension on append
    integer(i4)                            :: f_handle
    integer(i4)                            :: d_unlimit ! index of unlimited dimension
    integer(i4)                            :: u_dimid   ! dimid of unlimited dimension
    integer(i4)                            :: u_len     ! length of unlimited dimension
    integer(i4), dimension(:), allocatable :: dims
    integer(i4), dimension(:), allocatable :: dimid     ! netcdf IDs of each dimension
    integer(i4), dimension(:), allocatable :: varid     ! dimension variables and var id
    integer(i4)                            :: i         ! loop indices
    integer(i4), dimension(:), allocatable :: dummy_count
    integer(i4), dimension(1)              :: dummy     ! dummy read
    logical                                :: openfile  ! tmp logical
    !
    ndim = size( dnames, 1 )
    ! consistency checks
    d_unlimit    = 0_i4
    if ( present( dim_unlimited ) ) d_unlimit = dim_unlimited
    if ( ( ndim .eq. ndim_const + 1 ) .and. ( d_unlimit .ne. ndim_const+1 ) ) then
       print *, '***ERROR one more dimension name specified than dimension of array, but the last one is not unlimited'
       stop '***ERROR see StdOut'
    end if
    if ( ndim .gt. ndim_const + 1 ) then
       print *, '***ERROR too many dimension name specified, should be atmost ndim_const + 1'
       stop '***ERROR see StdOut'
    end if
    if  ( ( ( ndim .eq. ndim_const ) .and. ( d_unlimit .gt. ndim_const ) ) .or. &
         ( d_unlimit .lt. 0_i4 ) ) then
       print*, '***ERROR unlimited dimension out of bounds, must be positive but not greater than number of given dimensions'
       print*, '***Dims: ', ndim, ndim_const, d_unlimit, ndim_const, d_unlimit
       stop '***ERROR see StdOut'
    end if
    !
    allocate( chunksizes(  ndim ) )
    allocate( start(       ndim ) )
    allocate( counter(     ndim ) )
    allocate( dims(        ndim ) )
    allocate( dimid(       ndim ) )
    allocate( varid(  1 +  ndim ) )
    allocate( dummy_count( ndim ) )
    ! initialize
    deflate      = 1
    if ( ndim .gt. ndim_const ) then
       chunksizes     = (/ size( arr, 1), size( arr, 2), 1 /)
       dims(1:ndim-1) = shape( arr )
       dims(ndim)     = 1
    else
       chunksizes   = (/ size( arr, 1), size( arr, 2) /)
       if ( d_unlimit .gt. 0 ) chunksizes( d_unlimit ) = 1 
       dims(1:ndim_const) = shape( arr )
    end if
    start(:)     = 1_i4
    counter(:)   = dims
    dummy_count  = 1_i4
    dummy        = nf90_fill_int
    ! open the netcdf file
    if (present(ncid)) then
       if (ncid < 0_i4) then
          openfile = .true.
       else
          openfile = .false.
          f_handle = ncid
       endif
    else
       openfile = .true.
    endif
    if (openfile) then
       create_loc = .false.
       if (present(create)) create_loc = create
       f_handle = open_netcdf(f_name, create=create_loc)
    endif
    ! check whether variable exists
    if ( nf90_noerr .eq. nf90_inq_varid( f_handle, v_name, varid(ndim+1)) ) then
       ! append
       call check(nf90_inquire_variable(f_handle, varid(ndim+1), ndims=idim, dimids=dimid))
       ! consistency checks
       if (idim .ne. ndim) stop "var2nc_2d_i4: number of variable dimensions /= number of file variable dimensions."
       ! check unlimited dimension
       call check(nf90_inquire( f_handle, unlimitedDimId = u_dimid ))
       if ( u_dimid .eq. -1 ) stop 'var2nc_2d_i4: cannot append, no unlimited dimension defined'
       ! check for unlimited dimension
       if ( dimid( d_unlimit ) .ne. u_dimid ) stop 'var2nc_2d_i4: unlimited dimension not specified correctly'
       if (present(nrec)) then
          start(d_unlimit) = nrec
       else
          ! get length of unlimited dimension
          call check(nf90_inquire_dimension( f_handle, u_dimid, len = u_len ) )
          ! adapt start, that is find last written chunk
          do i = u_len, 1, -1
             if ( dummy(1) /= nf90_fill_int ) exit
             start(d_unlimit) = i
             call check( nf90_get_var( f_handle, varid(ndim+1), dummy, start, dummy_count ) )
          end do
          start(d_unlimit) = start(d_unlimit) + 1
       endif
    else
       ! define dimensions
       do i = 1, ndim
          ! check whether dimension exists
          if ( nf90_noerr .ne. nf90_inq_dimid( f_handle, dnames(i), dimid(i)) ) then
             ! create dimension
             if ( i .eq. d_unlimit ) then
                ! define unlimited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), NF90_UNLIMITED, dimid(i) ))
             else
                ! define limited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), dims(i), dimid(i) ))
             end if
          end if
       end do
       ! define variable
       call check(nf90_def_var(f_handle, v_name, NF90_INT, dimid, varid(ndim+1),&
            chunksizes=chunksizes(:), shuffle=.true., deflate_level=deflate))
       ! add attributes
       if ( present( long_name   ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'long_name', long_name ) )
       if ( present( units      ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'units', units ) )
       if ( present( missing_value ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'missing_value', missing_value ) )
       if ( present( attributes ) ) then
          do i = 1, size( attributes, dim = 1 )
             if ( trim( attributes(i,1)) .eq. 'missing_value' ) then
                ! write number
                read(attributes(i,2),'(I6)') dummy(1)
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), dummy(1) ) )
             else
                ! write string
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), trim(attributes(i,2)) ) )
             end if
          end do
       end if
       ! end definition
       call check(nf90_enddef(f_handle))
    end if
    ! check dimensions before writing
    do i=1, ndim_const
       call check(nf90_inquire_dimension(f_handle, dimid(i), dummy_name, dims(i)))
       if (trim(dummy_name) .ne. dnames(i)) &
            stop "var2nc_2d_i4: dimension name problem."
       if ( (dims(i) .ne. size(arr,i)) .and. ( d_unlimit .ne. i ) ) &
            stop "var2nc_2d_i4: variable dimension /= file variable dimension."
    end do
    ! write variable
    call check( nf90_put_var(f_handle, varid(ndim+1), arr, start, counter) )
    ! close netcdf_dataset
    if (present(ncid)) then
       if (ncid < 0_i4) ncid = f_handle
    else
       call close_netcdf( f_handle )
    endif
    !
  end subroutine var2nc_2d_i4

  subroutine var2nc_2d_sp( f_name, arr, dnames, v_name, dim_unlimited, &
       long_name, units, missing_value, attributes, create, ncid, nrec )
    !
    implicit none
    !
    integer(i4), parameter :: ndim_const = 2
    integer(i4)            :: ndim
    ! input variables
    character(len=*),                           intent(in) :: f_name
    real(sp),         dimension(:,:),           intent(in) :: arr
    character(len=*),                           intent(in) :: v_name
    character(len=*), dimension(:),             intent(in) :: dnames
    ! optional
    integer(i4),                      optional, intent(in) :: dim_unlimited
    character(len=*),                 optional, intent(in) :: long_name
    character(len=*),                 optional, intent(in) :: units
    real(sp),                         optional, intent(in) :: missing_value
    character(256),   dimension(:,:), optional, intent(in) :: attributes
    logical,                          optional, intent(in) :: create
    integer(i4),                      optional, intent(inout) :: ncid
    integer(i4),                      optional, intent(in) :: nrec
    ! local variables
    logical                                :: create_loc
    character(256)                         :: dummy_name
    integer(i4)                            :: deflate
    integer(i4), dimension(:), allocatable :: chunksizes
    integer(i4), dimension(:), allocatable :: start     ! start array for write
    integer(i4), dimension(:), allocatable :: counter   ! length array for write
    integer(i4)                            :: idim      ! read dimension on append
    integer(i4)                            :: f_handle
    integer(i4)                            :: d_unlimit ! index of unlimited dimension
    integer(i4)                            :: u_dimid   ! dimid of unlimited dimension
    integer(i4)                            :: u_len     ! length of unlimited dimension
    integer(i4), dimension(:), allocatable :: dims
    integer(i4), dimension(:), allocatable :: dimid     ! netcdf IDs of each dimension
    integer(i4), dimension(:), allocatable :: varid     ! dimension variables and var id
    integer(i4)                            :: i         ! loop indices
    integer(i4), dimension(:), allocatable :: dummy_count
    real(sp),    dimension(1)              :: dummy     ! dummy read
    logical                                :: openfile  ! tmp logical
    !
    ndim = size( dnames, 1 )
    ! consistency checks
    d_unlimit    = 0_i4
    if ( present( dim_unlimited ) ) d_unlimit = dim_unlimited
    if ( ( ndim .eq. ndim_const + 1 ) .and. ( d_unlimit .ne. ndim_const+1 ) ) then
       print *, '***ERROR one more dimension name specified than dimension of array, but the last one is not unlimited'
       stop '***ERROR see StdOut'
    end if
    if ( ndim .gt. ndim_const + 1 ) then
       print *, '***ERROR too many dimension name specified, should be atmost ndim_const + 1'
       stop '***ERROR see StdOut'
    end if
    if  ( ( ( ndim .eq. ndim_const ) .and. ( d_unlimit .gt. ndim_const ) ) .or. &
          ( d_unlimit .lt. 0_i4 ) ) then
       print*, '***ERROR unlimited dimension out of bounds, must be positive but not greater than number of given dimensions'
       print*, '***Dims: ', ndim, ndim_const, d_unlimit, ndim_const, d_unlimit
       stop '***ERROR see StdOut'
    end if
    !
    allocate( chunksizes(  ndim ) )
    allocate( start(       ndim ) )
    allocate( counter(     ndim ) )
    allocate( dims(        ndim ) )
    allocate( dimid(       ndim ) )
    allocate( varid(  1 +  ndim ) )
    allocate( dummy_count( ndim ) )
    ! initialize
    deflate      = 1
    if ( ndim .gt. ndim_const ) then
       chunksizes     = (/ size( arr, 1), size( arr, 2), 1 /)
       dims(1:ndim-1) = shape( arr )
       dims(ndim)     = 1
    else
       chunksizes   = (/ size( arr, 1), size( arr, 2) /)
       if ( d_unlimit .gt. 0 ) chunksizes( d_unlimit ) = 1 
       dims(1:ndim_const) = shape( arr )
    end if
    start(:)     = 1_i4
    counter(:)   = dims
    dummy_count  = 1_i4
    dummy        = nf90_fill_float
    ! open the netcdf file
    if (present(ncid)) then
       if (ncid < 0_i4) then
          openfile = .true.
       else
          openfile = .false.
          f_handle = ncid
       endif
    else
       openfile = .true.
    endif
    if (openfile) then
       create_loc = .false.
       if (present(create)) create_loc = create
       f_handle = open_netcdf(f_name, create=create_loc)
    endif
    ! check whether variable exists
    if ( nf90_noerr .eq. nf90_inq_varid( f_handle, v_name, varid(ndim+1)) ) then
       ! append
       call check(nf90_inquire_variable(f_handle, varid(ndim+1), ndims=idim, dimids=dimid))
       ! consistency checks
       if (idim .ne. ndim) stop "var2nc_2d_sp: number of variable dimensions /= number of file variable dimensions."
       ! check unlimited dimension
       call check(nf90_inquire( f_handle, unlimitedDimId = u_dimid ))
       if ( u_dimid .eq. -1 ) stop 'var2nc_2d_sp: cannot append, no unlimited dimension defined'
       ! check for unlimited dimension
       if ( dimid( d_unlimit ) .ne. u_dimid ) stop 'var2nc_2d_sp: unlimited dimension not specified correctly'
       if (present(nrec)) then
          start(d_unlimit) = nrec
       else
          ! get length of unlimited dimension
          call check(nf90_inquire_dimension( f_handle, u_dimid, len = u_len ) )
          ! adapt start, that is find last written chunk
          do i = u_len, 1, -1
             if ( ne(dummy(1),nf90_fill_float) ) exit
             start(d_unlimit) = i
             call check( nf90_get_var( f_handle, varid(ndim+1), dummy, start, dummy_count ) )
          end do
          start(d_unlimit) = start(d_unlimit) + 1
       endif
    else
       ! define dimensions
       do i = 1, ndim
          ! check whether dimension exists
          if ( nf90_noerr .ne. nf90_inq_dimid( f_handle, dnames(i), dimid(i)) ) then
             ! create dimension
             if ( i .eq. d_unlimit ) then
                ! define unlimited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), NF90_UNLIMITED, dimid(i) ))
             else
                ! define limited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), dims(i), dimid(i) ))
             end if
          end if
       end do
       ! define variable
       call check(nf90_def_var(f_handle, v_name, NF90_FLOAT, dimid, varid(ndim+1),&
            chunksizes=chunksizes(:), shuffle=.true., deflate_level=deflate))
       ! add attributes
       if ( present( long_name   ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'long_name', long_name ) )
       if ( present( units      ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'units', units ) )
       if ( present( missing_value ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'missing_value', missing_value ) )
       if ( present( attributes ) ) then
          do i = 1, size( attributes, dim = 1 )
             if ( trim( attributes(i,1)) .eq. 'missing_value' ) then
                ! write number
                read(attributes(i,2),'(F10.2)') dummy(1)
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), dummy(1) ) )
             else
                ! write string
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), trim(attributes(i,2)) ) )
             end if
          end do
       end if
       ! end definition
       call check(nf90_enddef(f_handle))
    end if
    ! check dimensions before writing
    do i=1, ndim_const
       call check(nf90_inquire_dimension(f_handle, dimid(i), dummy_name, dims(i)))
       if (trim(dummy_name) .ne. dnames(i)) &
            stop "var2nc_2d_sp: dimension name problem."
       if ( (dims(i) .ne. size(arr,i)) .and. ( d_unlimit .ne. i ) ) &
            stop "var2nc_2d_sp: variable dimension /= file variable dimension."
    end do
    ! write variable
    call check( nf90_put_var(f_handle, varid(ndim+1), arr, start, counter) )
    ! close netcdf_dataset
    if (present(ncid)) then
       if (ncid < 0_i4) ncid = f_handle
    else
       call close_netcdf( f_handle )
    endif
    !
  end subroutine var2nc_2d_sp

  subroutine var2nc_2d_dp( f_name, arr, dnames, v_name, dim_unlimited, &
       long_name, units, missing_value, attributes, create, ncid, nrec )
    !
    implicit none
    !
    integer(i4), parameter :: ndim_const = 2
    integer(i4)            :: ndim
    ! input variables
    character(len=*),                           intent(in) :: f_name
    real(dp),         dimension(:,:),           intent(in) :: arr
    character(len=*),                           intent(in) :: v_name
    character(len=*), dimension(:),             intent(in) :: dnames
    ! optional
    integer(i4),                      optional, intent(in) :: dim_unlimited
    character(len=*),                 optional, intent(in) :: long_name
    character(len=*),                 optional, intent(in) :: units
    real(dp),                         optional, intent(in) :: missing_value
    character(256),   dimension(:,:), optional, intent(in) :: attributes
    logical,                          optional, intent(in) :: create
    integer(i4),                      optional, intent(inout) :: ncid
    integer(i4),                      optional, intent(in) :: nrec
    ! local variables
    logical                                :: create_loc
    character(256)                         :: dummy_name
    integer(i4)                            :: deflate
    integer(i4), dimension(:), allocatable :: chunksizes
    integer(i4), dimension(:), allocatable :: start     ! start array for write
    integer(i4), dimension(:), allocatable :: counter   ! length array for write
    integer(i4)                            :: idim      ! read dimension on append
    integer(i4)                            :: f_handle
    integer(i4)                            :: d_unlimit ! index of unlimited dimension
    integer(i4)                            :: u_dimid   ! dimid of unlimited dimension
    integer(i4)                            :: u_len     ! length of unlimited dimension
    integer(i4), dimension(:), allocatable :: dims
    integer(i4), dimension(:), allocatable :: dimid     ! netcdf IDs of each dimension
    integer(i4), dimension(:), allocatable :: varid     ! dimension variables and var id
    integer(i4)                            :: i         ! loop indices
    integer(i4), dimension(:), allocatable :: dummy_count
    real(dp),    dimension(1)              :: dummy     ! dummy read
    logical                                :: openfile  ! tmp logical
    !
    ndim = size( dnames, 1 )
    ! consistency checks
    d_unlimit    = 0_i4
    if ( present( dim_unlimited ) ) d_unlimit = dim_unlimited
    if ( ( ndim .eq. ndim_const + 1 ) .and. ( d_unlimit .ne. ndim_const+1 ) ) then
       print *, '***ERROR one more dimension name specified than dimension of array, but the last one is not unlimited'
       stop '***ERROR see StdOut'
    end if
    if ( ndim .gt. ndim_const + 1 ) then
       print *, '***ERROR too many dimension name specified, should be atmost ndim_const + 1'
       stop '***ERROR see StdOut'
    end if
    if  ( ( ( ndim .eq. ndim_const ) .and. ( d_unlimit .gt. ndim_const ) ) .or. &
          ( d_unlimit .lt. 0_i4 ) ) then
       print*, '***ERROR unlimited dimension out of bounds, must be positive but not greater than number of given dimensions'
       print*, '***Dims: ', ndim, ndim_const, d_unlimit, ndim_const, d_unlimit
       stop '***ERROR see StdOut'
    end if
    !
    allocate( chunksizes(  ndim ) )
    allocate( start(       ndim ) )
    allocate( counter(     ndim ) )
    allocate( dims(        ndim ) )
    allocate( dimid(       ndim ) )
    allocate( varid(  1 +  ndim ) )
    allocate( dummy_count( ndim ) )
    ! initialize
    deflate      = 1
    if ( ndim .gt. ndim_const ) then
       chunksizes     = (/ size( arr, 1), size( arr, 2), 1 /)
       dims(1:ndim-1) = shape( arr )
       dims(ndim)     = 1
    else
       chunksizes   = (/ size( arr, 1), size( arr, 2) /)
       if ( d_unlimit .gt. 0 ) chunksizes( d_unlimit ) = 1 
       dims(1:ndim_const) = shape( arr )
    end if
    start(:)     = 1_i4
    counter(:)   = dims
    dummy_count  = 1
    dummy        = nf90_fill_double
    ! open the netcdf file
    if (present(ncid)) then
       if (ncid < 0_i4) then
          openfile = .true.
       else
          openfile = .false.
          f_handle = ncid
       endif
    else
       openfile = .true.
    endif
    if (openfile) then
       create_loc = .false.
       if (present(create)) create_loc = create
       f_handle = open_netcdf(f_name, create=create_loc)
    endif
    ! check whether variable exists
    if ( nf90_noerr .eq. nf90_inq_varid( f_handle, v_name, varid(ndim+1)) ) then
       ! append
       call check(nf90_inquire_variable(f_handle, varid(ndim+1), ndims=idim, dimids=dimid))
       ! consistency checks
       if (idim .ne. ndim) stop "var2nc_2d_dp: number of variable dimensions /= number of file variable dimensions."
       ! check unlimited dimension
       call check(nf90_inquire( f_handle, unlimitedDimId = u_dimid ))
       if ( u_dimid .eq. -1 ) stop 'var2nc_2d_dp: cannot append, no unlimited dimension defined'
       ! check for unlimited dimension
       if ( dimid( d_unlimit ) .ne. u_dimid ) stop 'var2nc_2d_dp: unlimited dimension not specified correctly'
       if (present(nrec)) then
          start(d_unlimit) = nrec
       else
          ! get length of unlimited dimension
          call check(nf90_inquire_dimension( f_handle, u_dimid, len = u_len ) )
          ! adapt start, that is find last written chunk
          do i = u_len, 1, -1
             if ( ne(dummy(1),nf90_fill_double) ) exit
             start(d_unlimit) = i
             call check( nf90_get_var( f_handle, varid(ndim+1), dummy, start, dummy_count ) )
          end do
          start(d_unlimit) = start(d_unlimit) + 1
       endif
    else
       ! define dimensions
       do i = 1, ndim
          ! check whether dimension exists
          if ( nf90_noerr .ne. nf90_inq_dimid( f_handle, dnames(i), dimid(i)) ) then
             ! create dimension
             if ( i .eq. d_unlimit ) then
                ! define unlimited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), NF90_UNLIMITED, dimid(i) ))
             else
                ! define limited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), dims(i), dimid(i) ))
             end if
          end if
       end do
       ! define variable
       call check(nf90_def_var(f_handle, v_name, NF90_DOUBLE, dimid, varid(ndim+1), &
            chunksizes=chunksizes(:), shuffle=.true., deflate_level=deflate))
       ! add attributes
       if ( present( long_name   ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'long_name', long_name ) )
       if ( present( units      ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'units', units ) )
       if ( present( missing_value ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'missing_value',missing_value ) )
       if ( present( attributes ) ) then
          do i = 1, size( attributes, dim = 1 )
             if ( trim( attributes(i,1)) .eq. 'missing_value' ) then
                ! write number
                read(attributes(i,2),'(F10.2)') dummy(1)
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), dummy(1) ) )
             else
                ! write string
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), trim(attributes(i,2)) ) )
             end if
          end do
       end if
       ! end definition
       call check(nf90_enddef(f_handle))
    end if
    ! check dimensions before writing
    do i=1, ndim_const
       call check(nf90_inquire_dimension(f_handle, dimid(i), dummy_name, dims(i)))
       if (trim(dummy_name) .ne. dnames(i)) &
            stop "var2nc_2d_dp: dimension name problem."
       if ( (dims(i) .ne. size(arr,i)) .and. ( d_unlimit .ne. i ) ) &
            stop "var2nc_2d_dp: variable dimension /= file variable dimension."
    end do
    ! write variable
    call check( nf90_put_var(f_handle, varid(ndim+1), arr, start, counter) )
    ! close netcdf_dataset
    if (present(ncid)) then
       if (ncid < 0_i4) ncid = f_handle
    else
       call close_netcdf( f_handle )
    endif
    !
  end subroutine var2nc_2d_dp

  subroutine var2nc_3d_i4( f_name, arr, dnames, v_name, dim_unlimited, &
       long_name, units, missing_value, attributes, create, ncid, nrec )
    !
    implicit none
    !
    integer(i4), parameter :: ndim_const = 3
    integer(i4)            :: ndim
    ! input variables
    character(len=*),                           intent(in) :: f_name
    integer(i4),      dimension(:,:,:),         intent(in) :: arr
    character(len=*),                           intent(in) :: v_name
    character(len=*), dimension(:),             intent(in) :: dnames
    ! optional
    integer(i4),                      optional, intent(in) :: dim_unlimited
    character(len=*),                 optional, intent(in) :: long_name
    character(len=*),                 optional, intent(in) :: units
    integer(i4),                      optional, intent(in) :: missing_value
    character(256),   dimension(:,:), optional, intent(in) :: attributes
    logical,                          optional, intent(in) :: create
    integer(i4),                      optional, intent(inout) :: ncid
    integer(i4),                      optional, intent(in) :: nrec
    ! local variables
    logical                                :: create_loc
    character(256)                         :: dummy_name
    integer(i4)                            :: deflate
    integer(i4), dimension(:), allocatable :: chunksizes
    integer(i4), dimension(:), allocatable :: start     ! start array for write
    integer(i4), dimension(:), allocatable :: counter   ! length array for write
    integer(i4)                            :: idim      ! read dimension on append
    integer(i4)                            :: f_handle
    integer(i4)                            :: d_unlimit ! index of unlimited dimension
    integer(i4)                            :: u_dimid   ! dimid of unlimited dimension
    integer(i4)                            :: u_len     ! length of unlimited dimension
    integer(i4), dimension(:), allocatable :: dims
    integer(i4), dimension(:), allocatable :: dimid     ! netcdf IDs of each dimension
    integer(i4), dimension(:), allocatable :: varid     ! dimension variables and var id
    integer(i4)                            :: i         ! loop indices
    integer(i4), dimension(:), allocatable :: dummy_count
    integer(i4), dimension(1)              :: dummy     ! dummy read
    logical                                :: openfile  ! tmp logical
    !
    ndim = size( dnames, 1 )
    ! consistency checks
    d_unlimit    = 0_i4
    if ( present( dim_unlimited ) ) d_unlimit = dim_unlimited
    if ( ( ndim .eq. ndim_const + 1 ) .and. ( d_unlimit .ne. ndim_const+1 ) ) then
       print *, '***ERROR one more dimension name specified than dimension of array, but the last one is not unlimited'
       stop '***ERROR see StdOut'
    end if
    if ( ndim .gt. ndim_const + 1 ) then
       print *, '***ERROR too many dimension name specified, should be atmost ndim_const + 1'
       stop '***ERROR see StdOut'
    end if
    if  ( ( ( ndim .eq. ndim_const ) .and. ( d_unlimit .gt. ndim_const ) ) .or. &
          ( d_unlimit .lt. 0_i4 ) ) then
       print*, '***ERROR unlimited dimension out of bounds, must be positive but not greater than number of given dimensions'
       print*, '***Dims: ', ndim, ndim_const, d_unlimit, ndim_const, d_unlimit
       stop '***ERROR see StdOut'
    end if
    !
    allocate( chunksizes(  ndim ) )
    allocate( start(       ndim ) )
    allocate( counter(     ndim ) )
    allocate( dims(        ndim ) )
    allocate( dimid(       ndim ) )
    allocate( varid(  1 +  ndim ) )
    allocate( dummy_count( ndim ) )
    ! initialize
    deflate      = 1
    if ( ndim .gt. ndim_const ) then
       chunksizes     = (/ size( arr, 1), size( arr, 2), size( arr, 3), 1 /)
       dims(1:ndim-1) = shape( arr )
       dims(ndim)     = 1
    else
       chunksizes   = (/ size( arr, 1), size( arr, 2), size( arr, 3) /)
       if ( d_unlimit .gt. 0 ) chunksizes( d_unlimit ) = 1 
       dims(1:ndim_const) = shape( arr )
    end if
    start(:)     = 1_i4
    counter(:)   = dims
    dummy_count  = 1_i4
    dummy        = nf90_fill_int
    ! open the netcdf file
    if (present(ncid)) then
       if (ncid < 0_i4) then
          openfile = .true.
       else
          openfile = .false.
          f_handle = ncid
       endif
    else
       openfile = .true.
    endif
    if (openfile) then
       create_loc = .false.
       if (present(create)) create_loc = create
       f_handle = open_netcdf(f_name, create=create_loc)
    endif
    ! check whether variable exists
    if ( nf90_noerr .eq. nf90_inq_varid( f_handle, v_name, varid(ndim+1)) ) then
       ! append
       call check(nf90_inquire_variable(f_handle, varid(ndim+1), ndims=idim, dimids=dimid))
       ! consistency checks
       if (idim .ne. ndim) stop "var2nc_3d_i4: number of variable dimensions /= number of file variable dimensions."
       ! check unlimited dimension
       call check(nf90_inquire( f_handle, unlimitedDimId = u_dimid ))
       if ( u_dimid .eq. -1 ) stop 'var2nc_3d_i4: cannot append, no unlimited dimension defined'
       ! check for unlimited dimension
       if ( dimid( d_unlimit ) .ne. u_dimid ) stop 'var2nc_3d_i4: unlimited dimension not specified correctly'
       if (present(nrec)) then
          start(d_unlimit) = nrec
       else
          ! get length of unlimited dimension
          call check(nf90_inquire_dimension( f_handle, u_dimid, len = u_len ) )
          ! adapt start, that is find last written chunk
          do i = u_len, 1, -1
             if ( dummy(1) /= nf90_fill_int ) exit
             start(d_unlimit) = i
             call check( nf90_get_var( f_handle, varid(ndim+1), dummy, start, dummy_count ) )
          end do
          start(d_unlimit) = start(d_unlimit) + 1
       endif
    else
       ! define dimensions
       do i = 1, ndim
          ! check whether dimension exists
          if ( nf90_noerr .ne. nf90_inq_dimid( f_handle, dnames(i), dimid(i)) ) then
             ! create dimension
             if ( i .eq. d_unlimit ) then
                ! define unlimited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), NF90_UNLIMITED, dimid(i) ))
             else
                ! define limited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), dims(i), dimid(i) ))
             end if
          end if
       end do
       ! define variable
       call check(nf90_def_var(f_handle, v_name, NF90_INT, dimid, varid(ndim+1), &
            chunksizes=chunksizes(:), shuffle=.true., deflate_level=deflate))
       ! add attributes
       if ( present( long_name   ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'long_name', long_name ) )
       if ( present( units      ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'units', units ) )
       if ( present( missing_value ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'missing_value', missing_value ) )
       if ( present( attributes ) ) then
          do i = 1, size( attributes, dim = 1 )
             if ( trim( attributes(i,1)) .eq. 'missing_value' ) then
                ! write number
                read(attributes(i,2),'(I6)') dummy(1)
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), dummy(1) ) )
             else
                ! write string
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), trim(attributes(i,2)) ) )
             end if
          end do
       end if
       ! end definition
       call check(nf90_enddef(f_handle))
    end if
    ! check dimensions before writing
    do i=1, ndim_const
       call check(nf90_inquire_dimension(f_handle, dimid(i), dummy_name, dims(i)))
       if (trim(dummy_name) .ne. dnames(i)) &
            stop "var2nc_3d_i4: dimension name problem."
       if ( (dims(i) .ne. size(arr,i)) .and. ( d_unlimit .ne. i ) ) &
            stop "var2nc_3d_i4: variable dimension /= file variable dimension."
    end do
    ! write variable
    call check( nf90_put_var(f_handle, varid(ndim+1), arr, start, counter) )
    ! close netcdf_dataset
    if (present(ncid)) then
       if (ncid < 0_i4) ncid = f_handle
    else
       call close_netcdf( f_handle )
    endif
    !
  end subroutine var2nc_3d_i4

  subroutine var2nc_3d_sp( f_name, arr, dnames, v_name, dim_unlimited, &
       long_name, units, missing_value, attributes, create, ncid, nrec )
    !
    implicit none
    !
    integer(i4), parameter :: ndim_const = 3
    integer(i4)            :: ndim
    ! input variables
    character(len=*),                           intent(in) :: f_name
    real(sp),         dimension(:,:,:),         intent(in) :: arr
    character(len=*),                           intent(in) :: v_name
    character(len=*), dimension(:),             intent(in) :: dnames
    ! optional
    integer(i4),                      optional, intent(in) :: dim_unlimited
    character(len=*),                 optional, intent(in) :: long_name
    character(len=*),                 optional, intent(in) :: units
    real(sp),                         optional, intent(in) :: missing_value
    character(256),   dimension(:,:), optional, intent(in) :: attributes
    logical,                          optional, intent(in) :: create
    integer(i4),                      optional, intent(inout) :: ncid
    integer(i4),                      optional, intent(in) :: nrec
    ! local variables
    logical                                :: create_loc
    character(256)                         :: dummy_name
    integer(i4)                            :: deflate
    integer(i4), dimension(:), allocatable :: chunksizes
    integer(i4), dimension(:), allocatable :: start     ! start array for write
    integer(i4), dimension(:), allocatable :: counter   ! length array for write
    integer(i4)                            :: idim      ! read dimension on append
    integer(i4)                            :: f_handle
    integer(i4)                            :: d_unlimit ! index of unlimited dimension
    integer(i4)                            :: u_dimid   ! dimid of unlimited dimension
    integer(i4)                            :: u_len     ! length of unlimited dimension
    integer(i4), dimension(:), allocatable :: dims
    integer(i4), dimension(:), allocatable :: dimid     ! netcdf IDs of each dimension
    integer(i4), dimension(:), allocatable :: varid     ! dimension variables and var id
    integer(i4)                            :: i         ! loop indices
    integer(i4), dimension(:), allocatable :: dummy_count
    real(sp),    dimension(1)              :: dummy     ! dummy read
    logical                                :: openfile  ! tmp logical
    !
    ndim = size( dnames, 1 )
    ! consistency checks
    d_unlimit    = 0_i4
    if ( present( dim_unlimited ) ) d_unlimit = dim_unlimited
    if ( ( ndim .eq. ndim_const + 1 ) .and. ( d_unlimit .ne. ndim_const+1 ) ) then
       print *, '***ERROR one more dimension name specified than dimension of array, but the last one is not unlimited'
       stop '***ERROR see StdOut'
    end if
    if ( ndim .gt. ndim_const + 1 ) then
       print *, '***ERROR too many dimension name specified, should be atmost ndim_const + 1'
       stop '***ERROR see StdOut'
    end if
    if  ( ( ( ndim .eq. ndim_const ) .and. ( d_unlimit .gt. ndim_const ) ) .or. &
          ( d_unlimit .lt. 0_i4 ) ) then
       print*, '***ERROR unlimited dimension out of bounds, must be positive but not greater than number of given dimensions'
       print*, '***Dims: ', ndim, ndim_const, d_unlimit, ndim_const, d_unlimit
       stop '***ERROR see StdOut'
    end if
    !
    allocate( chunksizes(  ndim ) )
    allocate( start(       ndim ) )
    allocate( counter(     ndim ) )
    allocate( dims(        ndim ) )
    allocate( dimid(       ndim ) )
    allocate( varid(  1 +  ndim ) )
    allocate( dummy_count( ndim ) )
    ! initialize
    deflate      = 1
    ! set chunk sizes and dimension names
    if ( ndim .gt. ndim_const ) then
       chunksizes     = (/ size( arr, 1), size( arr, 2), size( arr, 3), 1 /)
       dims(1:ndim-1) = shape( arr )
       dims(ndim)     = 1
    else
       chunksizes   = (/ size( arr, 1), size( arr, 2), size( arr, 3) /)
       if ( d_unlimit .gt. 0 ) chunksizes( d_unlimit ) = 1 
       dims(1:ndim_const) = shape( arr )
    end if
    start(:)     = 1_i4
    counter(:)   = dims
    dummy_count  = 1_i4
    dummy        = nf90_fill_float
    ! open the netcdf file
    if (present(ncid)) then
       if (ncid < 0_i4) then
          openfile = .true.
       else
          openfile = .false.
          f_handle = ncid
       endif
    else
       openfile = .true.
    endif
    if (openfile) then
       create_loc = .false.
       if (present(create)) create_loc = create
       f_handle = open_netcdf(f_name, create=create_loc)
    endif
    ! check whether variable exists
    if ( nf90_noerr .eq. nf90_inq_varid( f_handle, v_name, varid(ndim+1)) ) then
       ! append
       call check(nf90_inquire_variable(f_handle, varid(ndim+1), ndims=idim, dimids=dimid))
       ! consistency checks
       if (idim .ne. ndim) stop "var2nc_3d_sp: number of variable dimensions /= number of file variable dimensions."
       ! check unlimited dimension
       call check(nf90_inquire( f_handle, unlimitedDimId = u_dimid ))
       if ( u_dimid .eq. -1 ) stop 'var2nc_3d_sp: cannot append, no unlimited dimension defined'
       ! check for unlimited dimension
       if ( dimid( d_unlimit ) .ne. u_dimid ) stop 'var2nc_3d_sp: unlimited dimension not specified correctly'
       if (present(nrec)) then
          start(d_unlimit) = nrec
       else
          ! get length of unlimited dimension
          call check(nf90_inquire_dimension( f_handle, u_dimid, len = u_len ) )
          ! adapt start, that is find last written chunk
          do i = u_len, 1, -1
             if ( ne(dummy(1),nf90_fill_float) ) exit
             start(d_unlimit) = i
             call check( nf90_get_var( f_handle, varid(ndim+1), dummy, start, dummy_count ) )
          end do
          start(d_unlimit) = start(d_unlimit) + 1
       endif
    else
       ! define dimensions
       do i = 1, ndim
          ! check whether dimension exists
          if ( nf90_noerr .ne. nf90_inq_dimid( f_handle, dnames(i), dimid(i)) ) then
             ! create dimension
             if ( i .eq. d_unlimit ) then
                ! define unlimited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), NF90_UNLIMITED, dimid(i) ))
             else
                ! define limited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), dims(i), dimid(i) ))
             end if
          end if
       end do
       ! define variable
       call check(nf90_def_var(f_handle, v_name, NF90_FLOAT, dimid, varid(ndim+1), &
            chunksizes=chunksizes(:), shuffle=.true., deflate_level=deflate))
       ! add attributes
       if ( present( long_name   ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'long_name', long_name ) )
       if ( present( units      ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'units', units ) )
       if ( present( missing_value ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'missing_value', missing_value ) )
       if ( present( attributes ) ) then
          do i = 1, size( attributes, dim = 1 )
             if ( trim( attributes(i,1)) .eq. 'missing_value' ) then
                ! write number
                read(attributes(i,2),'(F10.2)') dummy(1)
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), dummy(1) ) )
             else
                ! write string
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), trim(attributes(i,2)) ) )
             end if
          end do
       end if
       ! end definition
       call check(nf90_enddef(f_handle))
    end if
    ! check dimensions before writing
    do i=1, ndim_const
       call check(nf90_inquire_dimension(f_handle, dimid(i), dummy_name, dims(i)))
       if (trim(dummy_name) .ne. dnames(i)) &
            stop "var2nc_3d_sp: dimension name problem."
       if ( (dims(i) .ne. size(arr,i)) .and. ( d_unlimit .ne. i ) ) &
            stop "var2nc_3d_sp: variable dimension /= file variable dimension."
    end do
    ! write variable
    call check( nf90_put_var(f_handle, varid(ndim+1), arr, start, counter) )
    ! close netcdf_dataset
    if (present(ncid)) then
       if (ncid < 0_i4) ncid = f_handle
    else
       call close_netcdf( f_handle )
    endif
    !
  end subroutine var2nc_3d_sp

  subroutine var2nc_3d_dp( f_name, arr, dnames, v_name, dim_unlimited, &
       long_name, units, missing_value, attributes, create, ncid, nrec )
    !
    implicit none
    !
    integer(i4), parameter :: ndim_const = 3
    integer(i4)            :: ndim
    ! input variables
    character(len=*),                           intent(in) :: f_name
    real(dp),         dimension(:,:,:),         intent(in) :: arr
    character(len=*),                           intent(in) :: v_name
    character(len=*), dimension(:),             intent(in) :: dnames
    ! optional
    integer(i4),                      optional, intent(in) :: dim_unlimited
    character(len=*),                 optional, intent(in) :: long_name
    character(len=*),                 optional, intent(in) :: units
    real(dp),                         optional, intent(in) :: missing_value
    character(256),   dimension(:,:), optional, intent(in) :: attributes
    logical,                          optional, intent(in) :: create
    integer(i4),                      optional, intent(inout) :: ncid
    integer(i4),                      optional, intent(in) :: nrec
    ! local variables
    logical                                :: create_loc
    character(256)                         :: dummy_name
    integer(i4)                            :: deflate
    integer(i4), dimension(:), allocatable :: chunksizes
    integer(i4), dimension(:), allocatable :: start     ! start array for write
    integer(i4), dimension(:), allocatable :: counter   ! length array for write
    integer(i4)                            :: idim      ! read dimension on append
    integer(i4)                            :: f_handle
    integer(i4)                            :: d_unlimit ! index of unlimited dimension
    integer(i4)                            :: u_dimid   ! dimid of unlimited dimension
    integer(i4)                            :: u_len     ! length of unlimited dimension
    integer(i4), dimension(:), allocatable :: dims
    integer(i4), dimension(:), allocatable :: dimid     ! netcdf IDs of each dimension
    integer(i4), dimension(:), allocatable :: varid     ! dimension variables and var id
    integer(i4)                            :: i         ! loop indices
    integer(i4), dimension(:), allocatable :: dummy_count
    real(dp),    dimension(1)              :: dummy     ! dummy read
    logical                                :: openfile  ! tmp logical
    !
    ndim = size( dnames, 1 )
    ! consistency checks
    d_unlimit    = 0_i4
    if ( present( dim_unlimited ) ) d_unlimit = dim_unlimited
    if ( ( ndim .eq. ndim_const + 1 ) .and. ( d_unlimit .ne. ndim_const+1 ) ) then
       print *, '***ERROR one more dimension name specified than dimension of array, but the last one is not unlimited'
       stop '***ERROR see StdOut'
    end if
    if ( ndim .gt. ndim_const + 1 ) then
       print *, '***ERROR too many dimension name specified, should be atmost ndim_const + 1'
       stop '***ERROR see StdOut'
    end if
    if  ( ( ( ndim .eq. ndim_const ) .and. ( d_unlimit .gt. ndim_const ) ) .or. &
          ( d_unlimit .lt. 0_i4 ) ) then
       print*, '***ERROR unlimited dimension out of bounds, must be positive but not greater than number of given dimensions'
       print*, '***Dims: ', ndim, ndim_const, d_unlimit, ndim_const, d_unlimit
       stop '***ERROR see StdOut'
    end if
    !
    allocate( chunksizes(  ndim ) )
    allocate( start(       ndim ) )
    allocate( counter(     ndim ) )
    allocate( dims(        ndim ) )
    allocate( dimid(       ndim ) )
    allocate( varid(  1 +  ndim ) )
    allocate( dummy_count( ndim ) )
    ! initialize
    deflate      = 1
    if ( ndim .gt. ndim_const ) then
       chunksizes     = (/ size( arr, 1), size( arr, 2), size( arr, 3), 1 /)
       dims(1:ndim-1) = shape( arr )
       dims(ndim)     = 1
    else
       chunksizes   = (/ size( arr, 1), size( arr, 2), size( arr, 3) /)
       if ( d_unlimit .gt. 0 ) chunksizes( d_unlimit ) = 1 
       dims(1:ndim_const) = shape( arr )
    end if
    start(:)     = 1_i4
    counter(:)   = dims
    dummy_count  = 1
    dummy        = nf90_fill_double
    ! open the netcdf file
    if (present(ncid)) then
       if (ncid < 0_i4) then
          openfile = .true.
       else
          openfile = .false.
          f_handle = ncid
       endif
    else
       openfile = .true.
    endif
    if (openfile) then
       create_loc = .false.
       if (present(create)) create_loc = create
       f_handle = open_netcdf(f_name, create=create_loc)
    endif
    ! check whether variable exists
    if ( nf90_noerr .eq. nf90_inq_varid( f_handle, v_name, varid(ndim+1)) ) then
       ! append
       call check(nf90_inquire_variable(f_handle, varid(ndim+1), ndims=idim, dimids=dimid))
       ! consistency checks
       if (idim .ne. ndim) stop "var2nc_3d_dp: number of variable dimensions /= number of file variable dimensions."
       ! check unlimited dimension
       call check(nf90_inquire( f_handle, unlimitedDimId = u_dimid ))
       if ( u_dimid .eq. -1 ) stop 'var2nc_3d_dp: cannot append, no unlimited dimension defined'
       ! check for unlimited dimension
       if ( dimid( d_unlimit ) .ne. u_dimid ) stop 'var2nc_3d_dp: unlimited dimension not specified correctly'
       if (present(nrec)) then
          start(d_unlimit) = nrec
       else
          ! get length of unlimited dimension
          call check(nf90_inquire_dimension( f_handle, u_dimid, len = u_len ) )
          ! adapt start, that is find last written chunk
          do i = u_len, 1, -1
             if ( ne(dummy(1),nf90_fill_double) ) exit
             start(d_unlimit) = i
             call check( nf90_get_var( f_handle, varid(ndim+1), dummy, start, dummy_count ) )
          end do
          start(d_unlimit) = start(d_unlimit) + 1
       endif
    else
       ! define dimensions
       do i = 1, ndim
          ! check whether dimension exists
          if ( nf90_noerr .ne. nf90_inq_dimid( f_handle, dnames(i), dimid(i)) ) then
             ! create dimension
             if ( i .eq. d_unlimit ) then
                ! define unlimited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), NF90_UNLIMITED, dimid(i) ))
             else
                ! define limited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), dims(i), dimid(i) ))
             end if
          end if
       end do
       ! define variable
       call check(nf90_def_var(f_handle, v_name, NF90_DOUBLE, dimid, varid(ndim+1), &
            chunksizes=chunksizes(:), shuffle=.true., deflate_level=deflate))
       ! add attributes
       if ( present( long_name   ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'long_name', long_name ) )
       if ( present( units      ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'units', units ) )
       if ( present( missing_value ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'missing_value',missing_value ) )
       if ( present( attributes ) ) then
          do i = 1, size( attributes, dim = 1 )
             if ( trim( attributes(i,1)) .eq. 'missing_value' ) then
                ! write number
                read(attributes(i,2),'(F10.2)') dummy(1)
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), dummy(1) ) )
             else
                ! write string
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), trim(attributes(i,2)) ) )
             end if
          end do
       end if
       ! end definition
       call check(nf90_enddef(f_handle))
    end if
    ! check dimensions before writing
    do i=1, ndim_const
       call check(nf90_inquire_dimension(f_handle, dimid(i), dummy_name, dims(i)))
       if (trim(dummy_name) .ne. dnames(i)) &
            stop "var2nc_3d_dp: dimension name problem."
       if ( (dims(i) .ne. size(arr,i)) .and. ( d_unlimit .ne. i ) ) &
            stop "var2nc_3d_dp: variable dimension /= file variable dimension."
    end do
    ! write variable
    call check( nf90_put_var(f_handle, varid(ndim+1), arr, start, counter) )
    ! close netcdf_dataset
    if (present(ncid)) then
       if (ncid < 0_i4) ncid = f_handle
    else
       call close_netcdf( f_handle )
    endif
    !
  end subroutine var2nc_3d_dp

  subroutine var2nc_4d_i4( f_name, arr, dnames, v_name, dim_unlimited, &
       long_name, units, missing_value, attributes, create, ncid, nrec )
    !
    implicit none
    !
    integer(i4), parameter :: ndim_const = 4
    integer(i4)            :: ndim
    ! input variables
    character(len=*),                           intent(in) :: f_name
    integer(i4),      dimension(:,:,:,:),       intent(in) :: arr
    character(len=*),                           intent(in) :: v_name
    character(len=*), dimension(:),             intent(in) :: dnames
    ! optional
    integer(i4),                      optional, intent(in) :: dim_unlimited
    character(len=*),                 optional, intent(in) :: long_name
    character(len=*),                 optional, intent(in) :: units
    integer(i4),                      optional, intent(in) :: missing_value
    character(256),   dimension(:,:), optional, intent(in) :: attributes
    logical,                          optional, intent(in) :: create
    integer(i4),                      optional, intent(inout) :: ncid
    integer(i4),                      optional, intent(in) :: nrec
    ! local variables
    logical                                :: create_loc
    character(256)                         :: dummy_name
    integer(i4)                            :: deflate
    integer(i4), dimension(:), allocatable :: chunksizes
    integer(i4), dimension(:), allocatable :: start     ! start array for write
    integer(i4), dimension(:), allocatable :: counter   ! length array for write
    integer(i4)                            :: idim      ! read dimension on append
    integer(i4)                            :: f_handle
    integer(i4)                            :: d_unlimit ! index of unlimited dimension
    integer(i4)                            :: u_dimid   ! dimid of unlimited dimension
    integer(i4)                            :: u_len     ! length of unlimited dimension
    integer(i4), dimension(:), allocatable :: dims
    integer(i4), dimension(:), allocatable :: dimid     ! netcdf IDs of each dimension
    integer(i4), dimension(:), allocatable :: varid     ! dimension variables and var id
    integer(i4)                            :: i         ! loop indices
    integer(i4), dimension(:), allocatable :: dummy_count
    integer(i4), dimension(1)              :: dummy     ! dummy read
    logical                                :: openfile  ! tmp logical
    !
    ndim = size( dnames, 1 )
    ! consistency checks
    d_unlimit    = 0_i4
    if ( present( dim_unlimited ) ) d_unlimit = dim_unlimited
    if ( ( ndim .eq. ndim_const + 1 ) .and. ( d_unlimit .ne. ndim_const+1 ) ) then
       print *, '***ERROR one more dimension name specified than dimension of array, but the last one is not unlimited'
       stop '***ERROR see StdOut'
    end if
    if ( ndim .gt. ndim_const + 1 ) then
       print *, '***ERROR too many dimension name specified, should be atmost ndim_const + 1'
       stop '***ERROR see StdOut'
    end if
    if  ( ( ( ndim .eq. ndim_const ) .and. ( d_unlimit .gt. ndim_const ) ) .or. &
          ( d_unlimit .lt. 0_i4 ) ) then
       print*, '***ERROR unlimited dimension out of bounds, must be positive but not greater than number of given dimensions'
       print*, '***Dims: ', ndim, ndim_const, d_unlimit, ndim_const, d_unlimit
       stop '***ERROR see StdOut'
    end if
    !
    allocate( chunksizes(  ndim ) )
    allocate( start(       ndim ) )
    allocate( counter(     ndim ) )
    allocate( dims(        ndim ) )
    allocate( dimid(       ndim ) )
    allocate( varid(  1 +  ndim ) )
    allocate( dummy_count( ndim ) )
    ! initialize
    deflate      = 1
    if ( ndim .gt. ndim_const ) then
       chunksizes     = (/ size( arr, 1), size( arr, 2), &
            size( arr, 3), size( arr, 4 ),  1 /)
       dims(1:ndim-1) = shape( arr )
       dims(ndim)     = 1
    else
       chunksizes   = (/ size( arr, 1), size( arr, 2), &
            size( arr, 3), size( arr, 4) /)
       if ( d_unlimit .gt. 0 ) chunksizes( d_unlimit ) = 1 
       dims(1:ndim_const) = shape( arr )
    end if
    start(:)     = 1_i4
    counter(:)   = dims
    dummy_count  = 1_i4
    dummy        = nf90_fill_int
    ! open the netcdf file
    if (present(ncid)) then
       if (ncid < 0_i4) then
          openfile = .true.
       else
          openfile = .false.
          f_handle = ncid
       endif
    else
       openfile = .true.
    endif
    if (openfile) then
       create_loc = .false.
       if (present(create)) create_loc = create
       f_handle = open_netcdf(f_name, create=create_loc)
    endif
    ! check whether variable exists
    if ( nf90_noerr .eq. nf90_inq_varid( f_handle, v_name, varid(ndim+1)) ) then
       ! append
       call check(nf90_inquire_variable(f_handle, varid(ndim+1), ndims=idim, dimids=dimid))
       ! consistency checks
       if (idim .ne. ndim) stop "var2nc_4d_i4: number of variable dimensions /= number of file variable dimensions."
       ! check unlimited dimension
       call check(nf90_inquire( f_handle, unlimitedDimId = u_dimid ))
       if ( u_dimid .eq. -1 ) stop 'var2nc_4d_i4: cannot append, no unlimited dimension defined'
       ! check for unlimited dimension
       if ( dimid( d_unlimit ) .ne. u_dimid ) stop 'var2nc_4d_sp: unlimited dimension not specified correctly'
       if (present(nrec)) then
          start(d_unlimit) = nrec
       else
          ! get length of unlimited dimension
          call check(nf90_inquire_dimension( f_handle, u_dimid, len = u_len ) )
          ! adapt start, that is find last written chunk
          do i = u_len, 1, -1
             if ( dummy(1) /= nf90_fill_int ) exit
             start(d_unlimit) = i
             call check( nf90_get_var( f_handle, varid(ndim+1), dummy, start, dummy_count ) )
          end do
          start(d_unlimit) = start(d_unlimit) + 1
       endif
    else
       ! define dimensions
       do i = 1, ndim
          ! check whether dimension exists
          if ( nf90_noerr .ne. nf90_inq_dimid( f_handle, dnames(i), dimid(i)) ) then
             ! create dimension
             if ( i .eq. d_unlimit ) then
                ! define unlimited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), NF90_UNLIMITED, dimid(i) ))
             else
                ! define limited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), dims(i), dimid(i) ))
             end if
          end if
       end do
       ! define variable
       call check(nf90_def_var(f_handle, v_name, NF90_INT, dimid, varid(ndim+1), &
            chunksizes=chunksizes(:), shuffle=.true., deflate_level=deflate))
       ! add attributes
       if ( present( long_name   ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'long_name', long_name ) )
       if ( present( units      ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'units', units ) )
       if ( present( missing_value ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'missing_value', missing_value ) )
       if ( present( attributes ) ) then
          do i = 1, size( attributes, dim = 1 )
             if ( trim( attributes(i,1)) .eq. 'missing_value' ) then
                ! write number
                read(attributes(i,2),'(I6)') dummy(1)
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), dummy(1) ) )
             else
                ! write string
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), trim(attributes(i,2)) ) )
             end if
          end do
       end if
       ! end definition
       call check(nf90_enddef(f_handle))
    end if
    ! check dimensions before writing
    do i=1, ndim_const
       call check(nf90_inquire_dimension(f_handle, dimid(i), dummy_name, dims(i)))
       if (trim(dummy_name) .ne. dnames(i)) &
            stop "var2nc_4d_i4: dimension name problem."
       if ( (dims(i) .ne. size(arr,i)) .and. ( d_unlimit .ne. i ) ) &
            stop "var2nc_4d_i4: variable dimension /= file variable dimension."
    end do
    ! write variable
    call check( nf90_put_var(f_handle, varid(ndim+1), arr, start, counter) )
    ! close netcdf_dataset
    if (present(ncid)) then
       if (ncid < 0_i4) ncid = f_handle
    else
       call close_netcdf( f_handle )
    endif
    !
  end subroutine var2nc_4d_i4

  subroutine var2nc_4d_sp( f_name, arr, dnames, v_name, dim_unlimited, &
       long_name, units, missing_value, attributes, create, ncid, nrec )
    !
    implicit none
    !
    integer(i4), parameter :: ndim_const = 4
    integer(i4)            :: ndim
    ! input variables
    character(len=*),                           intent(in) :: f_name
    real(sp),         dimension(:,:,:,:),       intent(in) :: arr
    character(len=*),                           intent(in) :: v_name
    character(len=*), dimension(:),             intent(in) :: dnames
    ! optional
    integer(i4),                      optional, intent(in) :: dim_unlimited
    character(len=*),                 optional, intent(in) :: long_name
    character(len=*),                 optional, intent(in) :: units
    real(sp),                         optional, intent(in) :: missing_value
    character(256),   dimension(:,:), optional, intent(in) :: attributes
    logical,                          optional, intent(in) :: create
    integer(i4),                      optional, intent(inout) :: ncid
    integer(i4),                      optional, intent(in) :: nrec
    ! local variables
    logical                                :: create_loc
    character(256)                         :: dummy_name
    integer(i4)                            :: deflate
    integer(i4), dimension(:), allocatable :: chunksizes
    integer(i4), dimension(:), allocatable :: start     ! start array for write
    integer(i4), dimension(:), allocatable :: counter   ! length array for write
    integer(i4)                            :: idim      ! read dimension on append
    integer(i4)                            :: f_handle
    integer(i4)                            :: d_unlimit ! index of unlimited dimension
    integer(i4)                            :: u_dimid   ! dimid of unlimited dimension
    integer(i4)                            :: u_len     ! length of unlimited dimension
    integer(i4), dimension(:), allocatable :: dims
    integer(i4), dimension(:), allocatable :: dimid     ! netcdf IDs of each dimension
    integer(i4), dimension(:), allocatable :: varid     ! dimension variables and var id
    integer(i4)                            :: i         ! loop indices
    integer(i4), dimension(:), allocatable :: dummy_count
    real(sp),    dimension(1)              :: dummy     ! dummy read
    logical                                :: openfile  ! tmp logical
    !
    ndim = size( dnames, 1 )
    ! consistency checks
    d_unlimit    = 0_i4
    if ( present( dim_unlimited ) ) d_unlimit = dim_unlimited
    if ( ( ndim .eq. ndim_const + 1 ) .and. ( d_unlimit .ne. ndim_const+1 ) ) then
       print *, '***ERROR one more dimension name specified than dimension of array, but the last one is not unlimited'
       stop '***ERROR see StdOut'
    end if
    if ( ndim .gt. ndim_const + 1 ) then
       print *, '***ERROR too many dimension name specified, should be atmost ndim_const + 1'
       stop '***ERROR see StdOut'
    end if
    if  ( ( ( ndim .eq. ndim_const ) .and. ( d_unlimit .gt. ndim_const ) ) .or. &
          ( d_unlimit .lt. 0_i4 ) ) then
       print*, '***ERROR unlimited dimension out of bounds, must be positive but not greater than number of given dimensions'
       print*, '***Dims: ', ndim, ndim_const, d_unlimit, ndim_const, d_unlimit
       stop '***ERROR see StdOut'
    end if
    !
    allocate( chunksizes(  ndim ) )
    allocate( start(       ndim ) )
    allocate( counter(     ndim ) )
    allocate( dims(        ndim ) )
    allocate( dimid(       ndim ) )
    allocate( varid(  1 +  ndim ) )
    allocate( dummy_count( ndim ) )
    ! initialize
    deflate      = 1
    if ( ndim .gt. ndim_const ) then
       chunksizes     = (/ size( arr, 1), size( arr, 2), &
            size( arr, 3), size( arr, 4 ),  1 /)
       dims(1:ndim-1) = shape( arr )
       dims(ndim)     = 1
    else
       chunksizes   = (/ size( arr, 1), size( arr, 2), &
            size( arr, 3), size( arr, 4) /)
       if ( d_unlimit .gt. 0 ) chunksizes( d_unlimit ) = 1 
       dims(1:ndim_const) = shape( arr )
    end if
    start(:)     = 1_i4
    counter(:)   = dims
    dummy_count  = 1_i4
    dummy        = nf90_fill_float
    ! open the netcdf file
    if (present(ncid)) then
       if (ncid < 0_i4) then
          openfile = .true.
       else
          openfile = .false.
          f_handle = ncid
       endif
    else
       openfile = .true.
    endif
    if (openfile) then
       create_loc = .false.
       if (present(create)) create_loc = create
       f_handle = open_netcdf(f_name, create=create_loc)
    endif
    ! check whether variable exists
    if ( nf90_noerr .eq. nf90_inq_varid( f_handle, v_name, varid(ndim+1)) ) then
       ! append
       call check(nf90_inquire_variable(f_handle, varid(ndim+1), ndims=idim, dimids=dimid))
       ! consistency checks
       if (idim .ne. ndim) stop "var2nc_4d_sp: number of variable dimensions /= number of file variable dimensions."
       ! check unlimited dimension
       call check(nf90_inquire( f_handle, unlimitedDimId = u_dimid ))
       if ( u_dimid .eq. -1 ) stop 'var2nc_4d_sp: cannot append, no unlimited dimension defined'
       ! check for unlimited dimension
       if ( dimid( d_unlimit ) .ne. u_dimid ) stop 'var2nc_4d_sp: unlimited dimension not specified correctly'
       if (present(nrec)) then
          start(d_unlimit) = nrec
       else
          ! get length of unlimited dimension
          call check(nf90_inquire_dimension( f_handle, u_dimid, len = u_len ) )
          ! adapt start, that is find last written chunk
          do i = u_len, 1, -1
             if ( ne(dummy(1),nf90_fill_float) ) exit
             start(d_unlimit) = i
             call check( nf90_get_var( f_handle, varid(ndim+1), dummy, start, dummy_count ) )
          end do
          start(d_unlimit) = start(d_unlimit) + 1
       endif
    else
       ! define dimensions
       do i = 1, ndim
          ! check whether dimension exists
          if ( nf90_noerr .ne. nf90_inq_dimid( f_handle, dnames(i), dimid(i)) ) then
             ! create dimension
             if ( i .eq. d_unlimit ) then
                ! define unlimited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), NF90_UNLIMITED, dimid(i) ))
             else
                ! define limited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), dims(i), dimid(i) ))
             end if
          end if
       end do
       ! define variable
       call check(nf90_def_var(f_handle, v_name, NF90_FLOAT, dimid, varid(ndim+1), &
            chunksizes=chunksizes(:), shuffle=.true., deflate_level=deflate))
       ! add attributes
       if ( present( long_name   ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'long_name', long_name ) )
       if ( present( units      ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'units', units ) )
       if ( present( missing_value ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'missing_value', missing_value ) )
       if ( present( attributes ) ) then
          do i = 1, size( attributes, dim = 1 )
             if ( trim( attributes(i,1)) .eq. 'missing_value' ) then
                ! write number
                read(attributes(i,2),'(F10.2)') dummy(1)
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), dummy(1) ) )
             else
                ! write string
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), trim(attributes(i,2)) ) )
             end if
          end do
       end if
       ! end definition
       call check(nf90_enddef(f_handle))
    end if
    ! check dimensions before writing
    do i=1, ndim_const
       call check(nf90_inquire_dimension(f_handle, dimid(i), dummy_name, dims(i)))
       if (trim(dummy_name) .ne. dnames(i)) &
            stop "var2nc_4d_sp: dimension name problem."
       if ( (dims(i) .ne. size(arr,i)) .and. ( d_unlimit .ne. i ) ) &
            stop "var2nc_4d_sp: variable dimension /= file variable dimension."
    end do
    ! write variable
    call check( nf90_put_var(f_handle, varid(ndim+1), arr, start, counter) )
    ! close netcdf_dataset
    if (present(ncid)) then
       if (ncid < 0_i4) ncid = f_handle
    else
       call close_netcdf( f_handle )
    endif
    !
  end subroutine var2nc_4d_sp

  subroutine var2nc_4d_dp( f_name, arr, dnames, v_name, dim_unlimited, &
       long_name, units, missing_value, attributes, create, ncid, nrec )
    !
    implicit none
    !
    integer(i4), parameter :: ndim_const = 4
    integer(i4)            :: ndim
    ! input variables
    character(len=*),                           intent(in) :: f_name
    real(dp),         dimension(:,:,:,:),       intent(in) :: arr
    character(len=*),                           intent(in) :: v_name
    character(len=*), dimension(:),             intent(in) :: dnames
    ! optional
    integer(i4),                      optional, intent(in) :: dim_unlimited
    character(len=*),                 optional, intent(in) :: long_name
    character(len=*),                 optional, intent(in) :: units
    real(dp),                         optional, intent(in) :: missing_value
    character(256),   dimension(:,:), optional, intent(in) :: attributes
    logical,                          optional, intent(in) :: create
    integer(i4),                      optional, intent(inout) :: ncid
    integer(i4),                      optional, intent(in) :: nrec
    ! local variables
    logical                                :: create_loc
    character(256)                         :: dummy_name
    integer(i4)                            :: deflate
    integer(i4), dimension(:), allocatable :: chunksizes
    integer(i4), dimension(:), allocatable :: start     ! start array for write
    integer(i4), dimension(:), allocatable :: counter   ! length array for write
    integer(i4)                            :: idim      ! read dimension on append
    integer(i4)                            :: f_handle
    integer(i4)                            :: d_unlimit ! index of unlimited dimension
    integer(i4)                            :: u_dimid   ! dimid of unlimited dimension
    integer(i4)                            :: u_len     ! length of unlimited dimension
    integer(i4), dimension(:), allocatable :: dims
    integer(i4), dimension(:), allocatable :: dimid     ! netcdf IDs of each dimension
    integer(i4), dimension(:), allocatable :: varid     ! dimension variables and var id
    integer(i4)                            :: i         ! loop indices
    integer(i4), dimension(:), allocatable :: dummy_count
    real(dp),    dimension(1)              :: dummy     ! dummy read
    logical                                :: openfile  ! tmp logical
    !
    ndim = size( dnames, 1 )
    ! consistency checks
    d_unlimit    = 0_i4
    if ( present( dim_unlimited ) ) d_unlimit = dim_unlimited
    if ( ( ndim .eq. ndim_const + 1 ) .and. ( d_unlimit .ne. ndim_const+1 ) ) then
       print *, '***ERROR one more dimension name specified than dimension of array, but the last one is not unlimited'
       stop '***ERROR see StdOut'
    end if
    if ( ndim .gt. ndim_const + 1 ) then
       print *, '***ERROR too many dimension name specified, should be atmost ndim_const + 1'
       stop '***ERROR see StdOut'
    end if
    if  ( ( ( ndim .eq. ndim_const ) .and. ( d_unlimit .gt. ndim_const ) ) .or. &
          ( d_unlimit .lt. 0_i4 ) ) then
       print*, '***ERROR unlimited dimension out of bounds, must be positive but not greater than number of given dimensions'
       print*, '***Dims: ', ndim, ndim_const, d_unlimit, ndim_const, d_unlimit
       stop '***ERROR see StdOut'
    end if
    !
    allocate( chunksizes(  ndim ) )
    allocate( start(       ndim ) )
    allocate( counter(     ndim ) )
    allocate( dims(        ndim ) )
    allocate( dimid(       ndim ) )
    allocate( varid(  1 +  ndim ) )
    allocate( dummy_count( ndim ) )
    ! Initialize
    deflate      = 1
    if ( ndim .gt. ndim_const ) then
       chunksizes     = (/ size( arr, 1), size( arr, 2), &
            size( arr, 3), size( arr, 4 ),  1 /)
       dims(1:ndim-1) = shape( arr )
       dims(ndim)     = 1
    else
       chunksizes   = (/ size( arr, 1), size( arr, 2), &
            size( arr, 3), size( arr, 4) /)
       if ( d_unlimit .gt. 0 ) chunksizes( d_unlimit ) = 1 
       dims(1:ndim_const) = shape( arr )
    end if
    start(:)     = 1_i4
    counter(:)   = dims
    dummy_count  = 1
    dummy        = nf90_fill_double
    ! open the netcdf file
    if (present(ncid)) then
       if (ncid < 0_i4) then
          openfile = .true.
       else
          openfile = .false.
          f_handle = ncid
       endif
    else
       openfile = .true.
    endif
    if (openfile) then
       create_loc = .false.
       if (present(create)) create_loc = create
       f_handle = open_netcdf(f_name, create=create_loc)
    endif
    ! check whether variable exists
    if ( nf90_noerr .eq. nf90_inq_varid( f_handle, v_name, varid(ndim+1)) ) then
       ! append
       call check(nf90_inquire_variable(f_handle, varid(ndim+1), ndims=idim, dimids=dimid))
       ! consistency checks
       if (idim .ne. ndim) stop "var2nc_4d_dp: number of variable dimensions /= number of file variable dimensions."
       ! check unlimited dimension
       call check(nf90_inquire( f_handle, unlimitedDimId = u_dimid ))
       if ( u_dimid .eq. -1 ) stop 'var2nc_4d_dp: cannot append, no unlimited dimension defined'
       ! check for unlimited dimension
       if ( dimid( d_unlimit ) .ne. u_dimid ) stop 'var2nc_4d_dp: unlimited dimension not specified correctly'
       if (present(nrec)) then
          start(d_unlimit) = nrec
       else
          ! get length of unlimited dimension
          call check(nf90_inquire_dimension( f_handle, u_dimid, len = u_len ) )
          ! adapt start, that is find last written chunk
          do i = u_len, 1, -1
             if ( ne(dummy(1),nf90_fill_double) ) exit
             start(d_unlimit) = i
             call check( nf90_get_var( f_handle, varid(ndim+1), dummy, start, dummy_count ) )
          end do
          start(d_unlimit) = start(d_unlimit) + 1
       endif
    else
       ! define dimensions
       do i = 1, ndim
          ! check whether dimension exists
          if ( nf90_noerr .ne. nf90_inq_dimid( f_handle, dnames(i), dimid(i)) ) then
             ! create dimension
             if ( i .eq. d_unlimit ) then
                ! define unlimited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), NF90_UNLIMITED, dimid(i) ))
             else
                ! define limited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), dims(i), dimid(i) ))
             end if
          end if
       end do
       ! define variable
       call check(nf90_def_var(f_handle, v_name, NF90_DOUBLE, dimid, varid(ndim+1), &
            chunksizes=chunksizes(:), shuffle=.true., deflate_level=deflate))
       ! add attributes
       if ( present( long_name   ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'long_name', long_name ) )
       if ( present( units      ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'units', units ) )
       if ( present( missing_value ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'missing_value',missing_value ) )
       if ( present( attributes ) ) then
          do i = 1, size( attributes, dim = 1 )
             if ( trim( attributes(i,1)) .eq. 'missing_value' ) then
                ! write number
                read(attributes(i,2),'(F10.2)') dummy(1)
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), dummy(1) ) )
             else
                ! write string
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), trim(attributes(i,2)) ) )
             end if
          end do
       end if
       ! end definition
       call check(nf90_enddef(f_handle))
    end if
    ! check dimensions before writing
    do i=1, ndim_const
       call check(nf90_inquire_dimension(f_handle, dimid(i), dummy_name, dims(i)))
       if (trim(dummy_name) .ne. dnames(i)) &
            stop "var2nc_4d_dp: dimension name problem."
       if ( (dims(i) .ne. size(arr,i)) .and. ( d_unlimit .ne. i ) ) &
            stop "var2nc_4d_dp: variable dimension /= file variable dimension."
    end do
    ! write variable
    call check( nf90_put_var(f_handle, varid(ndim+1), arr, start, counter) )
    ! close netcdf_dataset
    if (present(ncid)) then
       if (ncid < 0_i4) ncid = f_handle
    else
       call close_netcdf( f_handle )
    endif
    !
  end subroutine var2nc_4d_dp

  subroutine var2nc_5d_i4( f_name, arr, dnames, v_name, dim_unlimited, &
       long_name, units, missing_value, attributes, create, ncid, nrec )
    !
    implicit none
    !
    integer(i4), parameter :: ndim_const = 5
    integer(i4)            :: ndim
    ! input variables
    character(len=*),                           intent(in) :: f_name
    integer(i4),      dimension(:,:,:,:,:),     intent(in) :: arr
    character(len=*),                           intent(in) :: v_name
    character(len=*), dimension(:),             intent(in) :: dnames
    ! optional
    integer(i4),                      optional, intent(in) :: dim_unlimited
    character(len=*),                 optional, intent(in) :: long_name
    character(len=*),                 optional, intent(in) :: units
    integer(i4),                      optional, intent(in) :: missing_value
    character(256),   dimension(:,:), optional, intent(in) :: attributes
    logical,                          optional, intent(in) :: create
    integer(i4),                      optional, intent(inout) :: ncid
    integer(i4),                      optional, intent(in) :: nrec
    ! local variables
    logical                                :: create_loc
    character(256)                         :: dummy_name
    integer(i4)                            :: deflate
    integer(i4), dimension(:), allocatable :: chunksizes
    integer(i4), dimension(:), allocatable :: start     ! start array for write
    integer(i4), dimension(:), allocatable :: counter   ! length array for write
    integer(i4)                            :: idim      ! read dimension on append
    integer(i4)                            :: f_handle
    integer(i4)                            :: d_unlimit ! index of unlimited dimension
    integer(i4)                            :: u_dimid   ! dimid of unlimited dimension
    integer(i4)                            :: u_len     ! length of unlimited dimension
    integer(i4), dimension(:), allocatable :: dims
    integer(i4), dimension(:), allocatable :: dimid     ! netcdf IDs of each dimension
    integer(i4), dimension(:), allocatable :: varid     ! dimension variables and var id
    integer(i4)                            :: i         ! loop indices
    integer(i4), dimension(:), allocatable :: dummy_count
    integer(i4), dimension(1)              :: dummy     ! dummy read
    logical                                :: openfile  ! tmp logical
    !
    ndim = size( dnames, 1 )
    d_unlimit = 0_i4
    if (present(dim_unlimited)) d_unlimit = dim_unlimited
    ! consistency checks
    if ( ndim .gt. ndim_const ) then
       print *, '***ERROR more than five dimension names given'
       stop '***ERROR see StdOut'
    end if
    if  ( ( ( ndim .eq. ndim_const ) .and. ( d_unlimit .gt. ndim_const ) ) .or. &
          ( d_unlimit .lt. 0_i4 ) ) then
       print*, '***ERROR unlimited dimension out of bounds, must be positive but not greater than number of given dimensions'
       print*, '***Dims: ', ndim, ndim_const, d_unlimit, ndim_const, d_unlimit
       stop '***ERROR see StdOut'
    end if
    !
    allocate( chunksizes(  ndim ) )
    allocate( start(       ndim ) )
    allocate( counter(     ndim ) )
    allocate( dims(        ndim ) )
    allocate( dimid(       ndim ) )
    allocate( varid(  1 +  ndim ) )
    allocate( dummy_count( ndim ) )
    ! initialize
    deflate      = 1
    chunksizes   = (/ size( arr, 1), size( arr, 2), size( arr, 3), size( arr, 4), size( arr, 5) /)
    if ( d_unlimit .gt. 0 ) chunksizes( d_unlimit ) = 1 
    dims(1:ndim) = shape( arr )
    start(:)     = 1_i4
    counter(:)   = dims
    dummy_count  = 1_i4
    dummy        = nf90_fill_int
    d_unlimit    = 0_i4
    if ( present( dim_unlimited ) ) d_unlimit = dim_unlimited
    ! open the netcdf file
    if (present(ncid)) then
       if (ncid < 0_i4) then
          openfile = .true.
       else
          openfile = .false.
          f_handle = ncid
       endif
    else
       openfile = .true.
    endif
    if (openfile) then
       create_loc = .false.
       if (present(create)) create_loc = create
       f_handle = open_netcdf(f_name, create=create_loc)
    endif
    ! check whether variable exists
    if ( nf90_noerr .eq. nf90_inq_varid( f_handle, v_name, varid(ndim+1)) ) then
       ! append
       call check(nf90_inquire_variable(f_handle, varid(ndim+1), ndims=idim, dimids=dimid))
       ! consistency checks
       if (idim .ne. ndim) stop "var2nc_5d_i4: number of variable dimensions /= number of file variable dimensions."
       ! check unlimited dimension
       call check(nf90_inquire( f_handle, unlimitedDimId = u_dimid ))
       if ( u_dimid .eq. -1 ) stop 'var2nc_5d_i4: cannot append, no unlimited dimension defined'
       ! check for unlimited dimension
       if ( dimid( d_unlimit ) .ne. u_dimid ) stop 'var2nc_5d_sp: unlimited dimension not specified correctly'
       if (present(nrec)) then
          start(d_unlimit) = nrec
       else
          ! get length of unlimited dimension
          call check(nf90_inquire_dimension( f_handle, u_dimid, len = u_len ) )
          ! adapt start, that is find last written chunk
          do i = u_len, 1, -1
             if ( dummy(1) /= nf90_fill_int ) exit
             start(d_unlimit) = i
             call check( nf90_get_var( f_handle, varid(ndim+1), dummy, start, dummy_count ) )
          end do
          start(d_unlimit) = start(d_unlimit) + 1
       endif
    else
       ! define dimensions
       do i = 1, ndim
          ! check whether dimension exists
          if ( nf90_noerr .ne. nf90_inq_dimid( f_handle, dnames(i), dimid(i)) ) then
             ! create dimension
             if ( i .eq. d_unlimit ) then
                ! define unlimited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), NF90_UNLIMITED, dimid(i) ))
             else
                ! define limited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), dims(i), dimid(i) ))
             end if
          end if
       end do
       ! define variable
       call check(nf90_def_var(f_handle, v_name, NF90_INT, dimid, varid(ndim+1), &
            chunksizes=chunksizes(:), shuffle=.true., deflate_level=deflate))
       ! add attributes
       if ( present( long_name   ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'long_name', long_name ) )
       if ( present( units      ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'units', units ) )
       if ( present( missing_value ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'missing_value', missing_value ) )
       if ( present( attributes ) ) then
          do i = 1, size( attributes, dim = 1 )
             if ( trim( attributes(i,1)) .eq. 'missing_value' ) then
                ! write number
                read(attributes(i,2),'(I6)') dummy(1)
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), dummy(1) ) )
             else
                ! write string
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), trim(attributes(i,2)) ) )
             end if
          end do
       end if
       ! end definition
       call check(nf90_enddef(f_handle))
    end if
    ! check dimensions before writing
    do i=1, ndim_const
       call check(nf90_inquire_dimension(f_handle, dimid(i), dummy_name, dims(i)))
       if (trim(dummy_name) .ne. dnames(i)) &
            stop "var2nc_5d_i4: dimension name problem."
       if ( (dims(i) .ne. size(arr,i)) .and. ( d_unlimit .ne. i ) ) &
            stop "var2nc_5d_i4: variable dimension /= file variable dimension."
    end do
    ! write variable
    call check( nf90_put_var(f_handle, varid(ndim+1), arr, start, counter) )
    ! close netcdf_dataset
    if (present(ncid)) then
       if (ncid < 0_i4) ncid = f_handle
    else
       call close_netcdf( f_handle )
    endif
    !
  end subroutine var2nc_5d_i4

  subroutine var2nc_5d_sp( f_name, arr, dnames, v_name, dim_unlimited, &
       long_name, units, missing_value, attributes, create, ncid, nrec )
    !
    implicit none
    !
    integer(i4), parameter :: ndim_const = 5
    integer(i4)            :: ndim
    ! input variables
    character(len=*),                           intent(in) :: f_name
    real(sp),         dimension(:,:,:,:,:),     intent(in) :: arr
    character(len=*),                           intent(in) :: v_name
    character(len=*), dimension(:),             intent(in) :: dnames
    ! optional
    integer(i4),                      optional, intent(in) :: dim_unlimited
    character(len=*),                 optional, intent(in) :: long_name
    character(len=*),                 optional, intent(in) :: units
    real(sp),                         optional, intent(in) :: missing_value
    character(256),   dimension(:,:), optional, intent(in) :: attributes
    logical,                          optional, intent(in) :: create
    integer(i4),                      optional, intent(inout) :: ncid
    integer(i4),                      optional, intent(in) :: nrec
    ! local variables
    logical                                :: create_loc
    character(256)                         :: dummy_name
    integer(i4)                            :: deflate
    integer(i4), dimension(:), allocatable :: chunksizes
    integer(i4), dimension(:), allocatable :: start     ! start array for write
    integer(i4), dimension(:), allocatable :: counter   ! length array for write
    integer(i4)                            :: idim      ! read dimension on append
    integer(i4)                            :: f_handle
    integer(i4)                            :: d_unlimit ! index of unlimited dimension
    integer(i4)                            :: u_dimid   ! dimid of unlimited dimension
    integer(i4)                            :: u_len     ! length of unlimited dimension
    integer(i4), dimension(:), allocatable :: dims
    integer(i4), dimension(:), allocatable :: dimid     ! netcdf IDs of each dimension
    integer(i4), dimension(:), allocatable :: varid     ! dimension variables and var id
    integer(i4)                            :: i         ! loop indices
    integer(i4), dimension(:), allocatable :: dummy_count
    real(sp),    dimension(1)              :: dummy     ! dummy read
    logical                                :: openfile  ! tmp logical
    !
    ndim = size( dnames, 1 )
    d_unlimit = 0_i4
    if (present(dim_unlimited)) d_unlimit = dim_unlimited
    ! consistency checks
    if ( ndim .gt. ndim_const ) then
       print *, '***ERROR more than five dimension names given'
       stop '***ERROR see StdOut'
    end if
    if  ( ( ( ndim .eq. ndim_const ) .and. ( d_unlimit .gt. ndim_const ) ) .or. &
          ( d_unlimit .lt. 0_i4 ) ) then
       print*, '***ERROR unlimited dimension out of bounds, must be positive but not greater than number of given dimensions'
       print*, '***Dims: ', ndim, ndim_const, d_unlimit, ndim_const, d_unlimit
       stop '***ERROR see StdOut'
    end if
    !
    allocate( chunksizes(  ndim ) )
    allocate( start(       ndim ) )
    allocate( counter(     ndim ) )
    allocate( dims(        ndim ) )
    allocate( dimid(       ndim ) )
    allocate( varid(  1 +  ndim ) )
    allocate( dummy_count( ndim ) )
    ! initialize
    deflate      = 1
    chunksizes   = (/ size( arr, 1), size( arr, 2), size( arr, 3), size( arr, 4), size( arr, 5) /)
    if ( d_unlimit .gt. 0 ) chunksizes( d_unlimit ) = 1 
    dims(1:ndim) = shape( arr )
    start(:)     = 1_i4
    counter(:)   = dims
    dummy_count  = 1_i4
    dummy        = nf90_fill_float
    d_unlimit    = 0_i4
    if ( present( dim_unlimited ) ) d_unlimit = dim_unlimited
    ! open the netcdf file
    if (present(ncid)) then
       if (ncid < 0_i4) then
          openfile = .true.
       else
          openfile = .false.
          f_handle = ncid
       endif
    else
       openfile = .true.
    endif
    if (openfile) then
       create_loc = .false.
       if (present(create)) create_loc = create
       f_handle = open_netcdf(f_name, create=create_loc)
    endif
    ! check whether variable exists
    if ( nf90_noerr .eq. nf90_inq_varid( f_handle, v_name, varid(ndim+1)) ) then
       ! append
       call check(nf90_inquire_variable(f_handle, varid(ndim+1), ndims=idim, dimids=dimid))
       ! consistency checks
       if (idim .ne. ndim) stop "var2nc_5d_sp: number of variable dimensions /= number of file variable dimensions."
       ! check unlimited dimension
       call check(nf90_inquire( f_handle, unlimitedDimId = u_dimid ))
       if ( u_dimid .eq. -1 ) stop 'var2nc_5d_sp: cannot append, no unlimited dimension defined'
       ! check for unlimited dimension
       if ( dimid( d_unlimit ) .ne. u_dimid ) stop 'var2nc_5d_sp: unlimited dimension not specified correctly'
       if (present(nrec)) then
          start(d_unlimit) = nrec
       else
          ! get length of unlimited dimension
          call check(nf90_inquire_dimension( f_handle, u_dimid, len = u_len ) )
          ! adapt start, that is find last written chunk
          do i = u_len, 1, -1
             if ( ne(dummy(1),nf90_fill_float) ) exit
             start(d_unlimit) = i
             call check( nf90_get_var( f_handle, varid(ndim+1), dummy, start, dummy_count ) )
          end do
          start(d_unlimit) = start(d_unlimit) + 1
       endif
    else
       ! define dimensions
       do i = 1, ndim
          ! check whether dimension exists
          if ( nf90_noerr .ne. nf90_inq_dimid( f_handle, dnames(i), dimid(i)) ) then
             ! create dimension
             if ( i .eq. d_unlimit ) then
                ! define unlimited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), NF90_UNLIMITED, dimid(i) ))
             else
                ! define limited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), dims(i), dimid(i) ))
             end if
          end if
       end do
       ! define variable
       call check(nf90_def_var(f_handle, v_name, NF90_FLOAT, dimid, varid(ndim+1), &
            chunksizes=chunksizes(:), shuffle=.true., deflate_level=deflate))
       ! add attributes
       if ( present( long_name   ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'long_name', long_name ) )
       if ( present( units      ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'units', units ) )
       if ( present( missing_value ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'missing_value', missing_value ) )
       if ( present( attributes ) ) then
          do i = 1, size( attributes, dim = 1 )
             if ( trim( attributes(i,1)) .eq. 'missing_value' ) then
                ! write number
                read(attributes(i,2),'(F10.2)') dummy(1)
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), dummy(1) ) )
             else
                ! write string
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), trim(attributes(i,2)) ) )
             end if
          end do
       end if
       ! end definition
       call check(nf90_enddef(f_handle))
    end if
    ! check dimensions before writing
    do i=1, ndim_const
       call check(nf90_inquire_dimension(f_handle, dimid(i), dummy_name, dims(i)))
       if (trim(dummy_name) .ne. dnames(i)) &
            stop "var2nc_5d_sp: dimension name problem."
       if ( (dims(i) .ne. size(arr,i)) .and. ( d_unlimit .ne. i ) ) &
            stop "var2nc_5d_sp: variable dimension /= file variable dimension."
    end do
    ! write variable
    call check( nf90_put_var(f_handle, varid(ndim+1), arr, start, counter) )
    ! close netcdf_dataset
    if (present(ncid)) then
       if (ncid < 0_i4) ncid = f_handle
    else
       call close_netcdf( f_handle )
    endif
    !
  end subroutine var2nc_5d_sp

  subroutine var2nc_5d_dp( f_name, arr, dnames, v_name, dim_unlimited, &
       long_name, units, missing_value, attributes, create, ncid, nrec )
    !
    implicit none
    !
    integer(i4), parameter :: ndim_const = 5
    integer(i4)            :: ndim
    ! input variables
    character(len=*),                           intent(in) :: f_name
    real(dp),         dimension(:,:,:,:,:),     intent(in) :: arr
    character(len=*),                           intent(in) :: v_name
    character(len=*), dimension(:),             intent(in) :: dnames
    ! optional
    integer(i4),                      optional, intent(in) :: dim_unlimited
    character(len=*),                 optional, intent(in) :: long_name
    character(len=*),                 optional, intent(in) :: units
    real(dp),                         optional, intent(in) :: missing_value
    character(256),   dimension(:,:), optional, intent(in) :: attributes
    logical,                          optional, intent(in) :: create
    integer(i4),                      optional, intent(inout) :: ncid
    integer(i4),                      optional, intent(in) :: nrec
    ! local variables
    logical                                :: create_loc
    character(256)                         :: dummy_name
    integer(i4)                            :: deflate
    integer(i4), dimension(:), allocatable :: chunksizes
    integer(i4), dimension(:), allocatable :: start     ! start array for write
    integer(i4), dimension(:), allocatable :: counter   ! length array for write
    integer(i4)                            :: idim      ! read dimension on append
    integer(i4)                            :: f_handle
    integer(i4)                            :: d_unlimit ! index of unlimited dimension
    integer(i4)                            :: u_dimid   ! dimid of unlimited dimension
    integer(i4)                            :: u_len     ! length of unlimited dimension
    integer(i4), dimension(:), allocatable :: dims
    integer(i4), dimension(:), allocatable :: dimid     ! netcdf IDs of each dimension
    integer(i4), dimension(:), allocatable :: varid     ! dimension variables and var id
    integer(i4)                            :: i         ! loop indices
    integer(i4), dimension(:), allocatable :: dummy_count
    real(dp),    dimension(1)              :: dummy     ! dummy read
    logical                                :: openfile  ! tmp logical
    !
    ndim = size( dnames, 1 )
    d_unlimit = 0_i4
    if (present(dim_unlimited)) d_unlimit = dim_unlimited
    ! consistency checks
    if ( ndim .gt. ndim_const ) then
       print *, '***ERROR more than five dimension names given'
       stop '***ERROR see StdOut'
    end if
    if  ( ( ( ndim .eq. ndim_const ) .and. ( d_unlimit .gt. ndim_const ) ) .or. &
          ( d_unlimit .lt. 0_i4 ) ) then
       print*, '***ERROR unlimited dimension out of bounds, must be positive but not greater than number of given dimensions'
       print*, '***Dims: ', ndim, ndim_const, d_unlimit, ndim_const, d_unlimit
       stop '***ERROR see StdOut'
    end if
    !
    allocate( chunksizes(  ndim ) )
    allocate( start(       ndim ) )
    allocate( counter(     ndim ) )
    allocate( dims(        ndim ) )
    allocate( dimid(       ndim ) )
    allocate( varid(  1 +  ndim ) )
    allocate( dummy_count( ndim ) )
    ! initialize
    deflate      = 1
    chunksizes   = (/ size( arr, 1), size( arr, 2), size( arr, 3), size( arr, 4), size( arr, 5) /)
    if ( d_unlimit .gt. 0 ) chunksizes( d_unlimit ) = 1 
    dims(1:ndim) = shape( arr )
    start(:)     = 1_i4
    counter(:)   = dims
    dummy_count  = 1
    dummy        = nf90_fill_double
    d_unlimit    = 0_i4
    if ( present( dim_unlimited ) ) d_unlimit = dim_unlimited
    ! open the netcdf file
    if (present(ncid)) then
       if (ncid < 0_i4) then
          openfile = .true.
       else
          openfile = .false.
          f_handle = ncid
       endif
    else
       openfile = .true.
    endif
    if (openfile) then
       create_loc = .false.
       if (present(create)) create_loc = create
       f_handle = open_netcdf(f_name, create=create_loc)
    endif
    ! check whether variable exists
    if ( nf90_noerr .eq. nf90_inq_varid( f_handle, v_name, varid(ndim+1)) ) then
       ! append
       call check(nf90_inquire_variable(f_handle, varid(ndim+1), ndims=idim, dimids=dimid))
       ! consistency checks
       if (idim .ne. ndim) stop "var2nc_5d_dp: number of variable dimensions /= number of file variable dimensions."
       ! check unlimited dimension
       call check(nf90_inquire( f_handle, unlimitedDimId = u_dimid ))
       if ( u_dimid .eq. -1 ) stop 'var2nc_5d_dp: cannot append, no unlimited dimension defined'
       ! check for unlimited dimension
       if ( dimid( d_unlimit ) .ne. u_dimid ) stop 'var2nc_5d_dp: unlimited dimension not specified correctly'
       if (present(nrec)) then
          start(d_unlimit) = nrec
       else
          ! get length of unlimited dimension
          call check(nf90_inquire_dimension( f_handle, u_dimid, len = u_len ) )
          ! adapt start, that is find last written chunk
          do i = u_len, 1, -1
             if ( ne(dummy(1),nf90_fill_double) ) exit
             start(d_unlimit) = i
             call check( nf90_get_var( f_handle, varid(ndim+1), dummy, start, dummy_count ) )
          end do
          start(d_unlimit) = start(d_unlimit) + 1
       endif
    else
       ! define dimensions
       do i = 1, ndim
          ! check whether dimension exists
          if ( nf90_noerr .ne. nf90_inq_dimid( f_handle, dnames(i), dimid(i)) ) then
             ! create dimension
             if ( i .eq. d_unlimit ) then
                ! define unlimited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), NF90_UNLIMITED, dimid(i) ))
             else
                ! define limited dimension
                call check(nf90_def_dim(f_handle, trim(dnames(i)), dims(i), dimid(i) ))
             end if
          end if
       end do
       ! define variable
       call check(nf90_def_var(f_handle, v_name, NF90_DOUBLE, dimid, varid(ndim+1), &
            chunksizes=chunksizes(:), shuffle=.true., deflate_level=deflate))
       ! add attributes
       if ( present( long_name   ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'long_name', long_name ) )
       if ( present( units      ) ) call check(nf90_put_att (f_handle, varid(ndim+1), 'units', units ) )
       if ( present( missing_value ) ) call check(nf90_put_att(f_handle, varid(ndim+1), 'missing_value',missing_value ) )
       if ( present( attributes ) ) then
          do i = 1, size( attributes, dim = 1 )
             if ( trim( attributes(i,1)) .eq. 'missing_value' ) then
                ! write number
                read(attributes(i,2),'(F10.2)') dummy(1)
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), dummy(1) ) )
             else
                ! write string
                call check( nf90_put_att( f_handle, varid(ndim+1), &
                     trim(attributes(i,1)), trim(attributes(i,2)) ) )
             end if
          end do
       end if
       ! end definition
       call check(nf90_enddef(f_handle))
    end if
    ! check dimensions before writing
    do i=1, ndim_const
       call check(nf90_inquire_dimension(f_handle, dimid(i), dummy_name, dims(i)))
       if (trim(dummy_name) .ne. dnames(i)) &
            stop "var2nc_5d_dp: dimension name problem."
       if ( (dims(i) .ne. size(arr,i)) .and. ( d_unlimit .ne. i ) ) &
            stop "var2nc_5d_dp: variable dimension /= file variable dimension."
    end do
    ! write variable
    call check( nf90_put_var(f_handle, varid(ndim+1), arr, start, counter) )
    ! close netcdf_dataset
    if (present(ncid)) then
       if (ncid < 0_i4) ncid = f_handle
    else
       call close_netcdf( f_handle )
    endif
    !
  end subroutine var2nc_5d_dp

  ! ----------------------------------------------------------------------------

  ! NAME
  !     write_dynamic_netcdf

  ! PURPOSE
  !     This routine writes data, where one dimension has the unlimited attribute.
  !     Therefore, the number of the record which should be written has to be
  !     specified.

  ! CALLING SEQUENCE
  !     call write_dynamic_netcdf(nc, rec)

  ! INTENT(IN)
  !     integer(i4) :: nc - stream id of an open netcdf file where data should be written
  !                         can be obtained by an create_netcdf call

  ! INTENT(IN)
  !     integer(i4) :: rec - record id of record which will be written in the file

  ! RESTRICTIONS
  !     Writes only data, where the data pointers of the structure V are assigned
  !     and where one dimension has the unlimited attribute. Moreover only one
  !     record will be written.
  !     Writes only 1 to 4 dim arrays, integer, single or double precision.

  ! EXAMPLE
  !     see test_mo_ncwrite

  ! LITERATURE
  !     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90.html

  ! HISTORY
  !     Written  Luis Samaniego  Feb 2011
  !     Modified Stephan Thober  Dec 2011 - added comments and generalized
  !              Matthias Cuntz  Jan 2012 - Info
  !              Stephan Thober  Jan 2012 - iRec is not optional
  !              Matthias Cuntz  Mar 2013 - removed Info

  subroutine write_dynamic_netcdf(ncId, irec)

    implicit none

    ! netcdf related variables
    integer(i4), intent(in)           :: ncId
    integer(i4), intent(in)           :: iRec

    integer(i4)                       :: i
    ! NOTES: 1) netcdf file must be on *** data mode ***
    !        2) start and end of the data chuck is controled by
    !           V(:)%start and  V(:)%count

    ! set values for variables (one scalar or grid at a time)

    do i = 1, nVars
       if (.not. V(i)%unlimited) cycle
       if (.not. V(i)%wFlag) cycle
       V(i)%start (V(i)%nDims) = iRec
       select case (V(i)%xtype)
       case(NF90_BYTE)
          select case (V(i)%nDims)
          case (0)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G0_b, V(i)%start))
          case (1)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G1_b, V(i)%start, V(i)%count))
          case (2)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G2_b, V(i)%start, V(i)%count))
          case (3)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G3_b, V(i)%start, V(i)%count))
          case (4)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G4_b, V(i)%start, V(i)%count))
          end select
       case (NF90_INT)
          select case (V(i)%nDims-1)
          case (0)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G0_i, V(i)%start))
          case (1)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G1_i, V(i)%start, V(i)%count))
          case (2)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G2_i, V(i)%start, V(i)%count))
          case (3)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G3_i, V(i)%start, V(i)%count))
          case (4)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G4_i, V(i)%start, V(i)%count))
          end select
       case (NF90_FLOAT)
          select case (V(i)%nDims-1)
          case (0)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G0_f, V(i)%start))
          case (1)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G1_f, V(i)%start, V(i)%count))
          case (2)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G2_f, V(i)%start, V(i)%count))
          case (3)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G3_f, V(i)%start, V(i)%count))
          case (4)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G4_f, V(i)%start, V(i)%count))
          end select
       case (NF90_DOUBLE)
          select case (V(i)%nDims-1)
          case (0)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G0_d, V(i)%start))
          case (1)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G1_d, V(i)%start, V(i)%count))
          case (2)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G2_d, V(i)%start, V(i)%count))
          case (3)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G3_d, V(i)%start, V(i)%count))
          case (4)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G4_d, V(i)%start, V(i)%count))
          end select
       end select
    end do

  end subroutine write_dynamic_netcdf

  ! ----------------------------------------------------------------------------

  ! NAME
  !     write_static_netcdf

  ! PURPOSE
  !     This routines writes static data in the netcdf file that is data
  !     where no dimension has the unlimited attribute.

  ! CALLING SEQUENCE
  !     call write_static_netcdf(ncid)

  ! INTENT(IN)
  !     integer(i4) :: ncid - stream id of an open netcdf file where data should be written
  !                           can be obtained by an create_netcdf call

  ! RESTRICTIONS
  !     Writes only data, where the data pointers of the structure V are assigned.
  !     Writes only 1 to 4 dim arrays, integer, single or double precision.

  ! EXAMPLE
  !     see test_mo_ncwrite

  ! LITERATURE
  !     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90.html

  ! HISTORY
  !     Written  Luis Samaniego  Feb 2011
  !     Modified Stephan Thober  Dec 2011 - added comments and generalized
  !              Matthias Cuntz  Jan 2012 - Info
  !              Matthias Cuntz  Mar 2013 - removed Info

  subroutine write_static_netcdf(ncId)

    implicit none

    ! netcdf related variables
    integer(i4), intent(in)           :: ncId

    integer(i4)                       :: i

    ! write all static variables
    do i = 1, nVars
       if (V(i)%unlimited) cycle
       if (.not. V(i)%wFlag) cycle
       select case (V(i)%xtype)
       case(NF90_BYTE)
          select case (V(i)%nDims)
          case (0)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G0_b))
          case (1)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G1_b))
          case (2)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G2_b))
          case (3)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G3_b))
          case (4)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G4_b))
          end select
       case (NF90_INT)
          select case (V(i)%nDims)
          case (0)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G0_i))
          case (1)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G1_i))
          case (2)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G2_i))
          case (3)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G3_i))
          case (4)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G4_i))
          end select
       case (NF90_FLOAT)
          select case (V(i)%nDims)
          case (0)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G0_f))
          case (1)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G1_f))
          case (2)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G2_f))
          case (3)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G3_f))
          case (4)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G4_f))
          end select
       case (NF90_DOUBLE)
          select case (V(i)%nDims)
          case (0)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G0_d))
          case (1)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G1_d))
          case (2)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G2_d))
          case (3)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G3_d))
          case (4)
             call check(nf90_put_var(ncId,  V(i)%varId, V(i)%G4_d))
          end select
       end select
    end do

  end subroutine write_static_netcdf

  ! -----------------------------------------------------------------------------
  ! PRIVATE PART
  !

  ! -----------------------------------------------------------------------------

  ! private open netcdf function - returns file handle
  function open_netcdf( f_name, create )
    implicit none
    ! input variables
    character(len=*),           intent(in) :: f_name
    logical,                    intent(in) :: create ! flag indicates that file exists
    ! output variables
    integer(i4)   :: open_netcdf
    !
    if ( create ) then
       ! create file
       call check( nf90_create( trim(f_name), NF90_NETCDF4, open_netcdf ) )
    else
       ! open file
       call check( nf90_open( trim(f_name), NF90_WRITE, open_netcdf ) )
    end if
  end function open_netcdf


  ! -----------------------------------------------------------------------------

  !  private error checking routine
  subroutine check(status)

    implicit none

    integer(i4), intent(in) :: status

    if (status /= nf90_noerr) then
       write(*,*) 'mo_ncwrite.check error: ', trim(nf90_strerror(status))
       stop
    end if

  end subroutine check

  ! -----------------------------------------------------------------------------

end module mo_ncwrite
