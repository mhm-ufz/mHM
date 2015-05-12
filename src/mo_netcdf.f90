module mo_netcdf

  use mo_kind,         only: i4, i8, dp, sp
  use netcdf,          only: &
       nf90_open, nf90_close, NF90_OPEN, NF90_CREATE, NF90_WRITE,            &
       NF90_NETCDF4, NF90_NOERR, NF90_NOWRITE, nf90_def_dim,                 &
       NF90_UNLIMITED, NF90_SHORT, NF90_INT, NF90_INT64, NF90_FLOAT,         &
       NF90_DOUBLE, nf90_def_var, NF90_CHUNKED, NF90_CONTIGUOUS,             &
       NF90_ENDIAN_NATIVE, nf90_put_var,nf90_inquire, nf90_put_att,          &
       NF90_GLOBAL, nf90_strerror, nf90_inquire_dimension, nf90_inq_dimid,   &
       nf90_inquire_variable

  implicit none 

  type NcVariable

     integer(i4)   :: id
     integer(i4)   :: pid
     character(32) :: name
     
   contains

     procedure, private :: createVariableAttributeChar
     procedure, private :: createVariableAttributeDp
     procedure, private :: createVariableAttributeSp
     procedure, private :: createVariableAttributeI4
     procedure, private :: createVariableAttributeI8

     procedure, private :: putData1dI4
     procedure, private :: putData2dI4
     procedure, private :: putData3dI4
     procedure, private :: putData1dDp
     procedure, private :: putData2dDp
     procedure, private :: putData3dDp

     procedure, private :: preparePutArguments

     procedure, public  :: getNoDimensions
     
     generic, public :: putData => &
          putData1dI4, &
          putData2dI4, &
          putData3dI4, &
          putData1dDp, &
          putData2dDp, &
          putData3dDp
          
     generic, public :: createAttribute => &
          createVariableAttributeChar, &
          createVariableAttributeDp,   &
          createVariableAttributeSp,   &
          createVariableAttributeI4,   &
          createVariableAttributeI8
     
  end type NcVariable
  
  type NcDimension
     integer(i4)    :: id
     character(16)  :: name     
     integer(i4)    :: length
  end type NcDimension
  
  type NcDataset
     character(256)                       :: fname
     character(1)                         :: mode
     integer(i4)                          :: id
     
   contains
     
     procedure, private :: createGlobalAttributeChar
     procedure, private :: createGlobalAttributeDp
     procedure, private :: createGlobalAttributeSp
     procedure, private :: createGlobalAttributeI4
     procedure, private :: createGlobalAttributeI8
     
     procedure, private :: getDimensionByName
     procedure, private :: getDimensionById

     procedure, private :: createVariableWithTypes
     procedure, private :: createVariableWithNames
     procedure, private :: createVariableWithIds
     
     procedure, public  :: close
     procedure, public  :: createDimension

     generic,   public  :: createAttribute => &
          createGlobalAttributeChar, &
          createGlobalAttributeDp,   &
          createGlobalAttributeSp,   &
          createGlobalAttributeI4,   &
          createGlobalAttributeI8
     
     generic,   public  :: getDimension => &
          getDimensionById, &
          getDimensionByName

     generic,   public  :: createVariable => &
          createVariableWithNames, &
          createVariableWithTypes, &
          createVariableWithIds
     
  end type NcDataset

  interface NcDataset
     procedure initDataset
  end interface NcDataset
    
contains

  type(NcDataset)  function initDataset(fname,mode)
     ! mode:  
     !       r: read
     !       w: write/create
     !       a: alter

     character(*),intent(in)   :: fname
     character(1),intent(in)   :: mode
     integer(i4)     :: status
     type(NcDataset) :: nc
     
     ! try to open file
     select case(mode)
        case("w")
           status = nf90_create(trim(fname), NF90_NETCDF4, nc%id)
        case("r")
           status = nf90_open(trim(fname), NF90_NOWRITE, nc%id)           
        case("a")
           status = nf90_open(trim(fname), NF90_WRITE, nc%id)
        case default
           print*, "Mode argument must be in 'w','r','a' ! "
           stop 1
     end select        
     call check(status,"Failed to open file: " // fname)

     nc%fname = fname
     nc%mode  = mode
        
     initDataset = nc
     
   end function initDataset
   
   subroutine close(self)
     class(NcDataset) :: self          
     if (nf90_close(self%id) .ne. NF90_NOERR) then
        print*, "Failed to close file: '", self%fname, "'"
        stop 1
     end if

   end subroutine close
   
   subroutine check(status,msg)
     integer(i4) , intent(in) :: status
     character(*), intent(in) :: msg

     if (status .ne. NF90_NOERR) then
        print*, msg
        print*, nf90_strerror(status)
        stop 1
     end if
   end subroutine check

   function createDimension(self, name, length)

     ! length <= 0 -> unlimited 
     class(NcDataset), intent(in) :: self
     character(*)    , intent(in) :: name
     integer(i4)     , intent(in) :: length
     integer(i4)                  :: id,dimlength
     type(NcDimension)            :: createDimension
     
     if (length .le. 0) then
        dimlength = NF90_UNLIMITED
     else
        dimlength = length
     end if

     call check(nf90_def_dim(self%id, name, dimlength, id), &
          "Failed to create dimension: " // name)

     createDimension = NcDimension(id,trim(name),length)

   end function createDimension

   function createVariableWithIds(self, name, dtype, dimensions, contiguous, &
        chunksizes, deflate_level, shuffle, fletcher32, endianness, &
        cache_size, cache_nelems, cache_preemption)

     class(NcDataset),intent(in)                     :: self
     character(*),intent(in)                         :: name
     character(3),intent(in)                         :: dtype
     integer(i4),dimension(:),intent(in)             :: dimensions

     logical, optional, intent(in)                   :: contiguous,shuffle, fletcher32
     integer(i4), optional, intent(in)               :: endianness,deflate_level,cache_size, &
          cache_nelems, cache_preemption
     integer(i4), optional, dimension(:), intent(in) :: chunksizes
     
     type(NcVariable)                                :: createVariableWithIds
     integer(i4)                                     :: ncdtype,varid,status

     ncdtype = getDtype(dtype)
     status = nf90_def_var(self%id, name, ncdtype, invertArray(dimensions), varid, contiguous, &
          chunksizes, deflate_level, shuffle, fletcher32, endianness, &
          cache_size, cache_nelems, cache_preemption)
     call check(status, "Failed to create variable: " // name)
     
     createVariableWithIds = NcVariable(varid,self%id,name)

   end function createVariableWithIds
   
   function createVariableWithNames(self, name, dtype, dimensions, contiguous, &
        chunksizes, deflate_level, shuffle, fletcher32, endianness, &
        cache_size, cache_nelems, cache_preemption)
     
     class(NcDataset),intent(in)                     :: self
     character(*),intent(in)                         :: name
     character(3),intent(in)                         :: dtype
     character(*),dimension(:),intent(in)            :: dimensions

     logical, optional, intent(in)                   :: contiguous,shuffle, fletcher32
     integer(i4), optional, intent(in)               :: endianness,deflate_level,cache_size, &
          cache_nelems, cache_preemption
     integer(i4), optional, dimension(:), intent(in) :: chunksizes
     
     type(NcVariable)                                :: createVariableWithNames
     type(NcDimension)                               :: dim
     integer(i4),dimension(size(dimensions))         :: dimids
     integer(i4)                                     :: i

     do i = 1,size(dimensions)
        dim = self%getDimension(dimensions(i))
        dimids(i) = dim%id
     end do
     
     createVariableWithNames = createVariableWithIds(self, name, dtype, dimids, contiguous, &
        chunksizes, deflate_level, shuffle, fletcher32, endianness, &
        cache_size, cache_nelems, cache_preemption)
     
   end function createVariableWithNames

   function createVariableWithTypes(self, name, dtype, dimensions,contiguous, &
        chunksizes, deflate_level, shuffle, fletcher32, endianness, &
        cache_size, cache_nelems, cache_preemption)
     
     class(NcDataset),intent(in)                     :: self
     character(*),intent(in)                         :: name
     character(3),intent(in)                         :: dtype
     type(NcDimension),dimension(:),intent(in)       :: dimensions

     logical, optional, intent(in)                   :: contiguous,shuffle, fletcher32
     integer(i4), optional, intent(in)               :: &
          endianness,deflate_level,cache_size, cache_nelems, cache_preemption
     integer(i4), optional, dimension(size(dimensions)), intent(in) :: chunksizes
     
     type(NcVariable)                                :: createVariableWithTypes
     type(NcDimension)                               :: dim
     integer(i4),dimension(size(dimensions))         :: dimids
     integer(i4)                                     :: i
     
     do i = 1,size(dimensions)
        dim = dimensions(i)
        dimids(i) = dim%id
     end do

     createVariableWithTypes = createVariableWithIds(self, name, dtype, dimids, contiguous, &
        chunksizes, deflate_level, shuffle, fletcher32, endianness, &
        cache_size, cache_nelems, cache_preemption)
     
   end function createVariableWithTypes   

   
   function getDimensionById(self, id)
     class(NcDataset), intent(in)          :: self
     integer(i4)                           :: id
     type(NcDimension)                     :: getDimensionById
     character(32)                         :: msg

     write(msg,*) id
     call check(nf90_inquire_dimension(self%id,id,getDimensionById%name,getDimensionById%length), &
          "Could not inquire dimension: " // msg)
     getDimensionById%id = id     
   end function getDimensionById

   function getDimensionByName(self, name)
     class(NcDataset), intent(in)          :: self
     character(*)                          :: name
     type(NcDimension)                     :: getDimensionByName
     integer(i4)                           :: id

     call check(nf90_inq_dimid(self%id,name,id), &
          "Could not inquire dimension: " // name)
     getDimensionByName = self%getDimensionById(id)
   end function getDimensionByName

   function getNoDimensions(self)
     class(NcVariable), intent(in) :: self
     integer(i4)                   :: getNoDimensions
     
     call check(nf90_inquire_variable(self%pid,self%id,ndims=getNoDimensions),&
          "Could not inquire dimension: " // self%name)
   end function getNoDimensions

   
   ! The createGlobalAtrribute procedures
   subroutine createGlobalAttributeChar(self,name,data)     
     class(NcDataset),intent(in) :: self
     character(*), intent(in)    :: name
     character(*), intent(in)    :: data

     call check(nf90_put_att(self%id,NF90_GLOBAL,name,data), &
          "Failed to write attribute: " // name )     
   end subroutine createGlobalAttributeChar

   subroutine createGlobalAttributeDp(self,name,data)
     class(NcDataset),intent(in) :: self
     character(*), intent(in)    :: name
     real(dp), intent(in)        :: data

     call check(nf90_put_att(self%id,NF90_GLOBAL,name,data), &
          "Failed to write attribute: " // name )
   end subroutine createGlobalAttributeDp

   subroutine createGlobalAttributeSp(self,name,data)
     class(NcDataset),intent(in) :: self
     character(*), intent(in)    :: name
     real(sp), intent(in)        :: data

     call check(nf90_put_att(self%id,NF90_GLOBAL,name,data), &
          "Failed to write attribute: " // name )
   end subroutine createGlobalAttributeSp

   subroutine createGlobalAttributeI4(self,name,data)
     class(NcDataset),intent(in) :: self
     character(*), intent(in)    :: name
     integer(i4), intent(in)     :: data

     call check(nf90_put_att(self%id,NF90_GLOBAL,name,data), &
          "Failed to write attribute: " // name )
   end subroutine createGlobalAttributeI4

   subroutine createGlobalAttributeI8(self,name,data)
     class(NcDataset),intent(in) :: self
     character(*), intent(in)    :: name
     integer(i8), intent(in)     :: data

     call check(nf90_put_att(self%id,NF90_GLOBAL,name,data), &
          "Failed to write attribute: " // name )
   end subroutine createGlobalAttributeI8

   
   ! The createVariableAtrribute procedures 
   subroutine createVariableAttributeChar(self,name,data)
     class(NcVariable),intent(in) :: self
     character(*), intent(in)     :: name
     character(*), intent(in)     :: data

     call check(nf90_put_att(self%pid,self%id,name,data), &
          "Failed to write attribute: " // name )
   end subroutine createVariableAttributeChar

   subroutine createVariableAttributeDp(self,name,data)
     class(NcVariable),intent(in) :: self
     character(*), intent(in)     :: name
     real(dp), intent(in)         :: data

     call check(nf90_put_att(self%pid,self%id,name,data), &
          "Failed to write attribute: " // name )
   end subroutine createVariableAttributeDp

   subroutine createVariableAttributeSp(self,name,data)
     class(NcVariable),intent(in) :: self
     character(*), intent(in)     :: name
     real(sp), intent(in)         :: data

     call check(nf90_put_att(self%pid,self%id,name,data), &
          "Failed to write attribute: " // name )
   end subroutine createVariableAttributeSp

   subroutine createVariableAttributeI4(self,name,data)
     class(NcVariable),intent(in) :: self
     character(*), intent(in)     :: name
     integer(i4), intent(in)      :: data

     call check(nf90_put_att(self%pid,self%id,name,data), &
          "Failed to write attribute: " // name )
   end subroutine createVariableAttributeI4

   subroutine createVariableAttributeI8(self,name,data)
     class(NcVariable),intent(in) :: self
     character(*), intent(in)     :: name
     integer(i8), intent(in)      :: data

     call check(nf90_put_att(self%pid,self%id,name,data), &
          "Failed to write attribute: " // name )
   end subroutine createVariableAttributeI8

      
   ! putData i4
   subroutine putData1dI4(self,values,start,count,stride)
     class(NcVariable), intent(in)                   :: self
     integer(i4), dimension(:), intent(in)           :: values
     integer(i4), dimension(:), optional, intent(in) :: start, count, stride
     integer(i4), dimension(:),allocatable           :: pstart, pcount, pstride, pmap

     call self%preparePutArguments(shape(values), start, count, stride, pstart, pcount, pstride, pmap)
     call check(nf90_put_var(self%pid, self%id, values, pstart, pcount, pstride, pmap), &
          "Failed to write data into variable: " // trim(self%name))     
   end subroutine putData1dI4

   subroutine putData2dI4(self,values,start,count,stride)
     class(NcVariable), intent(in)                   :: self
     integer(i4), dimension(:,:), intent(in)         :: values
     integer(i4), dimension(:), optional, intent(in) :: start, count, stride
     integer(i4), dimension(:),allocatable           :: pstart, pcount, pstride, pmap

     call self%preparePutArguments(shape(values), start, count, stride, pstart, pcount, pstride, pmap)
     call check(nf90_put_var(self%pid, self%id, values, pstart, pcount, pstride, pmap), &
          "Failed to write data into variable: " // trim(self%name))          
   end subroutine putData2dI4

   subroutine putData3dI4(self,values,start,count,stride)
     class(NcVariable), intent(in)                   :: self
     integer(i4), dimension(:,:,:), intent(in)       :: values
     integer(i4), dimension(:), optional, intent(in) :: start, count, stride
     integer(i4), dimension(:),allocatable           :: pstart, pcount, pstride, pmap

     call self%preparePutArguments(shape(values), start, count, stride, pstart, pcount, pstride, pmap)
     call check(nf90_put_var(self%pid, self%id, values, pstart, pcount, pstride, pmap), &
          "Failed to write data into variable: " // trim(self%name))          
   end subroutine putData3dI4

   ! putData dp
   subroutine putData1dDp(self,values,start,count,stride)
     class(NcVariable), intent(in)                   :: self
     real(dp), dimension(:), intent(in)           :: values
     integer(i4), dimension(:), optional, intent(in) :: start, count, stride
     integer(i4), dimension(:),allocatable           :: pstart, pcount, pstride, pmap

     call self%preparePutArguments(shape(values), start, count, stride, pstart, pcount, pstride, pmap)
     call check(nf90_put_var(self%pid, self%id, values, pstart, pcount, pstride, pmap), &
          "Failed to write data into variable: " // trim(self%name))     
   end subroutine putData1dDp

   subroutine putData2dDp(self,values,start,count,stride)
     class(NcVariable), intent(in)                   :: self
     real(dp), dimension(:,:), intent(in)            :: values
     integer(i4), dimension(:), optional, intent(in) :: start, count, stride
     integer(i4), dimension(:),allocatable           :: pstart, pcount, pstride, pmap

     call self%preparePutArguments(shape(values), start, count, stride, pstart, pcount, pstride, pmap)
     call check(nf90_put_var(self%pid, self%id, values, pstart, pcount, pstride, pmap), &
          "Failed to write data into variable: " // trim(self%name))          
   end subroutine putData2dDp

   subroutine putData3dDp(self,values,start,count,stride)
     class(NcVariable), intent(in)                   :: self
     real(dp), dimension(:,:,:), intent(in)       :: values
     integer(i4), dimension(:), optional, intent(in) :: start, count, stride
     integer(i4), dimension(:),allocatable           :: pstart, pcount, pstride, pmap

     call self%preparePutArguments(shape(values), start, count, stride, pstart, pcount, pstride, pmap)
     call check(nf90_put_var(self%pid, self%id, values, pstart, pcount, pstride, pmap), &
          "Failed to write data into variable: " // trim(self%name))          
   end subroutine putData3dDp

   
   subroutine preparePutArguments(self,inshape,&
        instart,incount,instride, &
        outstart,outcount,outstride,outmap)
     class(NcVariable),intent(in)                    :: self
     integer(i4), dimension(:),intent(in)                         :: inshape
     integer(i4), dimension(:),intent(in), optional              :: instart,incount,instride
     integer(i4), dimension(:),allocatable, intent(out)          :: outstart,outcount,outstride,outmap
     integer(i4) :: datarank
     
     datarank = self%getNoDimensions()
     allocate(outstart(datarank),outcount(datarank),outstride(datarank),outmap(datarank))
     
     if (present(instart)) then
        outstart = invertArray(instart)
     else
        outstart = 1
     end if

     if (present(incount)) then
        outcount = invertArray(incount)
     else
        outcount = invertArray(inshape)
     end if
     
     if (present(instride)) then
        outstride = invertArray(instride)
     else
        outstride = 1
     end if

     ! That seems to be a hack ...
     outmap = invertArray(generateImap(invertArray(outcount)))
     
   end subroutine preparePutArguments
   
   
   function getDtype(dtypestr)
     integer(i4)          :: getDtype
     character(*)         :: dtypestr

     select case(dtypestr)
        case("f32")
           getDtype = NF90_FLOAT
        case("f64")
           getDtype = NF90_DOUBLE
        case("i16")
           getDtype = NF90_SHORT
        case("i32")
           getDtype = NF90_INT
        case("i64")
           getDtype = NF90_INT64
        case default
           print*, "Datatype not understood: ", dtypestr
           error stop
     end select

   end function getDtype

   function invertArray(inarray)
     integer(i4),dimension(:), intent(in) :: inarray
     integer(i4),dimension(size(inarray)) :: invertArray
     invertArray = inarray(ubound(inarray,1):lbound(inarray,1):-1)
   end function invertArray

   function generateImap(dshape)
     ! Generat imap. Needed as the dimension ordering in is not reversed (as id the default)
     ! in createVariable
     integer(i4), dimension(:), intent(in)  :: dshape
     integer(i4), dimension(size(dshape))   :: generateImap, tmp
     integer(i4)                            :: i     

     tmp = 1
     do i=1,size(dshape)-1
        tmp(i+1) = tmp(i)*dshape(i)
     end do
     ! print*, "imap: ",tmp
     generateImap = tmp

   end function generateImap

   
end module mo_netcdf

