!> \file mo_upscaling_operators.f90

!> \brief Module containing upscaling operators.

!> \details This module provides the routines for upscaling_operators.

!> \authors Giovanni Dalmasso, Rohini Kumar
!> \date Dec 2012

module mo_upscaling_operators

  ! This module contains the functions for upscaling grid L0_fineScale_2D_data.

  ! Written  Giovanni Dalmasso, Rohini Kumar, Dec 2012

  use mo_kind, only: i4, dp

  implicit none

  private

  public :: majority_statistics       ! upscale grid L0_fineScale_2D_data based on a majority statistics
  public :: L0_fractionalCover_in_Lx  ! fractional coverage of a given class of L0 fields in Lx field (Lx = L1 or L11)
  public :: upscale_arithmetic_mean   ! upscale grid L0_fineScale_2D_data based on a ARITHMETIC MEAN
  public :: upscale_harmonic_mean     ! upscale grid L0_fineScale_2D_data based on a HARMONIC MEAN
  public :: upscale_geometric_mean    ! upscale grid L0_fineScale_2D_data based on a GEOMETRIC MEAN

contains

  ! ----------------------------------------------------------------------------

  !      NAME
  !          majority_statistics

  !>         \brief majority statistics

  !>         \details upscale grid L0_fineScale_2D_data based on a majority statistics

  !      INTENT(IN)
  !>         \param[in] "integer(i4) :: nClass"                       - number of classes
  !>         \param[in] "integer(i4) :: L1_upper_rowId_cell(:)"       - upper row boundary (level-0) of a level-1 cell
  !>         \param[in] "integer(i4) :: L1_lower_rowId_cell(:)"       - lower row boundary (level-0) of a level-1 cell
  !>         \param[in] "integer(i4) :: L1_left_colonId_cell(:)"      - left colon boundary (level-0) of a level-1 cell
  !>         \param[in] "integer(i4) :: L1_right_colonId_cell(:)"     - right colon boundary (level-0) of a level-1 cell
  !>         \param[in] "integer(i4) :: L0_fineScale_2D_data(:,:)"    - high resolution data

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
  !>       \return     integer(i4) :: majority_statistics(:) &mdash; Upscaled variable based on majority.

  !     RESTRICTIONS
  !>       \note Input values must be floating points.

  !      EXAMPLE
  !          calling sequence
  !          variable = majority_statistics(   nClass,                 &
  !                                           L1_upper_rowId_cell,    &
  !                                           L1_lower_rowId_cell,    &
  !                                           L1_left_colonId_cell,   &
  !                                           L1_right_colonId_cell,  &
  !                                           L0_fineScale_2D_data,   &
  !                                         )
  !
  !     LITERATURE
  !         None

  !      HISTORY
  !>         \author Giovanni Dalmasso, Rohini Kumar
  !>         \date Dec 2012
  !          Written, Giovanni Dalmasso, Dec 2012

  function majority_statistics(   nClass, &   ! number of classes
       L1_upper_rowId_cell,               &   ! upper row boundary (level-0) of a level-1 cell
       L1_lower_rowId_cell,               &   ! lower row boundary (level-0) of a level-1 cell
       L1_left_colonId_cell,              &   ! left colon boundary (level-0) of a level-1 cell
       L1_right_colonId_cell,             &   ! right colon boundary (level-0) of a level-1 cell
       L0_fineScale_2D_data               &   ! high resolution data
       )

    implicit none
    ! input
    integer(i4),                    intent(in)              :: nClass
    integer(i4), dimension(:),      intent(in)              :: L1_upper_rowId_cell
    integer(i4), dimension(:),      intent(in)              :: L1_lower_rowId_cell
    integer(i4), dimension(:),      intent(in)              :: L1_left_colonId_cell
    integer(i4), dimension(:),      intent(in)              :: L1_right_colonId_cell
    integer(i4), dimension(:, :),   intent(in)              :: L0_fineScale_2D_data
    ! output
    integer(i4), dimension( size(L1_upper_rowId_cell,1) )   :: majority_statistics

    ! local variables
    integer(i4) :: L1_nCells
    integer(i4) :: iu, id, jl, jr
    integer(i4) :: nC
    integer(i4) :: max_val
    integer(i4) :: kk, ll

    L1_nCells = size(majority_statistics,1)

    do kk = 1, L1_nCells
       iu = L1_upper_rowId_cell(kk)
       id = L1_lower_rowId_cell(kk)
       jl = L1_left_colonId_cell(kk)
       jr = L1_right_colonId_cell(kk)

       max_val = -9999
       do ll = 1, nClass
          nC = count(L0_fineScale_2D_data(iu:id, jl:jr) == ll)
          if(nC > max_val) then
             majority_statistics(kk) = ll
             max_val = nC
          end if
       end do
    end do

  end function majority_statistics

  ! ------------------------------------------------------------------

  !      NAME
  !          L0_fractionalCover_in_Lx

  !>         \brief fractional coverage of a given class of L0 fields in Lx field (Lx = L1 or L11)

  !>         \details Fractional coverage of a given class of L0 fields in Lx field (Lx = L1 or L11). 
  !>                  For example, this routine can be used for calculating the karstic fraction.

  !     INTENT(IN)
  !>         \param[in] "integer(i4), dimension(:,:)  :: dataIn0"            input fields at finer scale
  !>         \param[in] "integer(i4)                  :: classId"            class id for which fraction has to be estimated
  !>         \param[in] "logical, dimension(:,:)      :: mask0"              finer scale L0 mask
  !>         \param[in] "integer(i4), dimension(:)    :: L0upBound_inLx"     row start at finer L0 scale
  !>         \param[in] "integer(i4), dimension(:)    :: L0downBound_inLx"   row end   at finer L0 scale
  !>         \param[in] "integer(i4), dimension(:)    :: L0leftBound_inLx"   col start at finer L0 scale
  !>         \param[in] "integer(i4), dimension(:)    :: L0rightBound_inLx"  col end   at finer L0 scale
  !>         \param[in] "integer(i4), dimension(:)    :: nTCells0_inLx"      total number of valid L0 cells in a given Lx cell

  !     INTENT(INOUT)
  !         None

  !      INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !          None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return     real(dp) :: L0_fractionalCover_in_Lx(:) &mdash; packed 1D fraction coverage (Lx) of given class id

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !      HISTORY
  !>         \author Rohini Kumar
  !>         \date Feb 2013

  ! ------------------------------------------------------------------

  function L0_fractionalCover_in_Lx( &
       dataIn0,                   &  ! input fields at finer scale
       classId,                   &  ! class id for which fraction has to be estimated
       mask0,                     &  ! finer scale L0 mask
       L0upBound_inLx,            &  ! row start at finer L0 scale
       L0downBound_inLx,          &  ! row end   at finer L0 scale
       L0leftBound_inLx,          &  ! col start at finer L0 scale
       L0rightBound_inLx,         &  ! col end   at finer L0 scale
       nTCells0_inLx              &  ! total number of valid L0 cells in a given Lx cell
       ) result(frac_cover_Lx)       ! packed 1D data at coarser scale (Lx = L1 or L11)

    use mo_mhm_constants, only: nodata_i4

    implicit none

    ! input
    integer(i4),dimension(:),   intent(in)   :: dataIn0
    integer(i4),                intent(in)   :: classId
    logical, dimension(:,:),    intent(in)   :: mask0
    integer(i4), dimension(:),  intent(in)   :: L0upBound_inLx
    integer(i4), dimension(:),  intent(in)   :: L0downBound_inLx
    integer(i4), dimension(:),  intent(in)   :: L0leftBound_inLx
    integer(i4), dimension(:),  intent(in)   :: L0rightBound_inLx
    integer(i4), dimension(:),  intent(in)   :: nTCells0_inLx

    ! return
    real(dp), dimension( size(L0upBound_inLx,1) )  :: frac_cover_Lx

    ! local variables
    integer(i4)                              :: kk, iu, id, jl, jr, nT
    integer(i4)                              :: nrows0, ncols0
    integer(i4), dimension(:,:), allocatable :: dummy_Matrix
    integer(i4), dimension(:,:), allocatable :: nodata_val
    integer(i4)                              :: nCells1

    ! estimate number of cells
    nCells1 =  size(L0upBound_inLx,1)

    ! get nrows and ncols
    nrows0 = size(mask0,1)
    ncols0 = size(mask0,2)

    !unpack input data from 1D to 2D
    allocate( dummy_Matrix(nrows0, ncols0) )
    allocate(   nodata_val(nrows0, ncols0) )
    nodata_val(:,:)   = nodata_i4
    dummy_Matrix(:,:) = unpack( dataIn0(:), mask0(:,:), nodata_val(:,:) )

    ! initalize return variable
    frac_cover_Lx(:) = 0.0_dp

    ! start calculation
    do kk = 1, nCells1
       iu = L0upBound_inLx(kk)
       id = L0downBound_inLx(kk)
       jl = L0leftBound_inLx(kk)
       jr = L0rightBound_inLx(kk)
       nT = nTCells0_inLx(kk)

       frac_cover_Lx(kk) = real( count(dummy_Matrix(iu:id, jl:jr) == classId), dp ) / real(nT, dp)

    end do

    ! free space
    deallocate(dummy_Matrix, nodata_val)

  end function L0_fractionalCover_in_Lx

  ! ----------------------------------------------------------------------------

  !      NAME
  !          upscale_arithmetic_mean

  !>         \brief aritmetic mean

  !>         \details upscaling of level-0 grid data to level-1 using aritmetic mean

  !      INTENT(IN)
  !>         \param[in] "integer(i4) :: nL0_cells_in_L1_cell(:)"  - number of level-0 cells within a level-1 cell
  !>         \param[in] "integer(i4) :: L1_upper_rowId_cell(:)"   - upper row boundary (level-0) of a level-1 cell
  !>         \param[in] "integer(i4) :: L1_lower_rowId_cell(:)"   - lower row boundary (level-0) of a level-1 cell
  !>         \param[in] "integer(i4) :: L1_left_colonId_cell(:)"  - left colon boundary (level-0) of a level-1 cell
  !>         \param[in] "integer(i4) :: L1_right_colonId_cell(:)" - right colon boundary (level-0) of a level-1 cell
  !>         \param[in] "integer(i4) :: L0_cellId(:,:)"           - cell ID at level-0
  !>         \param[in] "logical     :: mask0(:,:)"               - mask at Level 0
  !>         \param[in] "real(dp)    :: nodata_value"             - no data value
  !>         \param[in] "real(dp)    :: L0_fineScale_data(:,:)"- high resolution data

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
  !>       \return     real(dp) :: upscale_arithmetic_mean(:) &mdash; Upscaled variable from L0 to L1 using arithmetic mean

  !     RESTRICTIONS
  !         None

  !      EXAMPLE
  !          variable = upscale_arithmetic_mean( nL0_cells_in_L1_cell,                           &
  !                                              L1_upper_rowId_cell, L1_lower_rowId_cell,       &
  !                                              L1_left_colonId_cell, L1_right_colonId_cell,    &
  !                                              L0_cellId, mask0,                               &
  !                                              nodata_value, L0_fineScale_2D_data              &
  !                                             )

  !     LITERATURE
  !         None

  !      HISTORY
  !>         \author Giovanni Dalmasso, Rohini Kumar
  !>         \date Dec 2012
  !          Written,  Giovanni Dalmasso, Dec 2012
  !          Modified, Stephan Thober,    Feb 2013 - changed dimension of L0 input from 2d to 1d

  function upscale_arithmetic_mean(   &
       nL0_cells_in_L1_cell,   &   ! number of level-0 cells within a level-1 cell
       L1_upper_rowId_cell,    &   ! upper row boundary (level-0) of a level-1 cell
       L1_lower_rowId_cell,    &   ! lower row boundary (level-0) of a level-1 cell
       L1_left_colonId_cell,   &   ! left colon boundary (level-0) of a level-1 cell
       L1_right_colonId_cell,  &   ! right colon boundary (level-0) of a level-1 cell
       L0_cellId,              &   ! cell ID at level-0
       mask0,                  &   ! mask at level 0
       nodata_value,           &   ! no data value
       L0_fineScale_data       &   ! high resolution data
       )

    implicit none
    ! input
    integer(i4), dimension(:),      intent(in)              :: nL0_cells_in_L1_cell
    integer(i4), dimension(:),      intent(in)              :: L1_upper_rowId_cell
    integer(i4), dimension(:),      intent(in)              :: L1_lower_rowId_cell
    integer(i4), dimension(:),      intent(in)              :: L1_left_colonId_cell
    integer(i4), dimension(:),      intent(in)              :: L1_right_colonId_cell
    integer(i4), dimension(:),      intent(in)              :: L0_cellId
    logical,     dimension(:,:),    intent(in)              :: mask0
    real(dp),                       intent(in)              :: nodata_value
    real(dp),    dimension(:),      intent(in)              :: L0_fineScale_data
    ! output
    real(dp), dimension( size(nL0_cells_in_L1_cell,1) )     :: upscale_arithmetic_mean

    ! local variables
    integer(i4)                                             :: L1_nCells
    integer(i4)                                             :: iu, id, jl, jr
    integer(i4)                                             :: kk
    integer(i4), dimension(size(mask0,1),size(mask0,2))     :: nodata_2d
    integer(i4), dimension(size(mask0,1),size(mask0,2))     :: L0_cellId_2d
    real(dp),    dimension(size(mask0,1),size(mask0,2))     :: L0_fineScale_2D_data

    ! allocation and initialisation
    upscale_arithmetic_mean(:)      = 0.0_dp
    nodata_2d            = int(nodata_value, i4)
    L0_cellId_2d         = unpack( L0_cellId, mask0, nodata_2d )
    L0_fineScale_2D_data = unpack( L0_fineScale_data, mask0, nodata_value )

    L1_nCells = size(upscale_arithmetic_mean,1)

    do kk = 1, L1_nCells
       iu = L1_upper_rowId_cell(kk)
       id = L1_lower_rowId_cell(kk)
       jl = L1_left_colonId_cell(kk)
       jr = L1_right_colonId_cell(kk)
       upscale_arithmetic_mean(kk) = sum( L0_fineScale_2D_data(iu:id, jl:jr), L0_cellId_2d(iu:id, jl:jr) /= &
            int(nodata_value,i4) ) / real( nL0_cells_in_L1_cell(kk), dp)
    end do

  end function upscale_arithmetic_mean

  ! ----------------------------------------------------------------------------

  !      NAME
  !          upscale_harmonic_mean

  !>         \brief harmonic mean

  !>         \details upscaling of level-0 grid data to level-1 using harmonic mean

  !       INTENT(IN)
  !>         \param[in] "integer(i4) :: nL0_cells_in_L1_cell(:)"  - number of level-0 cells within a level-1 cell
  !>         \param[in] "integer(i4) :: L1_upper_rowId_cell(:)"   - upper row boundary (level-0) of a level-1 cell
  !>         \param[in] "integer(i4) :: L1_lower_rowId_cell(:)"   - lower row boundary (level-0) of a level-1 cell
  !>         \param[in] "integer(i4) :: L1_left_colonId_cell(:)"  - left colon boundary (level-0) of a level-1 cell
  !>         \param[in] "integer(i4) :: L1_right_colonId_cell(:)" - right colon boundary (level-0) of a level-1 cell
  !>         \param[in] "integer(i4) :: L0_cellId(:,:)"           - cell ID at level-0
  !>         \param[in] "logical     :: mask0(:,:)"               - mask at Level 0
  !>         \param[in] "real(dp)    :: nodata_value"             - no data value
  !>         \param[in] "real(dp)    :: L0_fineScale_2D_data(:,:)"- high resolution data

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
  !>       \return     real(dp) :: upscale_harmonic_mean(:) &mdash; Upscaled variable from L0 to L1 using harmonic mean

  !     RESTRICTIONS
  !>       \note Input values must be floating points.

  !      EXAMPLE
  !          calling sequence
  !          variable = upscale_harmonic_mean(   nL0_cells_in_L1_cell,                           &
  !                                              L1_upper_rowId_cell, L1_lower_rowId_cell,       &
  !                                              L1_left_colonId_cell, L1_right_colonId_cell,    &
  !                                              L0_cellId, mask0,                               &
  !                                              nodata_value, L0_fineScale_2D_data              &
  !                                             )

  !     LITERATURE
  !         None

  !      HISTORY
  !>         \author Giovanni Dalmasso, Rohini Kumar
  !>         \date Dec 2012
  !          Written, Giovanni Dalmasso, Dec 2012
  !          Modified,   Stephan Thober, Jan 2013 - change example calling sequence
  !                      Stephan Thober, Feb 2013 - added Level 0 mask

  function upscale_harmonic_mean( &
       nL0_cells_in_L1_cell,   &   ! number of level-0 cells within a level-1 cell
       L1_upper_rowId_cell,    &   ! upper row boundary (level-0) of a level-1 cell
       L1_lower_rowId_cell,    &   ! lower row boundary (level-0) of a level-1 cell
       L1_left_colonId_cell,   &   ! left colon boundary (level-0) of a level-1 cell
       L1_right_colonId_cell,  &   ! right colon boundary (level-0) of a level-1 cell
       L0_cellId,              &   ! cell ID at level-0
       mask0,                  &   ! mask at Level 0
       nodata_value,           &   ! no data value
       L0_fineScale_data       &   ! high resolution data
       )

    implicit none

    ! input
    integer(i4), dimension(:),      intent(in)              :: nL0_cells_in_L1_cell
    integer(i4), dimension(:),      intent(in)              :: L1_upper_rowId_cell
    integer(i4), dimension(:),      intent(in)              :: L1_lower_rowId_cell
    integer(i4), dimension(:),      intent(in)              :: L1_left_colonId_cell
    integer(i4), dimension(:),      intent(in)              :: L1_right_colonId_cell
    integer(i4), dimension(:),      intent(in)              :: L0_cellId
    logical,     dimension(:,:),    intent(in)              :: mask0
    real(dp),                       intent(in)              :: nodata_value
    real(dp),    dimension(:),      intent(in)              :: L0_fineScale_data
    ! output
    real(dp), dimension( size(nL0_cells_in_L1_cell,1) )     :: upscale_harmonic_mean

    ! local variables
    integer(i4)                                             :: L1_nCells
    integer(i4)                                             :: iu, id, jl, jr
    integer(i4)                                             :: kk
    integer(i4), dimension(size(mask0,1),size(mask0,2))     :: nodata_2d
    integer(i4), dimension(size(mask0,1),size(mask0,2))     :: L0_cellId_2d
    real(dp),    dimension(size(mask0,1),size(mask0,2))     :: L0_fineScale_2D_data

    ! allocation and initialisation
    upscale_harmonic_mean(:) = 0.0_dp
    nodata_2d            = int(nodata_value, i4)
    L0_cellId_2d         = unpack( L0_cellId, mask0, nodata_2d )
    L0_fineScale_2D_data = unpack( L0_fineScale_data, mask0, nodata_value )

    L1_nCells = size(upscale_harmonic_mean,1)

    do kk = 1, L1_nCells
       iu = L1_upper_rowId_cell(kk)
       id = L1_lower_rowId_cell(kk)
       jl = L1_left_colonId_cell(kk)
       jr = L1_right_colonId_cell(kk)
       upscale_harmonic_mean(kk) = real( nL0_cells_in_L1_cell(kk), dp ) &
            / sum( 1.0_dp/L0_fineScale_2D_data(iu:id, jl:jr), L0_cellId_2d(iu:id, jl:jr) /= int(nodata_value,i4) )
    end do

  end function upscale_harmonic_mean

  ! ----------------------------------------------------------------------------

  !      NAME
  !          upscale_geometric_mean

  !>         \brief geometric mean

  !>         \details upscaling of level-0 grid data to level-1 using geometric mean

  !       INTENT(IN)
  !>         \param[in] "integer(i4) :: L1_upper_rowId_cell(:)"   - upper row boundary (level-0) of a level-1 cell
  !>         \param[in] "integer(i4) :: L1_lower_rowId_cell(:)"   - lower row boundary (level-0) of a level-1 cell
  !>         \param[in] "integer(i4) :: L1_left_colonId_cell(:)"  - left colon boundary (level-0) of a level-1 cell
  !>         \param[in] "integer(i4) :: L1_right_colonId_cell(:)" - right colon boundary (level-0) of a level-1 cell
  !>         \param[in] "logical     :: mask0(:,:)"               - mask at Level 0
  !>         \param[in] "real(dp)    :: nodata_value"             - no data value

  !     INTENT(INOUT)
  !>         \param[in,out] real(dp) :: L0_fineScale_2D_data(:,:)  - high resolution data

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       \return     real(dp) :: upscale_geometric_mean(:) &mdash; Upscaled variable from L0 to L1 using geometric mean

  !     RESTRICTIONS
  !         None

  !      EXAMPLE
  !          calling sequence
  !          variable = upscale_geometric_mean(  L1_upper_rowId_cell, L1_lower_rowId_cell,       &
  !                                              L1_left_colonId_cell, L1_right_colonId_cell,    &
  !                                              mask0, nodata_value, L0_fineScale_2D_data       &
  !                                             )

  !     LITERATURE
  !         None

  !      HISTORY
  !>         \author Giovanni Dalmasso, Rohini Kumar
  !>         \date Dec 2012
  !          Written, Giovanni Dalmasso, Dec 2012

  function upscale_geometric_mean( &
       L1_upper_rowId_cell,        &   ! upper row boundary (level-0) of a level-1 cell
       L1_lower_rowId_cell,        &   ! lower row boundary (level-0) of a level-1 cell
       L1_left_colonId_cell,       &   ! left colon boundary (level-0) of a level-1 cell
       L1_right_colonId_cell,      &   ! right colon boundary (level-0) of a level-1 cell
       mask0,                      &   ! mask at level 0
       nodata_value,               &   ! no data value
       L0_fineScale_data           &   ! high resolution data
       )

    implicit none

    ! input
    integer(i4), dimension(:),      intent(in)            :: L1_upper_rowId_cell
    integer(i4), dimension(:),      intent(in)            :: L1_lower_rowId_cell
    integer(i4), dimension(:),      intent(in)            :: L1_left_colonId_cell
    integer(i4), dimension(:),      intent(in)            :: L1_right_colonId_cell
    logical,     dimension(:,:),    intent(in)            :: mask0
    real(dp),                       intent(in)            :: nodata_value
    real(dp),    dimension(:),      intent(in)            :: L0_fineScale_data
    ! output
    real(dp),    dimension( size(L1_upper_rowId_cell,1) ) :: upscale_geometric_mean

    ! loca variables
    integer(i4)                                           :: L1_nCells
    integer(i4)                                           :: iu, id, jl, jr
    integer(i4)                                           :: kk
    integer(i4)                                           :: nCells_L0_in_L1
    real(dp),    dimension(size(mask0,1),size(mask0,2))   :: L0_fineScale_2D_data
    real(dp),    dimension(size(mask0,1),size(mask0,2))   :: nodata_2d

    ! allocation and initialisation
    upscale_geometric_mean(:)      = nodata_value
    nodata_2d                      = nodata_value
    L0_fineScale_2D_data           = unpack( L0_fineScale_data, mask0, nodata_2d )

    L1_nCells = size(upscale_geometric_mean,1)
    where(L0_fineScale_2D_data > 0.0_dp) L0_fineScale_2D_data = log(L0_fineScale_2D_data)

    do kk = 1, L1_nCells
       iu = L1_upper_rowId_cell(kk)
       id = L1_lower_rowId_cell(kk)
       jl = L1_left_colonId_cell(kk)
       jr = L1_right_colonId_cell(kk)
       nCells_L0_in_L1 = count(L0_fineScale_2D_data(iu:id, jl:jr) > 0.0_dp)

       if(nCells_L0_in_L1 > 0_i4) then
          upscale_geometric_mean(kk) = sum(L0_fineScale_2D_data(iu:id, jl:jr),    &
               L0_fineScale_2D_data(iu:id, jl:jr) > 0.0_dp ) / real(nCells_L0_in_L1, dp)
          upscale_geometric_mean(kk) = exp( upscale_geometric_mean(kk) )
       else
          upscale_geometric_mean(kk) = 0.0_dp
       end if
    end do

  end function upscale_geometric_mean


end module mo_upscaling_operators
