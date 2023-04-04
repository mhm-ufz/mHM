!> \file mo_mrm_net_startup.f90
!> \brief \copybrief mo_mrm_net_startup
!> \details \copydetails mo_mrm_net_startup

!> \brief Startup drainage network for mHM.
!> \details This module initializes the drainage network at L11 in mHM.
!! - Delineation of drainage network at level 11.
!! - Setting network topology (i.e. nodes and link).
!! - Determining routing order.
!! - Determining cell locations for network links.
!! - Find drainage outlet.
!! - Determine stream (links) features.
!!
!> \changelog
!! - Rohini Kumar May 2014
!!   - cell area calulation based on a regular lat-lon grid or on a regular X-Y coordinate system
!! - Robert Schweppe Jun 2018
!!   - refactoring and reformatting
!> \authors Luis Samaniego
!> \date Dec 2012
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mrm
module mo_mrm_net_startup
  use mo_kind, only : i4, dp
  use mo_message, only : message, error_message

  implicit none

  PUBLIC :: L11_L1_mapping
  PUBLIC :: L11_flow_direction
  PUBLIC :: L11_set_network_topology
  PUBLIC :: L11_routing_order
  PUBLIC :: L11_link_location
  PUBLIC :: L11_set_drain_outlet_gauges
  PUBLIC :: L11_stream_features
  PUBLIC :: L11_fraction_sealed_floodplain
  PUBLIC :: get_distance_two_lat_lon_points
  PUBLIC :: L11_flow_accumulation
  PUBLIC :: L11_calc_celerity
contains

  !    NAME
  !        L11_L1_mapping

  !    PURPOSE
  !>       \brief TODO: add description

  !>       \details TODO: add description

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iDomain" domain

  !    HISTORY
  !>       \authors Robert Schweppe

  !>       \date Jun 2018

  ! Modifications:

  subroutine L11_L1_mapping(iDomain)

    use mo_append, only : append
    use mo_common_constants, only : nodata_i4
    use mo_common_variables, only : level1
    use mo_mrm_global_variables, only : L11_L1_ID, L1_L11_ID, level11

    implicit none

    ! domain
    integer(i4), intent(in) :: iDomain

    integer(i4) :: nrows1, ncols1

    integer(i4) :: nrows11, ncols11

    integer(i4) :: kk

    integer(i4) :: icc, jcc

    integer(i4) :: iu, id, jl, jr

    ! mapping of L11 Id on L1
    integer(i4), dimension(:, :), allocatable :: L11Id_on_L1

    ! mapping of L1 Id on L11
    integer(i4), dimension(:, :), allocatable :: L1Id_on_L11

    ! dummy ID
    integer(i4), dimension(:, :), allocatable :: dummy_2d_id

    real(dp) :: cellFactorRbyH

    integer(i4) :: cellFactorRbyH_inv


    nrows1 = level1(iDomain)%nrows
    nrows11 = level11(iDomain)%nrows
    ncols1 = level1(iDomain)%ncols
    ncols11 = level11(iDomain)%ncols

    ! allocate variables for mapping L11 Ids and L1 Ids
    allocate (L11Id_on_L1  (nrows1, ncols1))
    allocate (L1Id_on_L11  (nrows11, ncols11))
    L11Id_on_L1(:, :) = nodata_i4
    L1Id_on_L11(:, :) = nodata_i4

    ! set cell factor for routing
    cellFactorRbyH = level11(iDomain)%cellsize / level1(iDomain)%cellsize

    ! set mapping
    ! create mapping between L11 and L1 for L11 resolution higher than L1 resolution
    if (cellFactorRbyH .lt. 1._dp) then
      allocate (dummy_2d_id  (nrows1, ncols1))
      dummy_2d_id = unpack(level1(iDomain)%Id, level1(iDomain)%mask, nodata_i4)
      cellFactorRbyH_inv = int(1. / cellFactorRbyH, i4)
      kk = 0
      do jcc = 1, ncols1
        do icc = 1, nrows1
          if(.not. level1(iDomain)%mask(icc, jcc)) cycle
          kk = kk + 1
          !
          iu = (icc - 1) * cellFactorRbyH_inv + 1
          id = min(icc * cellFactorRbyH_inv, nrows11)
          jl = (jcc - 1) * cellFactorRbyH_inv + 1
          jr = min(jcc * cellFactorRbyH_inv, ncols11)

          L1Id_on_L11(iu : id, jl : jr) = merge(dummy_2d_id(icc, jcc), nodata_i4, level11(iDomain)%mask(iu : id, jl : jr))
        end do
      end do
    else
      allocate (dummy_2d_id  (nrows11, ncols11))
      dummy_2d_id = unpack(level11(iDomain)%Id, level11(iDomain)%mask, nodata_i4)

      kk = 0
      do jcc = 1, ncols11
        do icc = 1, nrows11
          if(.not. level11(iDomain)%mask(icc, jcc)) cycle
          kk = kk + 1

          ! coord. of all corners L11 -> of finer scale level-1
          iu = (icc - 1) * nint(cellFactorRbyH, i4) + 1
          id = icc * nint(cellFactorRbyH, i4)
          jl = (jcc - 1) * nint(cellFactorRbyH, i4) + 1
          jr = jcc * nint(cellFactorRbyH, i4)

          ! constrain the range of up, down, left, and right boundaries
          if(iu < 1) iu = 1
          if(id > nrows1) id = nrows1
          if(jl < 1) jl = 1
          if(jr > ncols1) jr = ncols1

          ! Delimitation of level-11 cells on level-1 for L11 resolution lower than L1 resolution
          L11Id_on_L1(iu : id, jl : jr) = merge(dummy_2d_id(icc, jcc), nodata_i4, level1(iDomain)%mask(iu : id, jl : jr))
        end do
      end do
    end if

    ! L1 data sets
    call append(L1_L11_Id, pack (L11Id_on_L1(:, :), level1(iDomain)%mask))
    ! L11 data sets
    call append(L11_L1_Id, pack (L1Id_on_L11(:, :), level11(iDomain)%mask))
    ! free space
    deallocate(L11Id_on_L1, L1Id_on_L11, dummy_2d_id)

  end subroutine L11_L1_mapping
  ! --------------------------------------------------------------------------

  !    NAME
  !        L11_flow_direction

  !    PURPOSE
  !>       \brief Determine the flow direction of the upscaled river
  !>       network at level L11.

  !>       \details The hydrographs generated at each cell are routed
  !>       through the drainage network at level-11 towards their
  !>       outlets. The drainage network at level-11 is conceptualized as a
  !>       graph whose nodes are hypothetically located at the center of
  !>       each grid cell connected by links that represent the river
  !>       reaches. The flow direction of a link correspond to the
  !>       direction towards a neighboring cell in which the net flow
  !>       accumulation (outflows minus inflows) attains its maximum
  !>       value. The net flow accumulation across a cell's boundary at
  !>       level-11 is estimated based on flow direction and flow
  !>       accumulation obtained at level-0 (\ref fig_routing "Routing
  !>       Network"). Note: level-1 denotes the modeling level, whereas
  !>       level-L11 is at least as coarse as level-1. Experience has
  !>       shown that routing can be done at a coarser resolution as
  !>       level-1, hence the level-11 was introduced.
  !>       \image html  routing.png "Upscaling routing network from L0 to L1 (or L11)"
  !>       \anchor fig_routing \image latex routing.pdf "Upscaling routing network from L0 to L1 (or L11)" width=14cm
  !>       The left panel depicts a schematic derivation of a drainage
  !>       network at the level-11 based on level-0 flow direction and
  !>       flow accumulation. The dotted line circle denotes the point
  !>       with the highest flow accumulation within a grid cell. The
  !>       topology of a tipical drainage routing network at level-11 is
  !>       shown in the right panel. Gray color areas denote the flood
  !>       plains estimated in mo_init_mrm, where the network
  !>       upscaling is also carried out.
  !>       For the sake of simplicity, it is assumed that all runoff leaving
  !>       a given cell would exit through a major direction.
  !>       Note that multiple outlets can exist within the modelling domain.
  !>       If a variable is added or removed here, then it also has to
  !>       be added or removed in the subroutine L11_config_set in
  !>       module mo_restart and in the subroutine set_L11_config in module
  !>       mo_set_netcdf_restart
  !>       ADDITIONAL INFORMATION
  !>       L11_flow_direction

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iDomain" Domain Id

  !    HISTORY
  !>       \authors Luis Samaniego

  !>       \date Dec 2005

  ! Modifications:
  ! Luis Samaniego Jan 2013 - modular version
  ! Rohini Kumar   Apr 2014 - Case of L0 is same as L11 implemented
  ! Stephan Thober Aug 2015 - ported to mRM
  ! Stephan Thober Sep 2015 - create mapping between L11 and L1 if L11 resolution is higher than L1 resolution
  ! Stephan Thober May 2016 - introducing multiple outlets
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine L11_flow_direction(iDomain)

    use mo_append, only : append
    use mo_common_constants, only : nodata_i4
    use mo_common_types, only: Grid
    use mo_common_variables, only : domainMeta, level0
    use mo_mrm_global_variables, only : L0_draSC, L0_fAcc, L0_fDir, L0_l11_remap, L11_colOut, L11_fDir, &
                                        L11_nOutlets, L11_rowOut, domain_mrm, level11
    use mo_string_utils, only : num2str

    implicit none

    ! Domain Id
    integer(i4), intent(in) :: iDomain

    integer(i4) :: nCells0

    integer(i4) :: nrows0, ncols0

    integer(i4) :: s0, e0

    integer(i4) :: nrows11, ncols11

    ! =  ncells11
    integer(i4) :: nNodes

    integer(i4) :: ii, jj, kk, ic, jc

    integer(i4) :: iu, id

    integer(i4) :: jl, jr

    integer(i4) :: iRow, jCol

    integer(i4), dimension(:, :), allocatable :: iD0

    integer(i4), dimension(:, :), allocatable :: fDir0

    integer(i4), dimension(:, :), allocatable :: fAcc0

    integer(i4), dimension(:, :), allocatable :: fDir11

    ! northing cell loc. of the Outlet
    integer(i4), dimension(:), allocatable :: rowOut

    ! easting cell loc. of the Outlet
    integer(i4), dimension(:), allocatable :: colOut

    integer(i4), dimension(:, :), allocatable :: draSC0

    ! output location in L0
    integer(i4), dimension(:, :), allocatable :: oLoc

    integer(i4) :: side

    integer(i4) :: fAccMax, idMax

    ! Number of outlet found
    integer(i4) :: Noutlet

    ! Number of outlets before this Domain
    integer(i4) :: old_Noutlet

    ! flag whether outlet is found
    logical :: is_outlet

    type(Grid), pointer :: level0_iDomain => null()


    !--------------------------------------------------------
    ! STEPS:
    ! 1) Estimate each variable locally for a given Domain
    ! 2) Pad each variable to its corresponding global one
    !--------------------------------------------------------
    ! initialize
    Noutlet = 0_i4
    level0_iDomain => level0(domainMeta%L0DataFrom(iDomain))
    !------------------------------------------------------------------
    !                Set Flow Direction at Level 11
    !                       Searching order
    !                             jl    jr
    !                          iu  +  4  +
    !                              3     1
    !                          id  +  2  +
    !------------------------------------------------------------------

    ! flow direction at level-11

    ! allocate
    nrows0 = level0_iDomain%nrows
    ncols0 = level0_iDomain%ncols
    nCells0 = level0_iDomain%ncells
    nrows11 = level11(iDomain)%nrows
    ncols11 = level11(iDomain)%ncols
    nNodes = level11(iDomain)%ncells
    s0 = level0_iDomain%iStart
    e0 = level0_iDomain%iEnd

    allocate (iD0         (nrows0, ncols0))
    allocate (fAcc0       (nrows0, ncols0))
    allocate (fDir0       (nrows0, ncols0))
    allocate (draSC0      (nrows0, ncols0))
    allocate (fDir11      (nrows11, ncols11))
    allocate (rowOut      (nNodes))
    allocate (colOut      (nNodes))
    allocate (oLoc        (1, 2))

    ! initialize
    iD0(:, :) = nodata_i4
    fAcc0(:, :) = nodata_i4
    fDir0(:, :) = nodata_i4
    draSC0(:, :) = nodata_i4
    fDir11(:, :) = nodata_i4
    rowOut(:) = nodata_i4
    colOut(:) = nodata_i4
    oLoc(:, :) = nodata_i4


    ! get iD, fAcc, fDir at L0
    iD0(:, :) = UNPACK(level0_iDomain%Id, level0_iDomain%mask, nodata_i4)
    fAcc0(:, :) = UNPACK(L0_fAcc (s0 : e0), level0_iDomain%mask, nodata_i4)
    fDir0(:, :) = UNPACK(L0_fDir (s0 : e0), level0_iDomain%mask, nodata_i4)

    ! case where routing and input data scale is similar
    IF(nCells0 .EQ. nNodes) THEN
      oLoc(1, :) = maxloc(fAcc0, level0_iDomain%mask)
      kk = L0_L11_remap(iDomain)%lowres_id_on_highres(oLoc(1, 1), oLoc(1, 2))
      ! for a single node model run
      if(nCells0 .EQ. 1) then
        fDir11(1, 1) = fDir0(oLoc(1, 1), oLoc(1, 2))
      else
        fDir11(:, :) = fDir0(:, :)
      end if
      fDir11 (level11(iDomain)%CellCoor(kk, 1), level11(iDomain)%CellCoor(kk, 2)) = 0
      ! set location of main outlet in L11
      do kk = 1, nNodes
        ii = level11(iDomain)%CellCoor(kk, 1)
        jj = level11(iDomain)%CellCoor(kk, 2)
        rowOut(kk) = ii
        colOut(kk) = jj
      end do
      do kk = 1, ncells0
        ii = level0_iDomain%CellCoor(kk, 1)
        jj = level0_iDomain%CellCoor(kk, 2)
        draSC0(ii, jj) = kk
      end do

      ! case where routing and input data scale differs
    ELSE
      ! =======================================================================
      ! ST: find all cells whose downstream cells are outside the domain
      ! =======================================================================
      do ii = 1, nCells0
        iRow = level0_iDomain%CellCoor(ii, 1)
        jCol = level0_iDomain%CellCoor(ii, 2)
        call moveDownOneCell(fDir0(iRow, jCol), iRow, jCol)
        ! check whether new location is inside bound
        is_outlet = .False.
        if ((iRow .le. 0_i4) .or. (iRow .gt. nrows0) .or. &
            (jCol .le. 0_i4) .or. (jCol .gt. ncols0)) then
          is_outlet = .True.
        else
          if (fdir0(iRow, jCol) .le. 0) is_outlet = .True.
        end if
        !
        if (is_outlet) then
          Noutlet = Noutlet + 1_i4
          ! cell is an outlet
          if (Noutlet .eq. 1) then
            oLoc(1, :) = level0_iDomain%CellCoor(ii, :)
          else
            call append(oLoc, level0_iDomain%CellCoor(ii : ii, :))
          end if
          ! drain this cell into corresponding L11 cell
          kk = L0_L11_remap(iDomain)%lowres_id_on_highres(oLoc(Noutlet, 1), oLoc(Noutlet, 2))
          draSC0(oLoc(Noutlet, 1), oLoc(Noutlet, 2)) = kk
          ! check whether cell has maximum flow accumulation
          ! coord. of all corners
          iu = l0_l11_remap(iDomain)%upper_bound   (kk)
          id = l0_l11_remap(iDomain)%lower_bound (kk)
          jl = l0_l11_remap(iDomain)%left_bound (kk)
          jr = l0_l11_remap(iDomain)%right_bound(kk)
          if (maxval(facc0(iu : id, jl : jr)) .eq. facc0(oLoc(Noutlet, 1), oLoc(Noutlet, 2)))  then
            ! set location of outlet at L11
            rowOut(kk) = oLoc(Noutlet, 1)
            colOut(kk) = oLoc(Noutlet, 2)
            fdir11(level11(iDomain)%CellCoor(kk, 1), level11(iDomain)%CellCoor(kk, 2)) = 0
          end if
        end if
      end do

      ! finding cell L11 outlets -  using L0_fAcc

      do kk = 1, nNodes

        ! exclude outlet L11
        if (rowOut(kk) > 0) cycle

        ic = level11(iDomain)%CellCoor(kk, 1)
        jc = level11(iDomain)%CellCoor(kk, 2)

        ! coord. of all corners
        iu = l0_l11_remap(iDomain)%upper_bound (kk)
        id = l0_l11_remap(iDomain)%lower_bound (kk)
        jl = l0_l11_remap(iDomain)%left_bound (kk)
        jr = l0_l11_remap(iDomain)%right_bound(kk)

        fAccMax = -9
        idMax = 0
        side = -1
        ! searching on side 4
        do jj = jl, jr
          if ((fAcc0(iu, jj) > fAccMax)  .and. &
              (fDir0(iu, jj) ==  32 .or.  &
              fDir0(iu, jj) ==  64 .or.  &
              fDir0(iu, jj) == 128)) then
            fAccMax = fAcc0(iu, jj)
            idMax = id0(iu, jj)
            side = 4
          end if
        end do

        ! searching on side 1
        do ii = iu, id
          if ((fAcc0(ii, jr) > fAccMax)  .and. &
              (fDir0(ii, jr) ==   1 .or.  &
              fDir0(ii, jr) ==   2 .or.  &
              fDir0(ii, jr) == 128)) then
            fAccMax = fAcc0(ii, jr)
            idMax = id0(ii, jr)
            side = 1
          end if
        end do

        ! searching on side 2
        do jj = jl, jr
          if ((fAcc0(id, jj) > fAccMax)  .and. &
              (fDir0(id, jj) ==   2 .or.  &
              fDir0(id, jj) ==   4 .or.  &
              fDir0(id, jj) ==   8)) then
            fAccMax = fAcc0(id, jj)
            idMax = id0(id, jj)
            side = 2
          end if
        end do

        ! searching on side 3
        do ii = iu, id
          if ((fAcc0(ii, jl) > fAccMax)  .and. &
              (fDir0(ii, jl) ==   8 .or.  &
              fDir0(ii, jl) ==  16 .or.  &
              fDir0(ii, jl) ==  32)) then
            fAccMax = fAcc0(ii, jl)
            idMax = id0(ii, jl)
            side = 3
          end if
        end do

        ! set location of the cell-outlet (row, col) in L0
        ii = level0_iDomain%CellCoor(idMax, 1)
        jj = level0_iDomain%CellCoor(idMax, 2)
        rowOut(kk) = ii
        colOut(kk) = jj
        draSC0(ii, jj) = kk

        ! set fDir at L11
        if     (ii == iu .and.  jj == jl) then
          select case (fDir0(ii, jj))
          case (8, 16)
            fDir11(ic, jc) = 16
          case (32)
            fDir11(ic, jc) = 32
          case (64, 128)
            fDir11(ic, jc) = 64
          end select
        elseif (ii == iu .and.  jj == jr) then
          select case (fDir0(ii, jj))
          case (32, 64)
            fDir11(ic, jc) = 64
          case (128)
            fDir11(ic, jc) = 128
          case (1, 2)
            fDir11(ic, jc) = 1
          end select
        elseif (ii == id .and.  jj == jl) then
          select case (fDir0(ii, jj))
          case (2, 4)
            fDir11(ic, jc) = 4
          case (8)
            fDir11(ic, jc) = 8
          case (16, 32)
            fDir11(ic, jc) = 16
          end select
        elseif (ii == id .and.  jj == jr) then
          select case (fDir0(ii, jj))
          case (128, 1)
            fDir11(ic, jc) = 1
          case (2)
            fDir11(ic, jc) = 2
          case (4, 8)
            fDir11(ic, jc) = 4
          end select
        else
          ! cell on one side
          select case (side)
          case (1)
            fDir11(ic, jc) = 1
          case (2)
            fDir11(ic, jc) = 4
          case (3)
            fDir11(ic, jc) = 16
          case (4)
            fDir11(ic, jc) = 64
          case default
             call error_message('Error L11_flow_direction: side = -1')
          end select
        end if

      end do

    END IF
    !--------------------------------------------------------
    ! Start padding up local variables to global variables
    !--------------------------------------------------------

    ! allocate space for row and col Outlet
    allocate(domain_mrm(iDomain)%L0_rowOutlet(1))
    allocate(domain_mrm(iDomain)%L0_colOutlet(1))
    domain_mrm(iDomain)%L0_Noutlet = nodata_i4
    domain_mrm(iDomain)%L0_rowOutlet = nodata_i4
    domain_mrm(iDomain)%L0_colOutlet = nodata_i4

    ! L0 data sets
    ! check whether L0 data is shared
    !ToDo: check if change is correct
    if (iDomain .eq. 1) then
      call append(L0_draSC, PACK (draSC0(:, :), level0_iDomain%mask))
    else if (domainMeta%L0DataFrom(iDomain) == iDomain) then
      call append(L0_draSC, PACK (draSC0(:, :), level0_iDomain%mask))
    end if

    domain_mrm(iDomain)%L0_Noutlet = Noutlet
    ! set L0 outlet coordinates
    old_Noutlet = size(domain_mrm(iDomain)%L0_rowOutlet, dim = 1)
    if (Noutlet .le. old_Noutlet) then
      domain_mrm(iDomain)%L0_rowOutlet(: Noutlet) = oLoc(:, 1)
      domain_mrm(iDomain)%L0_colOutlet(: Noutlet) = oLoc(:, 2)
    else
      ! store up to size of old_Noutlet
      domain_mrm(iDomain)%L0_rowOutlet(: old_Noutlet) = oLoc(: old_Noutlet, 1)
      domain_mrm(iDomain)%L0_colOutlet(: old_Noutlet) = oLoc(: old_Noutlet, 2)
      ! enlarge rowOutlet and colOutlet in domain_mrm structure
      !TODO: do other domains also need to be enlarged accordingly???
      call append(domain_mrm(iDomain)%L0_rowOutlet, oLoc(old_Noutlet + 1 :, 1))
      call append(domain_mrm(iDomain)%L0_colOutlet, oLoc(old_Noutlet + 1 :, 2))
    end if

    ! L11 data sets
    call append(L11_nOutlets, count(fdir11 .eq. 0_i4))
    call append(L11_fDir, PACK (fDir11(:, :), level11(iDomain)%mask))
    call append(L11_rowOut, rowOut(:))
    call append(L11_colOut, colOut(:))

    ! communicate
    call message('')
    call message('    Domain: ' // num2str(iDomain, '(i3)'))
    call message('      Number of outlets found at Level 0:.. ' // num2str(Noutlet, '(i7)'))
    call message('      Number of outlets found at Level 11:. ' // num2str(count(fdir11 .eq. 0_i4), '(i7)'))

    ! free space
    deallocate(fDir0, fAcc0, fDir11, rowOut, colOut, draSC0)

  end subroutine L11_flow_direction

  ! ------------------------------------------------------------------

  !    NAME
  !        L11_set_network_topology

  !    PURPOSE
  !>       \brief Set network topology

  !>       \details Set network topology from and to node for all links
  !>       at level-11 (\ref fig_routing "Routing Network").
  !>       If a variable is added or removed here, then it also has to
  !>       be added or removed in the subroutine L11_config_set in
  !>       module mo_restart and in the subroutine set_L11_config in module
  !>       mo_set_netcdf_restart.

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iDomain" Domain Id

  !    HISTORY
  !>       \authors Luis Samaniego

  !>       \date Dec 2005

  ! Modifications:
  ! Luis Samaniego Jan 2013 - modular version
  ! Stephan Thober Aug 2015 - ported to mRM
  ! Stephan Thober May 2016 - moved calculation of sink here
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine L11_set_network_topology(iDomain)

    use mo_append, only : append
    use mo_common_constants, only : nodata_i4
    use mo_mrm_global_variables, only : L11_fDir, L11_fromN, L11_toN, level11

    implicit none

    ! Domain Id
    integer(i4), intent(in) :: iDomain

    integer(i4), dimension(:, :), allocatable :: fDir11

    integer(i4), dimension(:, :), allocatable :: dummy_2d_id

    integer(i4) :: jj, kk, ic, jc

    integer(i4) :: fn, tn

    integer(i4), dimension(:), allocatable :: nLinkFromN, nLinkToN


    !     Routing network vectors have nNodes size instead of nLinks to
    !     avoid the need of having two extra indices to identify a Domain.
    ! allocate
    allocate (nLinkFromN (level11(iDomain)%nCells))  ! valid from (1 : nLinks)
    allocate (nLinkToN   (level11(iDomain)%nCells))  ! "
    allocate (fDir11     (level11(iDomain)%nrows, level11(iDomain)%ncols))
    allocate (dummy_2d_id(level11(iDomain)%nrows, level11(iDomain)%ncols))
    dummy_2d_id = unpack(level11(iDomain)%Id, level11(iDomain)%mask, nodata_i4)


    ! initialize
    nLinkFromN(:) = nodata_i4
    nLinkToN(:) = nodata_i4
    fDir11(:, :) = nodata_i4

    ! get grids of L11
    fDir11(:, :) = UNPACK(L11_fDir (level11(iDomain)%iStart : level11(iDomain)%iEnd), level11(iDomain)%mask, nodata_i4)

    ! ------------------------------------------------------------------
    !  network topology
    ! ------------------------------------------------------------------

    jj = 0
    do kk = 1, level11(iDomain)%nCells
      ic = level11(iDomain)%CellCoor(kk, 1)
      jc = level11(iDomain)%CellCoor(kk, 2)
      fn = kk
      call moveDownOneCell(fDir11(ic, jc), ic, jc)
      tn = dummy_2d_id(ic, jc)
      if (fn == tn) cycle
      jj = jj + 1
      nLinkFromN(jj) = fn
      nLinkToN(jj) = tn
    end do

    !--------------------------------------------------------
    ! Start padding up local variables to global variables
    !--------------------------------------------------------

    ! L11 data sets
    call append(L11_fromN, nLinkFromN(:)) ! sinks are at the end
    call append(L11_toN, nLinkToN(:))

    ! free space
    deallocate (fDir11, nLinkFromN, nLinkToN)

  end subroutine L11_set_network_topology

  ! ------------------------------------------------------------------

  !    NAME
  !        L11_routing_order

  !    PURPOSE
  !>       \brief Find routing order, headwater cells and sink

  !>       \details Find routing order, headwater cells and sink.
  !>       If a variable is added or removed here, then it also has to
  !>       be added or removed in the subroutine L11_config_set in
  !>       module mo_restart and in the subroutine set_L11_config in module
  !>       mo_set_netcdf_restart

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iDomain" Domain Id

  !    HISTORY
  !>       \authors Luis Samaniego

  !>       \date Dec 2005

  ! Modifications:
  ! Luis Samaniego Jan 2013 - modular version
  ! Sa. Ku.        Jan 2015 - corrected initialization of nLinkSink
  ! Stephan Thober Aug 2015 - ported to mRM
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine L11_routing_order(iDomain)

    use mo_append, only : append
    use mo_common_constants, only : nodata_i4
    use mo_mrm_global_variables, only : L11_fDir, L11_fromN, L11_label, L11_nOutlets, L11_netPerm, L11_rOrder, L11_sink, L11_toN, &
                                        level11

    implicit none

    ! Domain Id
    integer(i4), intent(in) :: iDomain

    integer(i4) :: nLinks

    ! from node
    integer(i4), dimension(:), allocatable :: nLinkFromN

    ! to node
    integer(i4), dimension(:), allocatable :: nLinkToN

    ! network routing order
    integer(i4), dimension(:), allocatable :: nLinkROrder

    ! label Id [0='', 1=HeadWater, 2=Sink]
    integer(i4), dimension(:), allocatable :: nLinkLabel

    ! == .true. if sink node reached
    logical, dimension(:), allocatable :: nLinkSink

    ! routing order (permutation)
    integer(i4), dimension(:), allocatable :: netPerm

    integer(i4) :: ii, jj, kk

    logical :: flag


    nLinks = level11(iDomain)%nCells - L11_nOutlets(iDomain)
    !  Routing network vectors have nNodes size instead of nLinks to
    !  avoid the need of having two extra indices to identify a Domain.

    ! allocate
    allocate (nLinkFromN  (level11(iDomain)%nCells))  ! all vectors valid from (1 : nLinks)
    allocate (nLinkToN    (level11(iDomain)%nCells))
    allocate (nLinkROrder (level11(iDomain)%nCells))
    allocate (nLinkLabel  (level11(iDomain)%nCells))
    allocate (nLinkSink   (level11(iDomain)%nCells))
    allocate (netPerm     (level11(iDomain)%nCells))
    ! initialize
    nLinkFromN(:) = nodata_i4
    nLinkToN(:) = nodata_i4
    nLinkROrder(1 : nLinks) = 1
    nLinkROrder(level11(iDomain)%nCells) = nodata_i4
    nLinkLabel(1 : nLinks) = 0
    nLinkLabel(level11(iDomain)%nCells) = nodata_i4
    nLinkSink(:) = .FALSE.
    netPerm(:) = nodata_i4

    ! for a single node model run
    if(level11(iDomain)%nCells .GT. 1) then
      ! get network vectors of L11
      nLinkFromN(:) = L11_fromN (level11(iDomain)%iStart : level11(iDomain)%iEnd)
      nLinkToN(:) = L11_toN   (level11(iDomain)%iStart : level11(iDomain)%iEnd)

      loop1 : do ii = 1, nLinks
        loop2 : do jj = 1, nLinks
          if (jj == ii) cycle loop2
          if (nLinkFromN(ii) == nLinkToN(jj)) then
            nLinkROrder(ii) = -9
          end if
          if (nLinkROrder(ii) == -9) cycle loop1
        end do loop2
      end do loop1
      ! counting headwaters
      kk = 0
      do ii = 1, nLinks
        if (nLinkROrder(ii) == 1) then
          kk = kk + 1
          nLinkROrder(ii) = kk
          nLinkLabel(ii) = 1  ! 'Head Water'
        end if
      end do
      ! counting downstream
      do while (minval(nLinkROrder(1 : nLinks)) < 0)
        !!  print *, count(nLinkROrder .lt. 0), minval(nLinkROrder)
        loop3 : do ii = 1, nLinks
          if (.NOT. nLinkROrder(ii) == -9) cycle loop3
          flag = .TRUE.
          loop4 : do jj = 1, nLinks
            if (jj == ii .OR. nLinkFromN(ii)  /=  nLinkToN(jj)) then
              cycle loop4
            else if (.NOT. (nLinkFromN(ii)  == nLinkToN(jj)  .AND. nLinkROrder(jj) > 0)) then
              flag = .FALSE.
              exit loop4
            else
            end if
          end do loop4

          if (flag) then
            kk = kk + 1
            nLinkROrder(ii) = kk
          end if
        end do loop3
      end do

      ! identify sink cells
      do ii = 1, nLinks
        if (L11_fdir(level11(iDomain)%iStart + nLinkToN(ii) - 1_i4) .eq. 0_i4) nlinksink(ii) = .True.
      end do
      where(nlinksink) nLinkLabel = 2 !  'Sink'

      ! keep routing order
      do ii = 1, nLinks
        netPerm(nLinkROrder(ii)) = ii
      end do

      ! end of multi-node network design loop
    end if

    !--------------------------------------------------------
    ! Start padding up local variables to global variables
    !--------------------------------------------------------
    ! L11 network data sets
    call append(L11_rOrder, nLinkROrder(:))
    call append(L11_label, nLinkLabel(:))
    call append(L11_sink, nLinkSink(:))
    call append(L11_netPerm, netPerm(:))

    ! free space
    deallocate (nLinkFromN, nLinkToN, nLinkROrder, nLinkLabel, nLinkSink, netPerm)

  end subroutine L11_routing_order

  ! ------------------------------------------------------------------

  !    NAME
  !        L11_link_location

  !    PURPOSE
  !>       \brief Estimate the LO (row,col) location for each routing link at level L11

  !>       \details If a variable is added or removed here, then it also has to
  !>       be added or removed in the subroutine L11_config_set in
  !>       module mo_restart and in the subroutine set_L11_config in module
  !>       mo_set_netcdf_restart

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iDomain" Domain Id

  !    HISTORY
  !>       \authors Luis Samaniego

  !>       \date Dec 2005

  ! Modifications:
  ! Luis Samaniego Jan 2013 - modular version
  ! Stephan Thober Aug 2015 - ported to mRM
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine L11_link_location(iDomain)

    use mo_append, only : append
    use mo_common_constants, only : nodata_i4
    use mo_common_types, only: Grid
    use mo_common_variables, only : domainMeta, level0
    use mo_mrm_global_variables, only : L0_draSC, L0_fDir, L11_colOut, L11_fCol, L11_fRow, L11_fromN, &
                                        L11_nOutlets, L11_netPerm, L11_rowOut, L11_tCol, L11_tRow, domain_mrm, level11
    use mo_string_utils, only : num2str

    implicit none

    ! Domain Id
    integer(i4), intent(in) :: iDomain

    integer(i4) :: nLinks

    ! northing cell loc. of the Outlet
    integer(i4), dimension(:), allocatable :: rowOut

    ! easting cell loc. of the Outlet
    integer(i4), dimension(:), allocatable :: colOut

    integer(i4), dimension(:), allocatable :: nLinkFromN

    integer(i4), dimension(:), allocatable :: netPerm

    integer(i4), dimension(:), allocatable :: nLinkFromRow

    integer(i4), dimension(:), allocatable :: nLinkFromCol

    integer(i4), dimension(:), allocatable :: nLinkToRow

    integer(i4), dimension(:), allocatable :: nLinkToCol

    integer(i4), dimension(:, :), allocatable :: fDir0

    integer(i4), dimension(:, :), allocatable :: draSC0

    integer(i4) :: ii, rr, kk, s0, e0

    integer(i4) :: iNode, iRow, jCol, prevRow, prevCol

    ! output location in L0
    integer(i4), dimension(:, :), allocatable :: oLoc

    ! number of outlets in Domain
    integer(i4) :: nOutlets

    ! flag for finding outlet
    logical :: is_outlet

    type(Grid), pointer :: level0_iDomain => null()


    level0_iDomain => level0(domainMeta%L0DataFrom(iDomain))
    s0 = level0_iDomain%iStart
    e0 = level0_iDomain%iEnd

    nOutlets = L11_nOutlets(iDomain)

    nLinks = level11(iDomain)%nCells - nOutlets

    !  Routing network vectors have level11(iDomain)%nCells size instead of nLinks to
    !  avoid the need of having two extra indices to identify a Domain.
    ! allocate
    allocate (rowOut        (level11(iDomain)%nCells))
    allocate (colOut        (level11(iDomain)%nCells))
    allocate (nLinkFromN    (level11(iDomain)%nCells))  ! all network vectors valid from (1 : nLinks)
    allocate (netPerm       (level11(iDomain)%nCells))
    allocate (nLinkFromRow  (level11(iDomain)%nCells))
    allocate (nLinkFromCol  (level11(iDomain)%nCells))
    allocate (nLinkToRow    (level11(iDomain)%nCells))
    allocate (nLinkToCol    (level11(iDomain)%nCells))
    allocate (fDir0         (level0_iDomain%nrows, level0_iDomain%ncols))
    allocate (draSC0        (level0_iDomain%nrows, level0_iDomain%ncols))

    ! initialize
    rowOut = nodata_i4
    colOut = nodata_i4
    nLinkFromN = nodata_i4
    netPerm = nodata_i4
    nLinkFromRow = nodata_i4
    nLinkFromCol = nodata_i4
    nLinkToRow = nodata_i4
    nLinkToCol = nodata_i4
    fDir0 = nodata_i4
    draSC0 = nodata_i4

    ! for a single node model run
    if(level11(iDomain)%nCells .GT. 1) then
      ! get fDir at L0
      fDir0(:, :) = UNPACK(L0_fDir  (s0 : e0), level0_iDomain%mask, nodata_i4)
      draSC0(:, :) = UNPACK(L0_draSC (s0 : e0), level0_iDomain%mask, nodata_i4)

      ! get network vectors of L11
      nLinkFromN(:) = L11_fromN   (level11(iDomain)%iStart : level11(iDomain)%iEnd)
      netPerm(:) = L11_netPerm (level11(iDomain)%iStart : level11(iDomain)%iEnd)
      rowOut(:) = L11_rowOut  (level11(iDomain)%iStart : level11(iDomain)%iEnd)
      colOut(:) = L11_colOut  (level11(iDomain)%iStart : level11(iDomain)%iEnd)

      ! finding main outlet (row, col) in L0
      allocate(oLoc(Noutlets, 2))
      oLoc(:, 1) = domain_mrm(iDomain)%L0_rowOutlet(: Noutlets)
      oLoc(:, 2) = domain_mrm(iDomain)%L0_colOutlet(: Noutlets)

      ! Location of the stream-joint cells  (row, col)
      do rr = 1, nLinks

        ii = netPerm(rr)
        iNode = nLinkFromN(ii)
        iRow = rowOut(iNode)
        jCol = colOut(iNode)
        call moveDownOneCell(fDir0(iRow, jcol), iRow, jcol)
        ! set "from" cell
        nLinkFromRow(ii) = iRow
        nLinkFromCol(ii) = jCol

        ! check whether this location is an outlet
        is_outlet = .False.
        do kk = 1, Noutlets
          if (iRow .eq. oLoc(kk, 1) .and. jCol .eq. oLoc(kk, 2)) is_outlet = .True.
        end do

        if (is_outlet) then

          nLinkToRow(ii) = iRow
          nLinkToCol(ii) = jCol

        else

          do while (.not. (draSC0(iRow, jCol) > 0))
            prevRow = iRow
            prevCol = jCol
            call moveDownOneCell(fDir0(iRow, jcol), iRow, jCol)
            ! check whether this location is an outlet and exit
            do kk = 1, Noutlets
              if (iRow .eq. oLoc(kk, 1) .and. jCol .eq. oLoc(kk, 2)) exit
            end do
            if (prevRow .eq. iRow .and. prevCol .eq. jCol) then
              call error_message('Something went wrong during L11_link_location, ', &
                  'movedownonecell got stuck in infinite loop at cell (', num2str(iRow), ' ', &
                  num2str(jCol))
            end if
          end do
          ! set "to" cell (when an outlet is reached)
          nLinkToRow(ii) = iRow
          nLinkToCol(ii) = jCol

        end if
      end do

      ! end of multi-node network design loop
    end if

    !--------------------------------------------------------
    ! Start padding up local variables to global variables
    !--------------------------------------------------------
    ! L11 network data sets
    call append(L11_fRow, nLinkFromRow(:))
    call append(L11_fCol, nLinkFromCol(:))
    call append(L11_tRow, nLinkToRow(:))
    call append(L11_tCol, nLinkToCol(:))

    ! free space
    deallocate (rowOut, colOut, nLinkFromN, netPerm, nLinkFromRow, &
            nLinkFromCol, nLinkToRow, nLinkToCol, fDir0, draSC0)

  end subroutine L11_link_location

  ! ------------------------------------------------------------------

  !    NAME
  !        L11_set_drain_outlet_gauges

  !    PURPOSE
  !>       \brief Draining cell identification and Set gauging node

  !>       \details Perform the following tasks:
  !>       - Draining cell identification (cell at L0 to draining cell outlet at L11).
  !>       - Set gauging nodes
  !>       If a variable is added or removed here, then it also has to
  !>       be added or removed in the subroutine L11_config_set in
  !>       module mo_restart and in the subroutine set_L11_config in module
  !>       mo_set_netcdf_restart

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iDomain" Domain Id

  !    HISTORY
  !>       \authors Luis Samaniego

  !>       \date Dec 2005

  ! Modifications:
  ! Luis Samaniego Jan 2013 - modular version
  ! Matthias Zink  Mar 2014 - bugfix, added inflow gauge
  ! Rohini Kumar   Apr 2014 - variable index is changed to index_gauge
  ! Stephan Thober Aug 2015 - ported to mRM
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine L11_set_drain_outlet_gauges(iDomain)

    use mo_append, only : append
    use mo_common_constants, only : nodata_i4
    use mo_common_types, only: Grid
    use mo_common_variables, only : domainMeta, level0
    use mo_mrm_global_variables, only : L0_InflowgaugeLoc, L0_draCell, L0_draSC, L0_fDir, L0_gaugeLoc, domain_mrm, &
                                        l0_l11_remap

    implicit none

    ! Domain Id
    integer(i4), intent(in) :: iDomain

    integer(i4), dimension(:, :), allocatable :: draSC0

    integer(i4), dimension(:, :), allocatable :: fDir0

    integer(i4), dimension(:, :), allocatable :: gaugeLoc0

    integer(i4), dimension(:, :), allocatable :: InflowGaugeLoc0

    integer(i4), dimension(:, :), allocatable :: draCell0

    integer(i4) :: ii, jj, kk, ll, s0, e0

    integer(i4) :: iSc

    integer(i4) :: iRow, jCol

    type(Grid), pointer :: level0_iDomain => null()


    level0_iDomain => level0(domainMeta%L0DataFrom(iDomain))
    s0 = level0_iDomain%iStart
    e0 = level0_iDomain%iEnd


    ! allocate
    allocate (draSC0          (level0_iDomain%nrows, level0_iDomain%ncols))
    allocate (fDir0           (level0_iDomain%nrows, level0_iDomain%ncols))
    allocate (gaugeLoc0       (level0_iDomain%nrows, level0_iDomain%ncols))
    allocate (InflowGaugeLoc0 (level0_iDomain%nrows, level0_iDomain%ncols))
    allocate (draCell0        (level0_iDomain%nrows, level0_iDomain%ncols))

    ! initialize
    draSC0(:, :) = nodata_i4
    fDir0(:, :) = nodata_i4
    gaugeLoc0(:, :) = nodata_i4
    InflowGaugeLoc0(:, :) = nodata_i4
    draCell0(:, :) = nodata_i4

    draSC0(:, :) = UNPACK(L0_draSC          (s0 : e0), &
            level0_iDomain%mask, nodata_i4)
    fDir0(:, :) = UNPACK(L0_fDir           (s0 : e0), &
            level0_iDomain%mask, nodata_i4)
    gaugeLoc0(:, :) = UNPACK(L0_gaugeLoc       (s0 : e0), &
            level0_iDomain%mask, nodata_i4)
    InflowGaugeLoc0(:, :) = UNPACK(L0_InflowgaugeLoc (s0 : e0), &
            level0_iDomain%mask, nodata_i4)

    do kk = 1, level0_iDomain%nCells
      ii = level0_iDomain%CellCoor(kk, 1)
      jj = level0_iDomain%CellCoor(kk, 2)
      iSc = draSC0(ii, jj)
      ! find drainage path
      iRow = ii
      jCol = jj
      do while (.NOT. iSC > 0)
        ! move downstream
        call moveDownOneCell(fDir0(iRow, jCol), iRow, jCol)
        iSC = draSC0(iRow, jCol)
      end do
      draCell0(ii, jj) = iSC

      ! find cell at L11 corresponding to gauges in Domain at L0 ! >> L11_on_L0 is Id of
      ! the routing cell at level-11
      if (gaugeLoc0(ii, jj) .NE. nodata_i4) then
        ! evaluation gauges
        do ll = 1, domain_mrm(iDomain)%nGauges
          ! search for gaugeID in L0 grid and save ID on L11
          if (domain_mrm(iDomain)%gaugeIdList(ll) .EQ. gaugeLoc0(ii, jj)) then
            domain_mrm(iDomain)%gaugeNodeList(ll) = L0_L11_remap(iDomain)%lowres_id_on_highres(ii, jj)
          end if
        end do
      end if

      if (InflowGaugeLoc0(ii, jj) .NE. nodata_i4) then
        ! inflow gauges
        do ll = 1, domain_mrm(iDomain)%nInflowGauges
          ! search for gaugeID in L0 grid and save ID on L11
          if (domain_mrm(iDomain)%InflowGaugeIdList(ll) .EQ. InflowGaugeLoc0(ii, jj)) &
              domain_mrm(iDomain)%InflowGaugeNodeList(ll) = L0_L11_remap(iDomain)%lowres_id_on_highres(ii, jj)
        end do
      end if
    end do

    !--------------------------------------------------------
    ! Start padding up local variables to global variables
    !--------------------------------------------------------
    ! L0 data sets
    ! check whether L0 data is shared
    ! ToDo: check if change is correct
    if (iDomain .eq. 1) then
      call append(L0_draCell, PACK(draCell0(:, :), level0_iDomain%mask))
    else if (domainMeta%L0DataFrom(iDomain) == iDomain) then
      call append(L0_draCell, PACK(draCell0(:, :), level0_iDomain%mask))
    end if

    ! free space
    deallocate (draSC0, fDir0, gaugeLoc0, draCell0)

  end subroutine  L11_set_drain_outlet_gauges

  ! ------------------------------------------------------------------

  !    NAME
  !        L11_stream_features

  !    PURPOSE
  !>       \brief Stream features (stream network and floodplain)

  !>       \details Stream features (stream network and floodplain)
  !>       If a variable is added or removed here, then it also has to
  !>       be added or removed in the subroutine L11_config_set in
  !>       module mo_restart and in the subroutine set_L11_config in module
  !>       mo_set_netcdf_restart

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iDomain" Domain Id

  !    HISTORY
  !>       \authors Luis Samaniego

  !>       \date Dec 2005

  ! Modifications:
  ! Luis Samaniego Jan 2013 - modular version
  ! R. Kumar       Oct 2013 - stack size increased from nNodes to 100
  ! Stephan Thober Aug 2015 - ported to mRM
  ! Stephan Thober Nov 2016 - only read flood plain area if processMatrix for routing equals 1
  ! Robert Schweppe Jun 2018 - refactoring and reformatting
  ! Stephan Thober Jul 2018 - introduced cut off Length at 40 percentile to neglect short paths in headwaters for adaptive timesteps
  ! Stephan Thober, Pallav Kumar Shrestha, Sep 2020 - bug fix in cut off Length at 40 percentile, neglecting links with -9999. that occur if multiple outlets are present

  subroutine L11_stream_features(iDomain)

    use mo_append, only : append
    use mo_common_constants, only : nodata_dp, nodata_i4
    use mo_common_types, only: Grid
    use mo_common_variables, only : domainMeta, L0_elev, iFlag_cordinate_sys, level0, processMatrix
    use mo_mrm_global_variables, only : L0_fDir, &
        L0_floodPlain, L0_streamNet, L11_aFloodPlain, L11_fCol, L11_fRow, L11_length, &
        L11_nOutlets, L11_netPerm, L11_slope, L11_tCol, L11_tRow, level11
    use mo_percentile, only: percentile

    implicit none

    ! Domain Id
    integer(i4), intent(in) :: iDomain

    integer(i4) :: nLinks

    integer(i4), dimension(:, :), allocatable :: iD0

    integer(i4), dimension(:, :), allocatable :: fDir0

    real(dp), dimension(:, :), allocatable :: elev0

    real(dp), dimension(:, :), allocatable :: cellarea0

    integer(i4), dimension(:, :), allocatable :: streamNet0

    integer(i4), dimension(:, :), allocatable :: floodPlain0

    ! routing order (permutation)
    integer(i4), dimension(:), allocatable :: netPerm

    integer(i4), dimension(:), allocatable :: nLinkFromRow

    integer(i4), dimension(:), allocatable :: nLinkFromCol

    integer(i4), dimension(:), allocatable :: nLinkToRow

    integer(i4), dimension(:), allocatable :: nLinkToCol

    real(dp), dimension(:), allocatable :: nLinkLength

    real(dp), dimension(:), allocatable :: nLinkAFloodPlain

    real(dp), dimension(:), allocatable :: nLinkSlope

    integer(i4) :: ii, rr, ns, s0, e0

    integer(i4) :: frow, fcol

    integer(i4) :: fId, tId

    integer(i4), dimension(:, :), allocatable :: stack, append_chunk

    integer(i4), dimension(:), allocatable :: dummy_1d

    real(dp) :: length

    integer(i4), dimension(:, :), allocatable :: nodata_i4_tmp

    real(dp), dimension(:, :), allocatable :: nodata_dp_tmp

    type(Grid), pointer :: level0_iDomain => null()


    level0_iDomain => level0(domainMeta%L0DataFrom(iDomain))
    s0 = level0_iDomain%iStart
    e0 = level0_iDomain%iEnd
    nLinks = level11(iDomain)%nCells - L11_nOutlets(iDomain)


    ! allocate
    allocate (iD0           (level0_iDomain%nrows, level0_iDomain%ncols))
    allocate (elev0         (level0_iDomain%nrows, level0_iDomain%ncols))
    allocate (fDir0         (level0_iDomain%nrows, level0_iDomain%ncols))
    allocate (cellarea0     (level0_iDomain%nrows, level0_iDomain%ncols))
    allocate (streamNet0    (level0_iDomain%nrows, level0_iDomain%ncols))
    allocate (floodPlain0   (level0_iDomain%nrows, level0_iDomain%ncols))

    !  Routing network vectors have nNodes size instead of nLinks to
    !  avoid the need of having two extra indices to identify a Domain.
    allocate (stack             (level11(iDomain)%nCells, 2)) ! >> stack(nNodes, 2)
    allocate (dummy_1d          (2))
    allocate (append_chunk      (8, 2))
    allocate (netPerm           (level11(iDomain)%nCells))
    allocate (nLinkFromRow      (level11(iDomain)%nCells))
    allocate (nLinkFromCol      (level11(iDomain)%nCells))
    allocate (nLinkToRow        (level11(iDomain)%nCells))
    allocate (nLinkToCol        (level11(iDomain)%nCells))
    allocate (nLinkLength       (level11(iDomain)%nCells))
    allocate (nLinkAFloodPlain  (level11(iDomain)%nCells))
    allocate (nLinkSlope        (level11(iDomain)%nCells))

    allocate (nodata_i4_tmp      (level0_iDomain%nrows, level0_iDomain%ncols))
    allocate (nodata_dp_tmp      (level0_iDomain%nrows, level0_iDomain%ncols))

    ! initialize
    iD0(:, :) = nodata_i4
    elev0(:, :) = nodata_dp
    fDir0(:, :) = nodata_i4
    cellarea0(:, :) = nodata_dp
    streamNet0(:, :) = nodata_i4
    floodPlain0(:, :) = nodata_i4

    stack(:, :) = nodata_i4
    append_chunk(:, :) = nodata_i4
    netPerm(:) = nodata_i4
    nLinkFromRow(:) = nodata_i4
    nLinkFromCol(:) = nodata_i4
    nLinkToRow(:) = nodata_i4
    nLinkToCol(:) = nodata_i4
    nLinkLength(:) = nodata_dp
    nLinkAFloodPlain(:) = nodata_dp
    nLinkSlope(:) = nodata_dp

    nodata_i4_tmp(:, :) = nodata_i4
    nodata_dp_tmp(:, :) = nodata_dp

    ! for a single node model run
    if(level11(iDomain)%nCells .GT. 1) then
      ! get L0 fields
      iD0(:, :) = UNPACK(level0_iDomain%Id, level0_iDomain%mask, nodata_i4_tmp)
      elev0(:, :) = UNPACK(L0_elev (s0 : e0), &
          level0_iDomain%mask, nodata_dp_tmp)
      fDir0(:, :) = UNPACK(L0_fDir (s0 : e0), &
          level0_iDomain%mask, nodata_i4_tmp)
      cellarea0(:, :) = UNPACK(level0_iDomain%CellArea, level0_iDomain%mask, nodata_dp_tmp)

      ! get network vectors of L11
      netPerm(:) = L11_netPerm (level11(iDomain)%iStart : level11(iDomain)%iEnd)
      nLinkFromRow(:) = L11_fRow    (level11(iDomain)%iStart : level11(iDomain)%iEnd)
      nLinkFromCol(:) = L11_fCol    (level11(iDomain)%iStart : level11(iDomain)%iEnd)
      nLinkToRow(:) = L11_tRow    (level11(iDomain)%iStart : level11(iDomain)%iEnd)
      nLinkToCol(:) = L11_tCol    (level11(iDomain)%iStart : level11(iDomain)%iEnd)

      ! Flood plains:  stream network delineation
      streamNet0(:, :) = nodata_i4
      floodPlain0(:, :) = nodata_i4

      do rr = 1, nLinks

        ii = netPerm(rr)
        frow = nLinkFromRow(ii)
        fcol = nLinkFromCol(ii)

        ! Init
        streamNet0(frow, fcol) = ii
        floodPlain0(frow, fcol) = ii
        stack = 0
        append_chunk = 0
        ns = 1
        stack(ns, 1) = frow
        stack(ns, 2) = fcol

        call cellLength(iDomain, fDir0(frow, fcol), fRow, fCol, iFlag_cordinate_sys, nLinkLength(ii))
        nLinkSlope(ii) = elev0(frow, fcol)

        fId = iD0(frow, fcol)
        tId = iD0(nLinkToRow(ii), nLinkToCol(ii))

        do while (.NOT. (fId == tId))
          ! Search flood plain from point(frow,fcol) upwards, keep co-ordinates in STACK
          do while (ns > 0)
            if (ns + 8 .gt. size(stack, 1)) then
              call append(stack, append_chunk)
            end if
            call moveUp(elev0, fDir0, frow, fcol, stack, ns)
            stack(1, 1) = 0
            stack(1, 2) = 0
            ! stack = cshift(stack, SHIFT = 1, DIM = 1)
            ! substitute cshift <<<
            dummy_1d = stack(1, :)
            stack(: size(stack, dim = 1) - 1, :) = stack(2 :, :)
            stack(size(stack, dim = 1), :) = dummy_1d
            ! substitute cshift >>>
            if (stack(1, 1) > 0 .and. stack(1, 2) > 0) floodPlain0(stack(1, 1), stack(1, 2)) = ii
            ns = count(stack > 0) / 2
          end do

          ! move downstream
          call moveDownOneCell(fDir0(frow, fcol), frow, fcol)
          streamNet0(frow, fcol) = ii
          floodPlain0(frow, fcol) = ii
          fId = iD0(frow, fcol)
          stack = 0
          ns = 1
          stack(ns, 1) = frow
          stack(ns, 2) = fcol
          call cellLength(iDomain, fDir0(fRow, fCol), fRow, fCol, iFlag_cordinate_sys, length)
          nLinkLength(ii) = nLinkLength(ii) + length

        end do

        ! stream bed slope
        nLinkSlope(ii) = (nLinkSlope(ii) - elev0(frow, fcol)) / nLinkLength(ii)

        if (nLinkSlope(ii) < 0.0001_dp) nLinkSlope(ii) = 0.0001_dp

        ! calculate area of floodplains (avoid overwriting)
        nLinkAFloodPlain(ii) = sum (cellarea0(:, :), mask = (floodPlain0(:, :) == ii))
        !  old > real( count( floodPlain0(:,:,) == i), dp ) * cellarea0

      end do

      ! end of multi-node network design loop
    end if

    ! cut off Length at 40 percentile to neglect short paths in headwaters
    if ((processMatrix(8, 1) .eq. 2) .or. (processMatrix(8, 1) .eq. 3)) then
      if (count(nLinkLength(:) .ge. 0._dp) .gt. 2) then
        length = percentile(pack(nLinkLength(:), nLinkLength(:) .ge. 0._dp), 40._dp)
        nLinkLength(:) = merge(nLinkLength(:), length, (nLinkLength(:) .gt. length))
      end if
    end if

    !--------------------------------------------------------
    ! Start padding up local variables to global variables
    !--------------------------------------------------------

    ! L0 data sets
    ! check whether L0 data is shared
    ! ToDo: check if change is correct
    if (iDomain .eq. 1) then
      call append(L0_streamNet, PACK (streamNet0(:, :), level0_iDomain%mask))
      call append(L0_floodPlain, PACK (floodPlain0(:, :), level0_iDomain%mask))
    else if (domainMeta%L0DataFrom(iDomain) == iDomain) then
      call append(L0_streamNet, PACK (streamNet0(:, :), level0_iDomain%mask))
      call append(L0_floodPlain, PACK (floodPlain0(:, :), level0_iDomain%mask))
    end if


    ! L11 network data sets
    call append(L11_length, nLinkLength(:))
    call append(L11_aFloodPlain, nLinkAFloodPlain(:))
    call append(L11_slope, nLinkSlope(:))

    ! free space
    deallocate (&
        iD0, elev0, fDir0, streamNet0, floodPlain0, &
        cellarea0, stack, netPerm, nLinkFromRow, nLinkFromCol, nLinkToRow, nLinkToCol, &
        nLinkLength, nLinkAFloodPlain, nLinkSlope, dummy_1d)
    deallocate(nodata_i4_tmp, nodata_dp_tmp)

  end subroutine L11_stream_features

  ! ------------------------------------------------------------------

  !    NAME
  !        L11_fraction_sealed_floodplain

  !    PURPOSE
  !>       \brief Fraction of the flood plain with impervious cover

  !>       \details Fraction of the flood plain with impervious cover (\ref fig_routing "Routing
  !>       Network"). This proportion is used to regionalize the Muskingum parameters.
  !>       Samaniego et al. \cite SB05 found out that this fraction is one of the statistically
  !>       significant predictor variables of peak discharge in mesoscale Domains.
  !>       If a variable is added or removed here, then it also has to
  !>       be added or removed in the subroutine L11_config_set in
  !>       module mo_restart and in the subroutine set_L11_config in module
  !>       mo_set_netcdf_restart

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: LCClassImp" Impervious land cover class Id, e.g. = 2 (old code)
  !>       \param[in] "logical :: do_init"

  !    HISTORY
  !>       \authors Luis Samaniego

  !>       \date Dec 2005

  ! Modifications:
  ! Luis Samaniego Jan 2013 - modular version
  ! Stephan Thober Aug 2015 - ported to mRM
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine L11_fraction_sealed_floodplain(LCClassImp, do_init)

    use mo_append, only : append
    use mo_common_constants, only : nodata_dp
    use mo_common_types, only: Grid
    use mo_common_variables, only : domainMeta, L0_LCover, level0, domainMeta, nLCoverScene
    use mo_mrm_global_variables, only : L0_floodPlain, L11_aFloodPlain, &
                                        L11_nLinkFracFPimp, L11_nOutlets, level11

    implicit none

    ! Impervious land cover class Id, e.g. = 2 (old code)
    integer(i4), intent(in) :: LCClassImp

    logical, intent(in) :: do_init

    integer(i4) :: nLinks

    real(dp), dimension(:), pointer :: nLinkAFloodPlain => null()

    real(dp), dimension(:,:), allocatable :: temp_array

    integer(i4) :: ii, iDomain, iiLC, s0, e0

    type(Grid), pointer :: level0_iDomain => null()


    ! initialization
    do iDomain = 1, domainMeta%nDomains
      allocate(temp_array(level11(iDomain)%nCells, nLCoverScene))
      temp_array = nodata_dp
      if (do_init) then
        level0_iDomain => level0(domainMeta%L0DataFrom(iDomain))

        s0 = level0_iDomain%iStart
        e0 = level0_iDomain%iEnd
        nLinks = level11(iDomain)%nCells + 1 - L11_nOutlets(iDomain)
        nLinkAFloodPlain => L11_aFloodPlain(level11(iDomain)%iStart : level11(iDomain)%iEnd)

        do iiLC = 1, nLCoverScene
          ! for a single node model run
          if(nLinks .GT. 0) then
            do ii = 1, nLinks
              temp_array(ii, iiLC) = sum(level0_iDomain%CellArea(:), &
                  mask = (L0_floodPlain(s0 : e0) == ii .and. L0_LCover(s0 : e0, iiLC) == LCClassImp)) &
                  / nLinkAFloodPlain(ii)
            end do
          end if
        end do
      end if
      call append(L11_nLinkFracFPimp, temp_array(:,:))
      deallocate(temp_array)
    end do

  end subroutine L11_fraction_sealed_floodplain

  ! ------------------------------------------------------------------
  !  MOVE UPSTREAM FROM-TO
  ! ------------------------------------------------------------------
  !    NAME
  !        moveUp

  !    PURPOSE
  !>       \brief TODO: add description

  !>       \details TODO: add description

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:, :) :: elev0"
  !>       \param[in] "integer(i4), dimension(:, :) :: fDir0"
  !>       \param[in] "integer(i4) :: fi, fj"                 co-ordinate of the stream bed
  !>       \param[in] "integer(i4) :: fi, fj"                 co-ordinate of the stream bed

  !    INTENT(INOUT)
  !>       \param[inout] "integer(i4), dimension(:, :) :: ss"
  !>       \param[inout] "integer(i4) :: nn"

  !    HISTORY
  !>       \authors Robert Schweppe

  !>       \date Jun 2018

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine moveUp(elev0, fDir0, fi, fj, ss, nn)

    use mo_mrm_constants, only : deltaH
    use mo_utils, only : ge, le

    implicit none

    real(dp), dimension(:, :), allocatable, intent(IN) :: elev0

    integer(i4), dimension(:, :), allocatable, intent(IN) :: fDir0

    ! co-ordinate of the stream bed
    integer(i4), intent(IN) :: fi, fj

    integer(i4), dimension(:, :), intent(INOUT) :: ss

    integer(i4), intent(INOUT) :: nn

    integer(i4) :: ii, jj, ip, im, jp, jm

    integer(i4) :: nrows, ncols


    ii = ss(1, 1)
    jj = ss(1, 2)
    ip = ii + 1
    im = ii - 1
    jp = jj + 1
    jm = jj - 1

    nrows = size(fDir0, 1)
    ncols = size(fDir0, 2)

    !E
    if   (jp                <= ncols) then
      if ((fdir0(ii, jp) == 16)                   .and. &
          (le((elev0(ii, jp) - elev0(fi, fj)), deltaH))  .and. &
          (ge((elev0(ii, jp) - elev0(fi, fj)), 0.0_dp))        &
          ) then
        nn = nn + 1
        ss(nn, 1) = ii
        ss(nn, 2) = jp
        !print *, i,jp
      end if
    end if

    !SE
    if ((ip                <= nrows) .and. &
        (jp                <= ncols)) then
      if ((fdir0(ip, jp) == 32)                   .and. &
          (le((elev0(ip, jp) - elev0(fi, fj)), deltaH))  .and. &
          (ge((elev0(ii, jp) - elev0(fi, fj)), 0.0_dp))        &
          ) then
        nn = nn + 1
        ss(nn, 1) = ip
        ss(nn, 2) = jp
        !print *, ip,jp
      end if
    end if

    !S
    if ((ip               <= nrows)  .and. &
        (jp               <= ncols)) then
      if ((fdir0(ip, jj) == 64)                 .and. &
          (le((elev0(ip, jj) - elev0(fi, fj)), deltaH))  .and. &
          (ge((elev0(ii, jp) - elev0(fi, fj)), 0.0_dp))        &
          ) then
        nn = nn + 1
        ss(nn, 1) = ip
        ss(nn, 2) = jj
        !print *, ip,j
      end if
    end if

    !SW
    if ((ip                <= nrows) .and. &
        (jp                <= ncols) .and. &
        (jm                >= 1)) then
      if ((fdir0(ip, jm) == 128)                 .and. &
          (le((elev0(ip, jm) - elev0(fi, fj)), deltaH))  .and. &
          (ge((elev0(ii, jp) - elev0(fi, fj)), 0.0_dp))        &
          ) then
        nn = nn + 1
        ss(nn, 1) = ip
        ss(nn, 2) = jm
        !print *, ip,jm
      end if
    end if

    !W
    if ((jm                 >= 1) .and. &
        (jp                 <= ncols)) then
      if ((fdir0(ii, jm)  == 1)                 .and. &
          (le((elev0(ii, jm) - elev0(fi, fj)), deltaH))  .and. &
          (ge((elev0(ii, jp) - elev0(fi, fj)), 0.0_dp))        &
          ) then
        nn = nn + 1
        ss(nn, 1) = ii
        ss(nn, 2) = jm
        !print *, i,jm
      end if
    end if

    !NW
    if ((im                >= 1) .and. &
        (jp                <= ncols) .and. &
        (jm                >= 1))  then
      if ((fdir0(im, jm) == 2)                 .and. &
          (le((elev0(im, jm) - elev0(fi, fj)), deltaH))  .and. &
          (ge((elev0(ii, jp) - elev0(fi, fj)), 0.0_dp))        &
          ) then
        nn = nn + 1
        ss(nn, 1) = im
        ss(nn, 2) = jm
        !print *, im,jm
      end if
    end if

    !N
    if ((im                >= 1) .and. &
        (jp                 <= ncols)) then
      if ((fdir0(im, jj)  == 4)                 .and. &
          (le((elev0(im, jj) - elev0(fi, fj)), deltaH))  .and. &
          (ge((elev0(ii, jp) - elev0(fi, fj)), 0.0_dp))        &
          ) then
        nn = nn + 1
        ss(nn, 1) = im
        ss(nn, 2) = jj
        !print *, im,j
      end if
    end if

    !NE
    if ((im                >= 1) .and. &
        (jp                <= ncols))  then
      if ((fdir0(im, jp) == 8)               .and. &
          (le((elev0(im, jp) - elev0(fi, fj)), deltaH))  .and. &
          (ge((elev0(ii, jp) - elev0(fi, fj)), 0.0_dp))        &
          ) then
        nn = nn + 1
        ss(nn, 1) = im
        ss(nn, 2) = jp
        !print *, im,jp
      end if
    end if

  end subroutine moveUp

  ! ------------------------------------------------------------------
  !  MOVE DOWNSTREAM
  ! ------------------------------------------------------------------
  !    NAME
  !        moveDownOneCell

  !    PURPOSE
  !>       \brief TODO: add description

  !>       \details TODO: add description

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: fDir"

  !    INTENT(INOUT)
  !>       \param[inout] "integer(i4) :: iRow, jCol"
  !>       \param[inout] "integer(i4) :: iRow, jCol"

  !    HISTORY
  !>       \authors Robert Schweppe

  !>       \date Jun 2018

  ! Modifications:

  subroutine moveDownOneCell(fDir, iRow, jCol)
    implicit none

    integer(i4), intent(IN) :: fDir

    integer(i4), intent(INOUT) :: iRow, jCol


    select case (fDir)
    case(1)   !E
      jCol = jCol + 1
    case(2)   !SE
      iRow = iRow + 1
      jCol = jCol + 1
    case(4)   !S
      iRow = iRow + 1
    case(8)   !SW
      iRow = iRow + 1
      jCol = jCol - 1
    case(16)  !W
      jCol = jCol - 1
    case(32)  !NW
      iRow = iRow - 1
      jCol = jCol - 1
    case(64)  !N
      iRow = iRow - 1
    case(128) !NE
      iRow = iRow - 1
      jCol = jCol + 1
    case default !sink
      ! do nothing
    end select

  end subroutine moveDownOneCell

  ! ------------------------------------------------------------------
  !  CELL LENGTH
  ! ------------------------------------------------------------------
  !    NAME
  !        cellLength

  !    PURPOSE
  !>       \brief TODO: add description

  !>       \details TODO: add description

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iDomain"
  !>       \param[in] "integer(i4) :: fDir"
  !>       \param[in] "integer(i4) :: iRow"
  !>       \param[in] "integer(i4) :: jCol"
  !>       \param[in] "integer(i4) :: iCoorSystem"

  !    INTENT(OUT)
  !>       \param[out] "real(dp) :: length"

  !    HISTORY
  !>       \authors Robert Schweppe

  !>       \date Jun 2018

  ! Modifications:

  subroutine cellLength(iDomain, fDir, iRow, jCol, iCoorSystem, length)

    use mo_common_types, only: Grid
    use mo_common_variables, only : domainMeta, level0
    use mo_constants, only : SQRT2_dp

    implicit none

    integer(i4), intent(IN) :: iDomain

    integer(i4), intent(IN) :: fDir

    integer(i4), intent(IN) :: iRow

    integer(i4), intent(IN) :: jCol

    integer(i4), intent(IN) :: iCoorSystem

    real(dp), intent(OUT) :: length

    integer(i4) :: iRow_to, jCol_to

    real(dp) :: lat_1, long_1, lat_2, long_2

    type(Grid), pointer :: level0_iDomain => null()


    level0_iDomain => level0(domainMeta%L0DataFrom(iDomain))

    ! regular X-Y cordinate system
    IF(iCoorSystem .EQ. 0) THEN

      select case (fDir)
      case(1, 4, 16, 64)       ! E, S, W, N
        length = 1.0_dp
      case(2, 8, 32, 128)      ! SE, SW, NW, NE
        length = SQRT2_dp
      end select
      length = length * level0_iDomain%cellsize

      ! regular lat-lon cordinate system
    ELSE IF(iCoorSystem .EQ. 1) THEN
      iRow_to = iRow
      jCol_to = jCol

      ! move in the direction of flow
      call moveDownOneCell(fDir, iRow_to, jCol_to)

      ! estimate lat-lon points
      lat_1 = level0_iDomain%yllcorner + real((level0_iDomain%ncols - jCol), dp) * level0_iDomain%cellsize + &
              0.5_dp * level0_iDomain%cellsize
      long_1 = level0_iDomain%xllcorner + real((iRow - 1), dp) * level0_iDomain%cellsize + &
              0.5_dp * level0_iDomain%cellsize

      lat_2 = level0_iDomain%yllcorner + real((level0_iDomain%ncols - jCol_to), dp) * level0_iDomain%cellsize + &
              0.5_dp * level0_iDomain%cellsize
      long_2 = level0_iDomain%xllcorner + real((iRow_to - 1), dp) * level0_iDomain%cellsize + &
              0.5_dp * level0_iDomain%cellsize
      ! get distance between two points
      call get_distance_two_lat_lon_points(lat_1, long_1, lat_2, long_2, length)

    END IF
    !
  end subroutine cellLength


  ! --------------------------------------------------------------------------

  !    NAME
  !        get_distance_two_lat_lon_points

  !    PURPOSE
  !>       \brief estimate distance in [m] between two points in a lat-lon

  !>       \details estimate distance in [m] between two points in a lat-lon
  !>       Code is based on one that is implemented in the VIC-3L model

  !    INTENT(IN)
  !>       \param[in] "real(dp) :: lat1, long1, lat2, long2" latitude  of point-1
  !>       \param[in] "real(dp) :: lat1, long1, lat2, long2" longitude of point-1
  !>       \param[in] "real(dp) :: lat1, long1, lat2, long2" latitude  of point-2
  !>       \param[in] "real(dp) :: lat1, long1, lat2, long2" longitude of point-2

  !    INTENT(OUT)
  !>       \param[out] "real(dp) :: distance_out" distance between two points [m]

  !    HISTORY
  !>       \authors Rohini Kumar

  !>       \date May 2014

  ! Modifications:
  ! Stephan Thober Aug 2015 - ported to mRM
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine get_distance_two_lat_lon_points(lat1, long1, lat2, long2, distance_out)

    use mo_constants, only : RadiusEarth_dp, TWOPI_dp

    implicit none

    ! longitude of point-2
    real(dp), intent(in) :: lat1, long1, lat2, long2

    ! distance between two points [m]
    real(dp), intent(out) :: distance_out

    real(dp) :: theta1

    real(dp) :: phi1

    real(dp) :: theta2

    real(dp) :: phi2

    real(dp) :: dtor

    real(dp) :: term1

    real(dp) :: term2

    real(dp) :: term3

    real(dp) :: temp


    dtor = TWOPI_dp / 360.0_dp
    theta1 = dtor * long1
    phi1 = dtor * lat1
    theta2 = dtor * long2
    phi2 = dtor * lat2

    term1 = cos(phi1) * cos(theta1) * cos(phi2) * cos(theta2)
    term2 = cos(phi1) * sin(theta1) * cos(phi2) * sin(theta2)
    term3 = sin(phi1) * sin(phi2)
    temp = term1 + term2 + term3
    if(temp .GT. 1.0_dp) temp = 1.0_dp

    distance_out = RadiusEarth_dp * acos(temp);

  end subroutine get_distance_two_lat_lon_points

  ! --------------------------------------------------------------------------

  !     NAME
  !         L11_flow_accumulation

  !     PURPOSE

  !>       \brief Calculates L11 flow accumulation per grid cell
  !>       \details Calculates L11 flow accumulation per grid cell using L11_fDir
  !>                and L11_cellarea. L11_flow_accumulation contains the recursiv subroutine
  !>                calculate_L11_flow_accumulation

  !     INTENT(IN)
  !>        iDomain

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         L11_fAcc, L11_LinkIn_fAcc

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Matthias Kelbling
  !>        \date   Aug 2017
  !     Modified
  !         Stephan Thober, Jun 2018 - refactored to fit MPR_extract

  ! --------------------------------------------------------------------------
  subroutine L11_flow_accumulation(iDomain)

    use mo_mrm_global_variables, only: &
        L11_fDir,          & !  IN: flow direction at L11 (standard notation)
        L11_fAcc             ! OUT: flow accumulation at L11 [km^2]
    use mo_mrm_global_variables, only : level11
    use mo_common_constants, only : nodata_i4, nodata_dp
    use mo_append, only : append

    implicit none

    ! local variables
    integer(i4), intent(in)                   :: iDomain
    real(dp),    dimension(:,:), allocatable  :: fAcc11             ! L11_fAcc array
    integer(i4)                               :: ii, jj             ! row and col index
    integer(i4)                               :: s11, e11           ! Vec. position iDomain - fAcc
    integer(i4)                               :: nrows11, ncols11   ! array size Domain
    integer(i4), dimension(:,:), allocatable  :: fDir11             ! L11_fDir array
    logical,     dimension(:,:), allocatable  :: mask11             ! Domain mask

    ! initialize Domain info
    nrows11 = level11(iDomain)%nrows
    ncols11 = level11(iDomain)%ncols
    s11 = level11(iDomain)%iStart
    e11 = level11(iDomain)%iEnd
    mask11 = level11(iDomain)%mask

    ! allocate arrays
    allocate(fDir11      (nrows11, ncols11))
    allocate(fAcc11      (nrows11, ncols11))

    ! initialize
    fDir11(:,:)     = nodata_i4

    ! get data
    fDir11(:,:)  = UNPACK( L11_fDir(s11:e11),  mask11, nodata_i4 )

    ! initialize fAcc11 with cell area
    fAcc11 = UNPACK( level11(iDomain)%cellarea * 1.e-6,  mask11, nodata_dp )

    ! For each sink call "calculate_L11_flow_accumulation"
    do jj=1, ncols11
      do ii=1, nrows11
        if (fDir11(ii,jj) .eq. 0) then
          call calculate_L11_flow_accumulation(fDir = fDir11, &
              fAcc = fAcc11, &
              ii = ii, &
              jj = jj, &
              nrow = nrows11, &
              ncol = ncols11)
        end if
      end do
    end do

    ! Append
    call append( L11_fAcc, pack(fAcc11(:,:),mask11))

    ! free space
    deallocate(fDir11, fAcc11, mask11)

  contains

    recursive subroutine calculate_L11_flow_accumulation(fDir, fAcc, ii, jj, nrow, ncol)

      implicit none

      integer(i4), intent(in)            :: fDir(:,:)      ! flow Direction
      real(dp), intent (inout)           :: fAcc(:,:)      ! flow accumulation
      integer(i4), intent(in)            :: ii, jj         ! row and col index
      integer(i4), intent(in)            :: nrow, ncol     ! number of rows,cols in array

      ! Scan order:
      !
      !    6 7 8
      !    5 x 1
      !    4 3 2
      !
      ! Each:
      ! 1. Check if cell is inflow cell to current grid
      ! 2. If yes: Call calculate_subroutine and add result

      if (jj+1 .le. ncol) then
        if (fDir(ii,jj+1) .eq. 16_i4) then
          call calculate_L11_flow_accumulation(fDir, fAcc, ii, jj+1, nrow, ncol)
          fAcc(ii,jj) = fAcc(ii,jj) + fAcc(ii,jj+1)
        end if
      end if

      if ((ii+1 .le. nrow) .and. (jj+1 .le. ncol)) then
        if (fDir(ii+1,jj+1) .eq. 32_i4) then
          call calculate_L11_flow_accumulation(fDir, fAcc, ii+1, jj+1, nrow, ncol)
          fAcc(ii,jj) = fAcc(ii,jj) + fAcc(ii+1,jj+1)
        end if
      end if

      if (ii+1 .le. nrow) then
        if (fDir(ii+1,jj) .eq. 64_i4) then
          call calculate_L11_flow_accumulation(fDir, fAcc, ii+1, jj, nrow, ncol)
          fAcc(ii,jj) = fAcc(ii,jj) + fAcc(ii+1,jj)
        end if
      end if

      if ((ii+1 .le. nrow) .and. (jj-1 .ge. 1)) then
        if (fDir(ii+1,jj-1) .eq. 128_i4) then
          call calculate_L11_flow_accumulation(fDir, fAcc, ii+1, jj-1, nrow, ncol)
          fAcc(ii,jj) = fAcc(ii,jj) + fAcc(ii+1,jj-1)
        end if
      end if

      if (jj-1 .ge. 1) then
        if (fDir(ii,jj-1) .eq. 1_i4) then
          call calculate_L11_flow_accumulation(fDir, fAcc, ii, jj-1, nrow, ncol)
          fAcc(ii,jj) = fAcc(ii,jj) + fAcc(ii,jj-1)
        end if
      end if

      if ((ii-1 .ge. 1) .and. (jj-1 .ge. 1)) then
        if (fDir(ii-1,jj-1) .eq. 2_i4) then
          call calculate_L11_flow_accumulation(fDir, fAcc, ii-1, jj-1, nrow, ncol)
          fAcc(ii,jj) = fAcc(ii,jj) + fAcc(ii-1,jj-1)
        end if
      end if

      if (ii-1 .ge. 1) then
        if (fDir(ii-1,jj) .eq. 4_i4) then
          call calculate_L11_flow_accumulation(fDir, fAcc, ii-1, jj, nrow, ncol)
          fAcc(ii,jj) = fAcc(ii,jj) + fAcc(ii-1,jj)
        end if
      end if

      if ((ii-1 .ge. 1) .and. (jj+1 .le. ncol)) then
        if (fDir11(ii-1,jj+1) .eq. 8_i4) then
          call calculate_L11_flow_accumulation(fDir, fAcc, ii-1, jj+1, nrow, ncol)
          fAcc(ii,jj) = fAcc(ii,jj) + fAcc(ii-1,jj+1)
        end if
      end if

      ! print *, ii, jj, nrow, ncol, fAcc(ii, jj)

    end subroutine calculate_L11_flow_accumulation

  end subroutine L11_flow_accumulation


  ! ------------------------------------------------------------------

  !     NAME
  !         L11_calc_celerity

  !     PURPOSE
  !>        \brief L11 celerity based on L0 elevation and L0 fAcc

  !>        \details L11 celerity based on L0 elevation and L0 fAcc

  !     INTENT(IN)
  !>        \param[in]

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
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author  Matthias Kelbling
  !>        \date    Oct 2017

  ! ------------------------------------------------------------------

  subroutine L11_calc_celerity(iDomain, param)
    use mo_common_constants, only: nodata_i4, nodata_dp
    use mo_mad, only: mad
    use mo_append, only: append
    use mo_mpr_global_variables, only: &
        L0_slope               ! IN:    slope [%]
    use mo_common_types, only: Grid
    use mo_common_variables, only: &
        domainMeta,          & ! IN:    for L0 Domain indexer
        level0                 ! IN:    level 0 grid
    use mo_mrm_global_variables, only: &
        L0_fDir,             & ! IN:    flow direction (standard notation) L0
        L0_fAcc,             & ! IN:    flow accumulation (number of cells)?
        L0_streamNet,        & ! IN:    stream Network at Level 0
        level11,             & ! IN:    level 11 grid
        L11_fRow,            & ! IN:    from row in L0 grid
        L11_fCol,            & ! IN:    from col in L0 grid
        L11_tRow,            & ! IN:    to row in L0 grid
        L11_tCol,            & ! IN:    to col in L0 grid
        L11_netPerm,         & ! IN:    routing order (permutation)
        L11_nOutlets,        & ! IN:    Number of Outlets/Sinks
        L11_celerity,        & ! INOUT: averaged celerity
        L0_celerity            ! INOUT

    implicit none

    integer(i4), intent(in)                  :: iDomain         ! Domain
    real(dp), dimension(:), intent(in)       :: param

    ! local
    integer(i4)                              :: nCells0
    integer(i4)                              :: nNodes
    integer(i4)                              :: nLinks
    integer(i4)                              :: nrows0, ncols0
    integer(i4)                              :: iStart0, iEnd0
    integer(i4)                              :: nrows11, ncols11
    integer(i4)                              :: iStart11, iEnd11
    logical,     dimension(:,:), allocatable :: mask0
    integer(i4), dimension(:,:), allocatable :: iD0
    integer(i4), dimension(:,:), allocatable :: fDir0
    integer(i4), dimension(:,:), allocatable :: fAcc0
    real(dp),    dimension(:,:), allocatable :: slope0
    real(dp),    dimension(:), allocatable :: slope_tmp
    real(dp),    dimension(:,:), allocatable :: cellarea0
    integer(i4), dimension(:),   allocatable :: netPerm         ! routing order (permutation)
    integer(i4), dimension(:),   allocatable :: nLinkFromRow
    integer(i4), dimension(:),   allocatable :: nLinkFromCol
    integer(i4), dimension(:),   allocatable :: nLinkToRow
    integer(i4), dimension(:),   allocatable :: nLinkToCol
    integer(i4)                              :: ii, rr, ns
    integer(i4)                              :: frow, fcol
    integer(i4)                              :: fId,  tId
    real(dp),    dimension(:),   allocatable :: stack, append_chunk ! Stacks celerity along the L0 river-path
    integer(i4), dimension(:),   allocatable :: dummy_1d

    real(dp)                                 :: L0_link_slope
    real(dp),    dimension(:),   allocatable :: celerity11
    real(dp),    dimension(:,:), allocatable :: celerity0

    integer(i4), dimension(:,:), allocatable :: nodata_i4_tmp
    real(dp),    dimension(:,:), allocatable :: nodata_dp_tmp
    logical,     dimension(:),   allocatable :: slopemask0

    type(Grid), pointer :: level0_iDomain

    ! level-0 information
    level0_iDomain => level0(domainMeta%L0DataFrom(iDomain))
    nrows0 = level0_iDomain%nrows
    ncols0 = level0_iDomain%ncols
    nCells0 = level0_iDomain%ncells
    iStart0 = level0_iDomain%iStart
    iEnd0 = level0_iDomain%iEnd
    mask0 = level0_iDomain%mask

    ! level-11 information
    iStart11 = level11(iDomain)%iStart
    iEnd11 = level11(iDomain)%iEnd
    nrows11 = level11(iDomain)%nrows
    ncols11 = level11(iDomain)%ncols
    nNodes = level11(iDomain)%ncells

    nLinks  = nNodes - L11_nOutlets(iDomain)

    ! allocate
    allocate ( iD0         ( nrows0, ncols0 ) )
    allocate ( slope0      ( nrows0, ncols0 ) )
    allocate ( fDir0       ( nrows0, ncols0 ) )
    allocate ( fAcc0       ( nrows0, ncols0 ) )
    allocate ( cellarea0   ( nrows0, ncols0 ) )
    allocate ( celerity0   ( nrows0, ncols0 ) )
    allocate ( slopemask0  ( ncells0 ) )

    !  Routing network vectors have nNodes size instead of nLinks to
    !  avoid the need of having two extra indices to identify a Domain.
    allocate ( stack       ( 1 ) )
    allocate ( append_chunk( 1 ) )
    allocate ( dummy_1d    ( 2 ))
    allocate ( netPerm     ( nNodes ) )
    allocate ( nLinkFromRow( nNodes ) )
    allocate ( nLinkFromCol( nNodes ) )
    allocate ( nLinkToRow  ( nNodes ) )
    allocate ( nLinkToCol  ( nNodes ) )
    allocate ( celerity11  ( nNodes ) )
    allocate ( slope_tmp   ( nNodes ) )

    allocate (nodata_i4_tmp( nrows0, ncols0 ) )
    allocate (nodata_dp_tmp( nrows0, ncols0 ) )

    ! initialize
    iD0(:,:)           = nodata_i4
    fDir0(:,:)         = nodata_i4
    fAcc0(:,:)         = nodata_i4
    cellarea0(:,:)     = nodata_dp
    slope0(:,:)        = nodata_dp

    stack(:)           = nodata_dp
    append_chunk(:)    = nodata_dp
    netPerm(:)         = nodata_i4
    nLinkFromRow(:)    = nodata_i4
    nLinkFromCol(:)    = nodata_i4
    nLinkToRow(:)      = nodata_i4
    nLinkToCol(:)      = nodata_i4
    celerity11(:)      = nodata_dp
    celerity0(:,:)     = nodata_dp
    slopemask0(:)      = .False.

    nodata_i4_tmp(:,:) = nodata_i4
    nodata_dp_tmp(:,:) = nodata_dp

    ! for a single node model run

    if(nNodes .GT. 1) then
      ! get L0 fields
      iD0(:,:) = UNPACK(level0_iDomain%Id(1:nCells0), mask0, nodata_i4_tmp)
      fDir0(:,:) = UNPACK(L0_fDir(iStart0:iEnd0), mask0, nodata_i4_tmp)
      fAcc0(:,:) = UNPACK(L0_fAcc(iStart0:iEnd0), mask0, nodata_i4_tmp)
      cellarea0(:,:) = UNPACK(level0_iDomain%cellarea(1:nCells0), mask0, nodata_dp_tmp)

      ! smoothing river slope
      slope_tmp = L0_slope(iStart0:iEnd0)
      where ( slope_tmp .lt. 0.1_dp ) slope_tmp = 0.1_dp

      slopemask0(:) = (L0_streamNet(iStart0:iEnd0) .ne. nodata_i4)

      ! smooth river cells if there is more than one cell
      if( count(slopemask0) .GT. 1) then
        slope_tmp = mad(arr = slope_tmp, z = 2.25_dp, mask = slopemask0, tout='u', mval=0.1_dp)
      end if
      slope0(:,:) = UNPACK(slope_tmp,  mask0, nodata_dp_tmp )

      ! get network vectors of L11
      netPerm(:)      = L11_netPerm ( iStart11 : iEnd11 )
      nLinkFromRow(:) = L11_fRow    ( iStart11 : iEnd11 )
      nLinkFromCol(:) = L11_fCol    ( iStart11 : iEnd11 )
      nLinkToRow(:)   = L11_tRow    ( iStart11 : iEnd11 )
      nLinkToCol(:)   = L11_tCol    ( iStart11 : iEnd11 )

      do rr = 1, nLinks

        ii   = netPerm(rr)
        frow = nLinkFromRow(ii)
        fcol = nLinkFromCol(ii)

        ! Init
        stack(:) = 0_dp
        ns = 1

        fId = iD0( frow, fcol )
        tId = iD0( nLinkToRow(ii) , nLinkToCol(ii) )
        do
          L0_link_slope = slope0(frow, fcol) / 100._dp

          ! celerity parametrization
          stack(ns) = param(1) * sqrt(L0_link_slope)

          celerity0(frow, fcol) = stack(ns)
          ns = ns + 1
          fId = iD0(frow, fcol)
          if( .NOT. (fID == tID)) then
            call append(stack, append_chunk)
          else
            exit
          end if
          ! move downstream
          call moveDownOneCell( fDir0(frow,fcol), frow, fcol )
        end do

        celerity11(ii) = size(stack) / sum(1/stack(:))
        deallocate(stack)
        allocate(stack(1))

      end do

    else

      ! There is only one cell, so no routing is taking place
      ! set dummy value of 1 m / s
      celerity11(:) =  1._dp

    end if

    ! Write celerity
    L11_celerity(iStart11:iEnd11) = celerity11(:)
    L0_celerity(iStart0:iEnd0) = PACK(celerity0(:,:), mask0)

    ! free space
    deallocate (&
        mask0, iD0, slope_tmp, slopemask0, &
        slope0, fDir0, cellarea0,   &
        stack, netPerm, nLinkFromRow, nLinkFromCol, nLinkToRow, nLinkToCol)

  end subroutine L11_calc_celerity

end module mo_mrm_net_startup
