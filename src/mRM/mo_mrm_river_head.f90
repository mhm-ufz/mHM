module mo_mrm_river_head
    use mo_common_variables,     only : level0
    use mo_mrm_global_variables, only : L0_L11_remap, L11_bankfull_runoff_in, &
                                        L0_slope, L0_channel_elevation
    use mo_kind,                 only : i4, dp
    use mo_append,               only : append


    implicit none

    private

    public :: init_river_head ! allocates memory

    public :: calc_channel_elevation ! calculates the channel elevation from the bakfull river discharge
    public :: calc_river_head ! calculates the river head


    contains

    subroutine init_river_head(iBasin, river_head)
        use mo_common_constants, only : nodata_dp
        integer(i4), intent(in) :: iBasin
        real(dp), dimension(:), allocatable, intent(inout) :: river_head
        real(dp), dimension(:), allocatable :: dummy_1D

        allocate(dummy_1D(level0(iBasin)%nCells))
        dummy_1D = nodata_dp
        call append(river_head, dummy_1D)
        deallocate(dummy_1D)
    end subroutine init_river_head


    subroutine calc_channel_elevation()
        use mo_common_constants, only : nodata_i4, nodata_dp
        use mo_common_variables, only : nBasins, L0_elev
        use mo_mrm_global_variables, only : L0_fDir, L0_channel_depth
        real(dp), dimension(:,:), allocatable :: channel_dpth
        real(dp), dimension(:,:), allocatable :: channel_elev
        real(dp), dimension(:,:), allocatable :: slope
        real(dp), dimension(:,:), allocatable :: elev0
        integer(i4), dimension(:,:), allocatable :: fDir0
        real(dp) n ! Manning's roughness coefficient
        integer(i4) :: nrows0, ncols0
        integer(i4) :: s0, e0
        integer(i4) i, j, k
        integer(i4) iBasin

        n = .045_dp ! m^-1/3 s from Sutanudjaja et al. 2011

        do iBasin = 1, nBasins
            nrows0 = level0(iBasin)%nrows
            ncols0 = level0(iBasin)%ncols
            s0 = level0(iBasin)%iStart
            e0 = level0(iBasin)%iEnd
            allocate(channel_dpth(nrows0, ncols0))
            allocate(channel_elev(nrows0, ncols0))
            allocate(elev0(nrows0, ncols0))
            allocate(fDir0(nrows0, ncols0))
            allocate(slope(nrows0, ncols0))
            channel_dpth(:,:) = nodata_dp
            channel_elev(:,:) = nodata_dp
            slope(:,:) = nodata_dp

            elev0(:,:) = unpack(L0_elev(s0:e0), level0(iBasin)%mask, nodata_dp)
            fDir0(:,:) = unpack(L0_fDir(s0:e0), level0(iBasin)%mask, nodata_i4)

            do k = 1, level0(iBasin)%nCells
                i = level0(iBasin)%CellCoor(k, 1)
                j = level0(iBasin)%CellCoor(k, 2)
                if (fDir0(i,j) /= nodata_i4) then
                    slope(i,j) = calc_slope(iBasin, elev0, fDir0, i, j)
                    channel_dpth(i,j) = &
    ((n * L11_bankfull_runoff_in(L0_L11_remap(iBasin)%lowres_id_on_highres(i,j))) &
    / (4.8 * slope(i,j)**.5))**.6
                    channel_elev(i,j) = elev0(i,j) - channel_dpth(i,j)
                end if
            end do

            call append(L0_channel_depth, pack(channel_dpth(:,:), &
                        level0(iBasin)%mask))
            call append(L0_channel_elevation, pack(channel_elev(:,:), &
                        level0(iBasin)%mask))
            call append(L0_slope, pack(slope(:,:), &
                        level0(iBasin)%mask))

            deallocate(channel_dpth)
            deallocate(channel_elev)
            deallocate(elev0)
            deallocate(fDir0)
            deallocate(slope)
        end do
    end subroutine calc_channel_elevation


    subroutine calc_river_head(iBasin, L11_Qmod, river_head)
        integer(i4), intent(in) :: iBasin
        real(dp), dimension(:), intent(in) :: L11_Qmod ! modelled discharge at each grid cell
        real(dp), dimension(:), allocatable, intent(inout) :: river_head
        !real(dp), dimension(:,:), allocatable :: channel_elev
        !real(dp), dimension(:,:), allocatable :: river_head_2d
        !real(dp), dimension(:,:), allocatable :: slope
        real(dp) :: n ! Manning's roughness coefficient
        integer(i4) :: nrows0, ncols0
        integer(i4) :: s0, e0
        integer(i4) i, j, k, L11_ind

        n = .045_dp ! m^-1/3 s from Sutanudjaja et al. 2011

        nrows0 = level0(iBasin)%nrows
        ncols0 = level0(iBasin)%ncols
        s0 = level0(iBasin)%iStart
        e0 = level0(iBasin)%iEnd

        do k = s0, e0
            i = level0(iBasin)%CellCoor(k, 1)
            j = level0(iBasin)%CellCoor(k, 2)
            if (i >= 0 .and. i < 99999 .and. j  >= 0 .and. j < 99999) then
                L11_ind = L0_L11_remap(iBasin)%lowres_id_on_highres(i,j)
                ! TODO L11_Qmid(L11_ind) causes IEEE_UNDERFLOW_FLAG IEEE_DENORMAL
                river_head(k) = L0_channel_elevation(k) + &
                                (n * L11_Qmod(L11_ind) / &
                                L11_bankfull_runoff_in(L11_ind) / L0_slope(k)**.5)**.6
            end if
        end do

    end subroutine calc_river_head


    function calc_slope(iBasin, elev0, fDir0, i, j) result(slope)
        use mo_common_variables, only: iFlag_cordinate_sys
        use mo_mrm_net_startup,      only: cellLength, moveDownOneCell
        integer(i4), intent(in) :: iBasin
        integer(i4), intent(in) :: i, j
        integer(i4), intent(in), dimension(:,:), allocatable :: fDir0
        real(dp), intent(in), dimension(:,:), allocatable :: elev0
        real(dp) :: slope, length
        integer(i4) :: i_down, j_down

        call cellLength(iBasin, fDir0(i, j), i, j, &
                        iFlag_cordinate_sys, length)
        i_down = i
        j_down = j
        call moveDownOneCell(fDir0(i, j), i_down, j_down)

        slope = (elev0(i,j) - elev0(i_down, j_down)) / length
        if(slope < 0.0001_dp .OR. slope > 20._dp .OR. isnan(slope)) then
           slope = 0.0001_dp
        end if
    end function calc_slope

end module mo_mrm_river_head
