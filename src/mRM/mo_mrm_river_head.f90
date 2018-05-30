module mo_mrm_river_head
    !use mo_mhm_global_variables, only: L0_slope
    use mo_mrm_global_variables, only: nBasins, L0_elev_mRM, &
        L11_bankfull_runoff_in, L0_channel_depth, L0_channel_elevation, &
        L0_streamNet, L0_L11_Id, L0_Id, L0_fDir, iFlag_cordinate_sys
    use mo_mrm_constants,        only: nodata_i4, nodata_dp
    use mo_mrm_tools,            only: get_basin_info_mrm
    use mo_mrm_net_startup,      only: cellLength, moveDownOneCell
    use mo_julian,               only: dec2date
    use mo_kind,                 only: i4, dp
    use mo_common_variables,     only: period
    use mo_message,              only: message
    use mo_mhm_constants,        only: DaySecs

    implicit none

    private

    public :: calc_channel_elevation ! calculates the channel elevation from the bakfull river discharge

    contains

    subroutine calc_channel_elevation()
        implicit none

        real(dp) n ! Manning's roughness coefficient
        real(dp), dimension(:,:), allocatable :: elev0
        integer(i4), dimension(:,:), allocatable :: iD0            
        integer(i4), dimension(:,:), allocatable :: fDir0
        logical, dimension(:,:), allocatable :: mask0
        real(dp),    dimension(:,:), allocatable :: nodata_dp_tmp
        integer(i4), dimension(:,:), allocatable :: nodata_i4_tmp
        real(dp) :: length, slope
        integer(i4) :: nCells0
        integer(i4) :: nrows0, ncols0
        integer(i4) :: iStart0, iEnd0
        integer(i4) i, j, i_down, j_down
        integer(i4) iBasin

        n = .045_dp ! m^-1/3 s from Sutanudjaja et al. 2011

        allocate(L0_channel_depth(nBasins, size(L0_streamNet)))
        allocate(L0_channel_elevation(nBasins, size(L0_streamNet)))

        do iBasin = 1, nBasins
            do i = 1, size(L0_streamNet)
                L0_channel_depth(iBasin, i) = nodata_dp
                L0_channel_elevation(iBasin, i) = nodata_dp
            end do
            ! level-0 information
            call get_basin_info_mrm(iBasin, 0, nrows0, ncols0, ncells=nCells0, &
                                    iStart=iStart0, iEnd=iEnd0, mask=mask0) 

            allocate(iD0(nrows0, ncols0))
            allocate(elev0(nrows0, ncols0))
            allocate(fDir0(nrows0, ncols0))
            allocate(nodata_dp_tmp(nrows0, ncols0))
            allocate(nodata_i4_tmp(nrows0, ncols0))

            nodata_dp_tmp(:,:) = nodata_dp
            nodata_i4_tmp(:,:) = nodata_i4

            elev0(:,:) = unpack(L0_elev_mRM(iStart0:iEnd0), mask0, nodata_dp_tmp)
            iD0(:,:) = unpack(L0_Id(iStart0:iEnd0), mask0, nodata_i4_tmp)
            fDir0(:,:) = unpack(L0_fDir(iStart0:iEnd0), mask0, nodata_i4_tmp)

            do j = 1, ncols0
                do i = 1, nrows0
                    if (mask0(i,j)) then
                        call cellLength(iBasin, fDir0(i, j), i, j, &
                                        iFlag_cordinate_sys, length)
                        call moveDownOneCell(fDir0(i, j), i_down, j_down)
                        slope = (elev0(i,j) - elev0(i_down, j_down)) / length
                        if (slope < 0.0001_dp .OR. slope > 20._dp .OR. isnan(slope)) then
                            slope = 0.0001_dp
                        end if
                        L0_channel_depth(iBasin, iD0(i,j)) = &
                            ((n * L11_bankfull_runoff_in(L0_L11_Id(iD0(i,j)))) &
                            / (4.8 * slope**.5))**.6
                        L0_channel_elevation(iBasin, iD0(i,j)) = &
                            elev0(i,j) - L0_channel_depth(iBasin, ID0(i,j))
                    end if
                end do
            end do

            deallocate(iD0)
            deallocate(elev0)
            deallocate(fDir0)
            deallocate(nodata_dp_tmp)
            deallocate(nodata_i4_tmp)
        end do
    end subroutine calc_channel_elevation

end module mo_mrm_river_head
