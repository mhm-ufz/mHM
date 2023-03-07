!> \file mo_mpr_neutrons.f90
!> \brief \copybrief mo_mpr_neutrons
!> \details \copydetails mo_mpr_neutrons

!> \brief   Multiscale parameter regionalization (MPR) for neutrons
!> \details This module contains all routines required for parametrizing neutrons processes.
!> \author Maren Kaluza
!> \date Dec 2017
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mpr
module mo_mpr_neutrons

  use mo_kind, only: i4, dp

  implicit none

  public :: mpr_neutrons

  private

contains
  ! ----------------------------------------------------------------------------

  !      NAME
  !         mpr_neutrons

  !>        \brief multiscale parameter regionalization for neutrons

  !>        \details calculates neutron variables on L0
  !>                 Global parameters needed (see mhm_parameter.nml):\n
  !>                    - param( 1) = Desilets_N0   \n
  !>                    - param( 2) = COSMIC_N0     \n
  !>                    - param( 3) = COSMIC_N1     \n
  !>                    - param( 4) = COSMIC_N2     \n
  !>                    - param( 5) = COSMIC_alpha0 \n
  !>                    - param( 6) = COSMIC_alpha1 \n
  !>                    - param( 7) = COSMIC_L30    \n
  !>                    - param( 8) = COSMIC_L31    \n
  !>                    - param( 9) = COSMIC_LW0    \n
  !>                    - param(10) = COSMIC_LW1    \n

  !      INTENT(IN)
  !>        \param[in] "real(dp)    :: param(10)"        - global parameters
  !>        \param[in] "integer(i4) :: is_present(:)"    - indicates whether soiltype is present
  !>        \param[in] "integer(i4) :: nHorizons(:)"     - Number of Horizons per soiltype2
  !>        \param[in] "integer(i4) :: nTillHorizons(:)" - Number of Tillage Horizons
  !>        \param[in] "integer(i4) :: LCover0(:)"       - land cover ids at level 0
  !>        \param[in] "real(dp)    :: clay(:,:)"        - clay content
  !>        \param[in] "real(dp)    :: DbM(:,:)"         - mineral Bulk density
  !>        \param[in] "real(dp)    :: Db(:,:)"          - Bulk density

  !     INTENT(INOUT)
  !         None

  !>      INTENT(OUT)
  !>        \param[out] "real(dp)   :: COSMIC_L3_till(:,:,:)" - COSMIC paramter L3 tillage layer
  !>        \param[out] "real(dp)   :: latWat_till(:,:,:)"    - lattice water content tillage layer
  !>        \param[out] "real(dp)   :: COSMIC_L3(:,:)"        - COSMIC paramter L3 tillage layer
  !>        \param[out] "real(dp)   :: latWat(:,:)"           - lattice water contente

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
  !>        \author Maren Kaluza
  !>        \date Dec 2017


  subroutine mpr_neutrons( process_case, & ! IN: process case
       param               , & ! IN:  global parameter set
       is_present          , & ! IN:  flag indicating presence of soil
       nHorizons           , & ! IN:  Number of Horizons of Soiltype
       nTillHorizons       , & ! IN:  Number of tillage Horizons
       LCover0             , & ! IN:  land cover ids at level 0
       clay                , & ! IN:  clay content
       DbM                 , & ! IN:  mineral Bulk density
       Db                  , & ! IN: Bulk density
       COSMIC_L3_till      , & ! OUT: COSMIC paramter L3 tillage layer
       latWat_till         , & ! OUT: lattice water content tillage layer
       COSMIC_L3           , & ! OUT: COSMIC paramter L3 tillage layer
       latWat                & ! OUT: lattice water contente
       )

    ! lots of lines copy-pasted from mo_mpr_soilmoist.f90
    use mo_message, only: error_message
    use mo_mpr_global_variables, only: iFlag_soilDB
    !$  use omp_lib

    implicit none

    ! Input --------------------------------------------------------------------
    integer(i4),                   intent(in)  :: process_case ! process case
    real(dp),    dimension(:),     intent(in)  :: param        ! global parameters   !! dim = 3 for case 1 and 9 for case 2
    integer(i4), dimension(:),     intent(in)  :: is_present   ! indicates whether soiltype is present
    integer(i4), dimension(:),     intent(in)  :: nHorizons    ! Number of Horizons per soiltype
    integer(i4), dimension(:),     intent(in)  :: nTillHorizons! Number of Tillage Horizons
    real(dp),    dimension(:,:),   intent(in)  :: DbM          ! mineral Bulk density
    real(dp),    dimension(:,:,:), intent(in)  :: Db           ! Bulk density
    integer(i4), dimension(:),     intent(in)  :: LCOVER0      ! land cover ids at level 0
    real(dp),    dimension(:,:),   intent(in)  :: clay         ! clay content

    ! Output -------------------------------------------------------------------
    real(dp),    dimension(:,:,:), intent(out) :: COSMIC_L3_till ! COSMIC parameter L3 tillage layer
    real(dp),    dimension(:,:,:), intent(out) :: latWat_till    ! lattice water content tillage layer
    real(dp),    dimension(:,:),   intent(out) :: COSMIC_L3      ! COSMIC parameter L3
    real(dp),    dimension(:,:),   intent(out) :: latWat         ! lattice water content

    ! Local variables
    integer(i4)                               :: i               ! loop index
    integer(i4)                               :: j               ! loop index
    integer(i4)                               :: l               ! loop index
    integer(i4)                               :: tmp_minSoilHorizon


    !min soil horizon
    tmp_minSoilHorizon = minval(nTillHorizons(:))

    ! with zero there will be problem with
    ! upscaling with harmonic mean for the COMSIC_L3
    ! in case of process_case .EQ. 1
    COSMIC_L3_till = 0.000001_dp
    COSMIC_L3      = 0.000001_dp
    latWat_till    = 0.000001_dp
    latWat         = 0.000001_dp

    ! select case according to a given soil database flag
    SELECT CASE(iFlag_soilDB)

       ! classical mHM soil database format
       CASE(0)
          do i = 1, size(is_present)
             if ( is_present(i) .lt. 1 ) cycle
             horizon: do j = 1, nHorizons(i)
                ! calculating other soil hydraulic properties
                ! tillage horizons
                if ( j .le. nTillHorizons(i) ) then
                   ! LC class
                   do L = 1, maxval( LCOVER0 )
                      if(process_case .EQ. 1) call latticeWater(param(2:3), clay(i,j), latWat_till(i,j,L))
                      if(process_case .EQ. 2) then
                         call calcL3(param(6:7), Db(i,j,L), COSMIC_L3_till(i,j,L))
                         call latticeWater(param(8:9), clay(i,j), latWat_till(i,j,L))
                      end if
                   end do
                ! deeper layers
                else
                   if(process_case .EQ. 1) call latticeWater(param(2:3), clay(i,j), latWat(i,j-tmp_minSoilHorizon))
                   if(process_case .EQ. 2) then
                      call calcL3(param(6:7), DbM(i,j), COSMIC_L3(i,j-tmp_minSoilHorizon))
                      call latticeWater(param(8:9), clay(i,j), latWat(i,j-tmp_minSoilHorizon))
                   end if
                end if
             end do horizon
          end do

       ! to handle multiple soil horizons with unique soil class
       CASE(1)
           do i = 1, size(is_present)
             if ( is_present(i) .lt. 1 ) cycle
             ! **** FOR THE TILLAGE TYPE OF SOIL *****
             ! there is actually no soil horizons/soil type in this case
             ! but we assign of j = 1 to use variables as defined in the classical option (iFlag_soil = 0)
             do j = 1, 1
                ! tillage horizons properties depending on the LC class
                do L = 1, maxval( LCOVER0 )
                   if(process_case .EQ. 1) call latticeWater(param(2:3), clay(i,j), latWat_till(i,j,L))
                   if(process_case .EQ. 2) then
                      call calcL3(param(6:7), Db(i,j,L), COSMIC_L3_till(i,j,L))
                      call latticeWater(param(8:9), clay(i,j), latWat_till(i,j,L))
                   end if
                end do

                ! *** FOR NON-TILLAGE TYPE OF SOILS ***
                ! note j = 1
                if(process_case .EQ. 1) call latticeWater(param(2:3), clay(i,j), latWat(i,j))
                if(process_case .EQ. 2) then
                   call calcL3(param(6:7), DbM(i,j), COSMIC_L3(i,j))
                   call latticeWater(param(8:9), clay(i,j), latWat(i,j))
                end if

             end do  ! >> HORIZON
          end do   ! >> SOIL TYPE

       CASE DEFAULT
          call error_message('***ERROR: iFlag_soilDB option given does not exist. Only 0 and 1 is taken at the moment.')
       END SELECT
       !

   end subroutine


  !! >> L3 parameter
  subroutine calcL3(param, bulkDensity, L3)
    ! param(1) = COSMIC_L30
    ! param(2) = COSMIC_L31
    implicit none
    real(dp), dimension(2),  intent(in)       :: param
    real(dp),                intent(in)       :: bulkDensity
    real(dp),                intent(inout)    :: L3

    L3 = bulkDensity*param(1) - param(2)
    if( bulkDensity .LT. 0.4_dp ) then ! bulkDensity<0.39 yields negative L3, bulkDensity=0.39 yields L3=0
       L3 = 1.0_dp                     ! Prevent division by zero later on; added by joost Iwema to COSMIC 1.13, Feb. 2017
    endif

  end subroutine calcL3


  !! >>>> lattice water
  subroutine latticeWater( param, clay, latWat )
    ! param(1) = COSMIC_LW0 or deslet_LW0
    ! param(2) = COSMIC_LW1 or deslet_LW0
    implicit none
    ! Input
    real(dp), dimension(2), intent(in)  :: param
    real(dp),               intent(in)  :: clay
    ! Output
    real(dp),               intent(out) :: latWat

    !Martin Schroen's dissertation
    latWat = ( param(1)*clay/100.0_dp + param(2) )

  end subroutine latticeWater

end module mo_mpr_neutrons
