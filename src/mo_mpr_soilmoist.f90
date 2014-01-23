!> \file mo_mpr_soilmoist.f90

!> \brief multiscale parameter regionalization for soil moisture

!> \details This module contains all routines required for parametrizing
!> soil moisture processes.

!> \author Stephan Thober, Rohini Kumar
!> \date Dec 2012

module mo_mpr_soilmoist
  

  use mo_kind, only: i4, dp

  implicit none

  public :: mpr_sm

  private

contains
  ! ----------------------------------------------------------------------------

  !      NAME
  !         mpr_sm

  !>        \brief multiscale parameter regionalization for soil moisture

  !>        \details This subroutine is a wrapper around all soil moisture
  !>        parameter routines. This subroutine requires 13 parameters. These
  !>        parameters have to correspond to the parameters in the original
  !>        parameter array at the following locations: 10-12, 13-18, 27-30.\n
  !>        Global parameters needed (see mhm_parameter.nml):\n
  !>           - param( 1) = orgMatterContent_forest     \n
  !>           - param( 2) = orgMatterContent_impervious \n
  !>           - param( 3) = orgMatterContent_pervious   \n
  !>           - param( 4) = PTF_lower66_5_constant      \n
  !>           - param( 5) = PTF_lower66_5_clay          \n
  !>           - param( 6) = PTF_lower66_5_Db            \n
  !>           - param( 7) = PTF_higher66_5_constant     \n
  !>           - param( 8) = PTF_higher66_5_clay         \n
  !>           - param( 9) = PTF_higher66_5_Db           \n
  !>           - param(10) = PTF_Ks_constant             \n
  !>           - param(11) = PTF_Ks_sand                 \n
  !>           - param(12) = PTF_Ks_clay                 \n
  !>           - param(13) = PTF_Ks_curveSlope           \n

  !      INTENT(IN)
  !>        \param[in] "real(dp)    :: param(13)"        - global parameters
  !>        \param[in] "real(dp)    :: nodata"           - no data value
  !>        \param[in] "integer(i4) :: is_present(:)"    - indicates whether soiltype is present
  !>        \param[in] "integer(i4) :: nHorizons(:)"     - Number of Horizons per soiltype2
  !>        \param[in] "integer(i4) :: nTillHorizons(:)" - Number of Tillage Horizons
  !>        \param[in] "real(dp)    :: sand(:,:)"        - sand content
  !>        \param[in] "real(dp)    :: clay(:,:)"        - clay content
  !>        \param[in] "real(dp)    :: DbM(:,:)"         - mineral Bulk density
  !>        \param[in] "integer(i4) :: L0_ID(:,:)"       - cell ids at level 0
  !>        \param[in] "integer(i4) :: L0_soilId(:,:)"   - soil ids at level 0
  !>        \param[in] "integer(i4) :: L0_LUC(:,:)"      - land cover ids at level 0

  !      INTENT(OUT)
  !>        \param[out] "real(dp)   :: thetaS_till(:,:,:)" - saturated soil moisture tillage layer
  !>        \param[out] "real(dp)   :: thetaFC_till(:,:,:)" - field capacity tillage layer
  !>        \param[out] "real(dp)   :: thetaPW_till(:,:,:)" - permanent wilting point tillage layer
  !>        \param[out] "real(dp)   :: thetaS(:,:)"     - saturated soil moisture
  !>        \param[out] "real(dp)   :: thetaFC(:,:)"    - field capacity
  !>        \param[out] "real(dp)   :: thetaPW(:,:)"    - permanent wilting point
  !>        \param[out] "real(dp)   :: Ks(:,:,:)"       - saturated hydraulic conductivity
  !>        \param[out] "real(dp)   :: Db(:,:,:)"       - Bulk density
  !>        \param[out] "real(dp)   :: L0_KsVar_H(:,:)" - relative variability of saturated
  !>                                                      hydraulic counductivity for Horizantal flow
  !>        \param[out] "real(dp)   :: L0_KsVar_V(:,:)" - relative variability of saturated
  !>                                                      hydraulic counductivity for Vertical flow
  !>        \param[out] "real(dp)   :: L0_SMs_tot(:,:)" - total saturated soil mositure content
  !>        \param[out] "real(dp)   :: L0_SMs_FC(:,:)"  - soil mositure deficit from field
  !>                                                      capacity w.r.t to saturation

  !      HISTORY
  !>        \author Stephan Thober, Rohini Kumar
  !>        \date Dec 2012
  !         Written,  Stephan Thober, Dec 2012
  !         Modified, Juliane Mai,    Oct 2013 - OLD parametrization
  !                                                --> param(1) = orgMatterContent_forest
  !                                                --> param(2) = orgMatterContent_impervious 
  !                                                --> param(3) = orgMatterContent_pervious
  !                                                --> param(4:13) = ...
  !                                             -------------------------------
  !                                             orgMatterContent_forest = orgMatterContent_perv + delta_1
  !                                             -------------------------------
  !                                             NEW parametrization
  !                                                --> param(1) = delta_1 
  !                                                --> param(2) = orgMatterContent_impervious 
  !                                                --> param(3) = orgMatterContent_pervious  
  !                                                --> param(4:13) = ...
  !         Modified, Matthias Zink,  Nov 2013 - documentation, inouts --> out 
  !                                              moved constants to mhm_constants

  subroutine mpr_sm( &
       ! Input -----------------------------------------------------------------
       param        , & ! global parameter set
       nodata       , & ! nodata value
       is_present   , & ! flag indicating presence of soil
       nHorizons    , & ! Number of Horizons of Soiltype
       nTillHorizons, & ! Number of tillage Horizons
       sand         , & ! sand content
       clay         , & ! clay content
       DbM          , & ! mineral Bulk density
       ID0          , & ! cell ids at level 0
       soilId0      , & ! soil ids at level 0
       LCover0      , & ! land cover ids at level 0
       ! Output ----------------------------------------------------------------
       thetaS_till  , & ! saturated soil moisture tillage layer
       thetaFC_till , & ! field capacity tillage layer
       thetaPW_till , & ! permanent wilting point tillage layer
       thetaS       , & ! saturated soil moisture
       thetaFC      , & ! field capacity
       thetaPW      , & ! permanent wilting point
       Ks           , & ! saturated hydraulic conductivity
       Db           , & ! Bulk density
       KsVar_H0     , & ! relative variability of saturated
                        ! hydraulic counductivity for Horizantal flow
       KsVar_V0     , & ! relative variability of saturated
                        ! hydraulic counductivity for Horizantal flow
       SMs_tot0     , & ! total saturated soil mositure content
       SMs_FC0    )     ! soil mositure deficit from field capacity
                        ! w.r.t to saturation

    use mo_mhm_constants, only: BulkDens_OrgMatter

    implicit none

    ! Input --------------------------------------------------------------------
    real(dp),    dimension(13),    intent(in)  :: param        ! global parameters
    real(dp),                      intent(in)  :: nodata       ! no data value

    integer(i4), dimension(:),     intent(in)  :: is_present   ! indicates whether soiltype is present
    integer(i4), dimension(:),     intent(in)  :: nHorizons    ! Number of Horizons per soiltype
    integer(i4), dimension(:),     intent(in)  :: nTillHorizons! Number of Tillage Horizons
    real(dp),    dimension(:,:),   intent(in)  :: sand         ! sand content
    real(dp),    dimension(:,:),   intent(in)  :: clay         ! clay content
    real(dp),    dimension(:,:),   intent(in)  :: DbM          ! mineral Bulk density
    integer(i4), dimension(:),     intent(in)  :: ID0        ! cell ids at level 0
    integer(i4), dimension(:),     intent(in)  :: soilId0    ! soil ids at level 0
    integer(i4), dimension(:),     intent(in)  :: LCOVER0       ! land cover ids at level 0


    ! Output -------------------------------------------------------------------
    real(dp),    dimension(:,:,:), intent(out) :: thetaS_till   ! saturated soil moisture tillage layer
    real(dp),    dimension(:,:,:), intent(out) :: thetaFC_till  ! field capacity tillage layer
    real(dp),    dimension(:,:,:), intent(out) :: thetaPW_till  ! permanent wilting point tillage layer
    real(dp),    dimension(:,:),   intent(out) :: thetaS        ! saturated soil moisture
    real(dp),    dimension(:,:),   intent(out) :: thetaFC       ! field capacity
    real(dp),    dimension(:,:),   intent(out) :: thetaPW       ! permanent wilting point
    real(dp),    dimension(:,:,:), intent(out) :: Ks            ! saturated hydraulic conductivity
    real(dp),    dimension(:,:,:), intent(out) :: Db            ! Bulk density
    real(dp),    dimension(:),     intent(out) :: KsVar_H0      ! relative variability of
                                                                  ! saturated hydraulic
                                                                  ! cound. for Horizantal flow
    real(dp),    dimension(:),     intent(out) :: KsVar_V0      ! relative variability of
                                                                  ! saturated hydraulic
                                                                  ! cound. for vertical flow
    real(dp),    dimension(:),     intent(out) :: SMs_tot0      ! total saturated soilmositure
                                                                  ! content
    real(dp),    dimension(:),     intent(out) :: SMs_FC0       ! soil mositure deficit from
                                                                  ! field cap. w.r.t to saturation

    ! Local variables
    integer(i4)                               :: i               ! loop index
    integer(i4)                               :: j               ! loop index
    integer(i4)                               :: l               ! loop index
    real(dp)                                  :: pM
    real(dp)                                  :: pOM
    real(dp)                                  :: Ks_tmp          ! temporal saturated hydr. cond
    real(dp)                                  :: Genu_Mual_n     ! van Genuchten shape param
    real(dp)                                  :: Genu_Mual_alpha ! van Genuchten shape param
    real(dp)                                  :: tmp_orgMatterContent_forest
    real(dp)                                  :: tmp_orgMatterContent_pervious
    real(dp)                                  :: tmp_orgMatterContent_impervious

    tmp_orgMatterContent_forest     = param(3) + param(1)
    tmp_orgMatterContent_impervious = param(2)
    tmp_orgMatterContent_pervious   = param(3)

    ! initializing soil hydraulic properties
    KsVar_H0 = merge( 0.0_dp, nodata, ID0 /= int(nodata,i4) )
    KsVar_V0 = merge( 0.0_dp, nodata, ID0 /= int(nodata,i4) )
    SMs_tot0 = merge( 0.0_dp, nodata, ID0 /= int(nodata,i4) )
    SMs_FC0  = merge( 0.0_dp, nodata, ID0 /= int(nodata,i4) )

    ! Initialization
    thetaS_till  = 0.0_dp
    thetaFC_till = 0.0_dp
    thetaPW_till = 0.0_dp
    thetaS       = 0.0_dp
    thetaFC      = 0.0_dp
    thetaPW      = 0.0_dp
    Ks           = 0.0_dp  
    Db           = 0.0_dp

    do i = 1, size(is_present)

       if ( is_present(i) < 1 ) cycle
       
       horizon: do j = 1, nHorizons(i)

          ! calculating vertical hydraulic conductivity 
          call hydro_cond( Ks_tmp, param(10:13), sand(i,j), clay(i,j) )
          Ks(i,j,:) = Ks_tmp

          ! calculating other soil hydraulic
          ! tillage horizons
          
          if ( j <= nTillHorizons(i) ) then
             
             ! LC class
             do L = 1, maxval( LCOVER0 )
                select case (L)
                case(1)               ! forest
                   pOM = tmp_orgMatterContent_forest
                case(2)               ! impervious
                   pOM = tmp_orgMatterContent_impervious !param(2)  
                case(3)               ! permeable   
                   pOM = tmp_orgMatterContent_pervious
                case default
                   stop 'Error mpr_sm: pOM used uninitialized.'
                end select
                pM = 100.0_dp - pOM

                ! bulk density acording to Rawl's (1982) paper 
                Db(i,j,L) = 100.0_dp / ( (pOM/BulkDens_OrgMatter) + (pM/DbM(i,j)) ) 

                ! Effect of organic matter content
                ! This is taken into account in a simplified form by using
                ! the ratio of(Bd / BdOM)
                Ks_tmp = Ks_tmp * ( DbM(i,j) / Db(i,j,L))
                Ks(i,j,L) =  Ks_tmp

                ! estimated SMs_till & van Genuchten's shape parameter (n)
                call Genuchten( thetaS_till(i,j,L), Genu_Mual_n, Genu_Mual_alpha, & 
                                param(4:9), sand(i,j), clay(i,j), Db(i,j,L)       )

                ! estimating field capacity
                call field_cap( thetaFC_till(i,j,L), Ks_tmp, thetaS_till(i,j,L), Genu_Mual_n )

                ! estimating permanent wilting point
                call PWP( Genu_Mual_n, Genu_Mual_alpha, thetaS_till(i,j,L), thetaPW_till(i,j,L) )


                ! Soil properties over the whole soil coloum depth
                KsVar_H0 = merge( KsVar_H0 + thetaS_till (i,j,L) * Ks_tmp, &
                           KsVar_H0, ( soilId0 == i .and. LCOVER0 == L ) )
                KsVar_V0 = merge( KsVar_V0 + thetaS_till (i,j,L) / Ks_tmp, &
                           KsVar_V0, ( soilId0 == i .and. LCOVER0 == L ) )
                SMs_FC0  = merge( SMs_FC0  + thetaFC_till(i,j,L),          &
                           SMs_FC0,  ( soilId0 == i .and. LCOVER0 == L ) )
                SMs_tot0 = merge( SMs_tot0 + thetaS_till (i,j,L),          &
                           SMs_tot0, ( soilId0 == i .and. LCOVER0 == L ) )
 
             end do
             
             ! deeper layers
          else
             
             ! estimate SMs & van Genuchten's shape parameter (n)
             call Genuchten( thetaS(i, j-minval(nTillHorizons(:))), Genu_Mual_n, Genu_Mual_alpha, &
                  param(4:9), sand(i,j), clay(i,j), DbM(i,j) )
             
             ! estimate field capacity
             call field_cap( thetaFC(i,j-minval(nTillHorizons(:))), &
                  Ks_tmp, thetaS(i,j-minval(nTillHorizons(:))), Genu_Mual_n )
             
             ! estimate permanent wilting point
             call PWP( Genu_Mual_n, Genu_Mual_alpha, thetaS(i, j-minval(nTillHorizons(:))), &
                  thetaPW(i, j-minval(nTillHorizons(:))) )
             
             ! soil_properties over the whole soil column
             KsVar_H0 = merge( KsVar_H0 + thetaS (i, j-minval(nTillHorizons(:))) * Ks_tmp, &
                        KsVar_H0, ( soilId0 == i ) )
             KsVar_V0 = merge( KsVar_V0 + thetaS (i, j-minval(nTillHorizons(:))) / Ks_tmp, &
                        KsVar_V0, ( soilId0 == i ) )
             SMs_FC0  = merge( SMs_FC0  + thetaFC(i, j-minval(nTillHorizons(:))),          &
                        SMs_FC0, ( soilId0 == i ) )
             SMs_tot0 = merge( SMs_tot0 + thetaS (i, j-minval(nTillHorizons(:))),          &
                        SMs_tot0, ( soilId0 == i ) )
             
          end if
         
       end do horizon

       ! ------------------------------------------------------------------
       ! DETERMINE RELATIVE VARIABILITIES OF 
       !   Ks FOR HORIZONTAL FLOW (KsVar_H)
       !               &
       !   Ks FOR VERTICAL FLOW (KsVar_V)
       ! ------------------------------------------------------------------

       ! soil moisture saturation deficit relative to the field capacity soil moisture
       SMs_FC0  = merge( (SMs_tot0 - SMs_FC0) / SMs_tot0, SMs_FC0,  (soilId0 == i) )

       ! Ks variability over the whole soil coloum depth for
       ! both horizontal and vertical flows 
       KsVar_H0 = merge(  KsVar_H0 / SMs_tot0, KsVar_H0, (soilId0 == i) )
       KsVar_V0 = merge(  SMs_tot0 / KsVar_V0, KsVar_V0, (soilId0 == i) )

       ! relative variabilities
       KsVar_H0 = merge(  KsVar_H0 / param(13), KsVar_H0, (soilId0 == i) )
       KsVar_V0 = merge(  KsVar_V0 / param(13), KsVar_V0, (soilId0 == i) )

    end do
    
  end subroutine mpr_sm

  ! ------------------------------------------------------------------

  !      NAME
  !         PWP

  !>        \brief Permanent Wilting point

  !>        \details This subroutine calculates the permanent wilting
  !>        point according to Zacharias et al. (2007, Soil Phy.) and
  !>        using van Genuchten 1980's equation. For the water retention curve at 
  !>        a matrix potential of -1500 kPa, it is assumed that thetaR = 0.

  !      INTENT(IN)
  !>        \param[in] "real(dp) :: Genu_Mual_n"     - Genuchten shape parameter
  !>        \param[in] "real(dp) :: Genu_Mual_alpha" - Genuchten shape parameter
  !>        \param[in] "real(dp) :: thetaS"          - saturated water content

  !      INTENT(OUT)
  !>        \param[out] "real(dp) :: thetaPWP"       - Permanent Wilting point

  !      LITERATURE
  !         Zacharias et al. 2007, Soil Phy.

  !      HISTORY
  !>        \author Stephan Thober, Rohini Kumar
  !>        \date Dec, 2012
  !         Written, Stephan Thober, Dec 2012
  !         Modified, Matthias Zink,  Nov 2013 - documentation, moved constants to mhm_constants
  ! ------------------------------------------------------------------

  elemental pure subroutine PWP( &
       ! Input variables
       Genu_Mual_n     , & ! Genuchten shape parameters   
       Genu_Mual_alpha , & ! Genuchten shape parameters
       thetaS          , & ! saturated water content
       ! Output
       thetaPWP )          ! Permanent wilting point

    use mo_mhm_constants , only: PWP_c, PWP_matPot_ThetaR ! constant for m, 
    !                                                      ! matrix potential of 1500 kPa, assumed as thetaR = 0

    implicit none

    ! Input
    real(dp), intent(in)  :: Genu_Mual_n     ! Genuchten shape parameters
    real(dp), intent(in)  :: Genu_Mual_alpha ! Genuchten shape parameters
    real(dp), intent(in)  :: thetaS          ! saturated water content

    ! Output
    real(dp), intent(out) :: thetaPWP        ! Permanent wilting point

    ! Local variable
    real(dp)              :: x
    real(dp)              :: Genu_Mual_m     ! Genuchten shape parameter

    Genu_Mual_m = PWP_c - (PWP_c / Genu_Mual_n)
    x = PWP_c + exp(Genu_Mual_n * log(Genu_Mual_alpha * PWP_matPot_ThetaR))
    x = exp(Genu_Mual_m * log(x))
    ! constrain
    if( x < 1.0_dp) x = 1.0_dp
    thetaPWP = thetaS / x

  end subroutine PWP

  ! ----------------------------------------------------------------------------

  !      NAME
  !         field_cap

  !>        \brief calculates the field capacity

  !>        \details estimate Field capacity; FC -- Flux based
  !>        approach (Twarakavi, et. al. 2009, WRR) \n
  !>        According to the 
  !>        above reference FC is defined as the soil water content at 
  !>        which the drainage from a profile ceases under natural
  !>        conditions. Since drainage from a soil profile in a simulation
  !>        never becomes zero, we assume that drainage ceases when the
  !>        bottom flux from the soil reaches a value that is equivalent to
  !>        the minimum amount of precipitation that could be recorded
  !>        (i.e. 0.01 cm/d == 1 mm/d). It is assumed that ThetaR = 0.0_dp

  !      INTENT(IN)
  !>        \param[in] "real(dp) :: Ks"          - saturated hydraulic conductivity
  !>        \param[in] "real(dp) :: thetaS"      - saturated water content
  !>        \param[in] "real(dp) :: Genu_Mual_n" - Genuchten shape parameter

  !      INTENT(OUT)
  !>        \param[out] "real(dp) :: thetaFC" - Field capacity

  !      LITERATURE
  !         Twarakavi, et. al. 2009, WRR

  !      HISTORY
  !>        \author Stephan Thober, Rohini Kumar
  !>        \date Dec 2012
  !         Written, Stephan Thober, Dec 2012
  !         Modified, Matthias Zink,  Nov 2013 - documentation, moved constants to mhm_constants


  elemental pure subroutine field_cap( thetaFC, & ! Output
       Ks, thetaS, Genu_Mual_n )   ! Input

    use mo_mhm_constants , only: field_cap_c1, field_cap_c2

    implicit none

    ! Input
    real(dp), intent(in)  :: Ks          ! saturated hydraulic conductivity
    real(dp), intent(in)  :: thetaS      ! saturated water content
    real(dp), intent(in)  :: Genu_Mual_n ! Genuchten shape parameter

    ! Output
    real(dp), intent(out) :: thetaFC     ! Field capacity

    ! Local variable
    real(dp)              :: x

    x = (field_cap_c1) * (field_cap_c2 + log10( Ks ))
    thetaFC = thetaS * exp( x * log(Genu_Mual_n) )

  end subroutine field_cap

  ! ----------------------------------------------------------------------------

  !      NAME
  !         Genuchten

  !>        \brief calculates the Genuchten shape parameter

  !>        \details estimate SMs_till & van Genuchten's shape parameter (n)
  !>        (Zacharias et al, 2007, soil Phy.)\n
  !>        Global parameters needed (see mhm_parameter.nml):\n
  !>           - param( 1) = PTF_lower66_5_constant  \n
  !>           - param( 2) = PTF_lower66_5_clay      \n
  !>           - param( 3) = PTF_lower66_5_Db        \n
  !>           - param( 4) = PTF_higher66_5_constant \n
  !>           - param( 5) = PTF_higher66_5_clay     \n
  !>           - param( 6) = PTF_higher66_5_Db       \n

  !      INTENT(IN)
  !>        \param[in] "real(dp) :: param(6)"         - given parameters
  !>        \param[in] "real(dp) :: sand"             - [%] sand content
  !>        \param[in] "real(dp) :: clay"             - [%] clay content
  !>        \param[in] "real(dp) :: Db"               - [10^3 kg/m3] bulk density

  !      INTENT(OUT)
  !>        \param[out] "real(dp) :: thetaS"          - saturated water content
  !>        \param[out] "real(dp) :: Genu_Mual_n"     - van Genuchten shape parameter
  !>        \param[out] "real(dp) :: Genu_Mual_alpha" - van Genuchten shape parameter

  !      LITERATURE
  !         Zacharias et al, 2007, soil Phy.

  !      HISTORY
  !>        \author Stephan Thober, Rohini Kumar
  !>        \date Dec 2012
  !         Written, Stephan Thober, Dec 2012

  subroutine Genuchten(thetaS, Genu_Mual_n, Genu_Mual_alpha, & ! Output variables
       param, sand, clay, Db )                                 ! Input variables

    use mo_mhm_constants, only:  vGenuchten_sandtresh,&        ! van Genuchten snad treshold
                                 vGenuchtenN_c1 , &            ! constants for van Genuchten n 
                                 vGenuchtenN_c2 , &            ! constants for van Genuchten n 
                                 vGenuchtenN_c3 , &            ! constants for van Genuchten n 
                                 vGenuchtenN_c4 , &            ! constants for van Genuchten n 
                                 vGenuchtenN_c5 , &            ! constants for van Genuchten n 
                                 vGenuchtenN_c6 , &            ! constants for van Genuchten n 
                                 vGenuchtenN_c7 , &            ! constants for van Genuchten n 
                                 vGenuchtenN_c8 , &            ! constants for van Genuchten n 
                                 vGenuchtenN_c9 , &            ! constants for van Genuchten n 
                                 vGenuchtenN_c10, &            ! constants for van Genuchten n 
                                 vGenuchtenN_c11, &            ! constants for van Genuchten n 
                                 vGenuchtenN_c12, &            ! constants for van Genuchten n 
                                 vGenuchtenN_c13, &            ! constants for van Genuchten n 
                                 vGenuchtenN_c14, &            ! constants for van Genuchten n 
                                 vGenuchtenN_c15, &            ! constants for van Genuchten n 
                                 vGenuchtenN_c16, &            ! constants for van Genuchten n 
                                 vGenuchtenN_c17, &            ! constants for van Genuchten n 
                                 vGenuchtenN_c18               ! constants for van Genuchten n 

    implicit none

    ! Input
    real(dp), dimension(6), intent(in)  :: param           ! parameters
    real(dp),               intent(in)  :: sand            ! sand content
    real(dp),               intent(in)  :: clay            ! clay content
    real(dp),               intent(in)  :: Db              ! [10^3 kg/m3] bulk density

    ! Output
    real(dp),               intent(out) :: thetaS          ! saturated water content
    real(dp),               intent(out) :: Genu_Mual_n     ! van Genuchten shape parameter
    real(dp),               intent(out) :: Genu_Mual_alpha ! van Genuchten shape parameter

    ! Local variables
    real(dp)                            :: x               ! temporal variable

    ! estimate SMs_till & van Genuchten's parameters (alpha and n)
    if ( sand < vGenuchten_sandtresh ) then
       thetaS      =  param(1) + param(2) * clay + param(3) * Db
       Genu_Mual_n =  vGenuchtenN_c1  - vGenuchtenN_c2  * ( sand**(vGenuchtenN_c3) )     + &
                                        vGenuchtenN_c4  * ( clay**(vGenuchtenN_c5) ) 
       x           =  vGenuchtenN_c6  + vGenuchtenN_c7  * sand +   vGenuchtenN_c8 * clay - &
                                        vGenuchtenN_c9  * Db
    else
       thetaS      =  param(4) + param(5) * clay + param(6) * Db
       Genu_Mual_n =  vGenuchtenN_c10 + vGenuchtenN_c11 * (sand**(vGenuchtenN_c12) )    + &
                                        vGenuchtenN_c13 * (clay**(vGenuchtenN_c14) )
       x           = vGenuchtenN_c15  + vGenuchtenN_c16 * sand  + vGenuchtenN_c17 * clay - &
                                        vGenuchtenN_c18 * Db
    end if

    ! Maulaum alpha
    Genu_Mual_alpha = exp(x)

    if (thetaS < 0.0_dp) then
       write(*,*) 'JMJMJM-thetaS-BAD'
    end if

    ! Put constrains on theta_S
    if (thetaS < 0.0_dp) then
       thetaS = 0.0_dp
    end if

    if (Genu_Mual_n < 1.01000_dp) then
       write(*,*) 'JMJMJM-Genu_Mual_n-BAD'
    end if
    if (Genu_Mual_alpha < 0.00001_dp) then
       write(*,*) 'JMJMJM-Genu_Mual_alpha-BAD'
    end if


    ! Put constrains on  Genu_Mual_n
    ! according to (Zacharias et al, 2007, soil Phy.)
    if(Genu_Mual_n     < 1.01000_dp) Genu_Mual_n     = 1.01000_dp
    if(Genu_Mual_alpha < 0.00001_dp) Genu_Mual_alpha = 0.00001_dp
    
  end subroutine Genuchten

  ! ----------------------------------------------------------------------------

  !      NAME
  !         hydro_cond

  !>        \brief calculates the hydraulic conductivity Ks

  !>        \details By default save this value of Ks, particularly for the 
  !>        deeper layers where OM content plays relatively low or no role\n
  !>        Global parameters needed (see mhm_parameter.nml):\n
  !>           - param(1) = PTF_Ks_constant   \n
  !>           - param(2) = PTF_Ks_sand       \n
  !>           - param(3) = PTF_Ks_clay       \n
  !>           - param(4) = PTF_Ks_curveSlope \n

  !      INTENT(IN)
  !>        \param[in] "real(dp) :: param(4)" - given parameters
  !>        \param[in] "real(dp) :: sand"     - [%] sand content
  !>        \param[in] "real(dp) :: clay"     - [%] clay content

  !      INTENT(OUT)
  !>        \param[out] "real(dp) :: Ks"      - hydraulic conductivity

  !      HISTORY
  !>        \author Stephan Thober, Rohini Kumar
  !>        \date Dec 2012
  !         Written, Stephan Thober, Dec 2012
  !         Modified, Matthias Zink,  Nov 2013 - documentation, moved constants to mhm_constants

  subroutine hydro_cond( KS, param, sand, clay )

    use mo_mhm_constants, only: Ks_c 

    implicit none

    ! Input
    real(dp), dimension(4), intent(in)  :: param
    real(dp),               intent(in)  :: sand
    real(dp),               intent(in)  :: clay

    ! Output
    real(dp),               intent(out) :: KS

    ! Local variables
    real(dp)                            :: x ! temporal variable

    ! saturated vertical hydraulic conductivity, Ks (cm/d)
    ! (Cosby et. al.1984.; WRR)
    x  =  param(1) +  param(2) * sand -  param(3) * clay
    Ks =  param(4) * exp(X * log(Ks_c))


    if ( Ks < 1.10_dp ) then
       write(*,*) 'JMJMJM-Ks-BAD'
    end if

    ! minimum value of Ks = 1.1cm/d
    if(Ks < 1.10_dp) Ks = 1.10_dp

  end subroutine hydro_cond

end module mo_mpr_soilmoist
