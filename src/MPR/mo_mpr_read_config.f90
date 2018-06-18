!> \file mo_mpr_read_config.f90

!> \brief read mpr config

!> \details This module contains all mpr subroutines related to
!> reading the mpr configuration from file.

!> \authors Stephan Thober
!> \date Aug 2015
!         Modified, Robert Schweppe Dec 2017 - adapted for mpr


module mo_mpr_read_config

  use mo_kind, only : i4, dp

  implicit none

  public :: mpr_read_config

contains

  ! ------------------------------------------------------------------

  !     NAME
  !         read_mpr_config

  !     PURPOSE
  !>        \brief Read the general config of mpr
  !
  !>        \details Depending on the variable mrm_coupling_config, the
  !>        mRM config is either read from mrm.nml and parameters from
  !>        mrm_parameter.nml or copied from mHM.
  !
  !     INTENT(IN)
  !>        \param[in] "logical :: do_message" - flag for writing mHM standard messages
  !
  !     INTENT(INOUT)
  !         None
  !
  !     INTENT(OUT)
  !>        \param[out] "logical :: readLatLon" - flag for reading LatLon file
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !         None
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !       None
  !
  !     LITERATURE
  !       None

  !     HISTORY
  !>        \author Stephan Thober
  !>        \date Aug 2015
  !         Modified,
  !         Sep 2015, Stephan Thober - removed stop condition when routing resolution is smaller than hydrologic resolution
  !         Oct 2015, Stephan Thober - added NLoutputResults namelist, fileLatLon to directories_general namelist,
  !                                    and readLatLon flag
  !         Dec 2017, Robert Schweppe - adapted for MPR
  subroutine mpr_read_config(file_namelist, unamelist, file_namelist_param, unamelist_param)

    use mo_mpr_global_variables, only : &
            inputFormat_gridded_LAI, & ! format of gridded LAI data (nc only)
            timeStep_LAI_input, & ! time step of gridded LAI input
            iFlag_soilDB, &
            tillageDepth, & ! soil horizons info for mHM
            dirgridded_LAI, & ! directory where gridded LAI is located
            nGeoUnits, &
            HorizonDepth_mHM, nSoilHorizons_mHM, &
            fracSealed_cityArea ! land cover information
    use mo_mpr_constants, only : &
            maxGeoUnit, & ! maximum number of allowed geological classes
            maxNoSoilHorizons ! maximum number of allowed soil layers
    use mo_common_constants, only : &
            nColPars, & ! number of properties of the global variables
            maxNoBasins, & ! maximum number of allowed basins
            nodata_dp, &
            eps_dp
    use mo_common_variables, only : &
            nBasins, &
            global_parameters, &
            global_parameters_name, &
            processMatrix ! process configuration
    use mo_append, only : append
    use mo_utils, only : EQ
    use mo_common_functions, only : in_bound
    use mo_nml, only : open_nml, close_nml, position_nml
    use mo_string_utils, only : num2str
    use mo_message, only : message

    implicit none

    character(*), intent(in) :: file_namelist
    integer, intent(in) :: unamelist
    character(*), intent(in) :: file_namelist_param
    integer, intent(in) :: unamelist_param

    integer(i4) :: ii

    ! some dummy arrays for namelist read in (allocatables not allowed in namelists)
    real(dp), dimension(maxNoSoilHorizons) :: soil_Depth             ! depth of the single horizons
    character(256), dimension(maxNoBasins) :: dir_gridded_LAI        ! directory of gridded LAI data
    !                                                                         ! used when timeStep_LAI_input<0

    ! LOCAL variables
    ! PARAMETERS
    ! some dummy arrays for namelist read in (allocatables not allowed in namelists)
    character(256) :: dummy
    real(dp), dimension(5, nColPars) :: dummy_2d_dp ! space holder for routing parameters
    real(dp), dimension(1, nColPars) :: dummy_2d_dp_2 ! space holder for routing parameters
    ! interception
    real(dp), dimension(nColPars) :: canopyInterceptionFactor
    ! snow
    real(dp), dimension(nColPars) :: snowTreshholdTemperature
    real(dp), dimension(nColPars) :: degreeDayFactor_forest
    real(dp), dimension(nColPars) :: degreeDayFactor_impervious
    real(dp), dimension(nColPars) :: degreeDayFactor_pervious
    real(dp), dimension(nColPars) :: increaseDegreeDayFactorByPrecip
    real(dp), dimension(nColPars) :: maxDegreeDayFactor_forest
    real(dp), dimension(nColPars) :: maxDegreeDayFactor_impervious
    real(dp), dimension(nColPars) :: maxDegreeDayFactor_pervious
    ! soilmoisture
    real(dp), dimension(nColPars) :: orgMatterContent_forest
    real(dp), dimension(nColPars) :: orgMatterContent_impervious
    real(dp), dimension(nColPars) :: orgMatterContent_pervious
    real(dp), dimension(nColPars) :: PTF_lower66_5_constant
    real(dp), dimension(nColPars) :: PTF_lower66_5_clay
    real(dp), dimension(nColPars) :: PTF_lower66_5_Db
    real(dp), dimension(nColPars) :: PTF_higher66_5_constant
    real(dp), dimension(nColPars) :: PTF_higher66_5_clay
    real(dp), dimension(nColPars) :: PTF_higher66_5_Db
    real(dp), dimension(nColPars) :: infiltrationShapeFactor
    real(dp), dimension(nColPars) :: PTF_Ks_constant
    real(dp), dimension(nColPars) :: PTF_Ks_sand
    real(dp), dimension(nColPars) :: PTF_Ks_clay
    real(dp), dimension(nColPars) :: PTF_Ks_curveSlope
    real(dp), dimension(nColPars) :: rootFractionCoefficient_forest
    real(dp), dimension(nColPars) :: rootFractionCoefficient_impervious
    real(dp), dimension(nColPars) :: rootFractionCoefficient_pervious
    real(dp), dimension(nColPars) :: jarvis_sm_threshold_c1
    real(dp), dimension(nColPars) :: rootFractionCoefficient_sand
    real(dp), dimension(nColPars) :: rootFractionCoefficient_clay

    ! directRunoff
    real(dp), dimension(nColPars) :: imperviousStorageCapacity
    ! PET0
    real(dp), dimension(nColPars) :: PET_a_forest
    real(dp), dimension(nColPars) :: PET_a_impervious
    real(dp), dimension(nColPars) :: PET_a_pervious
    real(dp), dimension(nColPars) :: PET_b
    real(dp), dimension(nColPars) :: PET_c
    real(dp), dimension(nColPars) :: minCorrectionFactorPET
    real(dp), dimension(nColPars) :: maxCorrectionFactorPET
    real(dp), dimension(nColPars) :: aspectTresholdPET
    real(dp), dimension(nColPars) :: HargreavesSamaniCoeff
    real(dp), dimension(nColPars) :: PriestleyTaylorCoeff
    real(dp), dimension(nColPars) :: PriestleyTaylorLAIcorr
    real(dp), dimension(nColPars) :: canopyheigth_forest
    real(dp), dimension(nColPars) :: canopyheigth_impervious
    real(dp), dimension(nColPars) :: canopyheigth_pervious
    real(dp), dimension(nColPars) :: displacementheight_coeff
    real(dp), dimension(nColPars) :: roughnesslength_momentum_coeff
    real(dp), dimension(nColPars) :: roughnesslength_heat_coeff
    real(dp), dimension(nColPars) :: stomatal_resistance
    ! interflow
    real(dp), dimension(nColPars) :: interflowStorageCapacityFactor
    real(dp), dimension(nColPars) :: interflowRecession_slope
    real(dp), dimension(nColPars) :: fastInterflowRecession_forest
    real(dp), dimension(nColPars) :: slowInterflowRecession_Ks
    real(dp), dimension(nColPars) :: exponentSlowInterflow
    ! percolation
    real(dp), dimension(nColPars) :: rechargeCoefficient
    real(dp), dimension(nColPars) :: rechargeFactor_karstic
    real(dp), dimension(nColPars) :: gain_loss_GWreservoir_karstic
    ! routing moved to mRM
    ! geological parameters
    real(dp), dimension(maxGeoUnit, nColPars) :: GeoParam
    ! neutrons
    real(dp), dimension(nColPars) :: Desilets_N0
    real(dp), dimension(nColPars) :: COSMIC_N0
    real(dp), dimension(nColPars) :: COSMIC_N1
    real(dp), dimension(nColPars) :: COSMIC_N2
    real(dp), dimension(nColPars) :: COSMIC_alpha0
    real(dp), dimension(nColPars) :: COSMIC_alpha1
    real(dp), dimension(nColPars) :: COSMIC_L30
    real(dp), dimension(nColPars) :: COSMIC_L31

    ! namelist directories
    namelist /directories_MPR/ dir_gridded_LAI
    ! namelist soil database
    namelist /soildata/ iFlag_soilDB, tillageDepth, nSoilHorizons_mHM, soil_Depth
    ! namelist for LAI related data
    namelist /LAI_data_information/ inputFormat_gridded_LAI, timeStep_LAI_input
    ! namelist for land cover scenes
    namelist /LCover_MPR/ fracSealed_cityArea

    ! namelist parameters
    namelist /interception1/ canopyInterceptionFactor
    namelist /snow1/snowTreshholdTemperature, degreeDayFactor_forest, degreeDayFactor_impervious, &
            degreeDayFactor_pervious, increaseDegreeDayFactorByPrecip, maxDegreeDayFactor_forest, &
            maxDegreeDayFactor_impervious, maxDegreeDayFactor_pervious
    namelist /soilmoisture1/ orgMatterContent_forest, orgMatterContent_impervious, orgMatterContent_pervious, &
            PTF_lower66_5_constant, PTF_lower66_5_clay, PTF_lower66_5_Db, PTF_higher66_5_constant, &
            PTF_higher66_5_clay, PTF_higher66_5_Db, PTF_Ks_constant, &
            PTF_Ks_sand, PTF_Ks_clay, PTF_Ks_curveSlope, &
            rootFractionCoefficient_forest, rootFractionCoefficient_impervious, &
            rootFractionCoefficient_pervious, infiltrationShapeFactor
    namelist /soilmoisture2/ orgMatterContent_forest, orgMatterContent_impervious, orgMatterContent_pervious, &
            PTF_lower66_5_constant, PTF_lower66_5_clay, PTF_lower66_5_Db, PTF_higher66_5_constant, &
            PTF_higher66_5_clay, PTF_higher66_5_Db, PTF_Ks_constant, &
            PTF_Ks_sand, PTF_Ks_clay, PTF_Ks_curveSlope, &
            rootFractionCoefficient_forest, rootFractionCoefficient_impervious, &
            rootFractionCoefficient_pervious, infiltrationShapeFactor, jarvis_sm_threshold_c1
    namelist /soilmoisture3/ orgMatterContent_forest, orgMatterContent_impervious, orgMatterContent_pervious, &
            PTF_lower66_5_constant, PTF_lower66_5_clay, PTF_lower66_5_Db, PTF_higher66_5_constant, &
            PTF_higher66_5_clay, PTF_higher66_5_Db, PTF_Ks_constant, &
            PTF_Ks_sand, PTF_Ks_clay, PTF_Ks_curveSlope, &
            rootFractionCoefficient_forest, rootFractionCoefficient_impervious, &
            rootFractionCoefficient_pervious, infiltrationShapeFactor, jarvis_sm_threshold_c1, &
            rootFractionCoefficient_sand, rootFractionCoefficient_clay

    namelist /directRunoff1/ imperviousStorageCapacity
    ! PET is input, LAI driven correction
    namelist /PETminus1/  PET_a_forest, PET_a_impervious, PET_a_pervious, PET_b, PET_c
    ! PET is input, aspect driven correction
    namelist /PET0/  minCorrectionFactorPET, maxCorrectionFactorPET, aspectTresholdPET
    ! Hargreaves-Samani
    namelist /PET1/  minCorrectionFactorPET, maxCorrectionFactorPET, aspectTresholdPET, HargreavesSamaniCoeff
    ! Priestely-Taylor
    namelist /PET2/  PriestleyTaylorCoeff, PriestleyTaylorLAIcorr
    ! Penman-Monteith
    namelist /PET3/  canopyheigth_forest, canopyheigth_impervious, canopyheigth_pervious, displacementheight_coeff, &
            roughnesslength_momentum_coeff, roughnesslength_heat_coeff, stomatal_resistance
    namelist /interflow1/ interflowStorageCapacityFactor, interflowRecession_slope, fastInterflowRecession_forest, &
            slowInterflowRecession_Ks, exponentSlowInterflow
    namelist /percolation1/ rechargeCoefficient, rechargeFactor_karstic, gain_loss_GWreservoir_karstic
    namelist /neutrons1/ Desilets_N0, COSMIC_N0, COSMIC_N1, COSMIC_N2, COSMIC_alpha0, COSMIC_alpha1, COSMIC_L30, COSMIC_L31
    !
    namelist /geoparameter/ GeoParam

    !===============================================================
    ! INITIALIZATION
    !===============================================================
    soil_Depth = 0.0_dp
    dummy_2d_dp = nodata_dp
    dummy_2d_dp_2 = nodata_dp

    call open_nml(file_namelist, unamelist, quiet = .true.)

    !===============================================================
    !  Read namelist for LCover
    !===============================================================
    call position_nml('LCover_MPR', unamelist)
    read(unamelist, nml = LCover_MPR)

    !===============================================================
    ! Read soil layering information
    !===============================================================
    call position_nml('soildata', unamelist)
    read(unamelist, nml = soildata)

    allocate(HorizonDepth_mHM(nSoilHorizons_mHM))
    HorizonDepth_mHM(:) = 0.0_dp
    ! last layer is reset to 0 in MPR in case of iFlag_soilDB is 0
    HorizonDepth_mHM(1 : nSoilHorizons_mHM) = soil_Depth(1 : nSoilHorizons_mHM)

    ! counter checks -- soil horizons
    if (nSoilHorizons_mHM .GT. maxNoSoilHorizons) then
      call message()
      call message('***ERROR: Number of soil horizons is resticted to ', trim(num2str(maxNoSoilHorizons)), '!')
      stop
    end if

    ! the default is the HorizonDepths are all set up to last
    ! as is the default for option-1 where horizon specific information are taken into consideration
    if(iFlag_soilDB .eq. 0) then
      ! classical mhm soil database
      HorizonDepth_mHM(nSoilHorizons_mHM) = 0.0_dp
    else if(iFlag_soilDB .ne. 1) then
      call message()
      call message('***ERROR: iFlag_soilDB option given does not exist. Only 0 and 1 is taken at the moment.')
      stop
    end if

    ! some consistency checks for the specification of the tillage depth
    if(iFlag_soilDB .eq. 1) then
      if(count(abs(HorizonDepth_mHM(:) - tillageDepth) .lt. eps_dp)  .eq. 0) then
        call message()
        call message('***ERROR: Soil tillage depth must conform with one of the specified horizon (lower) depth.')
        stop
      end if
    end if

    !===============================================================
    ! Read LAI related information
    !===============================================================
    call position_nml('LAI_data_information', unamelist)
    read(unamelist, nml = LAI_data_information)

    if (timeStep_LAI_input .ne. 0) then
      !===============================================================
      !  Read namelist for main directories
      !===============================================================
      call position_nml('directories_MPR', unamelist)
      read(unamelist, nml = directories_MPR)

      allocate(dirgridded_LAI(nBasins))
      dirgridded_LAI = dir_gridded_LAI(1 : nBasins)

      if (timeStep_LAI_input .GT. 1) then
        call message()
        call message('***ERROR: option for selected timeStep_LAI_input not coded yet')
        stop
      end if
    end if

    call close_nml(unamelist)

    !===============================================================
    ! Read namelist global parameters
    !===============================================================
    call open_nml(file_namelist_param, unamelist_param, quiet = .true.)
    ! decide which parameters to read depending on specified processes

    ! Process 1 - interception
    select case (processMatrix(1, 1))
      ! 1 - maximum Interception
    case(1)
      call position_nml('interception1', unamelist_param)
      read(unamelist_param, nml = interception1)

      processMatrix(1, 2) = 1_i4
      processMatrix(1, 3) = 1_i4
      call append(global_parameters, reshape(canopyInterceptionFactor, (/1, nColPars/)))

      call append(global_parameters_name, (/  &
              'canopyInterceptionFactor'/))

      ! check if parameter are in range
      if (.not. in_bound(global_parameters)) then
        call message('***ERROR: parameter in namelist "interception1" out of bound in ', &
                trim(adjustl(file_namelist_param)))
        stop
      end if

    case DEFAULT
      call message()
      call message('***ERROR: Process description for process "interception" does not exist!')
      stop
    end select

    ! Process 2 - snow
    select case (processMatrix(2, 1))
      ! 1 - degree-day approach
    case(1)
      call position_nml('snow1', unamelist_param)
      read(unamelist_param, nml = snow1)

      processMatrix(2, 2) = 8_i4
      processMatrix(2, 3) = sum(processMatrix(1 : 2, 2))
      call append(global_parameters, reshape(snowTreshholdTemperature, (/1, nColPars/)))
      call append(global_parameters, reshape(degreeDayFactor_forest, (/1, nColPars/)))
      call append(global_parameters, reshape(degreeDayFactor_impervious, (/1, nColPars/)))
      call append(global_parameters, reshape(degreeDayFactor_pervious, (/1, nColPars/)))
      call append(global_parameters, reshape(increaseDegreeDayFactorByPrecip, (/1, nColPars/)))
      call append(global_parameters, reshape(maxDegreeDayFactor_forest, (/1, nColPars/)))
      call append(global_parameters, reshape(maxDegreeDayFactor_impervious, (/1, nColPars/)))
      call append(global_parameters, reshape(maxDegreeDayFactor_pervious, (/1, nColPars/)))

      call append(global_parameters_name, (/  &
              'snowTreshholdTemperature       ', &
                      'degreeDayFactor_forest         ', &
                      'degreeDayFactor_impervious     ', &
                      'degreeDayFactor_pervious       ', &
                      'increaseDegreeDayFactorByPrecip', &
                      'maxDegreeDayFactor_forest      ', &
                      'maxDegreeDayFactor_impervious  ', &
                      'maxDegreeDayFactor_pervious    '/))

      ! check if parameter are in range
      if (.not. in_bound(global_parameters)) then
        call message('***ERROR: parameter in namelist "snow1" out of bound in ', &
                trim(adjustl(file_namelist_param)))
        stop
      end if

    case DEFAULT
      call message()
      call message('***ERROR: Process description for process "snow" does not exist!')
      stop
    end select

    ! Process 3 - soilmoisture
    select case (processMatrix(3, 1))

      ! 1 - Feddes equation for PET reduction, bucket approach, Brooks-Corey like
    case(1)
      call position_nml('soilmoisture1', unamelist_param)
      read(unamelist_param, nml = soilmoisture1)
      processMatrix(3, 2) = 17_i4
      processMatrix(3, 3) = sum(processMatrix(1 : 3, 2))
      call append(global_parameters, reshape(orgMatterContent_forest, (/1, nColPars/)))
      call append(global_parameters, reshape(orgMatterContent_impervious, (/1, nColPars/)))
      call append(global_parameters, reshape(orgMatterContent_pervious, (/1, nColPars/)))
      call append(global_parameters, reshape(PTF_lower66_5_constant, (/1, nColPars/)))
      call append(global_parameters, reshape(PTF_lower66_5_clay, (/1, nColPars/)))
      call append(global_parameters, reshape(PTF_lower66_5_Db, (/1, nColPars/)))
      call append(global_parameters, reshape(PTF_higher66_5_constant, (/1, nColPars/)))
      call append(global_parameters, reshape(PTF_higher66_5_clay, (/1, nColPars/)))
      call append(global_parameters, reshape(PTF_higher66_5_Db, (/1, nColPars/)))
      call append(global_parameters, reshape(PTF_Ks_constant, (/1, nColPars/)))
      call append(global_parameters, reshape(PTF_Ks_sand, (/1, nColPars/)))
      call append(global_parameters, reshape(PTF_Ks_clay, (/1, nColPars/)))
      call append(global_parameters, reshape(PTF_Ks_curveSlope, (/1, nColPars/)))
      call append(global_parameters, reshape(rootFractionCoefficient_forest, (/1, nColPars/)))
      call append(global_parameters, reshape(rootFractionCoefficient_impervious, (/1, nColPars/)))
      call append(global_parameters, reshape(rootFractionCoefficient_pervious, (/1, nColPars/)))
      call append(global_parameters, reshape(infiltrationShapeFactor, (/1, nColPars/)))

      call append(global_parameters_name, (/     &
              'orgMatterContent_forest           ', &
                      'orgMatterContent_impervious       ', &
                      'orgMatterContent_pervious         ', &
                      'PTF_lower66_5_constant            ', &
                      'PTF_lower66_5_clay                ', &
                      'PTF_lower66_5_Db                  ', &
                      'PTF_higher66_5_constant           ', &
                      'PTF_higher66_5_clay               ', &
                      'PTF_higher66_5_Db                 ', &
                      'PTF_Ks_constant                   ', &
                      'PTF_Ks_sand                       ', &
                      'PTF_Ks_clay                       ', &
                      'PTF_Ks_curveSlope                 ', &
                      'rootFractionCoefficient_forest    ', &
                      'rootFractionCoefficient_impervious', &
                      'rootFractionCoefficient_pervious  ', &
                      'infiltrationShapeFactor           '/))

      ! check if parameter are in range
      if (.not. in_bound(global_parameters)) then
        call message('***ERROR: parameter in namelist "soilmoisture1" out of bound in ', &
                trim(adjustl(file_namelist_param)))
        stop
      end if

      ! 2- Jarvis equation for PET reduction, bucket approach, Brooks-Corey like
    case(2)
      call position_nml('soilmoisture2', unamelist_param)
      read(unamelist_param, nml = soilmoisture2)
      processMatrix(3, 2) = 18_i4
      processMatrix(3, 3) = sum(processMatrix(1 : 3, 2))
      call append(global_parameters, reshape(orgMatterContent_forest, (/1, nColPars/)))
      call append(global_parameters, reshape(orgMatterContent_impervious, (/1, nColPars/)))
      call append(global_parameters, reshape(orgMatterContent_pervious, (/1, nColPars/)))
      call append(global_parameters, reshape(PTF_lower66_5_constant, (/1, nColPars/)))
      call append(global_parameters, reshape(PTF_lower66_5_clay, (/1, nColPars/)))
      call append(global_parameters, reshape(PTF_lower66_5_Db, (/1, nColPars/)))
      call append(global_parameters, reshape(PTF_higher66_5_constant, (/1, nColPars/)))
      call append(global_parameters, reshape(PTF_higher66_5_clay, (/1, nColPars/)))
      call append(global_parameters, reshape(PTF_higher66_5_Db, (/1, nColPars/)))
      call append(global_parameters, reshape(PTF_Ks_constant, (/1, nColPars/)))
      call append(global_parameters, reshape(PTF_Ks_sand, (/1, nColPars/)))
      call append(global_parameters, reshape(PTF_Ks_clay, (/1, nColPars/)))
      call append(global_parameters, reshape(PTF_Ks_curveSlope, (/1, nColPars/)))
      call append(global_parameters, reshape(rootFractionCoefficient_forest, (/1, nColPars/)))
      call append(global_parameters, reshape(rootFractionCoefficient_impervious, (/1, nColPars/)))
      call append(global_parameters, reshape(rootFractionCoefficient_pervious, (/1, nColPars/)))
      call append(global_parameters, reshape(infiltrationShapeFactor, (/1, nColPars/)))
      call append(global_parameters, reshape(jarvis_sm_threshold_c1, (/1, nColPars/)))

      call append(global_parameters_name, (/     &
              'orgMatterContent_forest           ', &
                      'orgMatterContent_impervious       ', &
                      'orgMatterContent_pervious         ', &
                      'PTF_lower66_5_constant            ', &
                      'PTF_lower66_5_clay                ', &
                      'PTF_lower66_5_Db                  ', &
                      'PTF_higher66_5_constant           ', &
                      'PTF_higher66_5_clay               ', &
                      'PTF_higher66_5_Db                 ', &
                      'PTF_Ks_constant                   ', &
                      'PTF_Ks_sand                       ', &
                      'PTF_Ks_clay                       ', &
                      'PTF_Ks_curveSlope                 ', &
                      'rootFractionCoefficient_forest    ', &
                      'rootFractionCoefficient_impervious', &
                      'rootFractionCoefficient_pervious  ', &
                      'infiltrationShapeFactor           ', &
                      'jarvis_sm_threshold_c1            '/))


      ! 3- Jarvis equation for ET reduction and FC dependency on root fraction coefficient
    case(3)
      call position_nml('soilmoisture3', unamelist_param)
      read(unamelist_param, nml = soilmoisture3)
      processMatrix(3, 2) = 20_i4
      processMatrix(3, 3) = sum(processMatrix(1 : 3, 2))
      call append(global_parameters, reshape(orgMatterContent_forest, (/1, nColPars/)))
      call append(global_parameters, reshape(orgMatterContent_impervious, (/1, nColPars/)))
      call append(global_parameters, reshape(orgMatterContent_pervious, (/1, nColPars/)))
      call append(global_parameters, reshape(PTF_lower66_5_constant, (/1, nColPars/)))
      call append(global_parameters, reshape(PTF_lower66_5_clay, (/1, nColPars/)))
      call append(global_parameters, reshape(PTF_lower66_5_Db, (/1, nColPars/)))
      call append(global_parameters, reshape(PTF_higher66_5_constant, (/1, nColPars/)))
      call append(global_parameters, reshape(PTF_higher66_5_clay, (/1, nColPars/)))
      call append(global_parameters, reshape(PTF_higher66_5_Db, (/1, nColPars/)))
      call append(global_parameters, reshape(PTF_Ks_constant, (/1, nColPars/)))
      call append(global_parameters, reshape(PTF_Ks_sand, (/1, nColPars/)))
      call append(global_parameters, reshape(PTF_Ks_clay, (/1, nColPars/)))
      call append(global_parameters, reshape(PTF_Ks_curveSlope, (/1, nColPars/)))
      call append(global_parameters, reshape(rootFractionCoefficient_forest, (/1, nColPars/)))
      call append(global_parameters, reshape(rootFractionCoefficient_impervious, (/1, nColPars/)))
      call append(global_parameters, reshape(rootFractionCoefficient_pervious, (/1, nColPars/)))
      call append(global_parameters, reshape(infiltrationShapeFactor, (/1, nColPars/)))
      call append(global_parameters, reshape(rootFractionCoefficient_sand, (/1, nColPars/)))
      call append(global_parameters, reshape(rootFractionCoefficient_clay, (/1, nColPars/)))
      call append(global_parameters, reshape(jarvis_sm_threshold_c1, (/1, nColPars/)))

      call append(global_parameters_name, (/     &
              'orgMatterContent_forest           ', &
                      'orgMatterContent_impervious       ', &
                      'orgMatterContent_pervious         ', &
                      'PTF_lower66_5_constant            ', &
                      'PTF_lower66_5_clay                ', &
                      'PTF_lower66_5_Db                  ', &
                      'PTF_higher66_5_constant           ', &
                      'PTF_higher66_5_clay               ', &
                      'PTF_higher66_5_Db                 ', &
                      'PTF_Ks_constant                   ', &
                      'PTF_Ks_sand                       ', &
                      'PTF_Ks_clay                       ', &
                      'PTF_Ks_curveSlope                 ', &
                      'rootFractionCoefficient_forest    ', &
                      'rootFractionCoefficient_impervious', &
                      'rootFractionCoefficient_pervious  ', &
                      'infiltrationShapeFactor           ', &
                      'rootFractionCoefficient_sand      ', &
                      'rootFractionCoefficient_clay      ', &
                      'jarvis_sm_threshold_c1            '/))

      ! check if parameter are in range
      if (.not. in_bound(global_parameters)) then
        call message('***ERROR: parameter in namelist "soilmoisture1" out of bound in ', &
                trim(adjustl(file_namelist_param)))
        stop
      end if

    case DEFAULT
      call message()
      call message('***ERROR: Process description for process "soilmoisture" does not exist!')
      stop
    end select

    ! Process 4 - sealed area directRunoff
    select case (processMatrix(4, 1))
      ! 1 - bucket exceedance approach
    case(1)
      call position_nml('directRunoff1', unamelist_param)
      read(unamelist_param, nml = directRunoff1)
      processMatrix(4, 2) = 1_i4
      processMatrix(4, 3) = sum(processMatrix(1 : 4, 2))
      call append(global_parameters, reshape(imperviousStorageCapacity, (/1, nColPars/)))

      call append(global_parameters_name, (/'imperviousStorageCapacity'/))

      ! check if parameter are in range
      if (.not. in_bound(global_parameters)) then
        call message('***ERROR: parameter in namelist "directRunoff1" out of bound in ', &
                trim(adjustl(file_namelist_param)))
        stop
      end if

    case DEFAULT
      call message()
      call message('***ERROR: Process description for process "directRunoff" does not exist!')
      stop
    end select

    ! Process 5 - potential evapotranspiration (PET)
    select case (processMatrix(5, 1))
    case(-1) ! 0 - PET is input, correct PET by LAI
      call position_nml('PETminus1', unamelist_param)
      read(unamelist_param, nml = PETminus1)
      processMatrix(5, 2) = 5_i4
      processMatrix(5, 3) = sum(processMatrix(1 : 5, 2))
      call append(global_parameters, reshape(PET_a_forest, (/1, nColPars/)))
      call append(global_parameters, reshape(PET_a_impervious, (/1, nColPars/)))
      call append(global_parameters, reshape(PET_a_pervious, (/1, nColPars/)))
      call append(global_parameters, reshape(PET_b, (/1, nColPars/)))
      call append(global_parameters, reshape(PET_c, (/1, nColPars/)))

      call append(global_parameters_name, (/ &
              'PET_a_forest     ', &
                      'PET_a_impervious ', &
                      'PET_a_pervious   ', &
                      'PET_b            ', &
                      'PET_c            '/))

      ! check if parameter are in range
      if (.not. in_bound(global_parameters)) then
        call message('***ERROR: parameter in namelist "PETminus1" out of bound  n ', &
                trim(adjustl(file_namelist_param)))
        stop 1
      end if

    case(0) ! 0 - PET is input, correct PET by aspect
      call position_nml('PET0', unamelist_param)
      read(unamelist_param, nml = PET0)
      processMatrix(5, 2) = 3_i4
      processMatrix(5, 3) = sum(processMatrix(1 : 5, 2))
      call append(global_parameters, reshape(minCorrectionFactorPET, (/1, nColPars/)))
      call append(global_parameters, reshape(maxCorrectionFactorPET, (/1, nColPars/)))
      call append(global_parameters, reshape(aspectTresholdPET, (/1, nColPars/)))

      call append(global_parameters_name, (/ &
              'minCorrectionFactorPET ', &
                      'maxCorrectionFactorPET ', &
                      'aspectTresholdPET      '/))

      ! check if parameter are in range
      if (.not. in_bound(global_parameters)) then
        call message('***ERROR: parameter in namelist "PET0" out of bound in ', &
                trim(adjustl(file_namelist_param)))
        stop
      end if

    case(1) ! 1 - Hargreaves-Samani method (HarSam) - additional input needed: Tmin, Tmax
      call position_nml('PET1', unamelist_param)
      read(unamelist_param, nml = PET1)
      processMatrix(5, 2) = 4_i4
      processMatrix(5, 3) = sum(processMatrix(1 : 5, 2))
      call append(global_parameters, reshape(minCorrectionFactorPET, (/1, nColPars/)))
      call append(global_parameters, reshape(maxCorrectionFactorPET, (/1, nColPars/)))
      call append(global_parameters, reshape(aspectTresholdPET, (/1, nColPars/)))
      call append(global_parameters, reshape(HargreavesSamaniCoeff, (/1, nColPars/)))
      call append(global_parameters_name, (/ &
              'minCorrectionFactorPET', &
                      'maxCorrectionFactorPET', &
                      'aspectTresholdPET     ', &
                      'HargreavesSamaniCoeff '/))

      ! check if parameter are in range
      if (.not. in_bound(global_parameters)) then
        call message('***ERROR: parameter in namelist "PET1" out of bound in ', &
                trim(adjustl(file_namelist_param)))
        stop
      end if

    case(2) ! 2 - Priestley-Taylor method (PrieTay) - additional input needed: net_rad
      call position_nml('PET2', unamelist_param)
      read(unamelist_param, nml = PET2)
      processMatrix(5, 2) = 2_i4
      processMatrix(5, 3) = sum(processMatrix(1 : 5, 2))
      call append(global_parameters, reshape(PriestleyTaylorCoeff, (/1, nColPars/)))
      call append(global_parameters, reshape(PriestleyTaylorLAIcorr, (/1, nColPars/)))
      call append(global_parameters_name, (/ &
              'PriestleyTaylorCoeff  ', &
                      'PriestleyTaylorLAIcorr'/))

      ! check if parameter are in range
      if (.not. in_bound(global_parameters)) then
        call message('***ERROR: parameter in namelist "PET2" out of bound in ', &
                trim(adjustl(file_namelist_param)))
        stop
      end if

    case(3) ! 3 - Penman-Monteith method - additional input needed: net_rad, abs. vapour pressue, windspeed
      call position_nml('PET3', unamelist_param)
      read(unamelist_param, nml = PET3)
      processMatrix(5, 2) = 7_i4
      processMatrix(5, 3) = sum(processMatrix(1 : 5, 2))

      call append(global_parameters, reshape(canopyheigth_forest, (/1, nColPars/)))
      call append(global_parameters, reshape(canopyheigth_impervious, (/1, nColPars/)))
      call append(global_parameters, reshape(canopyheigth_pervious, (/1, nColPars/)))
      call append(global_parameters, reshape(displacementheight_coeff, (/1, nColPars/)))
      call append(global_parameters, reshape(roughnesslength_momentum_coeff, (/1, nColPars/)))
      call append(global_parameters, reshape(roughnesslength_heat_coeff, (/1, nColPars/)))
      call append(global_parameters, reshape(stomatal_resistance, (/1, nColPars/)))

      call append(global_parameters_name, (/ &
              'canopyheigth_forest           ', &
                      'canopyheigth_impervious       ', &
                      'canopyheigth_pervious         ', &
                      'displacementheight_coeff      ', &
                      'roughnesslength_momentum_coeff', &
                      'roughnesslength_heat_coeff    ', &
                      'stomatal_resistance           '/))

      ! check if parameter are in range
      if (.not. in_bound(global_parameters)) then
        call message('***ERROR: parameter in namelist "PET3" out of bound in ', &
                trim(adjustl(file_namelist_param)))
        stop
      end if

    case DEFAULT
      call message()
      call message('***ERROR: Process description for process "actualET" does not exist!')
      stop
    end select


    ! Process 6 - interflow
    select case (processMatrix(6, 1))
      ! 1 - parallel soil reservoir approach
    case(1)
      call position_nml('interflow1', unamelist_param)
      read(unamelist_param, nml = interflow1)
      processMatrix(6, 2) = 5_i4
      processMatrix(6, 3) = sum(processMatrix(1 : 6, 2))
      call append(global_parameters, reshape(interflowStorageCapacityFactor, (/1, nColPars/)))
      call append(global_parameters, reshape(interflowRecession_slope, (/1, nColPars/)))
      call append(global_parameters, reshape(fastInterflowRecession_forest, (/1, nColPars/)))
      call append(global_parameters, reshape(slowInterflowRecession_Ks, (/1, nColPars/)))
      call append(global_parameters, reshape(exponentSlowInterflow, (/1, nColPars/)))

      call append(global_parameters_name, (/ &
              'interflowStorageCapacityFactor', &
                      'interflowRecession_slope      ', &
                      'fastInterflowRecession_forest ', &
                      'slowInterflowRecession_Ks     ', &
                      'exponentSlowInterflow         '/))

      ! check if parameter are in range
      if (.not. in_bound(global_parameters)) then
        call message('***ERROR: parameter in namelist "interflow1" out of bound in ', &
                trim(adjustl(file_namelist_param)))
        stop
      end if

    case DEFAULT
      call message()
      call message('***ERROR: Process description for process "interflow" does not exist!')
      stop
    end select

    ! Process 7 - percolation
    select case (processMatrix(7, 1))
      ! 1 - GW layer is assumed as bucket
    case(1)
      call position_nml('percolation1', unamelist_param)
      read(unamelist_param, nml = percolation1)
      processMatrix(7, 2) = 3_i4
      processMatrix(7, 3) = sum(processMatrix(1 : 7, 2))
      call append(global_parameters, reshape(rechargeCoefficient, (/1, nColPars/)))
      call append(global_parameters, reshape(rechargeFactor_karstic, (/1, nColPars/)))
      call append(global_parameters, reshape(gain_loss_GWreservoir_karstic, (/1, nColPars/)))

      call append(global_parameters_name, (/ &
              'rechargeCoefficient          ', &
                      'rechargeFactor_karstic       ', &
                      'gain_loss_GWreservoir_karstic'/))

      ! check if parameter are in range
      if (.not. in_bound(global_parameters)) then
        call message('***ERROR: parameter in namelist "percolation1" out of bound in ', &
                trim(adjustl(file_namelist_param)))
        stop
      end if

    case DEFAULT
      call message()
      call message('***ERROR: Process description for process "percolation" does not exist!')
      stop
    end select

    ! Process 8 - routing
    select case (processMatrix(8, 1))
    case(0)
      ! 0 - deactivated
      call message()
      call message('***CAUTION: Routing is deativated! ')

      processMatrix(8, 2) = 0_i4
      processMatrix(8, 3) = sum(processMatrix(1 : 8, 2))
    case(1)
      ! parameter values and names are set in mRM
      ! 1 - Muskingum approach
      processMatrix(8, 2) = 5_i4
      processMatrix(8, 3) = sum(processMatrix(1 : 8, 2))
      call append(global_parameters, dummy_2d_dp)
      call append(global_parameters_name, (/'dummy', 'dummy', 'dummy', 'dummy', 'dummy'/))
    case(2)
      processMatrix(8, 2) = 1_i4
      processMatrix(8, 3) = sum(processMatrix(1 : 8, 2))
      call append(global_parameters, dummy_2d_dp_2)
      call append(global_parameters_name, (/'dummy'/))
    case DEFAULT
      call message()
      call message('***ERROR: Process description for process "routing" does not exist!')
      stop
    end select

    !===============================================================
    ! Geological formations
    !===============================================================
    dummy = dummy // ''   ! only to avoid warning

    ! Process 9 - geoparameter
    select case (processMatrix(9, 1))
    case(1)
      ! read in global parameters (NOT REGIONALIZED, i.e. these are <beta> and not <gamma>) for each geological formation used
      call position_nml('geoparameter', unamelist_param)
      GeoParam = nodata_dp
      read(unamelist_param, nml = geoparameter)

      ! search number of geological parameters
      do ii = 1, size(GeoParam, 1) ! no while loop to avoid risk of endless loop
        if (EQ(GeoParam(ii, 1), nodata_dp)) then
          nGeoUnits = ii - 1
          exit
        end if
      end do

      ! for geology parameters
      processMatrix(9, 2) = nGeoUnits
      processMatrix(9, 3) = sum(processMatrix(1 : 9, 2))

      call append(global_parameters, GeoParam(1 : nGeoUnits, :))

      ! create names
      do ii = 1, nGeoUnits
        dummy = 'GeoParam(' // trim(adjustl(num2str(ii))) // ',:)'
        call append(global_parameters_name, (/ trim(dummy) /))
      end do

      ! check if parameter are in range
      if (.not. in_bound(global_parameters)) then
        call message('***ERROR: parameter in namelist "geoparameter" out of bound in ', &
                trim(adjustl(file_namelist_param)))
        stop
      end if

    case DEFAULT
      call message()
      call message('***ERROR: Process description for process "geoparameter" does not exist!')
      stop
    end select

    ! Process 10 - neutrons
    !   0 - deactivated
    !   1 - inverse N0 based on Desilets et al. 2010
    !   2 - COSMIC forward operator by Shuttlworth et al. 2013
    if (processMatrix(10, 1) .gt. 0) then

      call position_nml('neutrons1', unamelist_param)
      read(unamelist_param, nml = neutrons1)

      processMatrix(10, 2) = 8_i4
      processMatrix(10, 3) = sum(processMatrix(1 : 10, 2))
      call append(global_parameters, reshape(Desilets_N0, (/1, nColPars/)))
      call append(global_parameters, reshape(COSMIC_N0, (/1, nColPars/)))
      call append(global_parameters, reshape(COSMIC_N1, (/1, nColPars/)))
      call append(global_parameters, reshape(COSMIC_N2, (/1, nColPars/)))
      call append(global_parameters, reshape(COSMIC_alpha0, (/1, nColPars/)))
      call append(global_parameters, reshape(COSMIC_alpha1, (/1, nColPars/)))
      call append(global_parameters, reshape(COSMIC_L30, (/1, nColPars/)))
      call append(global_parameters, reshape(COSMIC_L31, (/1, nColPars/)))

      call append(global_parameters_name, (/  &
              'Desilets_N0   ', &
                      'COSMIC_N0     ', &
                      'COSMIC_N1     ', &
                      'COSMIC_N2     ', &
                      'COSMIC_alpha0 ', &
                      'COSMIC_alpha1 ', &
                      'COSMIC_L30    ', &
                      'COSMIC_L31    '/))

      ! check if parameter are in range
      if (.not. in_bound(global_parameters)) then
        call message('***ERROR: parameter in namelist "neutrons1" out of bound in ', &
                trim(adjustl(file_namelist_param)))
        stop
      end if
    else
      call message(' INFO: Process (10, neutrons) is deactivated, so output will be suppressed.')
      ! this is done below, where nml_output is read
      processMatrix(10, 2) = 0_i4
      processMatrix(10, 3) = sum(processMatrix(1 : 10, 2))
    end if

    call close_nml(unamelist_param)

  end subroutine mpr_read_config

end module mo_mpr_read_config
