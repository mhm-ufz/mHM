!> \file mo_mpr_read_config.f90
!> \brief \copybrief mo_mpr_read_config
!> \details \copydetails mo_mpr_read_config

!> \brief read mpr config
!> \details This module contains all mpr subroutines related to reading the mpr configuration from file.
!> \changelog
!! - Robert Schweppe Dec 2017
!!   - adapted for MPR
!! - Robert Schweppe Jun 2018
!!   - refactoring and reformatting
!! - M. Cuneyd Demirel, Simon Stisen Jun 2020
!!   - added Feddes and FC dependency on root fraction coefficient processCase(3) = 4
!! - Rohini Kumar                    Oct 2021
!!   - Added Neutron count module to mHM integrate into develop branch (5.11.2)
!> \authors Stephan Thober
!> \date Aug 2015
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mpr
module mo_mpr_read_config

  use mo_kind, only : i4, dp

  implicit none

  public :: mpr_read_config

contains

  ! ------------------------------------------------------------------

  !    NAME
  !        mpr_read_config

  !    PURPOSE
  !>       \brief Read the general config of mpr

  !>       \details Depending on the variable mrm_coupling_config, the
  !>       mRM config is either read from mrm.nml and parameters from
  !>       mrm_parameter.nml or copied from mHM.

  !    INTENT(IN)
  !>       \param[in] "character(*) :: file_namelist"
  !>       \param[in] "character(*) :: file_namelist_param"

  !    HISTORY
  !>       \authors Stephan Thober

  !>       \date Aug 2015

  ! Modifications:
  ! Stephan Thober  Sep 2015 - removed stop condition when routing resolution is smaller than hydrologic resolution
  ! Stephan Thober  Oct 2015 - added NLoutputResults namelist, fileLatLon to directories_general namelist, and readLatLon flag
  ! Robert Schweppe Dec 2017 - adapted for MPR
  !  Rohini Kumar   Oct 2021 - Added Neutron count module to mHM integrate into develop branch (5.11.2)

  subroutine mpr_read_config(file_namelist, file_namelist_param)

    use mo_namelists, only : &
      nml_directories_mpr, &
      nml_soildata, &
      nml_lai_data_information, &
      nml_lcover_mpr, &
      nml_interception1, &
      nml_snow1, &
      nml_soilmoisture1, &
      nml_soilmoisture2, &
      nml_soilmoisture3, &
      nml_soilmoisture4, &
      nml_directrunoff1, &
      nml_petminus1, &
      nml_pet0, &
      nml_pet1, &
      nml_pet2, &
      nml_pet3, &
      nml_interflow1, &
      nml_percolation1, &
      nml_neutrons1, &
      nml_neutrons2, &
      nml_geoparameter
    use mo_append, only : append
    use mo_common_constants, only : eps_dp, maxNoDomains, nColPars, nodata_dp
    use mo_common_functions, only : in_bound
    use mo_common_variables, only : global_parameters, global_parameters_name, domainMeta, processMatrix
    use mo_message, only : message, error_message
    use mo_message, only : message, error_message
    use mo_mpr_constants, only : maxGeoUnit, &
                                 maxNoSoilHorizons
    use mo_mpr_global_variables, only : HorizonDepth_data, HorizonDepth_mHM, dirgridded_LAI, fracSealed_cityArea, &
                                        iFlag_soilDB, inputFormat_gridded_LAI, nGeoUnits, nSoilHorizons_data, &
                                        nSoilHorizons_mHM, soil_column_depth_weighting, tillageDepth, timeStep_LAI_input
    use mo_string_utils, only : num2str
    use mo_utils, only : EQ

    implicit none

    character(*), intent(in) :: file_namelist

    character(*), intent(in) :: file_namelist_param

    integer(i4) :: ii

    ! depth of the single horizons
    real(dp), dimension(maxNoSoilHorizons) :: soil_Depth

    ! lower depth of the horizons provided by the gridded soil input data
    real(dp), dimension(maxNoSoilHorizons) :: soil_data_depth

    ! directory of gridded LAI data
    ! used when timeStep_LAI_input<0
    character(256), dimension(maxNoDomains) :: dir_gridded_LAI

    character(256) :: dummy

    ! space holder for routing parameters
    real(dp), dimension(5, nColPars) :: dummy_2d_dp

    ! space holder for routing parameters
    real(dp), dimension(1, nColPars) :: dummy_2d_dp_2

    real(dp), dimension(nColPars) :: canopyInterceptionFactor

    real(dp), dimension(nColPars) :: snowTreshholdTemperature
    real(dp), dimension(nColPars) :: degreeDayFactor_forest
    real(dp), dimension(nColPars) :: degreeDayFactor_impervious
    real(dp), dimension(nColPars) :: degreeDayFactor_pervious
    real(dp), dimension(nColPars) :: increaseDegreeDayFactorByPrecip
    real(dp), dimension(nColPars) :: maxDegreeDayFactor_forest
    real(dp), dimension(nColPars) :: maxDegreeDayFactor_impervious
    real(dp), dimension(nColPars) :: maxDegreeDayFactor_pervious

    real(dp), dimension(nColPars) :: orgMatterContent_forest
    real(dp), dimension(nColPars) :: orgMatterContent_impervious
    real(dp), dimension(nColPars) :: orgMatterContent_pervious
    real(dp), dimension(nColPars) :: PTF_lower66_5_constant
    real(dp), dimension(nColPars) :: PTF_lower66_5_clay
    real(dp), dimension(nColPars) :: PTF_lower66_5_Db
    real(dp), dimension(nColPars) :: PTF_higher66_5_constant
    real(dp), dimension(nColPars) :: PTF_higher66_5_clay
    real(dp), dimension(nColPars) :: PTF_higher66_5_Db
    real(dp), dimension(nColPars) :: PTF_Ks_constant
    real(dp), dimension(nColPars) :: PTF_Ks_sand
    real(dp), dimension(nColPars) :: PTF_Ks_clay
    real(dp), dimension(nColPars) :: PTF_Ks_curveSlope
    real(dp), dimension(nColPars) :: rootFractionCoefficient_forest
    real(dp), dimension(nColPars) :: rootFractionCoefficient_impervious
    real(dp), dimension(nColPars) :: rootFractionCoefficient_pervious
    real(dp), dimension(nColPars) :: infiltrationShapeFactor
    real(dp), dimension(nColPars) :: jarvis_sm_threshold_c1

    real(dp), dimension(nColPars) :: FCmin_glob
    real(dp), dimension(nColPars) :: FCdelta_glob
    real(dp), dimension(nColPars) :: rootFractionCoefficient_sand
    real(dp), dimension(nColPars) :: rootFractionCoefficient_clay

    real(dp), dimension(nColPars) :: imperviousStorageCapacity

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

    real(dp), dimension(nColPars) :: interflowStorageCapacityFactor
    real(dp), dimension(nColPars) :: interflowRecession_slope
    real(dp), dimension(nColPars) :: fastInterflowRecession_forest
    real(dp), dimension(nColPars) :: slowInterflowRecession_Ks
    real(dp), dimension(nColPars) :: exponentSlowInterflow

    real(dp), dimension(nColPars) :: rechargeCoefficient
    real(dp), dimension(nColPars) :: rechargeFactor_karstic
    real(dp), dimension(nColPars) :: gain_loss_GWreservoir_karstic

    real(dp), dimension(maxGeoUnit, nColPars) :: GeoParam

    real(dp), dimension(nColPars) :: Desilets_N0
    real(dp), dimension(nColPars) :: Desilets_LW0
    real(dp), dimension(nColPars) :: Desilets_LW1


    real(dp), dimension(nColPars) :: COSMIC_N0
    real(dp), dimension(nColPars) :: COSMIC_N1
    real(dp), dimension(nColPars) :: COSMIC_N2
    real(dp), dimension(nColPars) :: COSMIC_alpha0
    real(dp), dimension(nColPars) :: COSMIC_alpha1
    real(dp), dimension(nColPars) :: COSMIC_L30
    real(dp), dimension(nColPars) :: COSMIC_L31
    real(dp), dimension(nColPars) :: COSMIC_LW0
    real(dp), dimension(nColPars) :: COSMIC_LW1

    integer(i4) :: iDomain, domainID

    !===============================================================
    ! INITIALIZATION
    !===============================================================
    dummy_2d_dp = nodata_dp
    dummy_2d_dp_2 = nodata_dp

    !===============================================================
    !  Read namelist for LCover
    !===============================================================
    call nml_lcover_mpr%read(file_namelist)
    fracSealed_cityArea = nml_lcover_mpr%fracSealed_cityArea

    !===============================================================
    ! Read soil layering information
    !===============================================================
    call nml_soildata%read(file_namelist)
    iFlag_soilDB = nml_soildata%iFlag_soilDB
    soil_column_depth_weighting = nml_soildata%soil_column_depth_weighting
    tillageDepth = nml_soildata%tillageDepth
    nSoilHorizons_mHM = nml_soildata%nSoilHorizons_mHM
    soil_Depth = nml_soildata%soil_Depth
    nSoilHorizons_data = nml_soildata%nSoilHorizons_data
    soil_data_depth = nml_soildata%soil_data_depth

    allocate(HorizonDepth_mHM(nSoilHorizons_mHM))
    HorizonDepth_mHM(:) = 0.0_dp
    ! last layer is reset to 0 in MPR in case of iFlag_soilDB is 0
    HorizonDepth_mHM(1 : nSoilHorizons_mHM) = soil_Depth(1 : nSoilHorizons_mHM)

    ! counter checks -- soil horizons
    if (nSoilHorizons_mHM .GT. maxNoSoilHorizons) &
      call error_message('***ERROR: Number of soil horizons is resticted to ', trim(num2str(maxNoSoilHorizons)), '!')

    ! the default is the HorizonDepths are all set up to last
    ! as is the default for option-1 where horizon specific information are taken into consideration
    if(iFlag_soilDB .eq. 0) then
      ! classical mhm soil database
      HorizonDepth_mHM(nSoilHorizons_mHM) = 0.0_dp
    else if(iFlag_soilDB .ne. 1) then
      call error_message('***ERROR: iFlag_soilDB option given does not exist. Only 0 and 1 is taken at the moment.')
      call error_message('***ERROR: iFlag_soilDB option given does not exist. Only 0 and 1 is taken at the moment.')
    end if

    ! some consistency checks for the specification of the tillage depth
    if(iFlag_soilDB .eq. 1) then
      if(count(abs(HorizonDepth_mHM(:) - tillageDepth) .lt. eps_dp)  .eq. 0) &
        call error_message('***ERROR: Soil tillage depth must conform with one of the specified horizon (lower) depth.')
    end if

    !---------------------------------------------------------------------------
    ! layering of the gridded soil input data (only relevant for iFlag_soilDB = 1)
    ! nSoilHorizons_data = 0 (default) means: the gridded soil input is provided
    ! on exactly the modelled horizons, i.e. no vertical re-layering is done
    !---------------------------------------------------------------------------
    if (iFlag_soilDB .ne. 1) then
      ! the classical soil database carries its own (soil type specific) layering
      nSoilHorizons_data = nSoilHorizons_mHM
      allocate(HorizonDepth_data(nSoilHorizons_data))
      HorizonDepth_data(:) = HorizonDepth_mHM(:)
    else
      if (nSoilHorizons_data .le. 0_i4) then
        nSoilHorizons_data = nSoilHorizons_mHM
        allocate(HorizonDepth_data(nSoilHorizons_data))
        HorizonDepth_data(:) = HorizonDepth_mHM(:)
      else
        if (nSoilHorizons_data .GT. maxNoSoilHorizons) &
          call error_message('***ERROR: Number of soil data horizons is resticted to ', &
                  trim(num2str(maxNoSoilHorizons)), '!')
        allocate(HorizonDepth_data(nSoilHorizons_data))
        HorizonDepth_data(:) = soil_data_depth(1 : nSoilHorizons_data)
        ! depths must be positive and strictly increasing downwards
        do ii = 1, nSoilHorizons_data
          if (HorizonDepth_data(ii) .le. 0.0_dp) &
            call error_message('***ERROR: soil_data_depth(', trim(num2str(ii)), ') must be greater than zero.')
          if (ii .gt. 1_i4) then
            if (HorizonDepth_data(ii) .le. HorizonDepth_data(ii - 1)) &
              call error_message('***ERROR: soil_data_depth must increase monotonically with the horizon number.')
          end if
        end do
        ! the input data must cover the whole modelled soil column
        if (HorizonDepth_data(nSoilHorizons_data) .lt. HorizonDepth_mHM(nSoilHorizons_mHM) - eps_dp) &
          call error_message('***ERROR: The soil input data (', &
                  trim(num2str(HorizonDepth_data(nSoilHorizons_data), '(F9.2)')), &
                  ' mm) do not cover the modelled soil column (', &
                  trim(num2str(HorizonDepth_mHM(nSoilHorizons_mHM), '(F9.2)')), ' mm).')
        call message()
        call message('Soil input data: ', trim(num2str(nSoilHorizons_data)), ' horizon(s) are read and')
        call message('                 aggregated onto ', trim(num2str(nSoilHorizons_mHM)), ' modelled horizon(s).')
      end if
    end if

    !===============================================================
    ! Read LAI related information
    !===============================================================
    call nml_lai_data_information%read(file_namelist)
    inputFormat_gridded_LAI = nml_lai_data_information%inputFormat_gridded_LAI
    timeStep_LAI_input = nml_lai_data_information%timeStep_LAI_input

    if (timeStep_LAI_input .ne. 0) then
      !===============================================================
      !  Read namelist for main directories
      !===============================================================
      call nml_directories_mpr%read(file_namelist)
      dir_gridded_LAI = nml_directories_mpr%dir_gridded_LAI

      allocate(dirgridded_LAI(domainMeta%nDomains))
      do iDomain = 1, domainMeta%nDomains
        domainID = domainMeta%indices(iDomain)
        dirgridded_LAI(iDomain) = dir_gridded_LAI(domainID)
      end do

      if (timeStep_LAI_input .GT. 1) &
        call error_message('***ERROR: option for selected timeStep_LAI_input not coded yet')
    end if

    !===============================================================
    ! Read namelist global parameters
    !===============================================================
    ! decide which parameters to read depending on specified processes

    ! Process 1 - interception
    select case (processMatrix(1, 1))
      ! 1 - maximum Interception
    case(1)
      call nml_interception1%read(file_namelist_param)
      canopyInterceptionFactor = nml_interception1%canopyInterceptionFactor

      processMatrix(1, 2) = 1_i4
      processMatrix(1, 3) = 1_i4
      call append(global_parameters, reshape(canopyInterceptionFactor, (/1, nColPars/)))

      call append(global_parameters_name, (/  &
              'canopyInterceptionFactor'/))

      ! check if parameter are in range
      if (.not. in_bound(global_parameters)) &
        call error_message('***ERROR: parameter in namelist "interception1" out of bound in ', trim(adjustl(file_namelist_param)))

    case DEFAULT
      call error_message('***ERROR: Process description for process "interception" does not exist!')
      stop
    end select

    ! Process 2 - snow
    select case (processMatrix(2, 1))
      ! 1 - degree-day approach
    case(1)
      call nml_snow1%read(file_namelist_param)
      snowTreshholdTemperature = nml_snow1%snowTreshholdTemperature
      degreeDayFactor_forest = nml_snow1%degreeDayFactor_forest
      degreeDayFactor_impervious = nml_snow1%degreeDayFactor_impervious
      degreeDayFactor_pervious = nml_snow1%degreeDayFactor_pervious
      increaseDegreeDayFactorByPrecip = nml_snow1%increaseDegreeDayFactorByPrecip
      maxDegreeDayFactor_forest = nml_snow1%maxDegreeDayFactor_forest
      maxDegreeDayFactor_impervious = nml_snow1%maxDegreeDayFactor_impervious
      maxDegreeDayFactor_pervious = nml_snow1%maxDegreeDayFactor_pervious

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
      if (.not. in_bound(global_parameters)) &
        call error_message('***ERROR: parameter in namelist "snow1" out of bound in ', trim(adjustl(file_namelist_param)))

    case DEFAULT
      call error_message('***ERROR: Process description for process "snow" does not exist!')
    end select

    ! Process 3 - soilmoisture
    select case (processMatrix(3, 1))

      ! 1 - Feddes equation for PET reduction, bucket approach, Brooks-Corey like
    case(1)
      call nml_soilmoisture1%read(file_namelist_param)
      orgMatterContent_forest = nml_soilmoisture1%orgMatterContent_forest
      orgMatterContent_impervious = nml_soilmoisture1%orgMatterContent_impervious
      orgMatterContent_pervious = nml_soilmoisture1%orgMatterContent_pervious
      PTF_lower66_5_constant = nml_soilmoisture1%PTF_lower66_5_constant
      PTF_lower66_5_clay = nml_soilmoisture1%PTF_lower66_5_clay
      PTF_lower66_5_Db = nml_soilmoisture1%PTF_lower66_5_Db
      PTF_higher66_5_constant = nml_soilmoisture1%PTF_higher66_5_constant
      PTF_higher66_5_clay = nml_soilmoisture1%PTF_higher66_5_clay
      PTF_higher66_5_Db = nml_soilmoisture1%PTF_higher66_5_Db
      PTF_Ks_constant = nml_soilmoisture1%PTF_Ks_constant
      PTF_Ks_sand = nml_soilmoisture1%PTF_Ks_sand
      PTF_Ks_clay = nml_soilmoisture1%PTF_Ks_clay
      PTF_Ks_curveSlope = nml_soilmoisture1%PTF_Ks_curveSlope
      rootFractionCoefficient_forest = nml_soilmoisture1%rootFractionCoefficient_forest
      rootFractionCoefficient_impervious = nml_soilmoisture1%rootFractionCoefficient_impervious
      rootFractionCoefficient_pervious = nml_soilmoisture1%rootFractionCoefficient_pervious
      infiltrationShapeFactor = nml_soilmoisture1%infiltrationShapeFactor

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
      if (.not. in_bound(global_parameters)) &
        call error_message('***ERROR: parameter in namelist "soilmoisture1" out of bound in ', trim(adjustl(file_namelist_param)))

      ! 2- Jarvis equation for PET reduction, bucket approach, Brooks-Corey like
    case(2)
      call nml_soilmoisture2%read(file_namelist_param)
      orgMatterContent_forest = nml_soilmoisture2%orgMatterContent_forest
      orgMatterContent_impervious = nml_soilmoisture2%orgMatterContent_impervious
      orgMatterContent_pervious = nml_soilmoisture2%orgMatterContent_pervious
      PTF_lower66_5_constant = nml_soilmoisture2%PTF_lower66_5_constant
      PTF_lower66_5_clay = nml_soilmoisture2%PTF_lower66_5_clay
      PTF_lower66_5_Db = nml_soilmoisture2%PTF_lower66_5_Db
      PTF_higher66_5_constant = nml_soilmoisture2%PTF_higher66_5_constant
      PTF_higher66_5_clay = nml_soilmoisture2%PTF_higher66_5_clay
      PTF_higher66_5_Db = nml_soilmoisture2%PTF_higher66_5_Db
      PTF_Ks_constant = nml_soilmoisture2%PTF_Ks_constant
      PTF_Ks_sand = nml_soilmoisture2%PTF_Ks_sand
      PTF_Ks_clay = nml_soilmoisture2%PTF_Ks_clay
      PTF_Ks_curveSlope = nml_soilmoisture2%PTF_Ks_curveSlope
      rootFractionCoefficient_forest = nml_soilmoisture2%rootFractionCoefficient_forest
      rootFractionCoefficient_impervious = nml_soilmoisture2%rootFractionCoefficient_impervious
      rootFractionCoefficient_pervious = nml_soilmoisture2%rootFractionCoefficient_pervious
      infiltrationShapeFactor = nml_soilmoisture2%infiltrationShapeFactor
      jarvis_sm_threshold_c1 = nml_soilmoisture2%jarvis_sm_threshold_c1

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

      ! check if parameter are in range
      if (.not. in_bound(global_parameters)) &
        call error_message('***ERROR: parameter in namelist "soilmoisture2" out of bound in ', trim(adjustl(file_namelist_param)))

      ! 3- Jarvis equation for ET reduction and FC dependency on root fraction coefficient
    case(3)
      call nml_soilmoisture3%read(file_namelist_param)
      orgMatterContent_forest = nml_soilmoisture3%orgMatterContent_forest
      orgMatterContent_impervious = nml_soilmoisture3%orgMatterContent_impervious
      orgMatterContent_pervious = nml_soilmoisture3%orgMatterContent_pervious
      PTF_lower66_5_constant = nml_soilmoisture3%PTF_lower66_5_constant
      PTF_lower66_5_clay = nml_soilmoisture3%PTF_lower66_5_clay
      PTF_lower66_5_Db = nml_soilmoisture3%PTF_lower66_5_Db
      PTF_higher66_5_constant = nml_soilmoisture3%PTF_higher66_5_constant
      PTF_higher66_5_clay = nml_soilmoisture3%PTF_higher66_5_clay
      PTF_higher66_5_Db = nml_soilmoisture3%PTF_higher66_5_Db
      PTF_Ks_constant = nml_soilmoisture3%PTF_Ks_constant
      PTF_Ks_sand = nml_soilmoisture3%PTF_Ks_sand
      PTF_Ks_clay = nml_soilmoisture3%PTF_Ks_clay
      PTF_Ks_curveSlope = nml_soilmoisture3%PTF_Ks_curveSlope
      rootFractionCoefficient_forest = nml_soilmoisture3%rootFractionCoefficient_forest
      rootFractionCoefficient_impervious = nml_soilmoisture3%rootFractionCoefficient_impervious
      rootFractionCoefficient_pervious = nml_soilmoisture3%rootFractionCoefficient_pervious
      infiltrationShapeFactor = nml_soilmoisture3%infiltrationShapeFactor
      rootFractionCoefficient_sand = nml_soilmoisture3%rootFractionCoefficient_sand
      rootFractionCoefficient_clay = nml_soilmoisture3%rootFractionCoefficient_clay
      FCmin_glob = nml_soilmoisture3%FCmin_glob
      FCdelta_glob = nml_soilmoisture3%FCdelta_glob
      jarvis_sm_threshold_c1 = nml_soilmoisture3%jarvis_sm_threshold_c1

      processMatrix(3, 2) = 22_i4
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
      call append(global_parameters, reshape(FCmin_glob, (/1, nColPars/)))
      call append(global_parameters, reshape(FCdelta_glob, (/1, nColPars/)))
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
              'FCmin_glob                        ', &
              'FCdelta_glob                      ', &
              'jarvis_sm_threshold_c1            '/))

      ! check if parameter are in range
      if (.not. in_bound(global_parameters)) &
        call error_message('***ERROR: parameter in namelist "soilmoisture3" out of bound in ', trim(adjustl(file_namelist_param)))

      ! 4- Feddes equation for ET reduction and FC dependency on root fraction coefficient
    case(4)
      call nml_soilmoisture4%read(file_namelist_param)
      orgMatterContent_forest = nml_soilmoisture4%orgMatterContent_forest
      orgMatterContent_impervious = nml_soilmoisture4%orgMatterContent_impervious
      orgMatterContent_pervious = nml_soilmoisture4%orgMatterContent_pervious
      PTF_lower66_5_constant = nml_soilmoisture4%PTF_lower66_5_constant
      PTF_lower66_5_clay = nml_soilmoisture4%PTF_lower66_5_clay
      PTF_lower66_5_Db = nml_soilmoisture4%PTF_lower66_5_Db
      PTF_higher66_5_constant = nml_soilmoisture4%PTF_higher66_5_constant
      PTF_higher66_5_clay = nml_soilmoisture4%PTF_higher66_5_clay
      PTF_higher66_5_Db = nml_soilmoisture4%PTF_higher66_5_Db
      PTF_Ks_constant = nml_soilmoisture4%PTF_Ks_constant
      PTF_Ks_sand = nml_soilmoisture4%PTF_Ks_sand
      PTF_Ks_clay = nml_soilmoisture4%PTF_Ks_clay
      PTF_Ks_curveSlope = nml_soilmoisture4%PTF_Ks_curveSlope
      rootFractionCoefficient_forest = nml_soilmoisture4%rootFractionCoefficient_forest
      rootFractionCoefficient_impervious = nml_soilmoisture4%rootFractionCoefficient_impervious
      rootFractionCoefficient_pervious = nml_soilmoisture4%rootFractionCoefficient_pervious
      infiltrationShapeFactor = nml_soilmoisture4%infiltrationShapeFactor
      rootFractionCoefficient_sand = nml_soilmoisture4%rootFractionCoefficient_sand
      rootFractionCoefficient_clay = nml_soilmoisture4%rootFractionCoefficient_clay
      FCmin_glob = nml_soilmoisture4%FCmin_glob
      FCdelta_glob = nml_soilmoisture4%FCdelta_glob

      processMatrix(3, 2) = 21_i4
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
      call append(global_parameters, reshape(FCmin_glob, (/1, nColPars/)))
      call append(global_parameters, reshape(FCdelta_glob, (/1, nColPars/)))

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
              'FCmin_glob                        ', &
              'FCdelta_glob                      '/))

      ! check if parameter are in range
      if (.not. in_bound(global_parameters)) &
        call error_message('***ERROR: parameter in namelist "soilmoisture4" out of bound in ', trim(adjustl(file_namelist_param)))

    case DEFAULT
      call error_message('***ERROR: Process description for process "soilmoisture" does not exist!')
      call error_message('***ERROR: Process description for process "soilmoisture" does not exist!')
    end select

    ! Process 4 - sealed area directRunoff
    select case (processMatrix(4, 1))
      ! 1 - bucket exceedance approach
    case(1)
      call nml_directrunoff1%read(file_namelist_param)
      imperviousStorageCapacity = nml_directrunoff1%imperviousStorageCapacity

      processMatrix(4, 2) = 1_i4
      processMatrix(4, 3) = sum(processMatrix(1 : 4, 2))
      call append(global_parameters, reshape(imperviousStorageCapacity, (/1, nColPars/)))

      call append(global_parameters_name, (/'imperviousStorageCapacity'/))

      ! check if parameter are in range
      if (.not. in_bound(global_parameters)) &
        call error_message('***ERROR: parameter in namelist "directRunoff1" out of bound in ', trim(adjustl(file_namelist_param)))

    case DEFAULT
      call error_message('***ERROR: Process description for process "directRunoff" does not exist!')
      call error_message('***ERROR: Process description for process "directRunoff" does not exist!')
    end select

    ! Process 5 - potential evapotranspiration (PET)
    select case (processMatrix(5, 1))
    case(-1) ! 0 - PET is input, correct PET by LAI
      call nml_petminus1%read(file_namelist_param)
      PET_a_forest = nml_petminus1%PET_a_forest
      PET_a_impervious = nml_petminus1%PET_a_impervious
      PET_a_pervious = nml_petminus1%PET_a_pervious
      PET_b = nml_petminus1%PET_b
      PET_c = nml_petminus1%PET_c

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
      if (.not. in_bound(global_parameters)) &
        call error_message('***ERROR: parameter in namelist "PETminus1" out of bound  n ', trim(adjustl(file_namelist_param)))

    case(0) ! 0 - PET is input, correct PET by aspect
      call nml_pet0%read(file_namelist_param)
      minCorrectionFactorPET = nml_pet0%minCorrectionFactorPET
      maxCorrectionFactorPET = nml_pet0%maxCorrectionFactorPET
      aspectTresholdPET = nml_pet0%aspectTresholdPET

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
      if (.not. in_bound(global_parameters)) &
        call error_message('***ERROR: parameter in namelist "PET0" out of bound in ', trim(adjustl(file_namelist_param)))

    case(1) ! 1 - Hargreaves-Samani method (HarSam) - additional input needed: Tmin, Tmax
      call nml_pet1%read(file_namelist_param)
      minCorrectionFactorPET = nml_pet1%minCorrectionFactorPET
      maxCorrectionFactorPET = nml_pet1%maxCorrectionFactorPET
      aspectTresholdPET = nml_pet1%aspectTresholdPET
      HargreavesSamaniCoeff = nml_pet1%HargreavesSamaniCoeff

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
      if (.not. in_bound(global_parameters)) &
        call error_message('***ERROR: parameter in namelist "PET1" out of bound in ', trim(adjustl(file_namelist_param)))

    case(2) ! 2 - Priestley-Taylor method (PrieTay) - additional input needed: net_rad
      call nml_pet2%read(file_namelist_param)
      PriestleyTaylorCoeff = nml_pet2%PriestleyTaylorCoeff
      PriestleyTaylorLAIcorr = nml_pet2%PriestleyTaylorLAIcorr

      processMatrix(5, 2) = 2_i4
      processMatrix(5, 3) = sum(processMatrix(1 : 5, 2))
      call append(global_parameters, reshape(PriestleyTaylorCoeff, (/1, nColPars/)))
      call append(global_parameters, reshape(PriestleyTaylorLAIcorr, (/1, nColPars/)))
      call append(global_parameters_name, (/ &
              'PriestleyTaylorCoeff  ', &
              'PriestleyTaylorLAIcorr'/))

      ! check if parameter are in range
      if (.not. in_bound(global_parameters)) &
        call error_message('***ERROR: parameter in namelist "PET2" out of bound in ', trim(adjustl(file_namelist_param)))

    case(3) ! 3 - Penman-Monteith method - additional input needed: net_rad, abs. vapour pressue, windspeed
      call nml_pet3%read(file_namelist_param)
      canopyheigth_forest = nml_pet3%canopyheigth_forest
      canopyheigth_impervious = nml_pet3%canopyheigth_impervious
      canopyheigth_pervious = nml_pet3%canopyheigth_pervious
      displacementheight_coeff = nml_pet3%displacementheight_coeff
      roughnesslength_momentum_coeff = nml_pet3%roughnesslength_momentum_coeff
      roughnesslength_heat_coeff = nml_pet3%roughnesslength_heat_coeff
      stomatal_resistance = nml_pet3%stomatal_resistance

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
      if (.not. in_bound(global_parameters)) &
        call error_message('***ERROR: parameter in namelist "PET3" out of bound in ', trim(adjustl(file_namelist_param)))

    case DEFAULT
      call error_message('***ERROR: Process description for process "actualET" does not exist!')
      call error_message('***ERROR: Process description for process "actualET" does not exist!')
    end select


    ! Process 6 - interflow
    select case (processMatrix(6, 1))
      ! 1 - parallel soil reservoir approach
    case(1)
      call nml_interflow1%read(file_namelist_param)
      interflowStorageCapacityFactor = nml_interflow1%interflowStorageCapacityFactor
      interflowRecession_slope = nml_interflow1%interflowRecession_slope
      fastInterflowRecession_forest = nml_interflow1%fastInterflowRecession_forest
      slowInterflowRecession_Ks = nml_interflow1%slowInterflowRecession_Ks
      exponentSlowInterflow = nml_interflow1%exponentSlowInterflow

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
      if (.not. in_bound(global_parameters)) &
        call error_message('***ERROR: parameter in namelist "interflow1" out of bound in ', trim(adjustl(file_namelist_param)))

    case DEFAULT
      call error_message('***ERROR: Process description for process "interflow" does not exist!')
      call error_message('***ERROR: Process description for process "interflow" does not exist!')
    end select

    ! Process 7 - percolation
    select case (processMatrix(7, 1))
      ! 1 - GW layer is assumed as bucket
    case(1)
      call nml_percolation1%read(file_namelist_param)
      rechargeCoefficient = nml_percolation1%rechargeCoefficient
      rechargeFactor_karstic = nml_percolation1%rechargeFactor_karstic
      gain_loss_GWreservoir_karstic = nml_percolation1%gain_loss_GWreservoir_karstic

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
      if (.not. in_bound(global_parameters)) &
        call error_message('***ERROR: parameter in namelist "percolation1" out of bound in ', trim(adjustl(file_namelist_param)))

    case DEFAULT
      call error_message('***ERROR: Process description for process "percolation" does not exist!')
      call error_message('***ERROR: Process description for process "percolation" does not exist!')
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
    case(3)
      processMatrix(8, 2) = 1_i4
      processMatrix(8, 3) = sum(processMatrix(1 : 8, 2))
      call append(global_parameters, dummy_2d_dp_2)
      call append(global_parameters_name, (/'dummy'/))
    case DEFAULT
      call error_message('***ERROR: Process description for process "routing" does not exist!')
    end select

    !===============================================================
    ! Geological formations
    !===============================================================
    dummy = dummy // ''   ! only to avoid warning

    ! Process 9 - geoparameter
    select case (processMatrix(9, 1))
    case(1)
      ! read in global parameters (NOT REGIONALIZED, i.e. these are <beta> and not <gamma>) for each geological formation used
      call nml_geoparameter%read(file_namelist_param)
      GeoParam = nml_geoparameter%GeoParam

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
      if (.not. in_bound(global_parameters)) &
        call error_message('***ERROR: parameter in namelist "geoparameter" out of bound in ', trim(adjustl(file_namelist_param)))

    case DEFAULT
      call error_message('***ERROR: Process description for process "geoparameter" does not exist!')
   end select

    !===============================================================
    ! NEUTRON COUNT
    !===============================================================
    ! Process 10 - neutrons
    !   0 - deactivated
    !   1 - inverse N0 based on Desilets et al. 2010
    !   2 - COSMIC forward operator by Shuttlworth et al. 2013
    select case (processMatrix(10, 1))
    case(0)
      ! 0 - deactivated
      call message()
      call message('***SELECTION: Neutron count routine is deativated! ')


    case(1)
      ! 1 - inverse N0 based on Desilets et al. 2010
      call nml_neutrons1%read(file_namelist_param)
      Desilets_N0 = nml_neutrons1%Desilets_N0
      Desilets_LW0 = nml_neutrons1%Desilets_LW0
      Desilets_LW1 = nml_neutrons1%Desilets_LW1

      processMatrix(10,2) = 3_i4
      processMatrix(10,3) = sum(processMatrix(1:10, 2))
      call append(global_parameters, reshape(Desilets_N0,  (/1, nColPars/)))
      call append(global_parameters, reshape(Desilets_LW0, (/1, nColPars/)))
      call append(global_parameters, reshape(Desilets_LW1, (/1, nColPars/)))

       call append(global_parameters_name, (/  &
               'Desilets_N0   ', &
               'Desilets_LW0  ', &
               'Desilets_LW1  '/))

      ! check if parameter are in range
      if (.not. in_bound(global_parameters)) &
        call error_message('***ERROR: parameter in namelist "neutrons1" out of bound in ', trim(adjustl(file_namelist_param)))

    case(2)
      ! 2 - COSMIC version
      call nml_neutrons2%read(file_namelist_param)
      COSMIC_N0 = nml_neutrons2%COSMIC_N0
      COSMIC_N1 = nml_neutrons2%COSMIC_N1
      COSMIC_N2 = nml_neutrons2%COSMIC_N2
      COSMIC_alpha0 = nml_neutrons2%COSMIC_alpha0
      COSMIC_alpha1 = nml_neutrons2%COSMIC_alpha1
      COSMIC_L30 = nml_neutrons2%COSMIC_L30
      COSMIC_L31 = nml_neutrons2%COSMIC_L31
      COSMIC_LW0 = nml_neutrons2%COSMIC_LW0
      COSMIC_LW1 = nml_neutrons2%COSMIC_LW1

      processMatrix(10,2) = 9_i4
      processMatrix(10,3) = sum(processMatrix(1:10, 2))
      call append(global_parameters, reshape(COSMIC_N0,     (/1, nColPars/)))
      call append(global_parameters, reshape(COSMIC_N1,     (/1, nColPars/)))
      call append(global_parameters, reshape(COSMIC_N2,     (/1, nColPars/)))
      call append(global_parameters, reshape(COSMIC_alpha0, (/1, nColPars/)))
      call append(global_parameters, reshape(COSMIC_alpha1, (/1, nColPars/)))
      call append(global_parameters, reshape(COSMIC_L30,    (/1, nColPars/)))
      call append(global_parameters, reshape(COSMIC_L31,    (/1, nColPars/)))
      call append(global_parameters, reshape(COSMIC_LW0,    (/1, nColPars/)))
      call append(global_parameters, reshape(COSMIC_LW1,    (/1, nColPars/)))

      call append(global_parameters_name, (/  &
              'COSMIC_N0     ', &
              'COSMIC_N1     ', &
              'COSMIC_N2     ', &
              'COSMIC_alpha0 ', &
              'COSMIC_alpha1 ', &
              'COSMIC_L30    ', &
              'COSMIC_L31    ', &
              'COSMIC_LW0    ', &
              'COSMIC_LW1    '/))
      ! check if parameter are in range
      if (.not. in_bound(global_parameters)) &
        call error_message('***ERROR: parameter in namelist "neutrons2" out of bound in ', trim(adjustl(file_namelist_param)))

    case DEFAULT
      call error_message('***ERROR: Process description for process "NEUTRON count" does not exist!')
    end select

  end subroutine mpr_read_config

end module mo_mpr_read_config
