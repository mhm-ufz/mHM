!> \file mo_mhm_read_config.f90
!> \brief   \copybrief mo_mhm_read_config
!> \details \copydetails mo_mhm_read_config

!> \brief Reading of main model configurations.
!> \details This routine reads the configurations of mHM including, input and
!!       output directories, module usage specification, simulation time periods,
!!       global parameters, ...
!> \authors Matthias Zink
!> \date Dec 2012
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mhm
MODULE mo_mhm_read_config

  USE mo_kind, ONLY : i4, dp
  use mo_message, only: message, error_message

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mhm_read_config ! read main directories

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        mhm_read_config

  !    PURPOSE
  !>       \brief Read main configurations for mHM

  !>       \details The main configurations in mHM are read from three files:
  !>       <ol>
  !>       <li> mhm.nml
  !>       <li> mhm_parameters.nml
  !>       <li> mhm_outputs.nml
  !>       </ol>
  !>       For details please refer to the above mentioned namelist files.

  !    INTENT(IN)
  !>       \param[in] "character(*) :: file_namelist"
  !>       \param[in] "integer :: unamelist"

  !    HISTORY
  !>       \authors Matthias Zink

  !>       \date Dec 2012

  ! Modifications:
  ! Luis Samaniego               Jan 2013 - messages Rohini Kumar
  ! Matthias Cuntz               Jan  2013 - namelist consolidation and positioning
  ! Matthias Zink                Jan  2013 - bug fix, added gaugeinfo reading
  ! Rohini Kumar                 Jun  2013 - added restart flags
  ! R. Kumar & S. Thober         Aug  2013 - code change to incorporate output timestep during
  !                                          writing of the netcdf file
  ! Rohini Kumar                 Aug  2013 - name changed from "inputFormat" to inputFormat_meteo_forcings
  ! Rohini Kumar                 Aug  2013 - added dirSoil_LUT and dirGeology_LUT, and changed
  !                                          in namelist made accordingly
  ! Rohini Kumar                 Aug  2013 - added new namelist for LAI related datasets, and changed in within
  !                                          the code made accordingly
  ! Matthias Zink                Aug  2013 - changed read in for land cover period
  ! Juliane Mai                  Oct  2013 - adding global_parameters_name
  ! Matthias Zink                Nov  2013 - edited documentation and included DEFAULT cases for ptocess Matrix
  ! Stephan Thober               Nov  2013 - added read of directories where latitude longitude fields are located
  ! Matthias Zink                Feb  2014 - added multiple options for PET process
  ! Matthias Zink                Mar  2014 - added inflow from upstream areas and gauge information as namelist
  ! Rohini Kumar                 May  2014 - added options for the model run coordinate system
  ! Stephan Thober               May  2014 - added switch for chunk read in
  ! Stephan Thober               Jun  2014 - added option for switching off mpr
  ! Matthias Cuntz & Juliane Mai Nov  2014 - LAI input from daily, monthly or yearly files
  ! Matthias Zink                Dec  2014 - adopted inflow gauges to ignore headwater cells
  ! Matthias Zink                Mar  2015 - added optional soil moisture read in for calibration
  ! Matthias Cuntz               Jul  2015 - removed adjustl from trim(adjustl()) of Geoparams for PGI compatibilty
  ! Stephan Thober               Aug  2015 - added read_config_routing and read_routing_params from mRM
  ! Oldrich Rakovec              Oct  2015 - added reading of the domain average TWS data
  ! Rohini Kumar                 Mar  2016 - options to handle different soil databases
  ! Stephan Thober               Nov  2016 - moved nProcesses and processMatrix to common variables
  ! Rohini Kumar                 Dec  2016 - option to handle monthly mean gridded fields of LAI
  ! M.Zink & M. Cuneyd Demirel   Mar  2017 - Added Jarvis soil water stress function at SM process(3)
  ! M.C. Demirel & Simon Stisen  Apr  2017 - Added FC dependency on root fraction coefficient (ET) at SM process(3)
  ! Robert Schweppe              Dec  2017 - switched from fractional julian day to integer
  ! Robert Schweppe              Jun  2018 - refactoring and reformatting

  subroutine mhm_read_config(file_namelist, unamelist)

    use mo_common_constants, only : maxNoDomains, nodata_i4
    use mo_common_mHM_mRM_read_config, only : common_check_resolution
    use mo_common_mhm_mrm_variables, only : opti_function, optimize
    use mo_common_variables, only : domainMeta, processMatrix
    use mo_file, only : file_defOutput, udefOutput
    use mo_global_variables, only : &
      L1_twsaObs, L1_etObs, L1_smObs, L1_neutronsObs, &
      evap_coeff, &
      nSoilHorizons_sm_input, outputFlxState, &
      timeStep_model_outputs, &
      output_deflate_level, output_double_precision, output_time_reference, &
      BFI_calc, BFI_obs
    use mo_mpr_constants, only : maxNoSoilHorizons
    use mo_mpr_global_variables, only : nSoilHorizons_mHM
    use mo_nml, only : close_nml, open_nml, position_nml
    use mo_string_utils, only : num2str

    implicit none

    character(*), intent(in) :: file_namelist

    integer, intent(in) :: unamelist

    integer(i4) :: iDomain, domainID

    ! soil moisture input
    character(256), dimension(maxNoDomains) :: dir_soil_moisture

    ! ground albedo neutron input
    character(256), dimension(maxNoDomains) :: dir_neutrons

    ! evapotranspiration input
    character(256), dimension(maxNoDomains) :: dir_evapotranspiration

    ! tws input
    character(256), dimension(maxNoDomains) :: dir_TWS

    integer(i4) :: timeStep_tws_input         ! time step of optional data: tws
    integer(i4) :: timeStep_et_input          ! time step of optional data: et
    integer(i4) :: timeStep_sm_input          ! time step of optional data: sm
    integer(i4) :: timeStep_neutrons_input    ! time step of optional data: neutrons


    ! define namelists
    ! optional data used for optimization
    namelist /optional_data/ &
            dir_soil_moisture, &
            nSoilHorizons_sm_input, &
            dir_neutrons, &
            dir_evapotranspiration, &
            dir_TWS, &
            timeStep_sm_input, &
            timeStep_neutrons_input, &
            timeStep_et_input, &
            timeStep_tws_input
    ! namelist for pan evaporation
    namelist /panEvapo/evap_coeff

    ! name list regarding output
    namelist /NLoutputResults/ &
            output_deflate_level, &
            output_double_precision, &
            output_time_reference, &
            timeStep_model_outputs, &
            outputFlxState
    ! namelist for baseflow index optimzation
    namelist /baseflow_config/ BFI_calc, BFI_obs

    !===============================================================
    !  Read namelist main directories
    !===============================================================
    call open_nml(file_namelist, unamelist, quiet = .true.)

    allocate(L1_twsaObs(domainMeta%nDomains))
    allocate(L1_etObs(domainMeta%nDomains))
    allocate(L1_smObs(domainMeta%nDomains))
    allocate(L1_neutronsObs(domainMeta%nDomains))
    ! observed baseflow indizes
    allocate(BFI_obs(domainMeta%nDomains))
    BFI_obs = -1.0_dp  ! negative value to flag missing values
    BFI_calc = .false.

    !===============================================================
    !  Read namelist of optional input data
    !===============================================================
    ! read optional optional data if necessary
    if (optimize) then
      select case (opti_function)
        case(10 : 13, 28)
          ! soil moisture
          call position_nml('optional_data', unamelist)
          read(unamelist, nml = optional_data)
          do iDomain = 1, domainMeta%nDomains
            domainID = domainMeta%indices(iDomain)
            L1_smObs(iDomain)%dir = dir_Soil_moisture(domainID)
            L1_smObs(iDomain)%timeStepInput = timeStep_sm_input
            L1_smObs(iDomain)%varname = 'sm'
          end do
          if (nSoilHorizons_sm_input .GT. nSoilHorizons_mHM) then
            call error_message('***ERROR: Number of soil horizons representative for input soil moisture exceeded', raise=.false.)
            call error_message('          defined number of soil horizions: ', adjustl(trim(num2str(maxNoSoilHorizons))), '!')
          end if
        case(17)
          ! neutrons
          call position_nml('optional_data', unamelist)
          read(unamelist, nml = optional_data)
          do iDomain = 1, domainMeta%nDomains
            domainID = domainMeta%indices(iDomain)
            L1_neutronsObs(iDomain)%dir = dir_neutrons(domainID)
            L1_neutronsObs(iDomain)%timeStepInput = timeStep_neutrons_input
            L1_neutronsObs(iDomain)%timeStepInput = -1 ! TODO: daily, hard-coded, to be flexibilized
            L1_neutronsObs(iDomain)%varname = 'neutrons'
          end do
        case(27, 29, 30)
          ! evapotranspiration
          call position_nml('optional_data', unamelist)
          read(unamelist, nml = optional_data)
          do iDomain = 1, domainMeta%nDomains
            domainID = domainMeta%indices(iDomain)
            L1_etObs(iDomain)%dir = dir_evapotranspiration(domainID)
            L1_etObs(iDomain)%timeStepInput = timeStep_et_input
            L1_etObs(iDomain)%varname = 'et'
          end do
        case(15)
          ! domain average TWS data
          call position_nml('optional_data', unamelist)
          read(unamelist, nml = optional_data)
          do iDomain = 1, domainMeta%nDomains
            domainID = domainMeta%indices(iDomain)
            L1_twsaObs(iDomain)%dir = dir_TWS(domainID)
            L1_twsaObs(iDomain)%timeStepInput = timeStep_tws_input
            L1_twsaObs(iDomain)%varname = 'twsa'
          end do
        case(33)
          ! evapotranspiration
          call position_nml('optional_data', unamelist)
          read(unamelist, nml = optional_data)
          do iDomain = 1, domainMeta%nDomains
            domainID = domainMeta%indices(iDomain)
            L1_etObs(iDomain)%dir = dir_evapotranspiration(domainID)
            L1_etObs(iDomain)%timeStepInput = timeStep_et_input
            L1_etObs(iDomain)%varname = 'et'
          end do

          ! domain average TWS data
          call position_nml('optional_data', unamelist)
          read(unamelist, nml = optional_data)
          do iDomain = 1, domainMeta%nDomains
            domainID = domainMeta%indices(iDomain)
            L1_twsaObs(iDomain)%dir = dir_TWS(domainID)
            L1_twsaObs(iDomain)%timeStepInput = timeStep_tws_input
            L1_twsaObs(iDomain)%varname = 'twsa'
          end do

        case(34)
          !baseflow index optimization
          call position_nml('baseflow_config', unamelist)
          read(unamelist, nml = baseflow_config)

      end select
    end if

    !===============================================================
    ! Read pan evaporation
    !===============================================================
    ! Evap. coef. for free-water surfaces
    call position_nml('panEvapo', unamelist)
    read(unamelist, nml = panEvapo)

    call common_check_resolution(.true., .false.)

    call close_nml(unamelist)

    !===============================================================
    ! Read output specifications for mHM
    !===============================================================
    call open_nml(file_defOutput, udefOutput, quiet = .true.)
    output_deflate_level = 6
    output_time_reference = 0
    output_double_precision = .true.
    outputFlxState = .FALSE.
    call position_nml('NLoutputResults', udefOutput)
    read(udefOutput, nml = NLoutputResults)
    call close_nml(udefOutput)

    call message('')
    call message('Following output will be written:')
    call message('  NetCDF deflate level: ', adjustl(trim(num2str(output_deflate_level))))
    if ( output_double_precision ) then
      call message('  NetCDF output precision: double')
    else
      call message('  NetCDF output precision: single')
    end if
    select case(output_time_reference)
      case(0)
        call message('    NetCDF output time reference point: start of time interval')
      case(1)
        call message('    NetCDF output time reference point: center of time interval')
      case(2)
        call message('    NetCDF output time reference point: end of time interval')
    end select
    call message('  STATES:')
    if (outputFlxState(1)) then
      call message('    interceptional storage                          (L1_inter) [mm]')
    end if
    if (outputFlxState(2)) then
      call message('    height of snowpack                           (L1_snowpack) [mm]')
    end if
    if (outputFlxState(3)) then
      call message('    soil water content in the single layers     (L1_soilMoist) [mm]')
    end if
    if (outputFlxState(4)) then
      call message('    volumetric soil moisture in the single layers              [mm/mm]')
    end if
    if (outputFlxState(5)) then
      call message('    mean volum. soil moisture averaged over all soil layers    [mm/mm]')
    end if
    if (outputFlxState(6)) then
      call message('    waterdepth in reservoir of sealed areas       (L1_sealSTW) [mm]')
    end if
    if (outputFlxState(7)) then
      call message('    waterdepth in reservoir of unsat. soil zone  (L1_unsatSTW) [mm]')
    end if
    if (outputFlxState(8)) then
      call message('    waterdepth in reservoir of sat. soil zone      (L1_satSTW) [mm]')
    end if
    if (processMatrix(10, 1) .eq. 0) outputFlxState(18) = .false. ! suppress output if process is off
    if (outputFlxState(18)) then
      call message('    ground albedo neutrons                       (L1_neutrons) [cph]')
    end if

    call message('  FLUXES:')
    if (outputFlxState(9)) then
      call message('    potential evapotranspiration PET                  (L1_pet) [mm/T]')
    end if
    if (outputFlxState(10)) then
      call message('    actual evapotranspiration aET               (L1_aETCanopy) [mm/T]')
    end if
    if (outputFlxState(11)) then
      call message('    total discharge generated per cell       (L1_total_runoff) [mm/T]')
    end if
    if (outputFlxState(12)) then
      call message('    direct runoff generated per cell           (L1_runoffSeal) [mm/T]')
    end if
    if (outputFlxState(13)) then
      call message('    fast interflow generated per cell          (L1_fastRunoff) [mm/T]')
    end if
    if (outputFlxState(14)) then
      call message('    slow interflow generated per cell          (L1_slowRunoff) [mm/T]')
    end if
    if (outputFlxState(15)) then
      call message('    baseflow generated per cell                  (L1_baseflow) [mm/T]')
    end if
    if (outputFlxState(16)) then
      call message('    groundwater recharge                           (L1_percol) [mm/T]')
    end if
    if (outputFlxState(17)) then
      call message('    infiltration                                (L1_infilSoil) [mm/T]')
    end if
    if (outputFlxState(19)) then
      call message('    actual evapotranspiration from soil layers    (L1_aETSoil) [mm/T]')
    end if
    if (outputFlxState(20)) then
      call message('    effective precipitation                     (L1_preEffect) [mm/T]')
    end if
    if (outputFlxState(21)) then
      call message('    snow melt                                        (L1_melt) [mm/T]')
    end if
    call message('')
    call message('FINISHED reading config')

    ! warning message
    if (any(outputFlxState) .and. optimize) then
      call message('WARNING: FLUXES and STATES netCDF will be not written since optimization flag is TRUE ')
    end if

  end subroutine mhm_read_config

END MODULE mo_mhm_read_config
