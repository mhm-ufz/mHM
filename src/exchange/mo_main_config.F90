!> \file    mo_main_config.f90
!> \brief   \copybrief mo_main_config
!> \details \copydetails mo_main_config

!> \brief   Module for a mHM process container.
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    Aug 2025
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_exchange
#include "logging.h"
module mo_main_config
  use mo_logging
  use mo_kind, only: i4, dp
  use mo_common_constants, only: nprocesses
  use mo_message, only: message, error_message
  use mo_string_utils, only: n2s => num2str
  use nml_helper, only: NML_OK
  ! mhm.nml
  use nml_config_project, only: nml_config_project_t
  use nml_config_processes, only: nml_config_processes_t
  ! parameters.nml
  use nml_interception1, only: nml_interception1_t
  use nml_snow1, only: nml_snow1_t
  use nml_soilmoisture1, only: nml_soilmoisture1_t
  use nml_soilmoisture2, only: nml_soilmoisture2_t
  use nml_soilmoisture3, only: nml_soilmoisture3_t
  use nml_soilmoisture4, only: nml_soilmoisture4_t
  use nml_directrunoff1, only: nml_directrunoff1_t
  use nml_PETm1, only: nml_PETm1_t
  use nml_PETm2, only: nml_PETm2_t
  use nml_PET1, only: nml_PET1_t
  use nml_PET2, only: nml_PET2_t
  use nml_PET3, only: nml_PET3_t
  use nml_interflow1, only: nml_interflow1_t
  use nml_percolation1, only: nml_percolation1_t
  use nml_routing1, only: nml_routing1_t
  use nml_routing2, only: nml_routing2_t
  use nml_routing3, only: nml_routing3_t
  use nml_neutrons1, only: nml_neutrons1_t
  use nml_neutrons2, only: nml_neutrons2_t
  use nml_rivertemp1, only: nml_rivertemp1_t
  use nml_geoparameter, only: nml_geoparameter_t

  !> \class   parameter_config_t
  !> \brief   Configuration for all parameters.
  type, public :: parameter_config_t
    type(nml_config_processes_t) :: processes !< configuration of the processes
    type(nml_interception1_t) :: interception1 !< interception1 configuration
    type(nml_snow1_t) :: snow1 !< snow1 configuration
    type(nml_soilmoisture1_t) :: soilmoisture1 !< soilmoisture1 configuration
    type(nml_soilmoisture2_t) :: soilmoisture2 !< soilmoisture2 configuration
    type(nml_soilmoisture3_t) :: soilmoisture3 !< soilmoisture3 configuration
    type(nml_soilmoisture4_t) :: soilmoisture4 !< soilmoisture4 configuration
    type(nml_directrunoff1_t) :: directrunoff1 !< directrunoff1 configuration
    type(nml_PETm1_t) :: PETm1 !< PETm1 configuration
    type(nml_PETm2_t) :: PETm2 !< PETm2 configuration
    type(nml_PET1_t) :: PET1 !< PET1 configuration
    type(nml_PET2_t) :: PET2 !< PET2 configuration
    type(nml_PET3_t) :: PET3 !< PET3 configuration
    type(nml_interflow1_t) :: interflow1 !< interflow1 configuration
    type(nml_percolation1_t) :: percolation1 !< percolation1 configuration
    type(nml_routing1_t) :: routing1 !< routing1 configuration
    type(nml_routing2_t) :: routing2 !< routing2 configuration
    type(nml_routing3_t) :: routing3 !< routing3 configuration
    type(nml_neutrons1_t) :: neutrons1 !< neutrons1 configuration
    type(nml_neutrons2_t) :: neutrons2 !< neutrons2 configuration
    type(nml_rivertemp1_t) :: rivertemp1 !< rivertemp1 configuration
    type(nml_geoparameter_t) :: geoparameter !< geoparameter configuration
  contains
    procedure :: set_dims => parameter_config_set_dims
    procedure :: read_parameter => parameter_config_read_parameter
  end type parameter_config_t

  !> \class   parameters_t
  !> \brief   Class for parameters and processes configuration.
  type, public :: parameters_t
    type(parameter_config_t) :: config !< configuration of the parameters
    !> Info about which process runs in which option and number of parameters necessary for this option
    !! - col1: process switch
    !! - col2: number of parameters for process case
    !! - col3: cumulated number of parameters in the array so far
    integer(i4), dimension(nprocesses, 3) :: process_matrix
    !> Matrix of global parameters definition (former: gamma)
    !! - col1: min
    !! - col2: max
    !! - col3: initial
    !! - col4: flag
    !! - col5: scaling
    real(dp), dimension(:, :), allocatable :: definition
    !> Array of global parameters names
    character(256), dimension(:), allocatable :: names
    !> Array of current parameter set
    real(dp), dimension(:), allocatable :: values
    integer(i4) :: nGeoUnits = 25_i4 !< Number of geological formations
  contains
    procedure :: set_dims => parameters_set_dims
    procedure :: configure => parameters_configure
    procedure :: initialize => parameters_initialize
    procedure :: set => parameters_set
    procedure :: get => parameters_get
    procedure :: get_process => parameters_get_process
    procedure :: print => parameters_print
    procedure :: mhm_active => parameters_mhm_active
    procedure :: meteo_active => parameters_meteo_active
    procedure :: mrm_active => parameters_mrm_active
  end type parameters_t

contains

  !> \brief Set runtime dimensions for generated parameter namelists.
  subroutine parameter_config_set_dims(self, n_geo_units)
    class(parameter_config_t), intent(inout) :: self
    integer(i4), intent(in) :: n_geo_units !< number of geological units in the global parameter set
    character(1024) :: errmsg
    integer :: status

    if (.not.self%geoparameter%is_configured) then
      status = self%geoparameter%set_dims(n_geo_units=n_geo_units, errmsg=errmsg)
      if (status /= NML_OK) then
        log_fatal(*) "Error setting geoparameter dimensions: ", trim(errmsg)
        error stop 1
      end if
    end if
  end subroutine parameter_config_set_dims

  !> \brief Set runtime dimensions for generated parameter namelists.
  subroutine parameters_set_dims(self, n_geo_units)
    class(parameters_t), intent(inout) :: self
    integer(i4), intent(in) :: n_geo_units !< number of geological units in the global parameter set

    self%nGeoUnits = n_geo_units
    call self%config%set_dims(n_geo_units)
  end subroutine parameters_set_dims

  !> \brief Check if mHM processes are active.
  logical function parameters_mhm_active(self)
    class(parameters_t), intent(in) :: self
    parameters_mhm_active = &
      (self%process_matrix(1, 1) /= 0_i4) .or. &
      (self%process_matrix(2, 1) /= 0_i4) .or. &
      (self%process_matrix(3, 1) /= 0_i4) .or. &
      (self%process_matrix(4, 1) /= 0_i4) .or. &
      (self%process_matrix(5, 1) /= 0_i4) .or. &
      (self%process_matrix(6, 1) /= 0_i4) .or. &
      (self%process_matrix(7, 1) /= 0_i4) .or. &
      (self%process_matrix(9, 1) /= 0_i4) .or. &
      (self%process_matrix(10, 1) /= 0_i4)
  end function parameters_mhm_active

  !> \brief Check if meteo processes are active.
  logical function parameters_meteo_active(self)
    class(parameters_t), intent(in) :: self
    parameters_meteo_active = self%mhm_active() .or. (self%process_matrix(4, 1) /= 0_i4) .or. (self%process_matrix(11, 1) /= 0_i4)
  end function parameters_meteo_active

  !> \brief Check if mRM processes are active.
  logical function parameters_mrm_active(self)
    class(parameters_t), intent(in) :: self
    parameters_mrm_active = (self%process_matrix(8, 1) /= 0_i4) .or. (self%process_matrix(11, 1) /= 0_i4)
  end function parameters_mrm_active

  !> \brief Read parameter configuration.
  subroutine parameter_config_read_parameter(self, file, processes)
    class(parameter_config_t), intent(inout) :: self
    character(*), intent(in) :: file !< file containing the parameter namelists
    type(nml_config_processes_t), intent(in) :: processes !< configured process switches
    character(1024) :: errmsg
    integer :: status
    log_info(*) "Read config parameter: ", file

    select case (processes%interception)
      case (1_i4)
        status = self%interception1%from_file(file=file, errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Error reading interception1 from: ", file, ", with error: ", trim(errmsg)
          error stop 1
        end if
    end select

    select case (processes%snow)
      case (1_i4)
        status = self%snow1%from_file(file=file, errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Error reading snow1 from: ", file, ", with error: ", trim(errmsg)
          error stop 1
        end if
    end select

    select case (processes%soil_moisture)
      case (1_i4)
        status = self%soilmoisture1%from_file(file=file, errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Error reading soilmoisture1 from: ", file, ", with error: ", trim(errmsg)
          error stop 1
        end if
      case (2_i4)
        status = self%soilmoisture2%from_file(file=file, errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Error reading soilmoisture2 from: ", file, ", with error: ", trim(errmsg)
          error stop 1
        end if
      case (3_i4)
        status = self%soilmoisture3%from_file(file=file, errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Error reading soilmoisture3 from: ", file, ", with error: ", trim(errmsg)
          error stop 1
        end if
      case (4_i4)
        status = self%soilmoisture4%from_file(file=file, errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Error reading soilmoisture4 from: ", file, ", with error: ", trim(errmsg)
          error stop 1
        end if
    end select

    select case (processes%direct_runoff)
      case (1_i4)
        status = self%directrunoff1%from_file(file=file, errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Error reading directrunoff1 from: ", file, ", with error: ", trim(errmsg)
          error stop 1
        end if
    end select

    select case (processes%pet)
      case (-2_i4)
        status = self%PETm2%from_file(file=file, errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Error reading PETm2 from: ", file, ", with error: ", trim(errmsg)
          error stop 1
        end if
      case (-1_i4)
        status = self%PETm1%from_file(file=file, errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Error reading PETm1 from: ", file, ", with error: ", trim(errmsg)
          error stop 1
        end if
      case (1_i4)
        status = self%PET1%from_file(file=file, errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Error reading PET1 from: ", file, ", with error: ", trim(errmsg)
          error stop 1
        end if
      case (2_i4)
        status = self%PET2%from_file(file=file, errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Error reading PET2 from: ", file, ", with error: ", trim(errmsg)
          error stop 1
        end if
      case (3_i4)
        status = self%PET3%from_file(file=file, errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Error reading PET3 from: ", file, ", with error: ", trim(errmsg)
          error stop 1
        end if
    end select

    select case (processes%interflow)
      case (1_i4)
        status = self%interflow1%from_file(file=file, errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Error reading interflow1 from: ", file, ", with error: ", trim(errmsg)
          error stop 1
        end if
    end select

    select case (processes%percolation)
      case (1_i4)
        status = self%percolation1%from_file(file=file, errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Error reading percolation1 from: ", file, ", with error: ", trim(errmsg)
          error stop 1
        end if
    end select

    select case (processes%baseflow)
      case (1_i4)
        status = self%geoparameter%from_file(file=file, errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Error reading geoparameter from: ", file, ", with error: ", trim(errmsg)
          error stop 1
        end if
    end select

    select case (processes%neutrons)
      case (1_i4)
        status = self%neutrons1%from_file(file=file, errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Error reading neutrons1 from: ", file, ", with error: ", trim(errmsg)
          error stop 1
        end if
      case (2_i4)
        status = self%neutrons2%from_file(file=file, errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Error reading neutrons2 from: ", file, ", with error: ", trim(errmsg)
          error stop 1
        end if
    end select

    select case (processes%routing)
      case (1_i4)
        status = self%routing1%from_file(file=file, errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Error reading routing1 from: ", file, ", with error: ", trim(errmsg)
          error stop 1
        end if
      case (2_i4)
        status = self%routing2%from_file(file=file, errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Error reading routing2 from: ", file, ", with error: ", trim(errmsg)
          error stop 1
        end if
      case (3_i4)
        status = self%routing3%from_file(file=file, errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Error reading routing3 from: ", file, ", with error: ", trim(errmsg)
          error stop 1
        end if
    end select

    select case (processes%temperature_routing)
      case (1_i4)
        status = self%rivertemp1%from_file(file=file, errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Error reading rivertemp1 from: ", file, ", with error: ", trim(errmsg)
          error stop 1
        end if
    end select

  end subroutine parameter_config_read_parameter

  !> \brief Configure process switches and process parameters.
  subroutine parameters_configure(self, main_file, para_file, root)
    use mo_os, only: path_join, path_normpath
    class(parameters_t), intent(inout) :: self
    character(*), intent(in), optional :: main_file !< file containing the processes namelists
    character(*), intent(in), optional :: para_file !< file containing the parameter namelists
    character(*), intent(in) :: root !< root directory for relative file paths
    character(1024) :: errmsg
    character(:), allocatable :: path
    integer :: status
    log_info(*) "Configure parameters"

    if (present(main_file)) then
      path = path_normpath(path_join(root, main_file))
      log_info(*) "Read config processes: ", path
      status = self%config%processes%from_file(file=path, errmsg=errmsg)
      if (status /= NML_OK) then
        log_fatal(*) "Error reading processes from: ", path, ", with error: ", trim(errmsg)
        error stop 1
      end if
    end if
    if (.not.self%config%processes%is_configured) then
      log_fatal(*) "Processes configuration not set."
      error stop 1
    end if
    status = self%config%processes%is_valid(errmsg=errmsg)
    if (status /= NML_OK) then
      log_fatal(*) "Processes configuration not valid: ", trim(errmsg)
      error stop
    end if

    if (present(para_file)) then
      path = path_normpath(path_join(root, para_file))
      call self%config%read_parameter(file=path, processes=self%config%processes)
    end if
    select case (self%config%processes%interception)
      case (1_i4)
        if (.not.self%config%interception1%is_configured) then
          log_fatal(*) "Interception parameters not set."
          error stop 1
        end if
        status = self%config%interception1%is_valid(errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Interception1 configuration not valid: ", trim(errmsg)
          error stop 1
        end if
    end select

    select case (self%config%processes%snow)
      case (1_i4)
        if (.not.self%config%snow1%is_configured) then
          log_fatal(*) "Snow parameters not set."
          error stop 1
        end if
        status = self%config%snow1%is_valid(errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Snow1 configuration not valid: ", trim(errmsg)
          error stop 1
        end if
    end select

    select case (self%config%processes%soil_moisture)
      case (1_i4)
        if (.not.self%config%soilmoisture1%is_configured) then
          log_fatal(*) "Soilmoisture parameters not set."
          error stop 1
        end if
        status = self%config%soilmoisture1%is_valid(errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Soilmoisture1 configuration not valid: ", trim(errmsg)
          error stop 1
        end if
      case (2_i4)
        if (.not.self%config%soilmoisture2%is_configured) then
          log_fatal(*) "Soilmoisture parameters not set."
          error stop 1
        end if
        status = self%config%soilmoisture2%is_valid(errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Soilmoisture2 configuration not valid: ", trim(errmsg)
          error stop 1
        end if
      case (3_i4)
        if (.not.self%config%soilmoisture3%is_configured) then
          log_fatal(*) "Soilmoisture parameters not set."
          error stop 1
        end if
        status = self%config%soilmoisture3%is_valid(errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Soilmoisture3 configuration not valid: ", trim(errmsg)
          error stop 1
        end if
      case (4_i4)
        if (.not.self%config%soilmoisture4%is_configured) then
          log_fatal(*) "Soilmoisture parameters not set."
          error stop 1
        end if
        status = self%config%soilmoisture4%is_valid(errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Soilmoisture4 configuration not valid: ", trim(errmsg)
          error stop 1
        end if
    end select

    select case (self%config%processes%direct_runoff)
      case (1_i4)
        if (.not.self%config%directrunoff1%is_configured) then
          log_fatal(*) "Direct runoff parameters not set."
          error stop 1
        end if
        status = self%config%directrunoff1%is_valid(errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Directrunoff1 configuration not valid: ", trim(errmsg)
          error stop 1
        end if
    end select

    select case (self%config%processes%pet)
      case (-2_i4)
        if (.not.self%config%PETm2%is_configured) then
          log_fatal(*) "PET parameters not set."
          error stop 1
        end if
        status = self%config%PETm2%is_valid(errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "PETm2 configuration not valid: ", trim(errmsg)
          error stop 1
        end if
      case (-1_i4)
        if (.not.self%config%PETm1%is_configured) then
          log_fatal(*) "PET parameters not set."
          error stop 1
        end if
        status = self%config%PETm1%is_valid(errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "PETm1 configuration not valid: ", trim(errmsg)
          error stop 1
        end if
      case (1_i4)
        if (.not.self%config%PET1%is_configured) then
          log_fatal(*) "PET parameters not set."
          error stop 1
        end if
        status = self%config%PET1%is_valid(errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "PET1 configuration not valid: ", trim(errmsg)
          error stop 1
        end if
      case (2_i4)
        if (.not.self%config%PET2%is_configured) then
          log_fatal(*) "PET parameters not set."
          error stop 1
        end if
        status = self%config%PET2%is_valid(errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "PET2 configuration not valid: ", trim(errmsg)
          error stop 1
        end if
      case (3_i4)
        if (.not.self%config%PET3%is_configured) then
          log_fatal(*) "PET parameters not set."
          error stop 1
        end if
        status = self%config%PET3%is_valid(errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "PET3 configuration not valid: ", trim(errmsg)
          error stop 1
        end if
    end select

    select case (self%config%processes%interflow)
      case (1_i4)
        if (.not.self%config%interflow1%is_configured) then
          log_fatal(*) "Interflow parameters not set."
          error stop 1
        end if
        status = self%config%interflow1%is_valid(errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Interflow1 configuration not valid: ", trim(errmsg)
          error stop 1
        end if
    end select

    select case (self%config%processes%percolation)
      case (1_i4)
        if (.not.self%config%percolation1%is_configured) then
          log_fatal(*) "Percolation parameters not set."
          error stop 1
        end if
        status = self%config%percolation1%is_valid(errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Percolation1 configuration not valid: ", trim(errmsg)
          error stop 1
        end if
    end select

    select case (self%config%processes%baseflow)
      case (1_i4)
        if (.not.self%config%geoparameter%is_configured) then
          log_fatal(*) "Geoparameter parameters not set."
          error stop 1
        end if
        status = self%config%geoparameter%is_valid(errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Geoparameter configuration not valid: ", trim(errmsg)
          error stop 1
        end if
    end select

    select case (self%config%processes%neutrons)
      case (1_i4)
        if (.not.self%config%neutrons1%is_configured) then
          log_fatal(*) "Neutrons parameters not set."
          error stop 1
        end if
        status = self%config%neutrons1%is_valid(errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Neutrons1 configuration not valid: ", trim(errmsg)
          error stop 1
        end if
      case (2_i4)
        if (.not.self%config%neutrons2%is_configured) then
          log_fatal(*) "Neutrons parameters not set."
          error stop 1
        end if
        status = self%config%neutrons2%is_valid(errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Neutrons2 configuration not valid: ", trim(errmsg)
          error stop 1
        end if
    end select

    select case (self%config%processes%routing)
      case (1_i4)
        if (.not.self%config%routing1%is_configured) then
          log_fatal(*) "Routing parameters not set."
          error stop 1
        end if
        status = self%config%routing1%is_valid(errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Routing1 configuration not valid: ", trim(errmsg)
          error stop 1
        end if
      case (2_i4)
        if (.not.self%config%routing2%is_configured) then
          log_fatal(*) "Routing parameters not set."
          error stop 1
        end if
        status = self%config%routing2%is_valid(errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Routing2 configuration not valid: ", trim(errmsg)
          error stop 1
        end if
      case (3_i4)
        if (.not.self%config%routing3%is_configured) then
          log_fatal(*) "Routing parameters not set."
          error stop 1
        end if
        status = self%config%routing3%is_valid(errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Routing3 configuration not valid: ", trim(errmsg)
          error stop 1
        end if
    end select

    select case (self%config%processes%temperature_routing)
      case (1_i4)
        if (.not.self%config%rivertemp1%is_configured) then
          log_fatal(*) "River temperature parameters not set."
          error stop 1
        end if
        status = self%config%rivertemp1%is_valid(errmsg=errmsg)
        if (status /= NML_OK) then
          log_fatal(*) "Rivertemp1 configuration not valid: ", trim(errmsg)
          error stop 1
        end if
    end select

    call self%initialize()
  end subroutine parameters_configure

  !> \brief Configure parameters
  !> \authors Stephan Thober
  !> \date Aug 2015
  ! Modifications:
  ! Stephan Thober  Sep 2015 - removed stop condition when routing resolution is smaller than hydrologic resolution
  ! Stephan Thober  Oct 2015 - added NLoutputResults namelist, fileLatLon to directories_general namelist, and readLatLon flag
  ! Robert Schweppe Dec 2017 - adapted for MPR
  !  Rohini Kumar   Oct 2021 - Added Neutron count module to mHM integrate into develop branch (5.11.2)
  subroutine parameters_initialize(self)

    use mo_append, only : append
    use mo_common_constants, only : nColPars
    use mo_utils, only : EQ

    implicit none
    class(parameters_t), intent(inout) :: self

    character(256) :: geo_name
    real(dp), allocatable, dimension(:, :) :: GeoParam
    integer(i4) :: ii, shp(2)
    character(1024) :: errmsg
    integer :: status

    shp = [1_i4, ncolpars] ! shape of a parameter matrix slice

    !===============================================================
    ! Read namelist global parameters
    !===============================================================
    ! decide which parameters to read depending on specified processes

    log_info(*) "Set processes and parameters"
    self%process_matrix = 0_i4
    self%process_matrix(1, 1) = self%config%processes%interception
    self%process_matrix(2, 1) = self%config%processes%snow
    self%process_matrix(3, 1) = self%config%processes%soil_moisture
    self%process_matrix(4, 1) = self%config%processes%direct_runoff
    self%process_matrix(5, 1) = self%config%processes%pet
    self%process_matrix(6, 1) = self%config%processes%interflow
    self%process_matrix(7, 1) = self%config%processes%percolation
    self%process_matrix(8, 1) = self%config%processes%routing
    self%process_matrix(9, 1) = self%config%processes%baseflow
    self%process_matrix(10, 1) = self%config%processes%neutrons
    self%process_matrix(11, 1) = self%config%processes%temperature_routing

    ! Process 1 - interception
    select case (self%process_matrix(1, 1))
      case(-1)
        ! -1 - pass through precipitation as throughfall
        log_debug(*) "SELECTION: interception pass-through is activated!"
        self%process_matrix(1, 3) = 0_i4
      case(0)
        ! 0 - no interception process selected
        log_debug(*) "SELECTION: interception process is deactivated!"
        self%process_matrix(1, 3) = 0_i4
      case(1)
        ! 1 - maximum Interception
        self%process_matrix(1, 2) = 1_i4
        self%process_matrix(1, 3) = 1_i4
        call append(self%definition, reshape(self%config%interception1%canopyInterceptionFactor, shp))

        call append(self%names, [  &
                'canopyInterceptionFactor'])

        ! check if parameter are in range
        if (.not. in_bound(self%definition)) then
          log_fatal(*) 'Parameter in namelist "interception1" out of bound.'
          error stop 1
        end if
    end select

    ! Process 2 - snow
    select case (self%process_matrix(2, 1))
      case(-1)
        ! -1 - pass through throughfall as effective precipitation
        log_debug(*) "SELECTION: snow pass-through is activated!"
        self%process_matrix(2, 3) = sum(self%process_matrix(1 : 2, 2))
      case(0)
        ! 0 - no snow process selected
        log_debug(*) "SELECTION: snow process is deactivated!"
        self%process_matrix(2, 3) = sum(self%process_matrix(1 : 2, 2))
      case(1)
        ! 1 - degree-day approach
        self%process_matrix(2, 2) = 8_i4
        self%process_matrix(2, 3) = sum(self%process_matrix(1 : 2, 2))
        call append(self%definition, reshape(self%config%snow1%snowTreshholdTemperature, shp))
        call append(self%definition, reshape(self%config%snow1%degreeDayFactor_forest, shp))
        call append(self%definition, reshape(self%config%snow1%degreeDayFactor_impervious, shp))
        call append(self%definition, reshape(self%config%snow1%degreeDayFactor_pervious, shp))
        call append(self%definition, reshape(self%config%snow1%increaseDegreeDayFactorByPrecip, shp))
        call append(self%definition, reshape(self%config%snow1%maxDegreeDayFactor_forest, shp))
        call append(self%definition, reshape(self%config%snow1%maxDegreeDayFactor_impervious, shp))
        call append(self%definition, reshape(self%config%snow1%maxDegreeDayFactor_pervious, shp))

        call append(self%names, [  &
                'snowTreshholdTemperature       ', &
                'degreeDayFactor_forest         ', &
                'degreeDayFactor_impervious     ', &
                'degreeDayFactor_pervious       ', &
                'increaseDegreeDayFactorByPrecip', &
                'maxDegreeDayFactor_forest      ', &
                'maxDegreeDayFactor_impervious  ', &
                'maxDegreeDayFactor_pervious    '])

        ! check if parameter are in range
        if (.not. in_bound(self%definition)) then
          log_fatal(*) 'Parameter in namelist "snow1" out of bound.'
          error stop 1
        end if
    end select

    ! Process 3 - soilmoisture
    select case (self%process_matrix(3, 1))
      case(0)
        ! 0 - no soil moisture process selected
        log_debug(*) "SELECTION: soil moisture is deactivated!"
        self%process_matrix(3, 3) = sum(self%process_matrix(1 : 3, 2))
      case(1)
        ! 1 - Feddes equation for PET reduction, bucket approach, Brooks-Corey like
        self%process_matrix(3, 2) = 17_i4
        self%process_matrix(3, 3) = sum(self%process_matrix(1 : 3, 2))
        call append(self%definition, reshape(self%config%soilmoisture1%orgMatterContent_forest, shp))
        call append(self%definition, reshape(self%config%soilmoisture1%orgMatterContent_impervious, shp))
        call append(self%definition, reshape(self%config%soilmoisture1%orgMatterContent_pervious, shp))
        call append(self%definition, reshape(self%config%soilmoisture1%PTF_lower66_5_constant, shp))
        call append(self%definition, reshape(self%config%soilmoisture1%PTF_lower66_5_clay, shp))
        call append(self%definition, reshape(self%config%soilmoisture1%PTF_lower66_5_Db, shp))
        call append(self%definition, reshape(self%config%soilmoisture1%PTF_higher66_5_constant, shp))
        call append(self%definition, reshape(self%config%soilmoisture1%PTF_higher66_5_clay, shp))
        call append(self%definition, reshape(self%config%soilmoisture1%PTF_higher66_5_Db, shp))
        call append(self%definition, reshape(self%config%soilmoisture1%PTF_Ks_constant, shp))
        call append(self%definition, reshape(self%config%soilmoisture1%PTF_Ks_sand, shp))
        call append(self%definition, reshape(self%config%soilmoisture1%PTF_Ks_clay, shp))
        call append(self%definition, reshape(self%config%soilmoisture1%PTF_Ks_curveSlope, shp))
        call append(self%definition, reshape(self%config%soilmoisture1%rootFractionCoefficient_forest, shp))
        call append(self%definition, reshape(self%config%soilmoisture1%rootFractionCoefficient_impervious, shp))
        call append(self%definition, reshape(self%config%soilmoisture1%rootFractionCoefficient_pervious, shp))
        call append(self%definition, reshape(self%config%soilmoisture1%infiltrationShapeFactor, shp))

        call append(self%names, [     &
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
                'infiltrationShapeFactor           '])

        ! check if parameter are in range
        if (.not. in_bound(self%definition)) then
          log_fatal(*) 'Parameter in namelist "soilmoisture1" out of bound.'
          error stop 1
        end if

        ! 2- Jarvis equation for PET reduction, bucket approach, Brooks-Corey like
      case(2)
        self%process_matrix(3, 2) = 18_i4
        self%process_matrix(3, 3) = sum(self%process_matrix(1 : 3, 2))
        call append(self%definition, reshape(self%config%soilmoisture2%orgMatterContent_forest, shp))
        call append(self%definition, reshape(self%config%soilmoisture2%orgMatterContent_impervious, shp))
        call append(self%definition, reshape(self%config%soilmoisture2%orgMatterContent_pervious, shp))
        call append(self%definition, reshape(self%config%soilmoisture2%PTF_lower66_5_constant, shp))
        call append(self%definition, reshape(self%config%soilmoisture2%PTF_lower66_5_clay, shp))
        call append(self%definition, reshape(self%config%soilmoisture2%PTF_lower66_5_Db, shp))
        call append(self%definition, reshape(self%config%soilmoisture2%PTF_higher66_5_constant, shp))
        call append(self%definition, reshape(self%config%soilmoisture2%PTF_higher66_5_clay, shp))
        call append(self%definition, reshape(self%config%soilmoisture2%PTF_higher66_5_Db, shp))
        call append(self%definition, reshape(self%config%soilmoisture2%PTF_Ks_constant, shp))
        call append(self%definition, reshape(self%config%soilmoisture2%PTF_Ks_sand, shp))
        call append(self%definition, reshape(self%config%soilmoisture2%PTF_Ks_clay, shp))
        call append(self%definition, reshape(self%config%soilmoisture2%PTF_Ks_curveSlope, shp))
        call append(self%definition, reshape(self%config%soilmoisture2%rootFractionCoefficient_forest, shp))
        call append(self%definition, reshape(self%config%soilmoisture2%rootFractionCoefficient_impervious, shp))
        call append(self%definition, reshape(self%config%soilmoisture2%rootFractionCoefficient_pervious, shp))
        call append(self%definition, reshape(self%config%soilmoisture2%infiltrationShapeFactor, shp))
        call append(self%definition, reshape(self%config%soilmoisture2%jarvis_sm_threshold_c1, shp))

        call append(self%names, [     &
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
                'jarvis_sm_threshold_c1            '])

        ! check if parameter are in range
        if (.not. in_bound(self%definition)) then
          log_fatal(*) 'Parameter in namelist "soilmoisture2" out of bound.'
          error stop 1
        end if

        ! 3- Jarvis equation for ET reduction and FC dependency on root fraction coefficient
      case(3)
        self%process_matrix(3, 2) = 22_i4
        self%process_matrix(3, 3) = sum(self%process_matrix(1 : 3, 2))
        call append(self%definition, reshape(self%config%soilmoisture3%orgMatterContent_forest, shp))
        call append(self%definition, reshape(self%config%soilmoisture3%orgMatterContent_impervious, shp))
        call append(self%definition, reshape(self%config%soilmoisture3%orgMatterContent_pervious, shp))
        call append(self%definition, reshape(self%config%soilmoisture3%PTF_lower66_5_constant, shp))
        call append(self%definition, reshape(self%config%soilmoisture3%PTF_lower66_5_clay, shp))
        call append(self%definition, reshape(self%config%soilmoisture3%PTF_lower66_5_Db, shp))
        call append(self%definition, reshape(self%config%soilmoisture3%PTF_higher66_5_constant, shp))
        call append(self%definition, reshape(self%config%soilmoisture3%PTF_higher66_5_clay, shp))
        call append(self%definition, reshape(self%config%soilmoisture3%PTF_higher66_5_Db, shp))
        call append(self%definition, reshape(self%config%soilmoisture3%PTF_Ks_constant, shp))
        call append(self%definition, reshape(self%config%soilmoisture3%PTF_Ks_sand, shp))
        call append(self%definition, reshape(self%config%soilmoisture3%PTF_Ks_clay, shp))
        call append(self%definition, reshape(self%config%soilmoisture3%PTF_Ks_curveSlope, shp))
        call append(self%definition, reshape(self%config%soilmoisture3%rootFractionCoefficient_forest, shp))
        call append(self%definition, reshape(self%config%soilmoisture3%rootFractionCoefficient_impervious, shp))
        call append(self%definition, reshape(self%config%soilmoisture3%rootFractionCoefficient_pervious, shp))
        call append(self%definition, reshape(self%config%soilmoisture3%infiltrationShapeFactor, shp))
        call append(self%definition, reshape(self%config%soilmoisture3%rootFractionCoefficient_sand, shp))
        call append(self%definition, reshape(self%config%soilmoisture3%rootFractionCoefficient_clay, shp))
        call append(self%definition, reshape(self%config%soilmoisture3%FCmin_glob, shp))
        call append(self%definition, reshape(self%config%soilmoisture3%FCdelta_glob, shp))
        call append(self%definition, reshape(self%config%soilmoisture3%jarvis_sm_threshold_c1, shp))

        call append(self%names, [     &
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
                'jarvis_sm_threshold_c1            '])

        ! check if parameter are in range
        if (.not. in_bound(self%definition)) then
          log_fatal(*) 'Parameter in namelist "soilmoisture3" out of bound.'
          error stop 1
        end if

        ! 4- Feddes equation for ET reduction and FC dependency on root fraction coefficient
      case(4)
        self%process_matrix(3, 2) = 21_i4
        self%process_matrix(3, 3) = sum(self%process_matrix(1 : 3, 2))
        call append(self%definition, reshape(self%config%soilmoisture4%orgMatterContent_forest, shp))
        call append(self%definition, reshape(self%config%soilmoisture4%orgMatterContent_impervious, shp))
        call append(self%definition, reshape(self%config%soilmoisture4%orgMatterContent_pervious, shp))
        call append(self%definition, reshape(self%config%soilmoisture4%PTF_lower66_5_constant, shp))
        call append(self%definition, reshape(self%config%soilmoisture4%PTF_lower66_5_clay, shp))
        call append(self%definition, reshape(self%config%soilmoisture4%PTF_lower66_5_Db, shp))
        call append(self%definition, reshape(self%config%soilmoisture4%PTF_higher66_5_constant, shp))
        call append(self%definition, reshape(self%config%soilmoisture4%PTF_higher66_5_clay, shp))
        call append(self%definition, reshape(self%config%soilmoisture4%PTF_higher66_5_Db, shp))
        call append(self%definition, reshape(self%config%soilmoisture4%PTF_Ks_constant, shp))
        call append(self%definition, reshape(self%config%soilmoisture4%PTF_Ks_sand, shp))
        call append(self%definition, reshape(self%config%soilmoisture4%PTF_Ks_clay, shp))
        call append(self%definition, reshape(self%config%soilmoisture4%PTF_Ks_curveSlope, shp))
        call append(self%definition, reshape(self%config%soilmoisture4%rootFractionCoefficient_forest, shp))
        call append(self%definition, reshape(self%config%soilmoisture4%rootFractionCoefficient_impervious, shp))
        call append(self%definition, reshape(self%config%soilmoisture4%rootFractionCoefficient_pervious, shp))
        call append(self%definition, reshape(self%config%soilmoisture4%infiltrationShapeFactor, shp))
        call append(self%definition, reshape(self%config%soilmoisture4%rootFractionCoefficient_sand, shp))
        call append(self%definition, reshape(self%config%soilmoisture4%rootFractionCoefficient_clay, shp))
        call append(self%definition, reshape(self%config%soilmoisture4%FCmin_glob, shp))
        call append(self%definition, reshape(self%config%soilmoisture4%FCdelta_glob, shp))

        call append(self%names, [     &
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
                'FCdelta_glob                      '])

        ! check if parameter are in range
        if (.not. in_bound(self%definition)) then
          log_fatal(*) 'Parameter in namelist "soilmoisture4" out of bound.'
          error stop 1
        end if
    end select

    ! Process 4 - sealed area directRunoff
    select case (self%process_matrix(4, 1))
      case(0)
        ! 0 - no direct runoff process selected
        log_debug(*) "SELECTION: direct runoff is deactivated!"
        self%process_matrix(4, 3) = sum(self%process_matrix(1 : 4, 2))
      case(1)
        ! 1 - bucket exceedance approach
        self%process_matrix(4, 2) = 1_i4
        self%process_matrix(4, 3) = sum(self%process_matrix(1 : 4, 2))
        call append(self%definition, reshape(self%config%directrunoff1%imperviousStorageCapacity, shp))

        call append(self%names, ['imperviousStorageCapacity'])

        ! check if parameter are in range
        if (.not. in_bound(self%definition)) then
          log_fatal(*) 'Parameter in namelist "directRunoff1" out of bound.'
          error stop 1
        end if
    end select

    ! Process 5 - potential evapotranspiration (PET)
    select case (self%process_matrix(5, 1))
      case(0) ! 0 - no PET process selected
        log_debug(*) "SELECTION: PET is deactivated!"
        self%process_matrix(5, 3) = sum(self%process_matrix(1 : 5, 2))
      case(-1) ! -1 - PET is input, correct PET by LAI
        self%process_matrix(5, 2) = 5_i4
        self%process_matrix(5, 3) = sum(self%process_matrix(1 : 5, 2))
        call append(self%definition, reshape(self%config%petm1%PET_a_forest, shp))
        call append(self%definition, reshape(self%config%petm1%PET_a_impervious, shp))
        call append(self%definition, reshape(self%config%petm1%PET_a_pervious, shp))
        call append(self%definition, reshape(self%config%petm1%PET_b, shp))
        call append(self%definition, reshape(self%config%petm1%PET_c, shp))

        call append(self%names, [ &
                'PET_a_forest     ', &
                'PET_a_impervious ', &
                'PET_a_pervious   ', &
                'PET_b            ', &
                'PET_c            '])

        ! check if parameter are in range
        if (.not. in_bound(self%definition)) then
          log_fatal(*) 'Parameter in namelist "PETm1" out of bound.'
          error stop 1
        end if

      case(-2) ! -2 - PET is input, correct PET by aspect
        self%process_matrix(5, 2) = 3_i4
        self%process_matrix(5, 3) = sum(self%process_matrix(1 : 5, 2))
        call append(self%definition, reshape(self%config%petm2%minCorrectionFactorPET, shp))
        call append(self%definition, reshape(self%config%petm2%maxCorrectionFactorPET, shp))
        call append(self%definition, reshape(self%config%petm2%aspectTresholdPET, shp))

        call append(self%names, [ &
                'minCorrectionFactorPET ', &
                'maxCorrectionFactorPET ', &
                'aspectTresholdPET      '])

        ! check if parameter are in range
        if (.not. in_bound(self%definition)) then
          log_fatal(*) 'Parameter in namelist "PETm2" out of bound.'
          error stop 1
        end if

      case(1) ! 1 - Hargreaves-Samani method (HarSam) - additional input needed: Tmin, Tmax
        self%process_matrix(5, 2) = 4_i4
        self%process_matrix(5, 3) = sum(self%process_matrix(1 : 5, 2))
        call append(self%definition, reshape(self%config%pet1%minCorrectionFactorPET, shp))
        call append(self%definition, reshape(self%config%pet1%maxCorrectionFactorPET, shp))
        call append(self%definition, reshape(self%config%pet1%aspectTresholdPET, shp))
        call append(self%definition, reshape(self%config%pet1%HargreavesSamaniCoeff, shp))
        call append(self%names, [ &
                'minCorrectionFactorPET', &
                'maxCorrectionFactorPET', &
                'aspectTresholdPET     ', &
                'HargreavesSamaniCoeff '])

        ! check if parameter are in range
        if (.not. in_bound(self%definition)) then
          log_fatal(*) 'Parameter in namelist "PET1" out of bound.'
          error stop 1
        end if

      case(2) ! 2 - Priestley-Taylor method (PrieTay) - additional input needed: net_rad
        self%process_matrix(5, 2) = 2_i4
        self%process_matrix(5, 3) = sum(self%process_matrix(1 : 5, 2))
        call append(self%definition, reshape(self%config%pet2%PriestleyTaylorCoeff, shp))
        call append(self%definition, reshape(self%config%pet2%PriestleyTaylorLAIcorr, shp))
        call append(self%names, [ &
                'PriestleyTaylorCoeff  ', &
                'PriestleyTaylorLAIcorr'])

        ! check if parameter are in range
        if (.not. in_bound(self%definition)) then
          log_fatal(*) 'Parameter in namelist "PET2" out of bound.'
          error stop 1
        end if

      case(3) ! 3 - Penman-Monteith method - additional input needed: net_rad, abs. vapour pressue, windspeed
        self%process_matrix(5, 2) = 7_i4
        self%process_matrix(5, 3) = sum(self%process_matrix(1 : 5, 2))

        call append(self%definition, reshape(self%config%pet3%canopyheigth_forest, shp))
        call append(self%definition, reshape(self%config%pet3%canopyheigth_impervious, shp))
        call append(self%definition, reshape(self%config%pet3%canopyheigth_pervious, shp))
        call append(self%definition, reshape(self%config%pet3%displacementheight_coeff, shp))
        call append(self%definition, reshape(self%config%pet3%roughnesslength_momentum_coeff, shp))
        call append(self%definition, reshape(self%config%pet3%roughnesslength_heat_coeff, shp))
        call append(self%definition, reshape(self%config%pet3%stomatal_resistance, shp))

        call append(self%names, [ &
                'canopyheigth_forest           ', &
                'canopyheigth_impervious       ', &
                'canopyheigth_pervious         ', &
                'displacementheight_coeff      ', &
                'roughnesslength_momentum_coeff', &
                'roughnesslength_heat_coeff    ', &
                'stomatal_resistance           '])

        ! check if parameter are in range
        if (.not. in_bound(self%definition)) then
          log_fatal(*) 'Parameter in namelist "PET3" out of bound.'
          error stop 1
        end if
    end select

    ! Process 6 - interflow
    select case (self%process_matrix(6, 1))
      case(0)
        ! 0 - deactivated
        log_debug(*) "SELECTION: Interflow is deactivated!"
        self%process_matrix(6, 3) = sum(self%process_matrix(1 : 6, 2))
      case(1)
        ! 1 - parallel soil reservoir approach
        self%process_matrix(6, 2) = 5_i4
        self%process_matrix(6, 3) = sum(self%process_matrix(1 : 6, 2))
        call append(self%definition, reshape(self%config%interflow1%interflowStorageCapacityFactor, shp))
        call append(self%definition, reshape(self%config%interflow1%interflowRecession_slope, shp))
        call append(self%definition, reshape(self%config%interflow1%fastInterflowRecession_forest, shp))
        call append(self%definition, reshape(self%config%interflow1%slowInterflowRecession_Ks, shp))
        call append(self%definition, reshape(self%config%interflow1%exponentSlowInterflow, shp))

        call append(self%names, [ &
                'interflowStorageCapacityFactor', &
                'interflowRecession_slope      ', &
                'fastInterflowRecession_forest ', &
                'slowInterflowRecession_Ks     ', &
                'exponentSlowInterflow         '])

        ! check if parameter are in range
        if (.not. in_bound(self%definition)) then
          log_fatal(*) 'Parameter in namelist "interflow1" out of bound.'
          error stop 1
        end if
    end select

    ! Process 7 - percolation
    select case (self%process_matrix(7, 1))
      case(0)
        ! 0 - deactivated
         log_debug(*) "SELECTION: Percolation is deativated!"
        self%process_matrix(7, 3) = sum(self%process_matrix(1 : 7, 2))
      case(1)
        ! 1 - GW layer is assumed as bucket
        self%process_matrix(7, 2) = 3_i4
        self%process_matrix(7, 3) = sum(self%process_matrix(1 : 7, 2))
        call append(self%definition, reshape(self%config%percolation1%rechargeCoefficient, shp))
        call append(self%definition, reshape(self%config%percolation1%rechargeFactor_karstic, shp))
        call append(self%definition, reshape(self%config%percolation1%gain_loss_GWreservoir_karstic, shp))

        call append(self%names, [ &
                'rechargeCoefficient          ', &
                'rechargeFactor_karstic       ', &
                'gain_loss_GWreservoir_karstic'])

        ! check if parameter are in range
        if (.not. in_bound(self%definition)) then
          log_fatal(*) 'Parameter in namelist "percolation1" out of bound.'
          error stop 1
        end if
    end select

    ! Process 8 - routing
    select case (self%process_matrix(8, 1))
      case(0)
        ! 0 - deactivated
         log_debug(*) "SELECTION: Routing is deactivated!"
        self%process_matrix(8, 3) = sum(self%process_matrix(1 : 8, 2))
      case(1)
        ! 1 - Muskingum approach
        self%process_matrix(8, 2) = 5_i4
        self%process_matrix(8, 3) = sum(self%process_matrix(1 : 8, 2))

        call append(self%definition, reshape(self%config%routing1%muskingumTravelTime_constant, shp))
        call append(self%definition, reshape(self%config%routing1%muskingumTravelTime_riverLength, shp))
        call append(self%definition, reshape(self%config%routing1%muskingumTravelTime_riverSlope, shp))
        call append(self%definition, reshape(self%config%routing1%muskingumTravelTime_impervious, shp))
        call append(self%definition, reshape(self%config%routing1%muskingumAttenuation_riverSlope, shp))

        call append(self%names, [&
                'muskingumTravelTime_constant   ', &
                'muskingumTravelTime_riverLength', &
                'muskingumTravelTime_riverSlope ', &
                'muskingumTravelTime_impervious ', &
                'muskingumAttenuation_riverSlope'])
        ! check if parameter are in range
        if (.not. in_bound(self%definition)) then
          log_fatal(*) 'Parameter in namelist "routing1" out of bound.'
          error stop 1
        end if
      case(2)
        self%process_matrix(8, 2) = 1_i4
        self%process_matrix(8, 3) = sum(self%process_matrix(1 : 8, 2))
        call append(self%definition, reshape(self%config%routing2%streamflow_celerity, shp))
        call append(self%names, ['streamflow_celerity'])
        ! check if parameter are in range
        if (.not. in_bound(self%definition)) then
          log_fatal(*) 'Parameter in namelist "routing2" out of bound.'
          error stop 1
        end if
      case(3)
        self%process_matrix(8, 2) = 1_i4
        self%process_matrix(8, 3) = sum(self%process_matrix(1 : 8, 2))
        call append(self%definition, reshape(self%config%routing3%slope_factor, shp))
        call append(self%names, ['slope_factor'])
        ! check if parameter are in range
        if (.not. in_bound(self%definition)) then
          log_fatal(*) 'Parameter in namelist "routing3" out of bound.'
          error stop 1
        end if
    end select

    ! Process 9 - geoparameter
    select case (self%process_matrix(9, 1))
      case(0)
        ! 0 - deactivated
        log_debug(*) "SELECTION: Baseflow is deactivated!"
        self%process_matrix(9, 3) = sum(self%process_matrix(1 : 9, 2))
      case(1)
        ! read in global parameters (NOT REGIONALIZED, i.e. these are <beta> and not <gamma>) for each geological formation used
        GeoParam = self%config%geoparameter%GeoParam
        if (size(GeoParam, 2) /= self%nGeoUnits) then
          log_fatal(*) 'GeoParam shape in namelist "geoparameter" does not match n_geo_units.'
          error stop 1
        end if

        ! for geology parameters
        self%process_matrix(9, 2) = self%nGeoUnits
        self%process_matrix(9, 3) = sum(self%process_matrix(1 : 9, 2))

        call append(self%definition, transpose(GeoParam(:, 1_i4:self%nGeoUnits)))

        ! create names
        do ii = 1, self%nGeoUnits
          geo_name = 'GeoParam(' // trim(adjustl(n2s(ii))) // ')'
          call append(self%names, [ trim(geo_name) ])
        end do

        ! check if parameter are in range
        if (.not. in_bound(self%definition)) then
          log_fatal(*) 'Parameter in namelist "geoparameter" out of bound.'
          error stop 1
        end if
    end select

    ! Process 10 - neutrons
    select case (self%process_matrix(10, 1))
      case(0)
        ! 0 - deactivated
        log_debug(*) "SELECTION: Neutron count routine is deativated!"
        self%process_matrix(10, 3) = sum(self%process_matrix(1 : 10, 2))
      case(1)
        ! 1 - inverse N0 based on Desilets et al. 2010
        self%process_matrix(10,2) = 3_i4
        self%process_matrix(10,3) = sum(self%process_matrix(1:10, 2))
        call append(self%definition, reshape(self%config%neutrons1%Desilets_N0,  shp))
        call append(self%definition, reshape(self%config%neutrons1%Desilets_LW0, shp))
        call append(self%definition, reshape(self%config%neutrons1%Desilets_LW1, shp))
        call append(self%names, [  &
                'Desilets_N0   ', &
                'Desilets_LW0  ', &
                'Desilets_LW1  '])
        ! check if parameter are in range
        if (.not. in_bound(self%definition)) then
          log_fatal(*) 'Parameter in namelist "neutrons1" out of bound.'
          error stop 1
        end if
      case(2)
        ! 2 - COSMIC forward operator by Shuttlworth et al. 2013
        self%process_matrix(10,2) = 9_i4
        self%process_matrix(10,3) = sum(self%process_matrix(1:10, 2))
        call append(self%definition, reshape(self%config%neutrons2%COSMIC_N0,     shp))
        call append(self%definition, reshape(self%config%neutrons2%COSMIC_N1,     shp))
        call append(self%definition, reshape(self%config%neutrons2%COSMIC_N2,     shp))
        call append(self%definition, reshape(self%config%neutrons2%COSMIC_alpha0, shp))
        call append(self%definition, reshape(self%config%neutrons2%COSMIC_alpha1, shp))
        call append(self%definition, reshape(self%config%neutrons2%COSMIC_L30,    shp))
        call append(self%definition, reshape(self%config%neutrons2%COSMIC_L31,    shp))
        call append(self%definition, reshape(self%config%neutrons2%COSMIC_LW0,    shp))
        call append(self%definition, reshape(self%config%neutrons2%COSMIC_LW1,    shp))
        call append(self%names, [  &
                'COSMIC_N0     ', &
                'COSMIC_N1     ', &
                'COSMIC_N2     ', &
                'COSMIC_alpha0 ', &
                'COSMIC_alpha1 ', &
                'COSMIC_L30    ', &
                'COSMIC_L31    ', &
                'COSMIC_LW0    ', &
                'COSMIC_LW1    '])
        ! check if parameter are in range
        if (.not. in_bound(self%definition)) then
          log_fatal(*) 'Parameter in namelist "neutrons2" out of bound.'
          error stop 1
        end if
    end select

    ! Process 11 - rivertemp
    select case (self%process_matrix(11, 1))
      case(0)
        ! 0 - deactivated
        log_debug(*) "SELECTION: River temperature routine is deativated!"
        self%process_matrix(11, 3) = sum(self%process_matrix(1 : 11, 2))
      case(1)
        ! 1 - rivertemp1
        self%process_matrix(11,2) = 4_i4
        self%process_matrix(11,3) = sum(self%process_matrix(1:11, 2))
        call append(self%definition, reshape(self%config%rivertemp1%albedo_water,       shp))
        call append(self%definition, reshape(self%config%rivertemp1%pt_a_water,         shp))
        call append(self%definition, reshape(self%config%rivertemp1%emissivity_water,   shp))
        call append(self%definition, reshape(self%config%rivertemp1%turb_heat_ex_coeff, shp))
        call append(self%names, [  &
                'albedo_water      ', &
                'pt_a_water        ', &
                'emissivity_water  ', &
                'turb_heat_ex_coeff'])
        ! check if parameter are in range
        if (.not. in_bound(self%definition)) then
          log_fatal(*) 'Parameter in namelist "rivertemp1" out of bound.'
          error stop 1
        end if
    end select

    ! set default values; pass-through-only configurations can legitimately have no parameters.
    if (allocated(self%definition)) then
      self%values = self%definition(:, 3)
    else
      allocate(self%values(0))
    end if

  end subroutine parameters_initialize

  !> set parameters from optimizer
  subroutine parameters_set(self, parameters)
    class(parameters_t), intent(inout) :: self
    real(dp), dimension(:), optional, intent(in) :: parameters !< parameter values
    if (present(parameters)) then
      if (size(parameters) /= size(self%definition, dim=1)) then
        log_fatal(*) "parameters%set: wrong number of parameters"
        error stop 1
      end if
      self%values = parameters
    else
      if (allocated(self%definition)) then
        self%values = self%definition(:, 3) ! default parameters from configuration
      else
        if (.not. allocated(self%values)) allocate(self%values(0))
      end if
    end if
  end subroutine parameters_set

  !> get parameter by name
  real(dp) function parameters_get(self, name)
    class(parameters_t), intent(in) :: self
    character(*), intent(in) ::  name !< parameter name
    integer(i4) :: i
    do i = 1_i4, size(self%names)
      if (trim(name) == trim(self%names(i))) then
        parameters_get = self%values(i)
        return
      end if
    end do
    log_fatal(*) "parameters%get: parameter '", trim(name), "' not present."
    error stop 1
  end function parameters_get

  !> get parameter array for selected process
  function parameters_get_process(self, id) result(array)
    class(parameters_t), intent(in) :: self
    integer(i4), intent(in) ::  id !< process id
    real(dp), allocatable :: array(:) !< resulting parameter array for selected process
    if (id < 1_i4 .or. id >= size(self%process_matrix, dim=1)) then
      log_fatal(*) "parameters%get_process: process id '", n2s(id), "' out of bounds."
      error stop 1
    end if
    if (self%process_matrix(id, 3_i4) == 0_i4) then
      allocate(array(0)) ! no parameters defined for this process
      return
    end if
    array = self%values(self%process_matrix(id, 3_i4) - self%process_matrix(id, 2_i4) + 1_i4 : self%process_matrix(id, 3_i4))
  end function parameters_get_process

  !> set parameters from optimizer
  subroutine parameters_print(self)
    class(parameters_t), intent(inout) :: self
    integer(i4) :: i
    do i = 1_i4, size(self%names)
      log_text(*) self%names(i)(1:35), n2s(self%values(i))
    end do
  end subroutine parameters_print

  !> check if all parameters are within their defined bounds
  logical function in_bound(params)
    real(dp), dimension(:, :), intent(in) :: params
    in_bound = .not.(any(params(:, 3) .lt. params(:, 1)) .or. any(params(:, 3) .gt. params(:, 2)))
  end function in_bound

end module mo_main_config
