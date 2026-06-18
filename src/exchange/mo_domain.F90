!> \file    mo_domain.F90
!> \copydoc mo_domain

!> \brief   Module for a domain container.
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    Aug 2025
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_exchange
#include "logging.h"
module mo_domain
  use mo_logging
  use mo_list, only: list
  use mo_kind, only: i4, dp
  use mo_exchange_type, only: exchange_t
  use mo_string_utils, only: n2s=>num2str
  use mo_datetime, only: datetime, timedelta
  ! config
  use mo_main_config, only: parameters_t
  ! containers
  use mo_input_container, only: input_t
  use mo_meteo_container, only: meteo_t
  use mo_mpr_container, only: mpr_t
  use mo_mhm_container, only: mhm_t
  use mo_mrm_container, only: mrm_t

  !> \class   domain_t
  !> \brief   Class for a single mHM domain.
  !> \details Always add the "target" attribute to instances of domain, or use allocated pointers (have implicit target attribute).
  type, public :: domain_t
    type(exchange_t) :: exchange !< the exchange container with all exchanged variables for this domain
    type(input_t) :: input !< the input container providing inputs from file or couplers
    type(meteo_t) :: meteo !< the meteorology container providing processed meteorological data
    type(mpr_t) :: mpr !< the MPR container providing parameter fields for the process containers
    type(mhm_t) :: mhm !< the mHM process container calculating vertical hydrological processes
    type(mrm_t) :: mrm !< the mRM process container for routing related processes
    logical :: is_finished = .false. !< whether the time-loop on this domain is finished
  contains
    procedure :: init => domain_init
    procedure :: configure => domain_configure
    procedure :: connect => domain_connect
    procedure :: initialize => domain_initialize
    procedure :: update => domain_update
    procedure :: finalize => domain_finalize
  end type domain_t

  !> \class   domain_list
  !> \brief   Class to hold a list of all domains as a linked list of pointers.
  type, extends(list) :: domain_list
    integer(i4) :: counter = 0_i4 !< internal counter for added domains
  contains
    procedure :: get_domain => domain_list_get
    procedure :: add_domain => domain_list_add
  end type domain_list

  type(domain_list), public :: domains
  integer(i4), allocatable, public :: selected_domains(:)

contains

  !> \brief Get pointer to desired domain from domain list.
  subroutine domain_list_get(self, key, domain)
    class(domain_list), intent(in) :: self
    integer(i4), intent(in) :: key !< domain key
    type(domain_t), pointer, intent(out) :: domain !< pointer to desired domain
    class(*), pointer :: p
    call self%get(key, p)
    if (associated(p)) then
      select type (p)
        type is (domain_t)
          domain => p
        class default
          log_fatal(*) "Domain '", key, "' not a domain type."
          error stop
      end select
    else
      log_fatal(*) "Domain '", key, "' not present."
      error stop
    end if
  end subroutine domain_list_get

  !> \brief Add a new domain to the domain list.
  function domain_list_add(self, key) result(new_key)
    class(domain_list), intent(inout) :: self
    integer(i4), optional, intent(in) :: key !< domain key
    integer(i4) :: new_key !< key of the added domain
    type(domain_t), save :: new_domain ! template for the new domain to be added (saved to avoid reallocation on each call)
    self%counter = self%counter + 1_i4
    if (present(key)) then
      new_key = key
    else
      new_key = self%counter
    end if
    call self%add_clone(new_key, new_domain)
  end function domain_list_add

  !> \brief Create a new domain.
  subroutine domain_init(self, meta_file, main_file, para_file, domain, cwd, run_n_domains, run_n_geo_units, run_max_layers, &
      read_domains_from_dirs)
    ! domain is always an item of a domain_list, which stores "allocated pointers" and these implicitly have the "target" attribute
    class(domain_t), intent(inout), target :: self ! needs "target" so components can safely point to "exchange"
    character(*), intent(in), optional :: meta_file !< file containing the metadata namelists
    character(*), intent(in), optional :: main_file !< file containing the main namelists
    character(*), intent(in), optional :: para_file !< file containing the parameter namelists
    integer(i4), intent(in), optional :: domain !< domain ID of the current domain in the configuration arrays (1 by default)
    character(len=*), intent(in), optional :: cwd !< current working directory to set relative paths
    integer(i4), intent(in), optional :: run_n_domains !< total number of domains in the run
    integer(i4), intent(in), optional :: run_n_geo_units !< total number of geological units in the global parameter set
    integer(i4), intent(in), optional :: run_max_layers !< maximum number of soil-layer entries in the run
    logical, intent(in), optional :: read_domains_from_dirs !< whether domains are read from separate directories
    log_info(*) "SETUP NEW DOMAIN"
    call self%exchange%init(meta_file, main_file, para_file, domain, cwd, run_n_domains, run_n_geo_units, run_max_layers, &
      read_domains_from_dirs)
    ! set exchange pointer in components
    self%input%exchange => self%exchange
    self%meteo%exchange => self%exchange
    self%mpr%exchange => self%exchange
    self%mhm%exchange => self%exchange
    self%mrm%exchange => self%exchange
  end subroutine domain_init

  !> \brief Configure the domain.
  subroutine domain_configure(self, file, out_file)
    ! domain is always an item of a domain_list, which stores "allocated pointers" and these implicitly have the "target" attribute
    class(domain_t), intent(inout), target :: self ! needs "target" so components can safely point to "exchange"
    character(*), intent(in), optional :: file !< file containing the namelists
    character(*), intent(in), optional :: out_file !< file containing the output namelists
    log_info(*) "CONFIGURE COMPONENTS"
    call self%input%configure(file)
    if (self%exchange%parameters%mhm_active()) call self%mpr%configure(file)
    if (self%exchange%parameters%meteo_active()) call self%meteo%configure(file)
    if (self%exchange%parameters%mhm_active()) call self%mhm%configure(file, out_file)
    if (self%exchange%parameters%mrm_active()) call self%mrm%configure(file, out_file)
  end subroutine domain_configure

  !> \brief Connect the domain components.
  !> \details Check for dependencies and connect exchanged arrays between components after configuration.
  subroutine domain_connect(self)
    ! domain is always an item of a domain_list, which stores "allocated pointers" and these implicitly have the "target" attribute
    class(domain_t), intent(inout), target :: self ! needs "target" so components can safely point to "exchange"
    log_info(*) "CONNECT COMPONENTS"
    call self%input%connect()
    if (self%exchange%parameters%mhm_active()) call self%mpr%connect()
    if (self%exchange%parameters%meteo_active()) call self%meteo%connect()
    if (self%exchange%parameters%mhm_active()) call self%mhm%connect()
    if (self%exchange%parameters%mrm_active()) call self%mrm%connect()
  end subroutine domain_connect

  !> \brief Initialize the domain and do the initial state calculations in the components.
  subroutine domain_initialize(self, parameters)
    ! domain is always an item of a domain_list, which stores "allocated pointers" and these implicitly have the "target" attribute
    class(domain_t), intent(inout), target :: self
    !> a set of global parameter (gamma) to run mHM, DIMENSION [no. of global_Parameters]
    real(dp), dimension(:), optional, intent(in) :: parameters
    log_info(*) "INITIALIZE DOMAIN"
    self%exchange%step_count = 0_i4
    self%exchange%time = self%exchange%start_time
    call self%exchange%parameters%set(parameters)
    call self%input%initialize()
    if (self%exchange%parameters%mhm_active()) call self%mpr%initialize()
    if (self%exchange%parameters%meteo_active()) call self%meteo%initialize()
    if (self%exchange%parameters%mhm_active()) call self%mhm%initialize()
    if (self%exchange%parameters%mrm_active()) call self%mrm%initialize()
  end subroutine domain_initialize

  !> \brief Update the domain for the current time step.
  subroutine domain_update(self)
    ! domain is always an item of a domain_list, which stores "allocated pointers" and these implicitly have the "target" attribute
    class(domain_t), intent(inout), target :: self
    call self%exchange%time%add(self%exchange%step)
    self%exchange%step_count = self%exchange%step_count + 1_i4
    call self%input%update()
    if (self%exchange%parameters%mhm_active()) call self%mpr%update()
    if (self%exchange%parameters%meteo_active()) call self%meteo%update()
    if (self%exchange%parameters%mhm_active()) call self%mhm%update()
    if (self%exchange%parameters%mrm_active()) call self%mrm%update()
    log_trace(*) "Time step: ", self%exchange%time%str()
    if (self%exchange%time%is_new_year()) then
      log_info(*) "Finished year: ", self%exchange%time%year - 1_i4
    end if
  end subroutine domain_update

  !> \brief Finalize the domain and its components after the simulation.
  subroutine domain_finalize(self)
    ! domain is always an item of a domain_list, which stores "allocated pointers" and these implicitly have the "target" attribute
    class(domain_t), intent(inout), target :: self
    log_info(*) "FINALIZE COMPONENTS"
    call self%input%finalize()
    if (self%exchange%parameters%mhm_active()) call self%mpr%finalize()
    if (self%exchange%parameters%meteo_active()) call self%meteo%finalize()
    if (self%exchange%parameters%mhm_active()) call self%mhm%finalize()
    if (self%exchange%parameters%mrm_active()) call self%mrm%finalize()
  end subroutine domain_finalize

end module mo_domain
