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
  use mo_kind, only: i4, dp
  use mo_exchange_type, only: exchange_t
  ! containers
  use mo_input_container, only: input_t
  use mo_meteo_container, only: meteo_t
  use mo_mpr_container, only: mpr_t
  use mo_mhm_container, only: mhm_t
  use mo_mrm_container, only: mrm_t

  !> \class   domain_t
  !> \brief   Class for a single mHM domain.
  !> \details Store domain instances with the "target" attribute so components can safely point to their exchange.
  type, public :: domain_t
    type(exchange_t) :: exchange !< the exchange container with all exchanged variables for this domain
    type(input_t) :: input !< the input container providing inputs from file or couplers
    type(meteo_t) :: meteo !< the meteorology container providing processed meteorological data
    type(mpr_t) :: mpr !< the MPR container providing parameter fields for the process containers
    type(mhm_t) :: mhm !< the mHM process container calculating vertical hydrological processes
    type(mrm_t) :: mrm !< the mRM process container for routing related processes
    logical :: is_finished = .false. !< whether the time-loop on this domain is finished
  contains
    procedure :: create => domain_create
    procedure :: set_dims => domain_set_dims
    procedure :: configure => domain_configure
    procedure :: connect => domain_connect
    procedure :: initialize => domain_initialize
    procedure :: update => domain_update
    procedure :: finalize => domain_finalize
  end type domain_t

contains

  !> \brief Create a new domain.
  subroutine domain_create(self, meta_file, domain, cwd)
    class(domain_t), intent(inout), target :: self ! needs "target" so components can safely point to "exchange"
    character(*), intent(in), optional :: meta_file !< file containing the run metadata namelists
    integer(i4), intent(in), optional :: domain !< domain ID of the current domain in the configuration arrays (1 by default)
    character(len=*), intent(in), optional :: cwd !< current working directory to set relative paths
    log_info(*) "CREATE DOMAIN"
    call self%exchange%create(meta_file=meta_file, domain=domain, cwd=cwd)
    ! set exchange pointer in components
    self%input%exchange => self%exchange
    self%meteo%exchange => self%exchange
    self%mpr%exchange => self%exchange
    self%mhm%exchange => self%exchange
    self%mrm%exchange => self%exchange
    call self%set_dims()
  end subroutine domain_create

  !> \brief Set runtime dimensions on all generated namelist configs owned by this domain.
  subroutine domain_set_dims(self)
    class(domain_t), intent(inout), target :: self

    log_info(*) "SET NAMELIST DIMENSIONS"
    call self%exchange%set_dims()
    call self%input%set_dims()
    call self%meteo%set_dims()
    call self%mpr%set_dims()
    call self%mhm%set_dims()
    call self%mrm%set_dims()
  end subroutine domain_set_dims

  !> \brief Configure the domain.
  subroutine domain_configure(self, main_file, para_file, out_file)
    class(domain_t), intent(inout), target :: self ! needs "target" so components can safely point to "exchange"
    character(*), intent(in), optional :: main_file !< file containing the main namelists
    character(*), intent(in), optional :: para_file !< file containing the parameter namelists
    character(*), intent(in), optional :: out_file !< file containing the output namelists
    character(:), allocatable :: domain_main_file

    log_info(*) "CONFIGURE COMPONENTS"
    call self%exchange%configure(main_file=main_file, para_file=para_file)
    if (self%exchange%from_dirs) then
      domain_main_file = trim(self%exchange%config%domain%domain_nmls(self%exchange%domain_id))
      call self%input%configure(domain_main_file)
      if (self%exchange%parameters%mhm_active()) call self%mpr%configure(domain_main_file)
      if (self%exchange%parameters%meteo_active()) call self%meteo%configure(domain_main_file)
      if (self%exchange%parameters%mhm_active()) call self%mhm%configure(domain_main_file, out_file)
      if (self%exchange%parameters%mrm_active()) call self%mrm%configure(domain_main_file, out_file)
    else
      call self%input%configure(main_file)
      if (self%exchange%parameters%mhm_active()) call self%mpr%configure(main_file)
      if (self%exchange%parameters%meteo_active()) call self%meteo%configure(main_file)
      if (self%exchange%parameters%mhm_active()) call self%mhm%configure(main_file, out_file)
      if (self%exchange%parameters%mrm_active()) call self%mrm%configure(main_file, out_file)
    end if
  end subroutine domain_configure

  !> \brief Connect the domain components.
  !> \details Check for dependencies and connect exchanged arrays between components after configuration.
  subroutine domain_connect(self)
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
    class(domain_t), intent(inout), target :: self
    log_info(*) "FINALIZE COMPONENTS"
    call self%input%finalize()
    if (self%exchange%parameters%mhm_active()) call self%mpr%finalize()
    if (self%exchange%parameters%meteo_active()) call self%meteo%finalize()
    if (self%exchange%parameters%mhm_active()) call self%mhm%finalize()
    if (self%exchange%parameters%mrm_active()) call self%mrm%finalize()
  end subroutine domain_finalize

end module mo_domain
