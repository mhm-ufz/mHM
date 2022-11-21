!> \file    wrapper.f90
!> \brief   Module to wrap mHM with f2py to control it with Python.
!> \authors Sebastian Mueller
!> \date    Nov 2022

!> \brief   Python wrapper module to control a mHM model.
module model
  implicit none
contains
  !> \brief Initialize a mHM model.
  subroutine init(namelist_mhm, namelist_mhm_param, namelist_mhm_output, namelist_mrm_output, cwd)
    use mo_mhm_interface, only: mhm_interface_init
    implicit none
    character(*), intent(in) :: namelist_mhm !< path to mHM configuration namelist
    character(*), intent(in) :: namelist_mhm_param !< path to mHM parameter namelist
    character(*), intent(in) :: namelist_mhm_output !< path to mHM output namelist
    character(*), intent(in) :: namelist_mrm_output !< path to mRM output namelist
    character(*), intent(in) :: cwd !< desired working directory
    call mhm_interface_init(namelist_mhm, namelist_mhm_param, namelist_mhm_output, namelist_mrm_output, cwd)
  end subroutine init

  !> \brief Execute a mHM model.
  subroutine run()
    use mo_mhm_interface, only: mhm_interface_run
    implicit none
    call mhm_interface_run
  end subroutine run

  !> \brief Execute a mHM model or an optimization depending on the configuration.
  subroutine run_or_optimize()
    use mo_common_mHM_mRM_variables, only: optimize
    use mo_mhm_interface, only: &
      mhm_interface_run, &
      mhm_interface_run_optimization
    implicit none
    ! RUN OR OPTIMIZE
    if (optimize) then
      call mhm_interface_run_optimization
    else
      ! single mhm run with current settings
      call mhm_interface_run
    end if
  end subroutine run_or_optimize

  !> \brief Get the mHM version.
  subroutine version(ver_string)
    use mo_file, only: mhm_version => version
    implicit none
    character(64), intent(out) :: ver_string
    ver_string = mhm_version
  end subroutine version

  !> \brief Finalize a mHM model.
  subroutine finalize()
    use mo_mhm_interface, only: mhm_interface_finalize
    implicit none
    call mhm_interface_finalize
  end subroutine finalize
end module model

!> \brief   Python wrapper module to control a mHM model run per time step.
module run
  implicit none
contains
  !> \brief Prepare a mHM model run.
  subroutine prepare()
    use mo_mhm_interface_run, only: mhm_interface_run_prepare
    implicit none
    call mhm_interface_run_prepare
  end subroutine prepare

  !> \brief Get the number of domains of the current mHM model run.
  subroutine get_ndomains(n)
    use mo_mhm_interface_run, only: mhm_interface_run_get_ndomains
    implicit none
    integer, intent(out) :: n !< number of domains
    call mhm_interface_run_get_ndomains(n)
  end subroutine get_ndomains

  !> \brief Prepare a certain domain of the current mHM model run.
  subroutine prepare_domain(domain)
    use mo_mhm_interface_run, only: mhm_interface_run_prepare_domain
    implicit none
    integer, intent(in) :: domain !< domain index (1 based and 1 by default)
    !f2py integer :: domain = 1
    call mhm_interface_run_prepare_domain(domain)
  end subroutine prepare_domain

  !> \brief Check if the current mHM model time loop is finished.
  subroutine finished(output)
    use mo_mhm_interface_run, only: mhm_interface_run_finished
    implicit none
    logical, intent(out) :: output !< whether the current time loop finished
    call mhm_interface_run_finished(output)
  end subroutine finished

  !> \brief Do one time-step on the current domain of the current mHM model run.
  subroutine do_time_step()
    use mo_mhm_interface_run, only: mhm_interface_run_do_time_step
    implicit none
    call mhm_interface_run_do_time_step
  end subroutine do_time_step

  !> \brief Write output for the current domain of the current mHM model run.
  subroutine write_output()
    use mo_mhm_interface_run, only: mhm_interface_run_write_output
    implicit none
    call mhm_interface_run_write_output
  end subroutine write_output

  !> \brief Finalize the current domain of the current mHM model run.
  subroutine finalize_domain()
    use mo_mhm_interface_run, only: mhm_interface_run_finalize_domain
    implicit none
    call mhm_interface_run_finalize_domain
  end subroutine finalize_domain

  !> \brief Finalize the current mHM model run.
  subroutine finalize()
    use mo_mhm_interface_run, only: mhm_interface_run_finalize
    implicit none
    call mhm_interface_run_finalize
  end subroutine finalize

  !> \brief Get the current time the current domain of the current mHM model run.
  subroutine current_time(year, month, day, hour)
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer, intent(out) :: year !< current year
    integer, intent(out) :: month !< current month
    integer, intent(out) :: day !< current day
    integer, intent(out) :: hour !< current hour
    year = run_cfg%domainDateTime%year
    month = run_cfg%domainDateTime%month
    day = run_cfg%domainDateTime%day
    hour = run_cfg%domainDateTime%hour
  end subroutine current_time
end module run

!> \brief   Python wrapper module to get internal variables of a mHM model run.
module get
  implicit none
contains
  !> \brief Get the shape of mHM model runoff output.
  subroutine runoff_shape(shp)
    use mo_mrm_global_variables, only: mRM_runoff
    implicit none
    integer, intent(out) :: shp(2) !< 2D shape of the runoff
    shp = shape(mRM_runoff)
  end subroutine runoff_shape

  !> \brief Get the mHM model runoff output.
  subroutine runoff(output, m, n)
    use mo_mrm_global_variables, only: mRM_runoff
    implicit none
    integer :: m !< number of time-steps
    integer :: n !< number of gauges
    real*8, intent(out) :: output(m, n) !< runoff
    output = mRM_runoff
  end subroutine runoff

  !> \brief Get number of unmasked celles on Level-0 of the mHM model.
  subroutine L0_domain_size(n, domain)
    use mo_common_variables, only : level0
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer, intent(out) :: n !< number of unmasked celles
    integer, intent(in) :: domain !< selected domain (0 by default for current domain)
    !f2py integer :: domain = 0
    integer :: iDomain, i
    i = domain
    if ( i == 0 ) i = run_cfg%selected_domain
    iDomain = run_cfg%get_domain_index(i)
    n = level0(iDomain)%nCells
  end subroutine L0_domain_size

  !> \brief Get the shape of Level-0 of the mHM model.
  subroutine L0_domain_shape(shp, domain)
    use mo_common_variables, only : level0
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer, intent(out) :: shp(2) !< shape of Level-0
    integer, intent(in) :: domain !< selected domain (0 by default for current domain)
    !f2py integer :: domain = 0
    integer :: iDomain, i
    i = domain
    if ( i == 0 ) i = run_cfg%selected_domain
    iDomain = run_cfg%get_domain_index(i)
    shp = shape(level0(iDomain)%mask)
  end subroutine L0_domain_shape

  !> \brief Get the mask of Level-0 of the mHM model.
  subroutine L0_domain_mask(mask, n, m, domain)
    use mo_common_variables, only : level0
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer :: m !< number of columns
    integer :: n !< number of rows
    logical, intent(out) :: mask(m, n) !< mask at Level-0
    integer, intent(in) :: domain !< selected domain (0 by default for current domain)
    !f2py integer :: domain = 0
    integer :: iDomain, i
    i = domain
    if ( i == 0 ) i = run_cfg%selected_domain
    iDomain = run_cfg%get_domain_index(i)
    mask = level0(iDomain)%mask
  end subroutine L0_domain_mask

  !> \brief Get the information of Level-0 of the mHM model.
  subroutine L0_domain_info(ncols, nrows, ncells, xll, yll, cell_size, no_data, domain)
    use mo_common_variables, only : level0
    use mo_common_constants, only : nodata_dp
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer, intent(out) :: ncols !< number of columns
    integer, intent(out) :: nrows !< number of rows
    integer, intent(out) :: ncells !< number of cells
    real*8, intent(out) :: xll !< x coordinate of lower left corner
    real*8, intent(out) :: yll !< y coordinate of lower left corner
    real*8, intent(out) :: cell_size !< cell-size
    real*8, intent(out) :: no_data !< no data value
    integer, intent(in) :: domain !< selected domain (0 by default for current domain)
    !f2py integer :: domain = 0
    integer :: iDomain, i
    i = domain
    if ( i == 0 ) i = run_cfg%selected_domain
    iDomain = run_cfg%get_domain_index(i)
    ncols = level0(iDomain)%ncols
    nrows = level0(iDomain)%nrows
    ncells = level0(iDomain)%nCells
    xll = level0(iDomain)%xllcorner
    yll = level0(iDomain)%yllcorner
    cell_size = level0(iDomain)%cellsize
    ! no_data = level0(iDomain)%nodata_value
    no_data = nodata_dp
  end subroutine L0_domain_info

  !> \brief Get a variable on Level-0 of the mHM model.
  subroutine L0_variable(output, n, name, idx)
    use mo_common_run_variables, only : run_cfg
    use mo_common_variables, only : level0, domainMeta
    use mo_mpr_global_variables, only : L0_gridded_LAI
    implicit none
    integer :: n !< size of the variable
    real*8, intent(out) :: output(n) !< the desired variable
    character(*), intent(in) :: name !< name to select the variable
    integer, intent(in) :: idx !< optional index if the variable has multiple layer (1 by default)
    !f2py integer :: idx = 1
    integer :: iDomain, s0, e0

    iDomain = run_cfg%get_domain_index(run_cfg%selected_domain)
    s0 = level0(domainMeta%L0DataFrom(iDomain))%iStart
    e0 = level0(domainMeta%L0DataFrom(iDomain))%iEnd
    select case(name)
      case("L0_GRIDDED_LAI")
        output = L0_gridded_LAI(s0 : e0, idx)
      case default
        print*, "unknown variable: " // name
        stop "get.L0_variable: unknown variable"
    end select
  end subroutine L0_variable

  !> \brief Get number of unmasked celles on Level-1 of the mHM model.
  subroutine L1_domain_size(n, domain)
    use mo_common_variables, only : level1
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer, intent(out) :: n !< number of unmasked celles
    integer, intent(in) :: domain !< selected domain (0 by default for current domain)
    !f2py integer :: domain = 0
    integer :: iDomain, i
    i = domain
    if ( i == 0 ) i = run_cfg%selected_domain
    iDomain = run_cfg%get_domain_index(i)
    n = level1(iDomain)%nCells
  end subroutine L1_domain_size

  !> \brief Get the shape of Level-1 of the mHM model.
  subroutine L1_domain_shape(shp, domain)
    use mo_common_variables, only : level1
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer, intent(out) :: shp(2) !< shape of Level-1
    integer, intent(in) :: domain !< selected domain (0 by default for current domain)
    !f2py integer :: domain = 0
    integer :: iDomain, i
    i = domain
    if ( i == 0 ) i = run_cfg%selected_domain
    iDomain = run_cfg%get_domain_index(i)
    shp = shape(level1(iDomain)%mask)
  end subroutine L1_domain_shape

  !> \brief Get the mask of Level-1 of the mHM model.
  subroutine L1_domain_mask(mask, n, m, domain)
    use mo_common_variables, only : level1
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer :: m !< number of columns
    integer :: n !< number of rows
    logical, intent(out) :: mask(m, n) !< mask at Level-1
    integer, intent(in) :: domain !< selected domain (0 by default for current domain)
    !f2py integer :: domain = 0
    integer :: iDomain, i
    i = domain
    if ( i == 0 ) i = run_cfg%selected_domain
    iDomain = run_cfg%get_domain_index(i)
    mask = level1(iDomain)%mask
  end subroutine L1_domain_mask

  !> \brief Get the information of Level-1 of the mHM model.
  subroutine L1_domain_info(ncols, nrows, ncells, xll, yll, cell_size, no_data, domain)
    use mo_common_variables, only : level1
    use mo_common_constants, only : nodata_dp
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer, intent(out) :: ncols !< number of columns
    integer, intent(out) :: nrows !< number of rows
    integer, intent(out) :: ncells !< number of cells
    real*8, intent(out) :: xll !< x coordinate of lower left corner
    real*8, intent(out) :: yll !< y coordinate of lower left corner
    real*8, intent(out) :: cell_size !< cell-size
    real*8, intent(out) :: no_data !< no data value
    integer, intent(in) :: domain !< selected domain (0 by default for current domain)
    !f2py integer :: domain = 0
    integer :: iDomain, i
    i = domain
    if ( i == 0 ) i = run_cfg%selected_domain
    iDomain = run_cfg%get_domain_index(i)
    ncols = level1(iDomain)%ncols
    nrows = level1(iDomain)%nrows
    ncells = level1(iDomain)%nCells
    xll = level1(iDomain)%xllcorner
    yll = level1(iDomain)%yllcorner
    cell_size = level1(iDomain)%cellsize
    ! no_data = level1(iDomain)%nodata_value
    no_data = nodata_dp
  end subroutine L1_domain_info

  !> \brief Get a variable on Level-1 of the mHM model.
  subroutine L1_variable(output, n, name, idx)
    use mo_common_run_variables, only : run_cfg
    use mo_mpr_global_variables, only : &
      L1_fSealed, & ! 1d
      L1_soilMoistSat ! 2d
    use mo_global_variables, only : &
      L1_inter, & ! 1d
      L1_snowPack, & ! 1d
      L1_soilMoist, & ! 2d
      L1_sealSTW, & ! 1d
      L1_unsatSTW, & ! 1d
      L1_satSTW, & ! 1d
      L1_neutrons, & ! 1d
      L1_pet_calc, & ! 1d
      L1_temp_calc, & ! 1d
      L1_prec_calc, & ! 1d
      L1_aETSoil, & ! 2d
      L1_aETCanopy, & ! 1d
      L1_aETSealed, & ! 1d
      L1_total_runoff, & ! 1d
      L1_runoffSeal, & ! 1d
      L1_fastRunoff, & ! 1d
      L1_slowRunoff, & ! 1d
      L1_baseflow, & ! 1d
      L1_percol, & ! 1d
      L1_infilSoil, & ! 2d
      L1_preEffect ! 1d
    implicit none
    integer :: n !< size of the variable
    real*8, intent(out) :: output(n) !< the desired variable
    character(*), intent(in) :: name !< name to select the variable
    integer, intent(in) :: idx !< optional index if the variable has multiple layer (1 by default)
    !f2py integer :: idx = 1
    select case(name)
      case("L1_FSEALED")
        output = L1_fSealed(run_cfg%s1 : run_cfg%e1, 1, run_cfg%domainDateTime%yId)
      case("L1_FNOTSEALED")
        output = run_cfg%L1_fNotSealed(run_cfg%s1 : run_cfg%e1, 1, run_cfg%domainDateTime%yId)
      case("L1_INTER")
        output = L1_inter(run_cfg%s1 : run_cfg%e1)
      case("L1_SNOWPACK")
        output = L1_snowPack(run_cfg%s1 : run_cfg%e1)
      case("L1_SEALSTW")
        output = L1_sealSTW(run_cfg%s1 : run_cfg%e1)
      case("L1_UNSATSTW")
        output = L1_unsatSTW(run_cfg%s1 : run_cfg%e1)
      case("L1_SATSTW")
        output = L1_satSTW(run_cfg%s1 : run_cfg%e1)
      case("L1_NEUTRONS")
        output = L1_neutrons(run_cfg%s1 : run_cfg%e1)
      case("L1_PET_CALC")
        output = L1_pet_calc(run_cfg%s1 : run_cfg%e1)
      case("L1_TEMP_CALC")
        output = L1_temp_calc(run_cfg%s1 : run_cfg%e1)
      case("L1_PREC_CALC")
        output = L1_prec_calc(run_cfg%s1 : run_cfg%e1)
      case("L1_AETCANOPY")
        output = L1_aETCanopy(run_cfg%s1 : run_cfg%e1)
      case("L1_AETSEALED")
        output = L1_aETSealed(run_cfg%s1 : run_cfg%e1)
      case("L1_TOTAL_RUNOFF")
        output = L1_total_runoff(run_cfg%s1 : run_cfg%e1)
      case("L1_RUNOFFSEAL")
        output = L1_runoffSeal(run_cfg%s1 : run_cfg%e1)
      case("L1_FASTRUNOFF")
        output = L1_fastRunoff(run_cfg%s1 : run_cfg%e1)
      case("L1_SLOWRUNOFF")
        output = L1_slowRunoff(run_cfg%s1 : run_cfg%e1)
      case("L1_BASEFLOW")
        output = L1_baseflow(run_cfg%s1 : run_cfg%e1)
      case("L1_PERCOL")
        output = L1_percol(run_cfg%s1 : run_cfg%e1)
      case("L1_PREEFFECT")
        output = L1_preEffect(run_cfg%s1 : run_cfg%e1)
      case("L1_SOILMOIST")
        output = L1_soilMoist(run_cfg%s1 : run_cfg%e1, idx)
      case("L1_SOILMOIST_VOL")
        output = L1_soilMoist(run_cfg%s1 : run_cfg%e1, idx) &
          / L1_soilMoistSat(run_cfg%s1 : run_cfg%e1, idx, run_cfg%domainDateTime%yId)
      case("L1_SOILMOIST_VOL_ALL")
        output = sum(L1_soilMoist(run_cfg%s1 : run_cfg%e1, :), dim = 2) &
          / sum(L1_soilMoistSat(run_cfg%s1 : run_cfg%e1, :, run_cfg%domainDateTime%yId), dim = 2)
      case("L1_SOILMOISTSAT")
        output = L1_soilMoistSat(run_cfg%s1 : run_cfg%e1, idx, run_cfg%domainDateTime%yId)
      case("L1_AETSOIL")
        output = L1_aETSoil(run_cfg%s1 : run_cfg%e1, idx)
      case("L1_INFILSOIL")
        output = L1_infilSoil(run_cfg%s1 : run_cfg%e1, idx)
      case default
        print*, "unknown variable: " // name
        stop "get.variable: unknown variable"
    end select
  end subroutine L1_variable

  !> \brief Get number of unmasked celles on Level-11 of the mHM model.
  subroutine L11_domain_size(n, domain)
    use mo_mrm_global_variables, only : level11
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer, intent(out) :: n !< number of unmasked celles
    integer, intent(in) :: domain !< selected domain (0 by default for current domain)
    !f2py integer :: domain = 0
    integer :: iDomain, i
    i = domain
    if ( i == 0 ) i = run_cfg%selected_domain
    iDomain = run_cfg%get_domain_index(i)
    n = level11(iDomain)%nCells
  end subroutine L11_domain_size

  !> \brief Get the shape of Level-11 of the mHM model.
  subroutine L11_domain_shape(shp, domain)
    use mo_mrm_global_variables, only : level11
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer, intent(out) :: shp(2) !< shape of Level-11
    integer, intent(in) :: domain !< selected domain (0 by default for current domain)
    !f2py integer :: domain = 0
    integer :: iDomain, i
    i = domain
    if ( i == 0 ) i = run_cfg%selected_domain
    iDomain = run_cfg%get_domain_index(i)
    shp = shape(level11(iDomain)%mask)
  end subroutine L11_domain_shape

  !> \brief Get the mask of Level-11 of the mHM model.
  subroutine L11_domain_mask(mask, n, m, domain)
    use mo_mrm_global_variables, only : level11
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer :: m !< number of columns
    integer :: n !< number of rows
    logical, intent(out) :: mask(m, n) !< mask at Level-11
    integer, intent(in) :: domain !< selected domain (0 by default for current domain)
    !f2py integer :: domain = 0
    integer :: iDomain, i
    i = domain
    if ( i == 0 ) i = run_cfg%selected_domain
    iDomain = run_cfg%get_domain_index(i)
    mask = level11(iDomain)%mask
  end subroutine L11_domain_mask

  !> \brief Get the information of Level-11 of the mHM model.
  subroutine L11_domain_info(ncols, nrows, ncells, xll, yll, cell_size, no_data, domain)
    use mo_mrm_global_variables, only : level11
    use mo_common_constants, only : nodata_dp
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer, intent(out) :: ncols !< number of columns
    integer, intent(out) :: nrows !< number of rows
    integer, intent(out) :: ncells !< number of cells
    real*8, intent(out) :: xll !< x coordinate of lower left corner
    real*8, intent(out) :: yll !< y coordinate of lower left corner
    real*8, intent(out) :: cell_size !< cell-size
    real*8, intent(out) :: no_data !< no data value
    integer, intent(in) :: domain !< selected domain (0 by default for current domain)
    !f2py integer :: domain = 0
    integer :: iDomain, i
    i = domain
    if ( i == 0 ) i = run_cfg%selected_domain
    iDomain = run_cfg%get_domain_index(i)
    ncols = level11(iDomain)%ncols
    nrows = level11(iDomain)%nrows
    ncells = level11(iDomain)%nCells
    xll = level11(iDomain)%xllcorner
    yll = level11(iDomain)%yllcorner
    cell_size = level11(iDomain)%cellsize
    ! no_data = level11(iDomain)%nodata_value
    no_data = nodata_dp
  end subroutine L11_domain_info

  !> \brief Get a variable on Level-11 of the mHM model.
  subroutine L11_variable(output, n, name, idx)
    use mo_common_run_variables, only : run_cfg
    use mo_mrm_global_variables, only : &
      L11_qMod, &
      L11_qOUT, &
      L11_qTIN, &
      L11_qTR
    implicit none
    integer :: n !< size of the variable
    real*8, intent(out) :: output(n) !< the desired variable
    character(*), intent(in) :: name !< name to select the variable
    integer, intent(in) :: idx !< optional index if the variable has multiple layer (1 by default)
    !f2py integer :: idx = 1
    select case(name)
      case("L11_QMOD")
        output = L11_qMod(run_cfg%s11 : run_cfg%e11)
      case("L11_QOUT")
        output = L11_qOUT(run_cfg%s11 : run_cfg%e11)
      case("L11_QTIN")
        output = L11_qTIN(run_cfg%s11 : run_cfg%e11, idx)
      case("L11_QTR")
        output = L11_qTR(run_cfg%s11 : run_cfg%e11, idx)
      case default
        print*, "unknown variable: " // name
        stop "get.L11_variable: unknown variable"
    end select
  end subroutine L11_variable

  !> \brief Get number of unmasked celles on Level-2 of the mHM model.
  subroutine L2_domain_size(n, domain)
    use mo_global_variables, only : level2
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer, intent(out) :: n !< number of unmasked celles
    integer, intent(in) :: domain !< selected domain (0 by default for current domain)
    !f2py integer :: domain = 0
    integer :: iDomain, i
    i = domain
    if ( i == 0 ) i = run_cfg%selected_domain
    iDomain = run_cfg%get_domain_index(i)
    n = level2(iDomain)%nCells
  end subroutine L2_domain_size

  !> \brief Get the shape of Level-2 of the mHM model.
  subroutine L2_domain_shape(shp, domain)
    use mo_global_variables, only : level2
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer, intent(out) :: shp(2) !< shape of Level-2
    integer, intent(in) :: domain !< selected domain (0 by default for current domain)
    !f2py integer :: domain = 0
    integer :: iDomain, i
    i = domain
    if ( i == 0 ) i = run_cfg%selected_domain
    iDomain = run_cfg%get_domain_index(i)
    shp = shape(level2(iDomain)%mask)
  end subroutine L2_domain_shape

  !> \brief Get the mask of Level-2 of the mHM model.
  subroutine L2_domain_mask(mask, n, m, domain)
    use mo_global_variables, only : level2
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer :: m !< number of columns
    integer :: n !< number of rows
    logical, intent(out) :: mask(m, n) !< mask at Level-2
    integer, intent(in) :: domain !< selected domain (0 by default for current domain)
    !f2py integer :: domain = 0
    integer :: iDomain, i
    i = domain
    if ( i == 0 ) i = run_cfg%selected_domain
    iDomain = run_cfg%get_domain_index(i)
    mask = level2(iDomain)%mask
  end subroutine L2_domain_mask

  !> \brief Get the information of Level-2 of the mHM model.
  subroutine L2_domain_info(ncols, nrows, ncells, xll, yll, cell_size, no_data, domain)
    use mo_global_variables, only : level2
    use mo_common_constants, only : nodata_dp
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer, intent(out) :: ncols !< number of columns
    integer, intent(out) :: nrows !< number of rows
    integer, intent(out) :: ncells !< number of cells
    real*8, intent(out) :: xll !< x coordinate of lower left corner
    real*8, intent(out) :: yll !< y coordinate of lower left corner
    real*8, intent(out) :: cell_size !< cell-size
    real*8, intent(out) :: no_data !< no data value
    integer, intent(in) :: domain !< selected domain (0 by default for current domain)
    !f2py integer :: domain = 0
    integer :: iDomain, i
    i = domain
    if ( i == 0 ) i = run_cfg%selected_domain
    iDomain = run_cfg%get_domain_index(i)
    ncols = level2(iDomain)%ncols
    nrows = level2(iDomain)%nrows
    ncells = level2(iDomain)%nCells
    xll = level2(iDomain)%xllcorner
    yll = level2(iDomain)%yllcorner
    cell_size = level2(iDomain)%cellsize
    ! no_data = level2(iDomain)%nodata_value
    no_data = nodata_dp
  end subroutine L2_domain_info

end module get

!> \brief   Python wrapper module to set internal variables of a mHM model run.
module set
  implicit none
contains
  !> \brief Set a variable on Level-0 of the mHM model.
  subroutine L0_variable(input, n, name, idx)
    use mo_common_run_variables, only : run_cfg
    use mo_common_variables, only : level0, domainMeta
    use mo_mpr_global_variables, only : L0_gridded_LAI
    implicit none
    integer :: n !< size of the variable
    real*8, intent(in) :: input(n) !< the variable value
    character(*), intent(in) :: name !< name to select the variable
    integer, intent(in) :: idx !< optional index if the variable has multiple layer (1 by default)
    !f2py integer :: idx = 1
    integer :: iDomain, s0, e0

    iDomain = run_cfg%get_domain_index(run_cfg%selected_domain)
    s0 = level0(domainMeta%L0DataFrom(iDomain))%iStart
    e0 = level0(domainMeta%L0DataFrom(iDomain))%iEnd
    select case(name)
      case("L0_GRIDDED_LAI")
        L0_gridded_LAI(s0 : e0, idx) = input
      case default
        print*, "unknown variable: " // name
        stop "set.L0_variable: unknown variable"
    end select
  end subroutine L0_variable
end module set
