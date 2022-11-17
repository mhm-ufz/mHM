module model
  implicit none
contains
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

  subroutine run()
    use mo_mhm_interface, only: mhm_interface_run
    implicit none
    call mhm_interface_run
  end subroutine run

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

  subroutine version(ver_string)
    use mo_file, only: mhm_version => version
    implicit none
    character(64), intent(out) :: ver_string
    ver_string = mhm_version
  end subroutine version

  subroutine finalize()
    use mo_mhm_interface, only: mhm_interface_finalize
    implicit none
    call mhm_interface_finalize
  end subroutine finalize
end module model

module run
  implicit none
contains
  subroutine prepare()
    use mo_mhm_interface_run, only: mhm_interface_run_prepare
    implicit none
    call mhm_interface_run_prepare
  end subroutine prepare

  subroutine get_ndomains(n)
    use mo_mhm_interface_run, only: mhm_interface_run_get_ndomains
    implicit none
    integer, intent(out) :: n
    call mhm_interface_run_get_ndomains(n)
  end subroutine get_ndomains

  subroutine prepare_domain(domain)
    use mo_mhm_interface_run, only: mhm_interface_run_prepare_domain
    implicit none
    integer, intent(in) :: domain
    !f2py integer :: domain = 1
    call mhm_interface_run_prepare_domain(domain)
  end subroutine prepare_domain

  subroutine finished(output)
    use mo_mhm_interface_run, only: mhm_interface_run_finished
    implicit none
    logical, intent(out) :: output
    call mhm_interface_run_finished(output)
  end subroutine finished

  subroutine do_time_step()
    use mo_mhm_interface_run, only: mhm_interface_run_do_time_step
    implicit none
    call mhm_interface_run_do_time_step
  end subroutine do_time_step

  subroutine write_output()
    use mo_mhm_interface_run, only: mhm_interface_run_write_output
    implicit none
    call mhm_interface_run_write_output
  end subroutine write_output

  subroutine finalize_domain()
    use mo_mhm_interface_run, only: mhm_interface_run_finalize_domain
    implicit none
    call mhm_interface_run_finalize_domain
  end subroutine finalize_domain

  subroutine finalize()
    use mo_mhm_interface_run, only: mhm_interface_run_finalize
    implicit none
    call mhm_interface_run_finalize
  end subroutine finalize

  subroutine current_time(year, month, day, hour)
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer, intent(out) :: year, month, day, hour
    year = run_cfg%domainDateTime%year
    month = run_cfg%domainDateTime%month
    day = run_cfg%domainDateTime%day
    hour = run_cfg%domainDateTime%hour
  end subroutine current_time
end module run

module get
  implicit none
contains
  subroutine runoff_shape(shp)
    use mo_mrm_global_variables, only: mRM_runoff
    implicit none
    integer, intent(out) :: shp(2)
    shp = shape(mRM_runoff)
  end subroutine runoff_shape

  subroutine runoff(output, m, n)
    use mo_mrm_global_variables, only: mRM_runoff
    implicit none
    integer :: m, n
    real*8, intent(out) :: output(m, n)
    output = mRM_runoff
  end subroutine runoff

  subroutine L0_domain_size(n, domain)
    use mo_common_variables, only : level0
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer, intent(out) :: n
    integer, intent(in) :: domain
    !f2py integer :: domain = 0
    integer :: iDomain, i
    i = domain
    if ( i == 0 ) i = run_cfg%selected_domain
    iDomain = run_cfg%get_domain_index(i)
    n = level0(iDomain)%nCells
  end subroutine L0_domain_size

  subroutine L0_domain_shape(shp, domain)
    use mo_common_variables, only : level0
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer, intent(out) :: shp(2)
    integer, intent(in) :: domain
    !f2py integer :: domain = 0
    integer :: iDomain, i
    i = domain
    if ( i == 0 ) i = run_cfg%selected_domain
    iDomain = run_cfg%get_domain_index(i)
    shp = shape(level0(iDomain)%mask)
  end subroutine L0_domain_shape

  subroutine L0_domain_mask(mask, n, m, domain)
    use mo_common_variables, only : level0
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer :: m, n
    logical, intent(out) :: mask(m, n)
    integer, intent(in) :: domain
    !f2py integer :: domain = 0
    integer :: iDomain, i
    i = domain
    if ( i == 0 ) i = run_cfg%selected_domain
    iDomain = run_cfg%get_domain_index(i)
    mask = level0(iDomain)%mask
  end subroutine L0_domain_mask

  subroutine L0_domain_info(ncols, nrows, ncells, xll, yll, cell_size, no_data, domain)
    use mo_common_variables, only : level0
    use mo_common_constants, only : nodata_dp
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer, intent(out) :: ncols, nrows, ncells
    real*8, intent(out) :: xll, yll, cell_size, no_data
    integer, intent(in) :: domain
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

  subroutine L0_variable(output, n, name, idx)
    use mo_common_run_variables, only : run_cfg
    use mo_common_variables, only : level0, domainMeta
    use mo_mpr_global_variables, only : L0_gridded_LAI
    implicit none
    integer :: n
    real*8, intent(out) :: output(n)
    character(*), intent(in) :: name
    integer, intent(in) :: idx
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

  subroutine L1_domain_size(n, domain)
    use mo_common_variables, only : level1
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer, intent(out) :: n
    integer, intent(in) :: domain
    !f2py integer :: domain = 0
    integer :: iDomain, i
    i = domain
    if ( i == 0 ) i = run_cfg%selected_domain
    iDomain = run_cfg%get_domain_index(i)
    n = level1(iDomain)%nCells
  end subroutine L1_domain_size

  subroutine L1_domain_shape(shp, domain)
    use mo_common_variables, only : level1
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer, intent(out) :: shp(2)
    integer, intent(in) :: domain
    !f2py integer :: domain = 0
    integer :: iDomain, i
    i = domain
    if ( i == 0 ) i = run_cfg%selected_domain
    iDomain = run_cfg%get_domain_index(i)
    shp = shape(level1(iDomain)%mask)
  end subroutine L1_domain_shape

  subroutine L1_domain_mask(mask, n, m, domain)
    use mo_common_variables, only : level1
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer :: m, n
    logical, intent(out) :: mask(m, n)
    integer, intent(in) :: domain
    !f2py integer :: domain = 0
    integer :: iDomain, i
    i = domain
    if ( i == 0 ) i = run_cfg%selected_domain
    iDomain = run_cfg%get_domain_index(i)
    mask = level1(iDomain)%mask
  end subroutine L1_domain_mask

  subroutine L1_domain_info(ncols, nrows, ncells, xll, yll, cell_size, no_data, domain)
    use mo_common_variables, only : level1
    use mo_common_constants, only : nodata_dp
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer, intent(out) :: ncols, nrows, ncells
    real*8, intent(out) :: xll, yll, cell_size, no_data
    integer, intent(in) :: domain
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
    integer :: n
    real*8, intent(out) :: output(n)
    character(*), intent(in) :: name
    integer, intent(in) :: idx
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

  subroutine L11_domain_size(n)
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer, intent(out) :: n
    n = run_cfg%e11 - run_cfg%s11 + 1
  end subroutine L11_domain_size

  subroutine L11_domain_shape(shp, domain)
    use mo_mrm_global_variables, only : level11
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer, intent(out) :: shp(2)
    integer, intent(in) :: domain
    !f2py integer :: domain = 0
    integer :: iDomain, i
    i = domain
    if ( i == 0 ) i = run_cfg%selected_domain
    iDomain = run_cfg%get_domain_index(i)
    shp = shape(level11(iDomain)%mask)
  end subroutine L11_domain_shape

  subroutine L11_domain_mask(mask, n, m, domain)
    use mo_mrm_global_variables, only : level11
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer :: m, n
    logical, intent(out) :: mask(m, n)
    integer, intent(in) :: domain
    !f2py integer :: domain = 0
    integer :: iDomain, i
    i = domain
    if ( i == 0 ) i = run_cfg%selected_domain
    iDomain = run_cfg%get_domain_index(i)
    mask = level11(iDomain)%mask
  end subroutine L11_domain_mask

  subroutine L11_domain_info(ncols, nrows, ncells, xll, yll, cell_size, no_data, domain)
    use mo_mrm_global_variables, only : level11
    use mo_common_constants, only : nodata_dp
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer, intent(out) :: ncols, nrows, ncells
    real*8, intent(out) :: xll, yll, cell_size, no_data
    integer, intent(in) :: domain
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

  subroutine L11_variable(output, n, name, idx)
    use mo_common_run_variables, only : run_cfg
    use mo_mrm_global_variables, only : &
      L11_qMod, &
      L11_qOUT, &
      L11_qTIN, &
      L11_qTR
    implicit none
    integer :: n
    real*8, intent(out) :: output(n)
    character(*), intent(in) :: name
    integer, intent(in) :: idx
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

  subroutine L2_domain_size(n, domain)
    use mo_global_variables, only : level2
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer, intent(out) :: n
    integer, intent(in) :: domain
    !f2py integer :: domain = 0
    integer :: iDomain, i
    i = domain
    if ( i == 0 ) i = run_cfg%selected_domain
    iDomain = run_cfg%get_domain_index(i)
    n = level2(iDomain)%nCells
  end subroutine L2_domain_size

  subroutine L2_domain_shape(shp, domain)
    use mo_global_variables, only : level2
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer, intent(out) :: shp(2)
    integer, intent(in) :: domain
    !f2py integer :: domain = 0
    integer :: iDomain, i
    i = domain
    if ( i == 0 ) i = run_cfg%selected_domain
    iDomain = run_cfg%get_domain_index(i)
    shp = shape(level2(iDomain)%mask)
  end subroutine L2_domain_shape

  subroutine L2_domain_mask(mask, n, m, domain)
    use mo_global_variables, only : level2
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer :: m, n
    logical, intent(out) :: mask(m, n)
    integer, intent(in) :: domain
    !f2py integer :: domain = 0
    integer :: iDomain, i
    i = domain
    if ( i == 0 ) i = run_cfg%selected_domain
    iDomain = run_cfg%get_domain_index(i)
    mask = level2(iDomain)%mask
  end subroutine L2_domain_mask

  subroutine L2_domain_info(ncols, nrows, ncells, xll, yll, cell_size, no_data, domain)
    use mo_global_variables, only : level2
    use mo_common_constants, only : nodata_dp
    use mo_common_run_variables, only : run_cfg
    implicit none
    integer, intent(out) :: ncols, nrows, ncells
    real*8, intent(out) :: xll, yll, cell_size, no_data
    integer, intent(in) :: domain
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

module set
  implicit none
contains
  subroutine L0_variable(input, n, name, idx)
    use mo_common_run_variables, only : run_cfg
    use mo_common_variables, only : level0, domainMeta
    use mo_mpr_global_variables, only : L0_gridded_LAI
    implicit none
    integer :: n
    real*8, intent(in) :: input(n)
    character(*), intent(in) :: name
    integer, intent(in) :: idx
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
