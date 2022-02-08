!> \file    mo_mhm_messages.f90
!> \brief   Module for mHM messages.
!> \details \copydetails mo_mhm_messages

!> \brief   Module for mHM messages.
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    Oct 2021
!> \details Write out messages of mHM (startup, checks, status, finish, ect.).
module mo_mhm_messages

    use mo_kind, only: i4
    use mo_message, only: message
    use mo_string_utils, only: num2str, separator

  implicit none

  private

  public :: startup_message
  public :: domain_dir_check_message
  public :: finish_message

contains

  !> \brief write startup message of mHM.
  subroutine startup_message()
    use mo_file, only: &
      version, &
      version_date, &
      file_namelist_mhm, &
      file_namelist_mhm_param, &
      file_defOutput
    use mo_os, only: path_isfile
    !$ use omp_lib, only: OMP_GET_NUM_THREADS
#ifdef NAG
    use f90_unix_dir, only: GETCWD
#endif

    implicit none

    character(4096) :: message_text
    integer(i4), dimension(8) :: datetime
    logical :: compiled_with_openmp = .false.
    character(len=255)  :: cur_work_dir
    logical :: compiled_with_mpi = .false.
    !$ integer(i4) :: n_threads

#ifdef MPI
    compiled_with_mpi = .true.
#endif

    ! check for working dir (optional argument to the executable)
    CALL getcwd(cur_work_dir)

    call message(separator)
    call message('              mHM-UFZ')
    call message()
    call message('    MULTISCALE HYDROLOGIC MODEL')
    call message('           Version: ', trim(version))
    call message('           Date:    ', trim(version_date))
    call message()
    call message('Originally by L. Samaniego & R. Kumar')
    call message(separator)

    call message()
    !$ compiled_with_openmp = .true.
    if (compiled_with_openmp) then
      call message('OpenMP used.')
    else
      call message('Openmp not used.')
    end if
    !$OMP PARALLEL
    !$ n_threads = OMP_GET_NUM_THREADS()
    !$OMP END PARALLEL
    !$ call message('Run with OpenMP with ', trim(num2str(n_threads)), ' threads.')

    call message()
    if (compiled_with_mpi) then
      call message('MPI used.')
    else
      call message('MPI not used.')
    end if

    call message()
    call date_and_time(values = datetime)
    message_text = trim(num2str(datetime(3), '(I2.2)')) // "." // trim(num2str(datetime(2), '(I2.2)')) &
      // "." // trim(num2str(datetime(1), '(I4.4)')) // " " // trim(num2str(datetime(5), '(I2.2)')) &
      // ":" // trim(num2str(datetime(6), '(I2.2)')) // ":" // trim(num2str(datetime(7), '(I2.2)'))
    call message('Start at ', trim(message_text), '.')
    call message('Working directory: ', trim(cur_work_dir))
    call message('Using namelists:')
    call message('     ', trim(file_namelist_mhm))
    call message('     ', trim(file_namelist_mhm_param))
    call message('     ', trim(file_defOutput))
    call message()

    call path_isfile(path = file_namelist_mhm, quiet_ = .true., throwError_ = .true.)
    call path_isfile(path = file_namelist_mhm_param, quiet_ = .true., throwError_ = .true.)
    call path_isfile(path = file_defOutput, quiet_ = .true., throwError_ = .true.)

  end subroutine startup_message

  !> \brief Check input directories for mHM.
  subroutine domain_dir_check_message()
    use mo_check, only: check_dir
    use mo_global_variables, only: &
      dirPrecipitation, &
      dirTemperature, &
      dirReferenceET, &
      dirMinTemperature, &
      dirMaxTemperature, &
      dirNetRadiation, &
      dirabsVapPressure, &
      dirwindspeed
    use mo_common_variables, only: &
      dirMorpho, &
      dirLCover, &
      dirOut, &
      domainMeta, &
      processMatrix

    implicit none

    integer(i4) :: domainID, iDomain

    call message()
    call message('# of domains:         ', trim(num2str(domainMeta%overallNumberOfDomains)))
    call message()
    call message('  Input data directories:')
    do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)
      call message('  --------------')
      call message('      DOMAIN                  ', num2str(domainID, '(I3)'))
      call message('  --------------')
      call check_dir(dirMorpho(iDomain), "Morphological directory:", .false., 4, 30)
      call check_dir(dirLCover(iDomain), "Land cover directory:", .false., 4, 30)
      call check_dir(dirPrecipitation(iDomain), "Precipitation directory:", .false., 4, 30)
      call check_dir(dirTemperature(iDomain), "Temperature directory:", .false., 4, 30)
      select case (processMatrix(5, 1))
        case(-1 : 0) ! PET is input
          call check_dir(dirReferenceET(iDomain), "PET directory:", .false., 4, 30)
        case(1) ! Hargreaves-Samani
          call check_dir(dirMinTemperature(iDomain), "Min. temperature directory:", .false., 4, 30)
          call check_dir(dirMaxTemperature(iDomain), "Max. temperature directory:", .false., 4, 30)
        case(2) ! Priestely-Taylor
          call check_dir(dirNetRadiation(iDomain), "Net radiation directory:", .false., 4, 30)
        case(3) ! Penman-Monteith
          call check_dir(dirNetRadiation(iDomain), "Net radiation directory:", .false., 4, 30)
          call check_dir(dirabsVapPressure(iDomain), "Abs. vap. press. directory:", .false., 4, 30)
          call check_dir(dirwindspeed(iDomain), "Windspeed directory:", .false., 4, 30)
      end select
      call check_dir(dirOut(iDomain), "Output directory:", .true., 4, 30)
      call message()
    end do
    call message()

  end subroutine domain_dir_check_message

  !> \brief Finish message for mHM.
  subroutine finish_message()
    use mo_common_variables, only: domainMeta
    use mo_common_mHM_mRM_variables, only: simPer, nTstepDay

    implicit none

    integer(i4), dimension(8) :: datetime
    integer(i4) :: nTimeSteps
    character(4096) :: message_text

    nTimeSteps = maxval(simPer(1 : domainMeta%nDomains)%julEnd - simPer(1 : domainMeta%nDomains)%julStart + 1) * nTstepDay
    call date_and_time(values = datetime)
    call message()
    message_text = 'Done ' // trim(num2str(nTimeSteps, '(I10)')) // " time steps."
    call message(trim(message_text))
    message_text = trim(num2str(datetime(3), '(I2.2)')) // "." // trim(num2str(datetime(2), '(I2.2)')) &
            // "." // trim(num2str(datetime(1), '(I4.4)')) // " " // trim(num2str(datetime(5), '(I2.2)')) &
            // ":" // trim(num2str(datetime(6), '(I2.2)')) // ":" // trim(num2str(datetime(7), '(I2.2)'))
    call message('Finished at ', trim(message_text), '.')
    call message()
    call message(separator)
    call message('mHM: Finished!')
    call message(separator)

  end subroutine finish_message

end module mo_mhm_messages
