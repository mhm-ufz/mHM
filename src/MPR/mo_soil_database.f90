!> \file mo_soil_database.f90
!> \brief \copybrief mo_soil_database
!> \details \copydetails mo_soil_database

!> \brief Generating soil database from input file.
!> \details This module provides the routines for generating the soil database for mHM from an ASCII input file.
!! One routine \e read_soil_LUT reads a soil LookUpTable, performs some consistency checks and returns an initial
!! soil database.
!! The second routine \e generate_soil_database calculates based on the initial one the proper soil database.
!> \authors Juliane Mai
!> \date Dec 2012
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mpr
MODULE mo_soil_database

  ! This module to provide a soil database for mHM.

  ! Written  Juliane Mai, Dec 2012

  use mo_kind, only : i4, dp
  use mo_message, only: message, error_message
  use mo_string_utils, only : num2str
  use mo_os, only : check_path_isfile

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_soil_LUT             ! Reads the soil LUT file
  PUBLIC :: generate_soil_database    ! Generates the soil database

  ! ------------------------------------------------------------------

  !> \class soil_database_t
  !> \brief Domain-owned soil LUT, horizon controls, and support weights.
  !> \details Soil array representation retained from Robert Schweppe's MPR state module (2017).
  type, public :: soil_database_t
    integer(i4) :: mode = 0_i4 !< Soil database layout selector.
    integer(i4) :: n_layers = 0_i4 !< Model soil horizon count.
    integer(i4) :: n_types = 0_i4 !< LUT soil type count.
    real(dp) :: tillage_depth = 0.0_dp !< Depth of tillage in mm.
    real(dp), allocatable :: horizon_depth(:) !< Model horizon endpoints in mm.
    ! dim1 =  nSoilType (e.g. i=1..72 for BUEK)
    ! dim2 =  the maximum of nHorizons
    ! dim3 =  land cover classes
    ! input data
    integer(i4), dimension(:), allocatable :: id                 !            Soil Id
    integer(i4), dimension(:), allocatable :: nHorizons          !            Number of horizons
    integer(i4), dimension(:), allocatable :: is_present         !            Wether this soil type is present in
    !                                                                !            this domain or not
    real(dp), dimension(:, :), allocatable :: UD                 ! [mm]       Upper Bound of depth
    real(dp), dimension(:, :), allocatable :: LD                 ! [mm]       Lower Bound of depth
    real(dp), dimension(:, :), allocatable :: clay               ! [%]        Clay content
    real(dp), dimension(:, :), allocatable :: sand               ! [%]        Sand content
    real(dp), dimension(:, :), allocatable :: DbM                ! [g/cm2]    Mineral Bulk density
    real(dp), dimension(:, :), allocatable :: depth              ! [mm]       Depth of the soil Horizon
    real(dp), dimension(:), allocatable :: RZdepth            ! [mm]       Total soil depth
    real(dp), dimension(:, :, :), allocatable :: Wd                 ! [1]        Weights of mHM Horizons according to
    !                                                                !            horizons provided in soil database
    integer(i4), dimension(:), allocatable :: nTillHorizons      ! [1]        Number of tillage horizons

    ! derived soil hydraulic properties
    real(dp), dimension(:, :, :), allocatable :: thetaS_Till        ! [1]        Saturated water content of soil horizons
    !                                                                !            tillage depth - f(OM, management)
    real(dp), dimension(:, :), allocatable :: thetaS             ! [1]        Saturated water content of soil horizons
    !                                                                !            after tillage depth
    real(dp), dimension(:, :, :), allocatable :: Db                 ! [g/cm2]    Bulk density, LUC dependent
    !                                                                !            = f( OM, management)
    real(dp), dimension(:, :, :), allocatable :: thetaFC_Till       ! [1]        Field capacity of tillage layers;
    !                                                                !            LUC dependent - f(OM, management)
    real(dp), dimension(:, :), allocatable :: thetaFC            ! [1]        Field capacity of deeper layers
    real(dp), dimension(:, :, :), allocatable :: thetaPW_Till       ! [1]        Permament wilting point of tillage layers;
    !                                                                !            LUC dependent - f(OM, management)
    real(dp), dimension(:, :), allocatable :: thetaPW            ! [1]        Permanent wilting point of deeper layers
    real(dp), dimension(:, :, :), allocatable :: Ks                 ! [cm/d]     Saturated hydaulic conductivity
  contains
    procedure :: destroy => soil_database_destroy
  end type soil_database_t

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        read_soil_LUT

  !    PURPOSE
  !>       \brief Reads the soil LUT file.

  !>       \details Reads the soil LookUpTable file and checks for consistency.

  !    INTENT(IN)
  !>       \param[in] "character(len = *) :: filename" filename of the soil LUT

  !    HISTORY
  !>       \authors Juliane Mai

  !>       \date Dec 2012

  ! Modifications:
  ! Luis Samaniego   Nov 2013 - transform relation op. == -> .eq. etc
  ! Rohini Kumar     Mar 2016 - new variables for handling different soil databases
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine read_soil_LUT(filename, db)

    use mo_constants, only : eps_dp, nodata_dp, nodata_i4
    use mo_mpr_constants, only : nLCover_class


    implicit none

    ! filename of the soil LUT
    character(len = *), intent(in) :: filename !< Soil lookup table path.
    type(soil_database_t), intent(inout) :: db !< Owning domain soil database.
    integer :: soil_unit

    character(len = 10) :: dummy

    integer(i4) :: ios

    integer(i4) :: ii, jj, kk

    integer(i4) :: nR, nH, max_i

    real(dp) :: up, down

    real(dp) :: cly, snd, bd

    real(dp) :: dmin

    integer(i4) :: maxNumberHorizons

    integer(i4) :: minNumberTillHorizons, maxNumberTillHorizons


    SELECT CASE (db%mode)
      ! classical mHM soil database
    CASE(0)
      ios = 0_i4
      !checking whether the file exists
      call check_path_isfile(path = filename, raise=.true.)
      open(newunit=soil_unit, file = filename, status = 'old', iostat = ios)
      read(soil_unit, *) dummy, db%n_types
      dummy = dummy // ''   ! only to avoid warning

      ! allocate space
      allocate(db%Id(db%n_types))
      allocate(db%nHorizons(db%n_types))
      allocate(db%nTillHorizons(db%n_types))
      allocate(db%RZdepth(db%n_types))
      ! initialize
      db%Id(:) = nodata_i4
      db%nHorizons(:) = nodata_i4
      db%nTillHorizons(:) = nodata_i4
      db%RZdepth(:) = nodata_dp

      ! check if max soil_type id is equal to db%n_types
      max_i = 0_i4
      ! initalise total rows to read
      nR = 0_i4
      read(soil_unit, *) dummy
      do while (.NOT. (ios .ne. 0))
        read(soil_unit, *, IOSTAT = ios) ii, jj, up, down, cly, snd, bd

        ! Checks
        if(up .ge. down) then
          call error_message('read_soil_LUT: ERROR occurred: Mixed horizon depths in soil type', &
                  num2str(ii, '(I3)'), ' and horizon no.', num2str(jj, '(I3)'))
        end if
        if(cly .lt. 0.0_dp .OR. cly .gt. 100.0_dp .OR. &
                snd .lt. 0.0_dp .OR. snd .gt. 100.0_dp .OR. &
                bd  .lt. 0.0_dp .OR. bd  .gt.   5.0_dp) then
          call error_message('read_soil_LUT: ERROR occurred: Inappropriate soil properties in soil type', &
                  num2str(ii, '(I3)'), ' and horizon no.', num2str(jj, '(I3)'))
        end if

        ! initalise soil id
        if ( ii > size(db%Id) ) call error_message( &
          "ERROR: nSoil_Types (", num2str(size(db%Id)), &
          ") in soil_classdefinition.txt seems to be to low! Tried to read: ", num2str(ii) &
        )
        db%Id(ii) = ii
        db%nHorizons(ii) = jj
        if( anint(down, dp) .gt. db%RZdepth(ii) ) db%RZdepth(ii) = anint(down, dp)

        nR = nR + 1_i4
        max_i = max(max_i, ii)

      end do

      if ( max_i < db%n_types ) call error_message( &
        "ERROR: nSoil_Types (", num2str(size(db%Id)), &
        ") in soil_classdefinition.txt seems to be to high! Highest read value: ", num2str(max_i) &
      )

      ! initalise minimum root zone depth among all soil types
      dMin = minval(db%RZdepth(:), db%RZdepth(:) .gt. 0.0_dp)

      ! check the tillage depth...(if possible adjust it..)
      if(db%tillage_depth .gt. dMin) then
        call error_message('read_soil_LUT: ERROR occurred: ', raise=.false.)
        call error_message('    Tillage depth is greater than overall minimum total soil depth ', raise=.false.)
        call error_message('    So tillage depth should be at least', num2str(dMin, '(F7.2)'), raise=.false.)
        call error_message('    Please adjust!')
      end if

      ! insert a new tillage soil layer, only in those soil types, in which it is not present
      rewind(soil_unit)
      read(soil_unit, *) dummy
      read(soil_unit, *) dummy

      ! Last row is read twice so, read only upto (nR - 1)
      do ii = 1, nR - 1
        read(soil_unit, *) jj, nH, up, down, cly, snd, bd
        if(db%tillage_depth .gt. anint(up, dp) .and. db%tillage_depth .lt. anint(down, dp)) then
          db%nHorizons(jj) = db%nHorizons(jj) + 1_i4
        end if

        ! identify upto which soil horizon does the tillage depth goes in
        if(db%tillage_depth .gt. anint(up, dp) .and. db%tillage_depth .le. anint(down, dp)) then
          db%nTillHorizons(jj) = nH
        end if
      end do

      maxNumberHorizons = maxval(db%nHorizons(:))

      ! the variables of soilType are all allocated with maximal number of Horizons although not for each of the
      ! db%n_types the array will be used
      ! loops for array(i,j) should be: j=1, nHorizons(i), otherwise a nodata_dp will appear
      allocate(db%UD(db%n_types, maxNumberHorizons))
      allocate(db%LD(db%n_types, maxNumberHorizons))
      allocate(db%depth(db%n_types, maxNumberHorizons))
      allocate(db%clay(db%n_types, maxNumberHorizons))
      allocate(db%sand(db%n_types, maxNumberHorizons))
      allocate(db%dbM(db%n_types, maxNumberHorizons))

      db%UD(:, :) = nodata_dp
      db%LD(:, :) = nodata_dp
      db%depth(:, :) = nodata_dp
      db%clay(:, :) = nodata_dp
      db%sand(:, :) = nodata_dp
      db%dbM(:, :) = nodata_dp

      ! allocate space for other derived soil hydraulic properties ...
      minNumberTillHorizons = minval(db%nTillHorizons(:))
      maxNumberTillHorizons = maxval(db%nTillHorizons(:))

      allocate(db%thetaS_till(db%n_types, 1 : maxNumberTillHorizons, nLCover_class))
      allocate(db%thetaS(db%n_types, minNumberTillHorizons + 1 : maxNumberHorizons))

      allocate(db%thetaFC_till(db%n_types, 1 : maxNumberTillHorizons, nLCover_class))
      allocate(db%thetaFC(db%n_types, minNumberTillHorizons + 1 : maxNumberHorizons))

      allocate(db%thetaPW_till(db%n_types, 1 : maxNumberTillHorizons, nLCover_class))
      allocate(db%thetaPW(db%n_types, minNumberTillHorizons + 1 : maxNumberHorizons))

      allocate(db%Db(db%n_types, 1 : maxNumberTillHorizons, nLCover_class))

      allocate(db%Ks(db%n_types, 1 : maxNumberHorizons, nLCover_class))

      ! Initialize with default nodata_dp
      db%thetaS_till(:, :, :) = nodata_dp
      db%thetaS(:, :) = nodata_dp

      db%thetaFC_till(:, :, :) = nodata_dp
      db%thetaFC(:, :) = nodata_dp

      db%thetaPW_till(:, :, :) = nodata_dp
      db%thetaPW(:, :) = nodata_dp

      db%Db(:, :, :) = nodata_dp
      db%Ks(:, :, :) = nodata_dp

      ! Read again soil properties  from the data base
      rewind(soil_unit)
      read(soil_unit, *) dummy
      read(soil_unit, *) dummy

      ! Last row is read twice so, read only upto (nR - 1)
      kk = 0_i4
      do ii = 1, nR - 1_i4
        read(soil_unit, *) jj, nH, up, down, cly, snd, bd

        ! to avoid numerical errors in PTF
        if(cly .lt. 1.0_dp) cly = 1.0_dp
        if(snd .lt. 1.0_dp) snd = 1.0_dp

        ! Physical consistency
        if((cly + snd) .gt. 100.0_dp) then
          cly = cly / (cly + snd)
          snd = snd / (cly + snd)
        end if

        ! check for an extra tillage horizon (if not exists create a layer)
        if(db%tillage_depth .gt. anint(up, dp) .and. db%tillage_depth .lt. anint(down, dp)) then

          db%UD(jj, nH) = anint(up, dp)
          db%LD(jj, nH) = db%tillage_depth
          db%depth(jj, nH) = db%LD(jj, nH) - db%UD(jj, nH)
          db%clay(jj, nH) = cly
          db%sand(jj, nH) = snd
          db%dbM(jj, nH) = bd

          ! initalise the upper depth to new value... and the increment counter
          up = db%tillage_depth
          kk = 1_i4
        end if

        ! increment nH by one once it encounter the upper loop..
        if(kk .eq. 1_i4) nH = nH + 1_i4

        db%UD(jj, nH) = anint(up, dp)
        db%LD(jj, nH) = anint(down, dp)
        db%depth(jj, nH) = db%LD(jj, nH) - db%UD(jj, nH)
        db%clay(jj, nH) = cly
        db%sand(jj, nH) = snd
        db%dbM(jj, nH) = bd

        ! check for number of soil horizons...
        if(nH .gt. db%nHorizons(jj)) then
          call error_message('read_soil_LUT: ERROR occurred: ', raise=.false.)
          call error_message('    There is something wrong in allocating horizons in soil data base.', raise=.false.)
          call error_message('    Please check in code !')
        end if

        ! initalise the increment counter to zero
        if(nH .eq. db%nHorizons(jj)) kk = 0_i4

      end do
      close(soil_unit)

      ! soil database for the horizon specific case
    CASE(1)

      allocate(db%RZdepth(1))
      allocate(db%UD(1, 1))
      allocate(db%LD(1, 1))
      allocate(db%depth(1, 1))
      allocate(db%thetaS_till(1, 1, 1))
      allocate(db%thetaS(1, 1))
      allocate(db%thetaFC_till(1, 1, 1))
      allocate(db%thetaFC(1, 1))
      allocate(db%thetaPW_till(1, 1, 1))
      allocate(db%thetaPW(1, 1))
      allocate(db%Db(1, 1, 1))
      allocate(db%Ks(1, 1, 1))



      !checking whether the file exists
      call check_path_isfile(path = filename, raise=.true.)
      open(newunit=soil_unit, file = filename, status = 'old', action = 'read')
      read(soil_unit, *) dummy, db%n_types
      dummy = dummy // ''   ! only to avoid warning
      allocate(db%Id  (db%n_types))
      allocate(db%clay(db%n_types, 1))
      allocate(db%sand(db%n_types, 1))
      allocate(db%dbM (db%n_types, 1))
      db%clay(:, :) = nodata_dp
      db%sand(:, :) = nodata_dp
      db%dbM (:, :) = nodata_dp
      read(soil_unit, *) dummy
      do kk = 1, db%n_types
        read(soil_unit, *, iostat=ios) jj, cly, snd, bd
        if ( ios /= 0 ) call error_message( &
          "ERROR: nSoil_Types (", num2str(db%n_types), ") in soil_classdefinition_iFlag_soilDB_1.txt ", &
          "seems to be higher than available soil types!" &
        )
        ! to avoid numerical errors in PTF
        if(cly .lt. 1.0_dp) cly = 1.0_dp
        if(snd .lt. 1.0_dp) snd = 1.0_dp
        db%Id(kk) = jj
        db%clay(kk, 1) = cly
        db%sand(kk, 1) = snd
        db%dbM (kk, 1) = bd
      end do
      close(soil_unit)

      ! assign up to which horizon layer a soil is treated as tillage layer
      !   since our horizon information is well defined for modeling too
      !   this information is uniform across all soils/modeling cells
      ! for compatibility with db%mode == 0, assign
      ! both nHorizons & tillage horizons
      allocate(db%nHorizons(1))
      allocate(db%nTillHorizons(1))
      db%nHorizons(:) = db%n_layers
      db%nTillHorizons(:) = -9
      do kk = 1, db%n_layers
        if(abs(db%horizon_depth(kk) - db%tillage_depth) .lt. eps_dp) then
          db%nTillHorizons(1) = kk
        end if
      end do

      ! check
      if(db%nTillHorizons(1) .eq. -9) then
        ! rarely could happen *** since this is checked in reading of horizons depths only
        ! but is checked here for double confirmation
        call error_message('***ERROR: specification of tillage depths is not confirming', raise=.false.)
        call error_message('          with given depths of soil horizons to be modeled.')
      else
        call message()
        call message('Tillage layers: the tillage horizons are modelled ')
        call message('                upto mHM layers: ', trim(num2str(db%nTillHorizons(1))))
      end if

    CASE DEFAULT
      call error_message('***ERROR: db%mode option given does not exist. Only 0 and 1 is taken at the moment.')

    END SELECT

  end subroutine read_soil_LUT

  ! ------------------------------------------------------------------

  !    NAME
  !        generate_soil_database

  !    PURPOSE
  !>       \brief Generates soil database.

  !>       \details Calculates the proper soil database using the initialized soil database from read_soil_LUT.

  !    HISTORY
  !>       \authors Juliane Mai

  !>       \date Dec 2012

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine generate_soil_database(db)

    use mo_constants, only : nodata_dp, nodata_i4


    implicit none

    type(soil_database_t), intent(inout) :: db !< Owning domain soil database.
    integer(i4) :: ii, jj, kk

    real(dp) :: dmin

    real(dp) :: dpth_f, dpth_t

    integer(i4) :: layer_f, layer_t

    real(dp), parameter :: small = 0.000001_dp

    ! [mm]       soil depth accuracy
    real(dp), parameter :: soil_dAccuracy = 0.5_dp


    SELECT CASE (db%mode)
      ! classical mHM soil database
    CASE(0)
      ! initalise minimum root zone depth among all soil types
      dMin = minval(db%RZdepth(:), db%RZdepth(:) > 0.0_dp)

      ! check
      if (db%n_layers > 1_i4) then
        if (db%horizon_depth(db%n_layers - 1) .ge. dMin) then
        call error_message('generate_soil_database: ERROR occurred: ', raise=.false.)
        call error_message('   The depth of soil Horizons provided for modelling is not appropriate', raise=.false.)
        call error_message('   The global minimum of total soil horizon depth among all soil type is ', num2str(dMin, '(F7.2)'), raise=.false.)
        call error_message('   Adjust your modeling soil horizon depth in this range', raise=.false.)
        call error_message('   OR Increase the soil depth in data base for only those soil types', raise=.false.)
        call error_message('   whose total depth is smaller than your given modeling depth.')
        end if
      end if

      ! allocate and initialise depth weights

      allocate(db%Wd(db%n_types, db%n_layers, maxval(db%nHorizons(:))))
      db%Wd(:, :, :) = 0.0_dp

      ! Process further to estimate weight of each horizons
      ! weightage according to soil depths
      do ii = 1, db%n_types
        db%Wd(ii, :, db%nHorizons(ii) + 1_i4 : maxval(db%nHorizons(:))) = nodata_dp

        ! initalise last horizon depth to model w.r.t to surface
        ! NOTE:: it is different for each soil
        db%horizon_depth(db%n_layers) = db%RZdepth(ii)

        ! Estimate soil properties for each modeling layers
        do jj = 1, db%n_layers

          ! modeling depth ( **from --> to ** )
          ! take into account the depth accuracy [0.5mm, defined in module..]
          dpth_f = 0.0_dp
          if(jj .ne. 1_i4) dpth_f = db%horizon_depth(jj - 1)
          dpth_t = db%horizon_depth(jj) - soil_dAccuracy

          ! identify to which layer of batabase this mHM  horizon lies
          layer_f = nodata_i4
          layer_t = nodata_i4
          do kk = 1, db%nHorizons(ii)
            if(dpth_f .ge. db%UD(ii, kk) .and. dpth_f .le. (db%LD(ii, kk) - soil_dAccuracy)) layer_f = kk
            if(dpth_t .ge. db%UD(ii, kk) .and. dpth_t .le. (db%LD(ii, kk) - soil_dAccuracy)) layer_t = kk
          end do

          ! Check
          if(layer_f .le. 0_i4 .or. layer_t .le. 0_i4) then
            call error_message('generate_soil_database: ERROR occurred: ', raise=.false.)
            call error_message('     Horizon depths to model do not lie in database for soil type', num2str(ii, '(I3)'), raise=.false.)
            call error_message('     Please check!')
          end if
          if(layer_f .gt. layer_t) then
            call error_message('generate_soil_database: ERROR occurred: ', raise=.false.)
            call error_message('     Something is wrong in assignment of modeling soil horizons or', raise=.false.)
            call error_message('     database of soil type ', num2str(ii, '(I3)'), raise=.false.)
            call error_message('     Please check!')
          end if

          ! iF modeling depth of a given horizon falls in a same soil layer
          if(layer_f .eq. layer_t) then
            db%Wd(ii, jj, layer_f) = 1.0_dp

            ! else estimate depth weightage...
          else

            ! for starting layer...
            db%Wd(ii, jj, layer_f) = db%LD(ii, layer_f) - dpth_f

            ! for ending layer
            db%Wd(ii, jj, layer_t) = (dpth_t + soil_dAccuracy) - db%UD(ii, layer_t)

            ! other intermediate layers weight, if exit
            if(layer_t - layer_f .gt. 1_i4) then
              do kk = layer_f + 1, layer_t - 1
                db%Wd(ii, jj, kk) = db%LD(ii, kk) - db%UD(ii, kk)
              end do
            end if

            ! Estimate depth weightage
            if(jj .ne. 1_i4) then
              db%Wd(ii, jj, 1 : db%nHorizons(ii)) = db%Wd(ii, jj, 1 : db%nHorizons(ii)) / &
                      (db%horizon_depth(jj) - db%horizon_depth(jj - 1_i4))
            else
              db%Wd(ii, jj, 1 : db%nHorizons(ii)) = db%Wd(ii, jj, 1 : db%nHorizons(ii)) / &
                      db%horizon_depth(jj)
            end if

            ! Check (small margin for numerical errors)
            if(sum(db%Wd(ii, jj, :), db%Wd(ii, jj, :) .gt. 0.0_dp) .le. 1.0_dp - small .or. &
                    sum(db%Wd(ii, jj, :), db%Wd(ii, jj, :) .gt. 0.0_dp) .ge. 1.0_dp + small) then
              call error_message('generate_soil_database: ERROR occurred: ', raise=.false.)
              call error_message('     Weight assigned for each soil horizons are not correct.', raise=.false.)
              call error_message('     Please check!')
            end if

          end if

        end do
      end do
      ! soil database for the horizon specific case
    CASE(1)
      ! right now nothing is done here
      ! *** reserved for future changes
      allocate(db%Wd(1,1,1))

    CASE DEFAULT
      call error_message('***ERROR: db%mode option given does not exist. Only 0 and 1 is taken at the moment.')
    END SELECT

  end subroutine generate_soil_database

  !> \brief Release all domain-owned soil data and reset its controls.
  subroutine soil_database_destroy(self)
    class(soil_database_t), intent(inout), target :: self !< Soil database to reset.
    type(soil_database_t) :: empty
    select type (self)
    type is (soil_database_t)
      self = empty
    end select
  end subroutine soil_database_destroy
END MODULE mo_soil_database
