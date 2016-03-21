!> \file mo_soil_database.f90

!> \brief Generating soil database from input file.

!> \details This module provides the routines for generating the soil database for mHM from an ASCII input file.\n
!> One routine \e read_soil_LUT reads a soil LookUpTable, performs some consistency checks and returns an initial soil database.
!> The second routine \e generate_soil_database calculates based on the initial one the proper soil database.

!> \authors Juliane Mai
!> \date Dec 2012

MODULE mo_soil_database

  ! This module to provide a soil database for mHM.

  ! Written  Juliane Mai, Dec 2012

  use mo_kind,         only: i4, dp
  use mo_message,      only: message
  use mo_string_utils, only: num2str

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_soil_LUT             ! Reads the soil LUT file
  PUBLIC :: generate_soil_database    ! Generates the soil database

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !      NAME
  !         read_soil_LUT

  !>        \brief Reads the soil LUT file.

  !>        \details Reads the soil LookUpTable file and checks for consistency.

  !     INDENT(IN)
  !>        \param[in] "character(len=*) :: filename"        filename of the soil LUT
  !>        \param[in] "integer(i4)      :: iFlag_option"    option for which kind of database to read
  
  
  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !>        \param[out] "type(soilType) :: soilDB"          initialized soil database

  !     INDENT(IN), OPTIONAL
  !         None

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !>       \note Soil LUT is an ASCII file.

  !     EXAMPLE
  !         call read_soil_LUT('soil_BUEK_1000_LUT.txt', soilDB)

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Dec 2012
  !         Modified: Luis Samaniego, Nov 2013, transform relation op. == -> .eq. etc
  !           Rohini Kumar,           Mar 2016 - new variables for handling different soil databases

  subroutine read_soil_LUT(filename, iFlag_option, soilDB)
    use mo_constants,        only: eps_dp                    ! epsilon(1.0) in double precision
    use mo_global_variables, only: soilType                  ! Type definition
    use mo_global_variables, only: tillageDepth              ! value is used
    use mo_global_variables, only: nSoilTypes                ! value is set
    use mo_global_variables, only: nSoilHorizons_mHM, HorizonDepth_mHM       ! value is used
    use mo_file,             only: usoil_database
    use mo_mhm_constants,    only: nLCover_class             ! number of land cover classes
    use mo_mhm_constants,    only: nodata_i4                 ! global nodata value (i4)
    use mo_mhm_constants,    only: nodata_dp                 ! global nodata value (dp)

    implicit none

    character(len=*), intent(in)  :: filename       ! filename of the soil LUT
    integer(i4), intent(in)       :: iFlag_option   ! option for which kind of soil database to read
    type(soilType),   intent(out) :: soilDB         ! initial soil database --> ready for further processing

    ! local variables
    character(len=10)        :: dummy
    integer(i4)              :: ios
    integer(i4)              :: ii, jj, kk
    integer(i4)              :: nR, nH
    real(dp)                 :: up, down
    real(dp)                 :: cly, snd, bd
    real(dp)                 :: dmin
    integer(i4)              :: maxNumberHorizons
    integer(i4)              :: minNumberTillHorizons, maxNumberTillHorizons

    SELECT CASE (iFlag_option)
       ! classical mHM soil database
       CASE(0)
          ios = 0_i4
          open(usoil_database, file=filename, status = 'old', iostat = ios)
          read(usoil_database, *) dummy, nSoilTypes
          dummy = dummy//''   ! only to avoid warning

          ! allocate space
          allocate(soilDB%Id(nSoilTypes))
          allocate(soilDB%nHorizons(nSoilTypes))
          allocate(soilDB%nTillHorizons(nSoilTypes))
          allocate(soilDB%RZdepth(nSoilTypes))
          ! initialize
          soilDB%Id(:)            = nodata_i4
          soilDB%nHorizons(:)     = nodata_i4
          soilDB%nTillHorizons(:) = nodata_i4
          soilDB%RZdepth(:)       = nodata_dp

          ! initalise total rows to read
          nR = 0_i4
          read(usoil_database, *) dummy
          do while (.NOT. ( ios .ne. 0 ) )
             read(usoil_database, *, IOSTAT=ios ) ii, jj, up, down, cly, snd, bd

             ! Checks
             if(up .ge. down) then
                call message('read_soil_LUT: ERROR occurred: Mixed horizon depths in soil type', &
                     num2str(ii,'(I3)'), ' and horizon no.', num2str(jj,'(I3)'))
                stop
             end if
             if(  cly .lt. 0.0_dp .OR. cly .gt. 100.0_dp .OR. &
                  snd .lt. 0.0_dp .OR. snd .gt. 100.0_dp .OR. &
                  bd  .lt. 0.0_dp .OR. bd  .gt.   5.0_dp      ) then
                call message('read_soil_LUT: ERROR occurred: Inappropriate soil properties in soil type', &
                     num2str(ii,'(I3)'), ' and horizon no.', num2str(jj,'(I3)'))
                stop
             end if

             ! initalise soil id
             soilDB%Id(ii)        = ii
             soilDB%nHorizons(ii) = jj
             if( anint(down, dp) .gt. soilDB%RZdepth(ii) ) soilDB%RZdepth(ii) = anint(down, dp)

             nR = nR + 1_i4

          end do

          ! initalise minimum root zone depth among all soil types
          dMin = minval(soilDB%RZdepth(:), soilDB%RZdepth(:) .gt. 0.0_dp)

          ! check the tillage depth...(if possible adjust it..)
          if(tillageDepth .gt. dMin) then
             call message('read_soil_LUT: ERROR occurred: ')
             call message('    Tillage depth is greater than overall minimum total soil depth ')
             call message('    So tillage depth should be at least', num2str(dMin,'(F7.2)'))
             call message('    Please adjust!')
             stop
          end if

          ! insert a new tillage soil layer, only in those soil types, in which it is not present
          rewind(usoil_database)
          read(usoil_database, *) dummy
          read(usoil_database, *) dummy

          ! Last row is read twice so, read only upto (nR - 1)
          do ii = 1, nR - 1
             read(usoil_database, *) jj, nH, up, down, cly, snd, bd
             if(tillageDepth .gt. anint(up, dp) .and. tillageDepth .lt. anint(down,dp) ) then
                soilDB%nHorizons(jj) = soilDB%nHorizons(jj) + 1_i4
             end if

             ! identify upto which soil horizon does the tillage depth goes in
             if(tillageDepth .gt. anint(up, dp) .and. tillageDepth .le. anint(down, dp) ) then
                soilDB%nTillHorizons(jj) = nH
             end if
          end do

          maxNumberHorizons = maxval(soilDB%nHorizons(:))

          ! the variables of soilType are all allocated with maximal number of Horizons although not for each of the
          ! nSoilTypes the array will be used
          ! loops for array(i,j) should be: j=1, nHorizons(i), otherwise a nodata_dp will appear
          allocate(soilDB%UD(    nSoilTypes, maxNumberHorizons ) )
          allocate(soilDB%LD(    nSoilTypes, maxNumberHorizons ) )
          allocate(soilDB%depth( nSoilTypes, maxNumberHorizons ) )
          allocate(soilDB%clay(  nSoilTypes, maxNumberHorizons ) )
          allocate(soilDB%sand(  nSoilTypes, maxNumberHorizons ) )
          allocate(soilDB%dbM(   nSoilTypes, maxNumberHorizons ) )

          soilDB%UD(:,:)        = nodata_dp
          soilDB%LD(:,:)        = nodata_dp
          soilDB%depth(:,:)     = nodata_dp
          soilDB%clay(:,:)      = nodata_dp
          soilDB%sand(:,:)      = nodata_dp
          soilDB%dbM(:,:)       = nodata_dp

          ! allocate space for other derived soil hydraulic properties ...
          minNumberTillHorizons = minval(soilDB%nTillHorizons(:))
          maxNumberTillHorizons = maxval(soilDB%nTillHorizons(:))

          allocate(soilDB%thetaS_till( nSoilTypes, 1 : maxNumberTillHorizons, nLCover_class ) )
          allocate(soilDB%thetaS( nSoilTypes, minNumberTillHorizons+1 : maxNumberHorizons ) )

          allocate(soilDB%thetaFC_till( nSoilTypes, 1 : maxNumberTillHorizons, nLCover_class ) )
          allocate(soilDB%thetaFC( nSoilTypes, minNumberTillHorizons+1 : maxNumberHorizons ) )

          allocate(soilDB%thetaPW_till( nSoilTypes, 1 : maxNumberTillHorizons, nLCover_class ) )
          allocate(soilDB%thetaPW( nSoilTypes, minNumberTillHorizons+1 : maxNumberHorizons ) )

          allocate(soilDB%Db( nSoilTypes, 1 : maxNumberTillHorizons, nLCover_class ) )

          allocate(soilDB%Ks( nSoilTypes, 1 : maxNumberHorizons, nLCover_class ) )

          ! Initialize with default nodata_dp
          soilDB%thetaS_till(:,:,:)  = nodata_dp
          soilDB%thetaS(:,:)         = nodata_dp

          soilDB%thetaFC_till(:,:,:) = nodata_dp
          soilDB%thetaFC(:,:)        = nodata_dp

          soilDB%thetaPW_till(:,:,:) = nodata_dp
          soilDB%thetaPW(:,:)        = nodata_dp

          soilDB%Db(:,:,:)           = nodata_dp
          soilDB%Ks(:,:,:)           = nodata_dp

          ! Read again soil properties  from the data base
          rewind(usoil_database)
          read(usoil_database, *) dummy
          read(usoil_database, *) dummy

          ! Last row is read twice so, read only upto (nR - 1)
          kk = 0_i4
          do ii = 1, nR - 1_i4
             read(usoil_database, *) jj, nH, up, down, cly, snd, bd

             ! to avoid numerical errors in PTF
             if(cly .lt. 1.0_dp) cly = 1.0_dp
             if(snd .lt. 1.0_dp) snd = 1.0_dp

             ! Physical consistency
             if((cly + snd) .gt. 100.0_dp) then
                cly = cly / (cly + snd)
                snd = snd / (cly + snd)
             end if

             ! check for an extra tillage horizon (if not exists create a layer)
             if(tillageDepth .gt. anint(up, dp) .and. tillageDepth .lt. anint(down, dp) ) then

                soilDB%UD(jj,nH)    = anint(up, dp)
                soilDB%LD(jj,nH)    = tillageDepth
                soilDB%depth(jj,nH) = soilDB%LD(jj,nH) - soilDB%UD(jj,nH)
                soilDB%clay(jj,nH)  = cly
                soilDB%sand(jj,nH)  = snd
                soilDB%dbM(jj,nH)   = bd

                ! initalise the upper depth to new value... and the increment counter
                up = tillageDepth
                kk = 1_i4
             end if

             ! increment nH by one once it encounter the upper loop..
             if(kk .eq. 1_i4) nH = nH + 1_i4

             soilDB%UD(jj,nH)    = anint(up, dp)
             soilDB%LD(jj,nH)    = anint(down, dp)
             soilDB%depth(jj,nH) = soilDB%LD(jj,nH) - soilDB%UD(jj,nH)
             soilDB%clay(jj,nH)  = cly
             soilDB%sand(jj,nH)  = snd
             soilDB%dbM(jj,nH)   = bd

             ! check for number of soil horizons...
             if(nH .gt. soilDB%nHorizons(jj) ) then
                call message('read_soil_LUT: ERROR occurred: ')
                call message('    There is something wrong in allocating horizons in soil data base.')
                call message('    Please check in code !')
                STOP
             end if

             ! initalise the increment counter to zero
             if(nH .eq. soilDB%nHorizons(jj) ) kk = 0_i4

          end do
          close(usoil_database)
          
       ! soil database for the horizon specific case
       CASE(1)
          open(usoil_database, file=filename, status = 'old', action = 'read')
          read(usoil_database, *) dummy, nSoilTypes
          dummy = dummy//''   ! only to avoid warning
          allocate( soilDB%Id  (nSoilTypes   ) )
          allocate( soilDB%clay(nSoilTypes, 1) )
          allocate( soilDB%sand(nSoilTypes, 1) )
          allocate( soilDB%dbM (nSoilTypes, 1) )
          soilDB%clay(:,:) = nodata_dp
          soilDB%sand(:,:) = nodata_dp
          soilDB%dbM (:,:) = nodata_dp
          read(usoil_database, *) dummy
          do kk = 1, nSoilTypes
             read(usoil_database, *) jj, cly, snd, bd
             ! to avoid numerical errors in PTF
             if(cly .lt. 1.0_dp) cly = 1.0_dp
             if(snd .lt. 1.0_dp) snd = 1.0_dp
             soilDB%Id(kk)   = jj
             soilDB%clay(kk,1) = cly
             soilDB%sand(kk,1) = snd
             soilDB%dbM (kk,1) = bd
          end do
          close(usoil_database)

          ! assign up to which horizon layer a soil is treated as tillage layer
          !   since our horizon information is well defined for modeling too
          !   this information is uniform across all soils/modeling cells
          ! for compatibility with iFlag_option == 0, assign
          ! both nHorizons & tillage horizons
          allocate( soilDB%nHorizons(1)     )
          allocate( soilDB%nTillHorizons(1) )
          soilDB%nHorizons(:)     = nSoilHorizons_mHM
          soilDB%nTillHorizons(:) = -9
          do kk = 1, nSoilHorizons_mHM
             if( abs( HorizonDepth_mHM(kk) - tillageDepth ) .lt. eps_dp  ) then
                soilDB%nTillHorizons(1) = kk
             end if
          end do
         
          ! check 
          if( soilDB%nTillHorizons(1) .eq. -9 ) then
             ! rarely could happen *** since this is checked in reading of horizons depths only
             ! but is checked here for double confirmation
             call message()
             call message('***ERROR: specification of tillage depths is not confirming')
             call message('          with given depths of soil horizons to be modeled.')
             stop
          else
             call message()
             call message('Tillage layers: the tillage horizons are modelled ')
             call message('                upto mHM layers: ', trim(num2str(soilDB%nTillHorizons(1))))
          end if
          
       CASE DEFAULT
          call message()
          call message('***ERROR: iFlag_soilDB option given does not exist. Only 0 and 1 is taken at the moment.')
          stop
          
       END SELECT
          
  end subroutine read_soil_LUT

  ! ------------------------------------------------------------------

  !      NAME
  !         generate_soil_database

  !>        \brief Generates soil database.

  !>        \details Calculates the proper soil database using the initialized soil database from read_soil_LUT.

  !     INDENT(IN)
  !         None

  !     INDENT(INOUT)
  !>        \param[in,out] "type(soilType) :: soilDB"          initialized/ proper soil database
  !>        \param[in] "integer(i4)        :: iFlag_option"    option for which kind of database to read
  !     INDENT(OUT)
  !         None

  !     INDENT(IN), OPTIONAL
  !         None

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         call generate_soil_database(soilDB)

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Dec 2012
  !         Modified, Rohini Kumar &
  !                    Juliane Mai, Apr, 2013 - wieght for no-contributing soil horizons are intialized with nodata value
  !                   Rohini Kumar, Mar, 2016 - new variables for handling different soil databases

  ! ------------------------------------------------------------------

  subroutine generate_soil_database(soilDB, iFlag_option)

    use mo_global_variables, only: soilType                                                     ! Type definition
    use mo_global_variables, only: nSoilHorizons_mHM, HorizonDepth_mHM, nSoilTypes                  ! value is used
    use mo_mhm_constants,    only: nodata_i4                 ! global nodata value (i4)
    use mo_mhm_constants,    only: nodata_dp                 ! global nodata value (dp)

    implicit none

    type(soilType),   intent(inout) :: soilDB
    integer(i4), intent(in)         :: iFlag_option   ! option for which kind of soil database to read

    ! local variables
    integer(i4)              :: ii, jj, kk
    real(dp)                 :: dmin
    real(dp)                 :: dpth_f,  dpth_t
    integer(i4)              :: layer_f, layer_t
    real(dp),    parameter   :: small= 0.000001_dp
    real(dp), parameter      :: soil_dAccuracy=0.5_dp ! [mm]       soil depth accuracy

    SELECT CASE (iFlag_option)
       ! classical mHM soil database
       CASE(0)
          ! initalise minimum root zone depth among all soil types
          dMin = minval(soilDB%RZdepth(:), soilDB%RZdepth(:) > 0.0_dp)

          ! check
          if (HorizonDepth_mHM(nSoilHorizons_mHM -1) .ge. dMin) then
             call message('generate_soil_database: ERROR occurred: ')
             call message('    The depth of soil Horizons provided for modelling is not appropriate')
             call message('    The global minimum of total soil horizon depth among all soil type is ', num2str(dMin,'(F7.2)'))
             call message('    Adjust your modeling soil horizon depth in this range')
             call message('    OR Increase the soil depth in data base for only those soil types')
             call message('    whose total depth is smaller than your given modeling depth.')
             STOP
          end if

          ! allocate and initalise depth weight
          allocate(soilDB%Wd(nSoilTypes, nSoilHorizons_mHM, maxval(soilDB%nHorizons(:)) ))
          soilDB%Wd(:,:,:) = 0.0_dp

          ! Process further to estimate weight of each horizons
          ! weightage according to soil depths
          do ii = 1, nSoilTypes
             soilDB%Wd(ii,:,soilDB%nHorizons(ii)+1_i4:maxval(soilDB%nHorizons(:))) = nodata_dp

             ! initalise last horizon depth to model w.r.t to surface
             ! NOTE:: it is different for each soil
             HorizonDepth_mHM(nSoilHorizons_mHM) = soilDB%RZdepth(ii)

             ! Estimate soil properties for each modeling layers
             do jj = 1, nSoilHorizons_mHM

                ! modeling depth ( **from --> to ** )
                ! take into account the depth accuracy [0.5mm, defined in module..]
                dpth_f = 0.0_dp
                if(jj .ne. 1_i4 ) dpth_f = HorizonDepth_mHM(jj-1)
                dpth_t = HorizonDepth_mHM(jj) - soil_dAccuracy

                ! identify to which layer of batabase this mHM  horizon lies
                layer_f = nodata_i4
                layer_t = nodata_i4
                do kk = 1, soilDB%nHorizons(ii)
                   if( dpth_f .ge. soilDB%UD(ii,kk) .and. dpth_f .le. (soilDB%LD(ii,kk) - soil_dAccuracy) ) layer_f = kk
                   if( dpth_t .ge. soilDB%UD(ii,kk) .and. dpth_t .le. (soilDB%LD(ii,kk) - soil_dAccuracy) ) layer_t = kk
                end do

                ! Check
                if(layer_f .le. 0_i4 .or. layer_t .le. 0_i4) then
                   call message('generate_soil_database: ERROR occurred: ')
                   call message('     Horizon depths to model do not lie in database for soil type', num2str(ii,'(I3)'))
                   call message('     Please check!')
                   STOP
                end if
                if(layer_f .gt. layer_t) then
                   call message('generate_soil_database: ERROR occurred: ')
                   call message('     Something is wrong in assignment of modeling soil horizons or')
                   call message('     database of soil type ', num2str(ii,'(I3)'))
                   call message('     Please check!')
                   STOP
                end if

                ! iF modeling depth of a given horizon falls in a same soil layer
                if(layer_f .eq. layer_t) then
                   soilDB%Wd(ii, jj, layer_f) = 1.0_dp

                   ! else estimate depth weightage...
                else

                   ! for starting layer...
                   soilDB%Wd(ii, jj, layer_f) = soilDB%LD(ii, layer_f) - dpth_f

                   ! for ending layer
                   soilDB%Wd(ii,jj, layer_t) = (dpth_t + soil_dAccuracy) - soilDB%UD(ii,layer_t)

                   ! other intermediate layers weight, if exit
                   if(layer_t - layer_f .gt. 1_i4) then
                      do kk = layer_f + 1, layer_t - 1
                         soilDB%Wd(ii, jj, kk) = soilDB%LD(ii, kk) - soilDB%UD(ii, kk)
                      end do
                   end if

                   ! Estimate depth weightage
                   if(jj .ne. 1_i4) then
                      soilDB%Wd(ii, jj, 1:soilDB%nHorizons(ii)) = soilDB%Wd(ii, jj,1:soilDB%nHorizons(ii)) / &
                           ( HorizonDepth_mHM(jj) -  HorizonDepth_mHM(jj-1_i4) )
                   else
                      soilDB%Wd(ii, jj, 1:soilDB%nHorizons(ii)) = soilDB%Wd(ii, jj,1:soilDB%nHorizons(ii)) / &
                           HorizonDepth_mHM(jj)
                   end if

                   ! Check (small margin for numerical errors)
                   if(  sum(soilDB%Wd(ii, jj, :), soilDB%Wd(ii, jj, :) .gt. 0.0_dp) .le. 1.0_dp-small .or. &
                        sum(soilDB%Wd(ii, jj, :), soilDB%Wd(ii, jj, :) .gt. 0.0_dp) .ge. 1.0_dp+small     ) then
                      call message('generate_soil_database: ERROR occurred: ')
                      call message('     Weight assigned for each soil horizons are not correct.')
                      call message('     Please check!')
                      STOP
                   end if

                end if

             end do
          end do
       ! soil database for the horizon specific case
       CASE(1)
          ! right now nothing is done here
          ! *** reserved for future changes
          return
          
       CASE DEFAULT
          call message()
          call message('***ERROR: iFlag_soilDB option given does not exist. Only 0 and 1 is taken at the moment.')
          stop

       END SELECT

  end subroutine generate_soil_database

END MODULE mo_soil_database
