!> \file mo_standard_score.f90

!> \brief Routines for calculating the normalization (anomaly)/standard score/z score and the 
!>        deseasonalized (standard score on monthly basis) values of a time series.

!> \details In environmental research often the centralization and standardization are estimated
!>          for characterizing the dynamics of a signal.

!> \author Matthias Zink
!> \date May 2015

MODULE mo_standard_score

  ! This module contains routines for the masked calculation of
  ! the standard_score of a time series (centralized and standarsized time series).

  ! Literature

  ! Written May 2015, Matthias Zink

  ! License
  ! -------
  ! This file is part of the UFZ Fortran library.

  ! The UFZ Fortran library is free software: you can redistribute it and/or modify
  ! it under the terms of the GNU Lesser General Public License as published by
  ! the Free Software Foundation, either version 3 of the License, or
  ! (at your option) any later version.

  ! The UFZ Fortran library is distributed in the hope that it will be useful,
  ! but WITHOUT ANY WARRANTY; without even the implied warranty of
  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  ! GNU Lesser General Public License for more details.

  ! You should have received a copy of the GNU Lesser General Public License
  ! along with the UFZ Fortran library (cf. gpl.txt and lgpl.txt).
  ! If not, see <http://www.gnu.org/licenses/>.

  ! Copyright 2014 Matthias Zink

  USE mo_kind, ONLY: i4, sp, dp

  IMPLICIT NONE

  PUBLIC :: standard_score                 ! standard score of a population
  PUBLIC :: classified_standard_score      ! standard score for classes of a population (e.g. classes=months)

  ! ------------------------------------------------------------------

  !     NAME
  !         standard_score

  !     PURPOSE
  !>         \brief    Calculates the standard score / normalization (anomaly) / z-score.
  !>         \details  In statistics, the standard score is the (signed) number of standard deviations an observation
  !>           or datum is above the mean. Thus, a positive standard score indicates a datum above the mean,        
  !>           while a negative standard score indicates a datum below the mean.
  !>           It is a dimensionless quantity obtained by subtracting the population mean from 
  !>           an individual raw score and then dividing the difference by the population standard deviation.
  !>           This conversion process is called standardizing or normalizing (however, "normalizing" can
  !>           refer to many types of ratios).\n
  !>           Standard scores are also called z-values, z-scores, normal scores, and standardized variables; the use
  !>           of "Z" is because the normal distribution is also known as the "Z distribution". They are most frequently
  !>           used to compare a sample to a standard normal deviate, though they can be defined without assumptions of
  !>           normality (Wikipedia, May 2015).
  !>
  !>          \f[ standard\_score = \frac{x - \mu_x}{\sigma_x} \f]
  !>           where \f$ \mu_x \f$ is the mean of a population \f$ x \f$ and \f$ \sigma_x \f$ its standard deviation.
  !>
  !>           If an optinal mask is given, the calculations are over those locations that correspond to true values in the mask.
  !>           data can be single or double precision. The result will have the same numerical precision.

  !     CALLING SEQUENCE
  !         out = standard_score(data, mask=mask)
  
  !     INDENT(IN)
  !>        \param[in] "real(sp/dp), dimension(:) :: data" data to calculate the standard score for

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !        None

  !     INDENT(IN), OPTIONAL
  !>        \param[in] "logical, dimension(:),optinal :: mask" indication which cells to use for calculation
  !>           If present, only those locations in mask having true values in mask are evaluated.

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None
  
  !     RETURN
  !>        \return real(sp/dp) :: standard_score &mdash; standard score / normalization (anomaly) / z-score
  
  !     RESTRICTIONS
  !         Input values must be floating points.

  !     EXAMPLE
  !         data = (/ 1., 2, 3., -999., 5., 6. /)
  !         out  = standard_score(data, mask=(data >= 0.))
  !         -> see also example in test directory
  
  !     LITERATURE
  !>         \note Richard J. Larsen and Morris L. Marx (2000) An Introduction to Mathematical Statistics and Its
  !>            Applications, Third Edition, ISBN 0-13-922303-7. p. 282. 

  !     HISTORY
  !>         \author Matthias Zink
  !>         \date   May 2015

  INTERFACE standard_score
     MODULE PROCEDURE standard_score_sp, standard_score_dp
  END INTERFACE standard_score

  ! ------------------------------------------------------------------

  !     NAME
  !         classified_standard_score

  !     PURPOSE
  !>         \brief    Calculates the  classified standard score (e.g. classes are months).
  !>         \details  In statistics, the standard score is the (signed) number of standard deviations an observation
  !>           or datum is above the mean. Thus, a positive standard score indicates a datum above the mean,        
  !>           while a negative standard score indicates a datum below the mean.
  !>           It is a dimensionless quantity obtained by subtracting the population mean from 
  !>           an individual raw score and then dividing the difference by the population standard deviation.
  !>           This conversion process is called standardizing or normalizing (however, "normalizing" can
  !>           refer to many types of ratios).\n
  !>           Standard scores are also called z-values, z-scores, normal scores, and standardized variables; the use
  !>           of "Z" is because the normal distribution is also known as the "Z distribution". They are most frequently
  !>           used to compare a sample to a standard normal deviate, though they can be defined without assumptions of
  !>           normality (Wikipedia, May 2015).\n
  !>           In this particular case the standard score is calculated for means and standard deviations derived from
  !>           classes of the time series. Such classes could be for example months. Thus, the output would be a
  !>           deseasonalized time series. 
  !>
  !>          \f[ classified\_standard\_score = \frac{x_i - \mu_{c_{x_i}}}{\sigma_{c_{x_i}}} \f]
  !>           where \f$ x_i \f$ is an element of class \f$ c_{x_i} \f$. \f$ x \f$ is a population, \f$ \mu_{c_{x_i}} \f$
  !>           is the mean of all members of a class \f$ c_{x_i} \f$ and \f$ \sigma_{c_{x_i}} \f$ its standard deviation.
  !>
  !>           If an optinal mask is given, the calculations are over those locations that correspond to true values in the mask.
  !>           data can be single or double precision. The result will have the same numerical precision.

  !     CALLING SEQUENCE
  !         out = classified_standard_score(data, mask=mask)
  
  !     INDENT(IN)
  !>        \param[in] "integer,     dimension(:) :: classes" classes to categorize data (e.g. months)
  !>        \param[in] "real(sp/dp), dimension(:) :: data"    data to calculate the standard score for

  
  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !        None

  !     INDENT(IN), OPTIONAL
  !>        \param[in] "logical, dimension(:), optional :: mask" indication which cells to use for calculation
  !>           If present, only those locations in mask having true values in mask are evaluated.

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !         None
  
  !     RETURN
  !>        \return real(sp/dp) :: classified_standard_score &mdash; classified standard score (e.g. deseasonalized
  !>                                                                                                 time series)
  
  !     RESTRICTIONS
  !         Input values must be floating points.

  !     EXAMPLE
  !         data    = (/ 1., 2, 3., -999., 5., 6. /)
  !         classes = (/ 1,  1, 1,     2,  2 , 2 /)
  !         out  = classified_standard_score(data, classes, mask=(data >= 0.))
  !         -> see also example in test directory
  
  !     LITERATURE
  !        None

  !     HISTORY
  !>         \author Matthias Zink
  !>         \date   May 2015

  INTERFACE classified_standard_score
     MODULE PROCEDURE classified_standard_score_sp, classified_standard_score_dp
  END INTERFACE classified_standard_score
  ! ------------------------------------------------------------------

  PRIVATE

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  FUNCTION standard_score_sp(data, mask)

    use mo_moment, only: average, stddev
    
    implicit none
    
    real(sp), dimension(:),             intent(in) :: data ! data arrau input
    logical,  dimension(:), optional,   intent(in) :: mask ! optional input
    real(sp), dimension(size(data, dim=1))         :: standard_score_sp

    ! local
    logical, dimension(size(data, dim=1))          :: maske
    
    ! check if optional mask matches shape of data
    if (present(mask)) then
       if (size(mask) .ne. size(data)) stop '***Error: standard_score_sp: size(mask) .ne. size(data)'
       maske = mask
    else
       maske(:) = .true.
    endif

    ! check if enough values (>1) are available
    if (count(maske) .LE. 2) stop '***Error: standard_score_sp: less than 2 elements avaiable'
    
    standard_score_sp = ( data(:) - average(data, mask=maske) ) / stddev(data,mask=maske)

  END FUNCTION standard_score_sp

  
  FUNCTION standard_score_dp(data, mask)

    use mo_moment, only: average, stddev
    
    implicit none
    
    real(dp), dimension(:),             intent(in) :: data ! data arrau input
    logical,  dimension(:), optional,   intent(in) :: mask ! optional input
    real(dp), dimension(size(data, dim=1))         :: standard_score_dp

    ! local
    logical, dimension(size(data, dim=1))          :: maske
    
    ! check if optional mask matches shape of data
    if (present(mask)) then
       if (size(mask) .ne. size(data)) stop '***Error: standard_score_dp: size(mask) .ne. size(data)'
       maske = mask
    else
       maske(:) = .true.
    endif

    ! check if enough values (>1) are available
    if (count(maske) .LE. 2) stop '***Error: standard_score_dp: less than 2 elements avaiable'
    
    standard_score_dp = ( data(:) - average(data, mask=maske) ) / stddev(data,mask=maske)

  END FUNCTION standard_score_dp

  ! ------------------------------------------------------------------

  FUNCTION classified_standard_score_sp(data, classes, mask)

    use mo_moment,    only: average, stddev
    use mo_orderpack, only: unista
    
    implicit none
    
    real(sp),    dimension(:),             intent(in) :: data    ! data array with input
    integer,     dimension(:),             intent(in) :: classes ! array indicateing classes
    logical,     dimension(:), optional,   intent(in) :: mask    ! array masking elements of data
    real(sp),    dimension(size(data, dim=1))         :: classified_standard_score_sp

    ! local
    integer(i4)                                       :: iclass, ielem        ! loop variable
    integer(i4)                                       :: number_of_classes    ! number of unique classes in vector
    ! classes
    integer(i4), dimension(size(data, dim=1))         :: unique_classes       ! vector of uniqe classes
    real(sp)                                          :: class_mean           ! mean of class
    real(sp)                                          :: class_stddev         ! standard deviation of class
    logical,     dimension(size(data, dim=1))         :: maske                ! data mask
    logical,     dimension(size(data, dim=1))         :: mask_class_maske     ! combined mask for current class and
    ! maske

    ! check if optional mask matches shape of data
    if (present(mask)) then
       if (size(mask) .ne. size(data)) stop '***Error: classified_standard_score_sp: size(mask) .ne. size(data)'
       maske = mask
    else
       maske(:) = .true.
    endif

    ! check if enough values (>1) are available
    if (count(maske) .LE. 2) stop '***Error: classified_standard_score_sp: less than 2 elements avaiable'

    ! initialization
    classified_standard_score_sp = 0.0_sp
    
    ! write classes to new array for getting unique array elements
    unique_classes = classes
    call unista(unique_classes, number_of_classes) ! (unique arry elements in the 1:number_of_classes 
    !                                              ! indexes of array unique_classes)
    
    ! loop over classes
    do iclass = 1, number_of_classes
       ! calculate mean and standard deviation for class
       mask_class_maske = (maske .AND. (classes==unique_classes(iclass)))
       class_mean   = average(data, mask=mask_class_maske)
       class_stddev =  stddev(data, mask=mask_class_maske)
       ! loop over array elements
       do ielem = 1, size(data, dim=1)
          if (.NOT. mask_class_maske(ielem)) cycle ! skip masked values and other classes
          classified_standard_score_sp(ielem) = ( data(ielem) - class_mean ) / class_stddev
       end do
    end do
    
  END FUNCTION classified_standard_score_sp

  
  FUNCTION classified_standard_score_dp(data, classes, mask)

    use mo_moment,    only: average, stddev
    use mo_orderpack, only: unista
    
    implicit none
    
    real(dp),    dimension(:),             intent(in) :: data    ! data array with input
    integer,     dimension(:),             intent(in) :: classes ! array indicateing classes
    logical,     dimension(:), optional,   intent(in) :: mask    ! array masking elements of data
    real(dp),    dimension(size(data, dim=1))         :: classified_standard_score_dp

    ! local
    integer(i4)                                       :: iclass, ielem        ! loop variable
    integer(i4)                                       :: number_of_classes    ! number of unique classes in vector classes
    integer(i4), dimension(size(data, dim=1))         :: unique_classes       ! vector of uniqe classes
    real(dp)                                          :: class_mean           ! mean of class
    real(dp)                                          :: class_stddev         ! standard deviation of class
    logical,     dimension(size(data, dim=1))         :: maske                ! data mask
    logical,     dimension(size(data, dim=1))         :: mask_class_maske     ! combined mask for current class and maske

    ! check if optional mask matches shape of data
    if (present(mask)) then
       if (size(mask) .ne. size(data)) stop '***Error: classified_standard_score_dp: size(mask) .ne. size(data)'
       maske = mask
    else
       maske(:) = .true.
    endif

    ! check if enough values (>1) are available
    if (count(maske) .LE. 2) stop '***Error: classified_standard_score_dp: less than 2 elements avaiable'

    ! initialization
    classified_standard_score_dp = 0.0_dp
    
    ! write classes to new array for getting unique array elements
    unique_classes = classes
    call unista(unique_classes, number_of_classes) ! (unique arry elements in the 1:number_of_classes 
    !                                              ! indexes of array unique_classes)
    
    ! loop over classes
    do iclass = 1, number_of_classes
       ! calculate mean and standard deviation for class
       mask_class_maske = (maske .AND. (classes==unique_classes(iclass)))
       class_mean   = average(data, mask=mask_class_maske)
       class_stddev =  stddev(data, mask=mask_class_maske)
       ! loop over array elements
       do ielem = 1, size(data, dim=1)
          if (.NOT. mask_class_maske(ielem)) cycle ! skip masked values and other classes
          classified_standard_score_dp(ielem) = ( data(ielem) - class_mean ) / class_stddev
       end do
    end do
    
  END FUNCTION classified_standard_score_dp
  
END MODULE mo_standard_score
