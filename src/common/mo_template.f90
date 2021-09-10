!>       \file mo_template.f90

!>       \brief Template for future module developments.

!>       \details This module serves as a template for future model developments.
!>       It shows the module structure, the coding style, and documentation.
!>       Please read the \ref style "Coding and documentation style" guide.

!>       \authors Matthias Cuntz, Christoph Schneider

!>       \date Dec 2012

! Modifications:

MODULE mo_template

  ! This module is a template for the UFZ CHS mesoscale hydrologic model mHM.

  ! Written  Matthias Cuntz, Nov 2011
  ! Modified Matthias Cuntz, Nov 2011 - add private
  !                          Nov 2011 - add public
  !          Matthias Cuntz, Christoph Schneider, Dec 2012 - adapted for mHM and doxygen

  ! Always use the number precisions of mo_kind
  ! All USE statements have the only clause
  USE mo_kind, ONLY : i4, sp, dp

  ! Of course
  IMPLICIT NONE

  ! Explicitly make public only the routines, parameters, etc. that shall be provided
  ! Sort alphabetically and give 1-line descriptions.
  PUBLIC :: circum  ! Circumference of circle
  PUBLIC :: mean    ! Average
  PUBLIC :: PI_dp   ! Constant Pi in double precision
  PUBLIC :: PI_sp   ! Constant Pi in single precision

  ! Documentation of routines is in front of the routines.
  ! But documentation of module interface routines must be in front of the interface definition.

  ! Interfaces for single and double precision routines; sort alphabetically
  ! Interfaces have to be in a public section for doxygen
  ! ------------------------------------------------------------------

  !    NAME
  !        mean

  !    PURPOSE
  !>       \brief The average.

  !>       \details Calculates the average value of a vector, i.e. the first moment of a series of numbers:
  !>       \f[ \bar{x} = \frac{1}{N} \sum_{i=1}^N x_i \f]
  !>       If an optinal mask is given, the mean is only over those locations that correspond
  !>       to true values in the mask.
  !>       x can be single or double precision. The result will have the same numerical precision.
  !>       ADDITIONAL INFORMATION
  !>       vec = (/ 1., 2, 3., -999., 5., 6. /)
  !>       m   = mean(vec, mask=(vec >= 0.))
  !>       -> see also example in test directory
  !>       Sokal RR & Rohlf FJ - Biometry: the principle and practice of statistics in biological research,
  !>       Freeman & Co., ISBN 0-7167-2411-1
  !>       Press WH, Teukolsky SA, Vetterling WT, & Flannery BP - Numerical Recipes in Fortran 90 -
  !>       The Art of Parallel Scientific Computing, 2nd Edition, Volume 2 of Fortran Numerical Recipes,
  !>       Cambridge University Press, UK, 1996

  !    RETURN
  !>       \return real(sp/dp) :: mean &mdash; \f$ \bar{x} \f$ average of all elements in vec

  !    HISTORY
  !>       \authors Matthias Cuntz

  !>       \date Nov 2011

  ! Modifications:
  ! Matthias Cuntz Nov 2011 - include mask
  ! Matthias Cuntz Nov 2011 - test size(mask) == size(dat)


  Interface mean
    MODULE PROCEDURE mean_sp, mean_dp
  END INTERFACE mean
  ! Make everything private by default, after the interfaces
  PRIVATE

  ! Public parameters
  !> Constant Pi in double precision
  REAL(dp), PARAMETER :: PI_dp = 3.141592653589793238462643383279502884197_dp
  !> Constant Pi in single precision
  REAL(sp), PARAMETER :: PI_sp = 3.141592653589793238462643383279502884197_sp

  ! Private global parameters (not used, only for demonstration)
  INTEGER(i4), PARAMETER :: iTest = 1

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        circum

  !    PURPOSE
  !>       \brief Circumference of a circle

  !>       \details Calculates the circumference of a circle
  !>       \f[ c = 2 \pi r \f]

  !    INTENT(IN)
  !>       \param[in] "real(dp) :: radius" Radius

  !    RETURN
  !>       \return real(dp) :: circum &mdash; circumference of circle.

  !    HISTORY
  !>       \authors Matthias Cuntz

  !>       \date Dec 2012

  ! Modifications:

  elemental pure function circum(radius)
    implicit none

    ! Radius
    real(dp), intent(in) :: radius

    real(dp) :: circum


    circum = 2.0_dp * pi_dp * radius

  end function circum

  ! ------------------------------------------------------------------

  FUNCTION mean_dp(dat, mask)
    use mo_message, only: error_message

    REAL(dp), DIMENSION(:), INTENT(IN) :: dat

    LOGICAL, DIMENSION(:), OPTIONAL, INTENT(IN) :: mask

    REAL(dp) :: mean_dp

    REAL(dp) :: n

    LOGICAL, DIMENSION(size(dat)) :: maske


    if (present(mask)) then
      if (size(mask) /= size(dat)) call error_message('Error mean_dp: size(mask) /= size(dat)')
      maske = mask
      n = real(count(maske), dp)
    else
      maske(:) = .true.
      n = real(size(dat), dp)
    end if
    if (n <= (1.0_dp + tiny(1.0_dp))) call error_message('mean_dp: n must be at least 2')

    ! Mean
    mean_dp = sum(dat(:), mask = maske) / n

  END FUNCTION mean_dp


  FUNCTION mean_sp(dat, mask)
    use mo_message, only: error_message

    REAL(sp), DIMENSION(:), INTENT(IN) :: dat

    LOGICAL, DIMENSION(:), OPTIONAL, INTENT(IN) :: mask

    REAL(sp) :: mean_sp

    REAL(sp) :: n

    LOGICAL, DIMENSION(size(dat)) :: maske


    if (present(mask)) then
      if (size(mask) /= size(dat)) call error_message('Error mean_sp: size(mask) /= size(dat)')
      maske = mask
      n = real(count(maske), sp)
    else
      maske(:) = .true.
      n = real(size(dat), sp)
    end if
    if (n <= (1.0_sp + tiny(1.0_sp))) call error_message('mean_sp: n must be at least 2')

    ! Mean
    mean_sp = sum(dat(:), mask = maske) / n

  END FUNCTION mean_sp

END MODULE mo_template
