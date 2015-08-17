!> \file mo_spatialsimilarity.f90

!> \brief Routines for bias insensitive comparison of spatial patterns.

!> \details These routines are based on the idea that spatial similarity can be assessed by comparing
!>          the magnitude of neighboring pixels (e.g. is the neighboring pixel larger or smaller).             

!> \author Matthias Zink
!> \date Mar 2013

MODULE mo_spatialsimilarity

  ! This module contains routines for the masked calculation of
  ! spatial similarity measures like NDV

  ! Written Nov 2012, Matthias Zink

  ! License
  ! -------
  ! This file is part of the UFZ Fortran library.

  ! It is NOT released under the GNU Lesser General Public License, yet. 
  ! If you use this routine, please contact Matthias Zink or Juliane Mai. 

  ! Copyright 2012-2015 Matthias Zink and Juliane Mai

  USE mo_kind,   ONLY: i4, sp, dp

  IMPLICIT NONE

  PUBLIC :: NNDV       ! number of neighboring dominating values
  PUBLIC :: PD         ! patter dissimilarity measure

  ! ------------------------------------------------------------------

  !     NAME
  !         NNDV - number of neighboring dominating values

  !     PURPOSE
  !>         \brief    Calculates the number of neighboring dominating values, a measure for spatial dissimilarity.
  !>         \details
  !>             NNDV             = 1 - sum(abs(dominating_neighbors(mat1) - dominating_neighbors(mat2))) / count(mask)
  !>             dominating_neighbors(mat1) = comparison if pixel is larger than its neighbouring values
  !>
  !>            An array element value is compared with its 8 neighbouring cells to check if these cells are larger
  !>            than the array element value. The result is a 3x3 matrix in which larger cells are indicated by a 
  !>            true value. For comparison this is done with both input arrays. The resulting array is the sum of
  !>            substraction of the 3x3 matrices for each of the both arrays. The resulting matrix is afterwards
  !>            normalized to its available neighbors. 
  !>            Furthermore an average over the entire field is calculated. The valid interval of
  !>            the values for NNDV is [0..1]. In which 1 indicates full agreement and 0 full dismatching.  
  !>             <pre>
  !>            EXAMPLE:\n
  !>            mat1 =  | 12 17  1 | , mat2 = |  7  9 12 | 
  !>                    |  4 10 11 |          | 12 11 11 | 
  !>                    | 15  2 20 |          |  5 13  7 | 
  !>            booleans determined for every grid cell following fortran array scrolling 
  !>            i.e. (/col1_row1, col1_row2, col1_row3, col2_row1, .. ,col3_row3/),(/3,3/)
  !>
  !>            comp1 = | FFF FFF FTF, FFF FFF FFF, FTT FFT FFF |
  !>                    | FFF TFT TTF, TFT TFF FTT, TFF FFT FFF |
  !>                    | FFF FFF FFF, TTF TFF TTF, FFF FFF FFF |
  !>
  !>            comp2 = | FFF FFT FTT, FFT FFT FTT, FFF FFF FFF |
  !>                    | FFF FFF FFT, FTF FFT TFF, FFT TFF FFF |
  !>                    | FFF TFF TTF, FFF FFF FFF, TTF TFF FFF |
  !>
  !>           NNDVMatrix =
  !>           abs( count(comp1) - count(comp2) ) = | 1-3, 0-4, 3-0 | = | 2, 4, 3 |
  !>                                                | 4-1, 5-3, 2-2 |   | 3, 2, 0 |
  !>                                                | 0-3, 5-0, 0-3 |   | 3, 5, 3 |
  !>
   !>                                    DISSIMILAR / VALID NEIGH CELLS
  !>           NNDVMatrix / VALID NEIGH CELLS = | 2, 4, 3 | / | 3, 5, 3 |
  !>                                            | 3, 2, 0 |   | 5, 8, 5 |   
  !>                                            | 3, 5, 3 |   | 3, 5, 3 |   
  !>
  !>                                          = | 0.66, 0.80, 1.00 |
  !>                                            | 0.60, 0.25, 0.00 |
  !>                                            | 1.0,  1.00, 1.00 |
  !>
  !>           NNDV = 1 - sum(NNDVMatrix) / count(mask) = 1 - (6.31666666 / 9) = 0.2981
  !>            </pre>
  !>
  !>           If an optinal mask is given, the calculations are over those locations that correspond to true values in the mask.
  !>           mat1 and mat2 can be single or double precision. The result will have the same numerical precision.

  !     CALLING SEQUENCE
  !         out = NNDV(mat1, mat2, mask=mask, valid=valid)
  
  !     INDENT(IN)
  !>        \param[in] "real(sp/dp), dimension(:,:) :: mat1" 2D-array with input numbers
  !>        \param[in] "real(sp/dp), dimension(:,:) :: mat2" 2D-array with input numbers

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !        None

  !     INDENT(IN), OPTIONAL
  !>        \param[in] "logical,dimension(:,:),optinal :: mask" 2D-array of logical values with size(mat1/mat2).
  !>           If present, only those locations in mat1/mat2 having true values in mask are evaluated.

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !>        \param[out] "logical              :: valid"   indicates if the function could determine a valid value
  !>                                                      result can be unvalid if entire mask is .false. for ex.
  !>                                                      in this case PatternDissim is set to 0 (worst case)
  
  !     RETURN
  !>        \return real(sp/dp) :: NNDV &mdash; Number of neighboring dominating values
  
  !     RESTRICTIONS
  !         Input values must be floating points.

  !     EXAMPLE
  !         mat1 = reshape(/ 1., 2, 3., -999., 5., 6. /, (/3,3/))
  !         mat2 = reshape(/ 1., 2, 3., -999., 5., 6. /, (/3,3/))
  !         out  = NNDV(mat1, mat2, mask=(mat1 >= 0. .and. mat2 >= 0.))
  !         -> see also example in test directory
  
  !     LITERATURE
  !>         \note routine based on algorithm by Luis Samaniego 2009

  !     HISTORY
  !>         \author Matthias Zink
  !>         \date   Nov 2012
  !          update  May 2015 created documentation
  INTERFACE NNDV                  
     MODULE PROCEDURE NNDV_sp, NNDV_dp
  END INTERFACE NNDV

  ! ------------------------------------------------------------------

  !     NAME
  !         PD

  !     PURPOSE
  !>         \brief Calculates pattern dissimilarity (PD) measure
  !>         \details
  !>             PD             = 1 - sum(dissimilarity(mat1, mat2)) / count(mask)
  !>             dissimilarity(mat1, mat2) = comparison if pixel is larger than its neighbouring values
  !>
  !>            An array element value is compared with its 8 neighbouring cells to check if these cells are larger
  !>            than the array element value. The result is a 3x3 matrix in which larger cells are indicated by a 
  !>            true value. For comparison this is done with both input arrays. The resulting array is the sum of
  !>            xor values of the 3x3 matrices for each of the both arrays.  This means only neighbourhood comparisons
  !>            which are different in the 2 matrices are counted. This resulting matrix is afterwards normalized to its
  !>            available neighbors. Furthermore an average over the entire field is calculated. The valid interval of
  !>            the values for PD is [0..1]. In which 1 indicates full agreement and 0 full dismatching.
  !>
  !>             <pre>
 !>            EXAMPLE:\n
  !>            mat1 =  | 12 17  1 | , mat2 = |  7  9 12 | 
  !>                    |  4 10 11 |          | 12 11 11 | 
  !>                    | 15  2 20 |          |  5 13  7 |
  !>
  !>            booleans determined for every grid cell following fortran array scrolling 
  !>            i.e. (/col1_row1, col1_row2, col1_row3, col2_row1, .. ,col3_row3/),(/3,3/)
  !>
  !>            comp1 = | FFF FFF FTF, FFF FFF FFF, FTT FFT FFF |
  !>                    | FFF TFT TTF, TFT TFF FTT, TFF FFT FFF |
  !>                    | FFF FFF FFF, TTF TFF TTF, FFF FFF FFF |
  !>
  !>            comp2 = | FFF FFT FTT, FFT FFT FTT, FFF FFF FFF |
  !>                    | FFF FFF FFT, FTF FFT TFF, FFT TFF FFF |
  !>                    | FFF TFF TTF, FFF FFF FFF, TTF TFF FFF |
  !>
  !>
  !>            xor=neq = | FFF FFT FFT, FFT FFT FTT, FTT FFT FFF |
  !>                      | FFF TFT TTT, TTT TFT TTT, TFT TFT FFF |
  !>                      | FFF TFF TTF, TTF TFF TTF, TTF TFF FFF |               
  !>  
  !>                        DISSIMILAR / VALID NEIGH CELLS
  !>            PDMatrix = | 2, 4, 3 | / | 3, 5, 3 | = | 0.66, 0.80, 1.00 |
  !>                       | 5, 8, 4 |   | 5, 8, 5 |   | 1.00, 1.00, 0.80 |
  !>                       | 3, 5, 3 |   | 3, 5, 3 |   | 1.00, 1.00, 1.00 |
  !>  
  !>           PD = 1 - sum(PDMatrix) / count(mask) = 1 - (8.2666666 / 9) = 0.08148
  !>            </pre>
  !>
  !>           If an optinal mask is given, the calculations are over those locations that correspond to true values in the mask.
  !>           mat1 and mat2 can be single or double precision. The result will have the same numerical precision.

  !     CALLING SEQUENCE
  !         out = PD(mat1, mat2, mask=mask, valid=valid)
  
  !     INDENT(IN)
  !>        \param[in] "real(sp/dp), dimension(:,:) :: mat1" 2D-array with input numbers
  !>        \param[in] "real(sp/dp), dimension(:,:) :: mat2" 2D-array with input numbers

  !     INDENT(INOUT)
  !         None

  !     INDENT(OUT)
  !         None

  !     INDENT(IN), OPTIONAL
  !>        \param[in] "logical,dimension(:,:),optinal :: mask" 2D-array of logical values with size(mat1/mat2)
  !>           If present, only those locations in mat1/mat2 having true values in mask are evaluated.

  !     INDENT(INOUT), OPTIONAL
  !         None

  !     INDENT(OUT), OPTIONAL
  !>        \param[out] "logical,optinal :: valid"   indicates if the function could determine a valid value
  !>                                                      result can be unvalid if entire mask is .false. for ex.
  !>                                                      in this case PD is set to 0 (worst case)

  !     RETURN
  !>        \return real(sp/dp) :: PD &mdash; pattern dissimilarity measure 
  
  !     RESTRICTIONS
  !         Input values must be floating points.

  !     EXAMPLE
  !         mat1 = reshape(/ 1., 2, 3., -999., 5., 6. /, (/3,3/))
  !         mat2 = reshape(/ 1., 2, 3., -999., 5., 6. /, (/3,3/))
  !         out  = PD(mat1, mat2, mask=(mat1 >= 0. .and. mat2 >= 0.))
  !         -> see also example in test directory

  !     LITERATURE
  !          None         

  !     HISTORY
  !>         \author Matthias Zink and Juliane Mai
  !>         \date   Jan 2013
  INTERFACE PD                  
     MODULE PROCEDURE PD_sp, PD_dp
  END INTERFACE PD

  ! ------------------------------------------------------------------

CONTAINS
  
  FUNCTION NNDV_sp(mat1, mat2, mask, valid)

    IMPLICIT NONE

    REAL(sp),    DIMENSION(:,:),                                            INTENT(IN)  :: mat1, mat2
    LOGICAL,     DIMENSION(:,:),                                  OPTIONAL, INTENT(IN)  :: mask
    LOGICAL,                                                      OPTIONAL, INTENT(OUT) :: valid
    REAL(sp)                                                                            :: NNDV_sp

    INTEGER(i4)                                                                         :: iCo, iRo
    INTEGER(i4)                                                                         :: noValidPixels
    INTEGER(i4), DIMENSION(size(shape(mat1)) )                                          :: shapemask
    INTEGER(i4), DIMENSION(size(mat1, dim=1), size(mat1, dim=2))                        :: validcount
    REAL(sp),    DIMENSION(size(mat1, dim=1)+2_i4, size(mat1, dim=2)+2_i4)              :: bufferedMat1, bufferedMat2
    REAL(sp)   , DIMENSION(size(mat1, dim=1), size(mat1, dim=2))                        :: NNDVMatrix              
    LOGICAL,     DIMENSION(size(mat1, dim=1)+2_i4, size(mat1, dim=2)+2_i4)              :: maske
    LOGICAL,     DIMENSION(size(mat1, dim=1), size(mat1, dim=2))                        :: validmaske

    ! check if input has all the same dimensions
    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(mat1)
    end if
    !
    if (any(shape(mat1) .NE. shape(mat2)))  &
         stop 'NNDV_sp: shapes of input matrix 1 and input matrix 2 are not matching'
    if (any(shape(mat1) .NE. shapemask))    &
         stop 'NNDV_sp: shapes of input matrices and mask are not matching'
    !
    ! additional 2 rows and 2 cols added for checking the border cells without crashing the search agorithm
    ! so the search windows can always cover 9 cells even if it is checking a border cell
    ! buffer rows are initialized as false values within mask, i.e. they are not considered for the calculation
    ! of the criteria
    !
    ! initialize mask with default=.false.
    maske = .false.
    if (present(mask)) then
       maske(2:(size(maske,dim=1)-1_i4), 2:(size(maske,dim=2)-1_i4)) = mask
    else
       maske(2:(size(maske,dim=1)-1_i4), 2:(size(maske,dim=2)-1_i4)) = .true.
    end if
    !
    ! initialize bufferedMat1 & bufferedMat2 
    bufferedMat1 = 0.0_sp
    bufferedMat2 = 0.0_sp
    bufferedMat1(2:(size(maske,dim=1)-1), 2:(size(maske,dim=2)-1)) = mat1
    bufferedMat2(2:(size(maske,dim=1)-1), 2:(size(maske,dim=2)-1)) = mat2
    !
    ! initialize NNDV
    NNDVMatrix = 0.0_sp
    !
    NNDV_sp = 0.0_sp
    do iCo = 2_i4, size(bufferedMat1, dim=2) - 1
       do iRo = 2_i4, size(bufferedMat1, dim=1) - 1
          if (.NOT. maske(iRo,iCo)) cycle
          NNDVMatrix(iRo-1_i4,iCo-1_i4) =   &
               real(    &
                  abs(  &
                     count((bufferedMat1(iRo-1:iRo+1 , iCo-1:iCo+1) - bufferedMat1(iRo,iCo) > epsilon(0.0_sp)) .AND. &
                           (maske(iRo-1:iRo+1 , iCo-1:iCo+1))) -                                                     &
                     count((bufferedMat2(iRo-1:iRo+1 , iCo-1:iCo+1) - bufferedMat2(iRo,iCo) > epsilon(0.0_sp)) .AND. &
                           (maske(iRo-1:iRo+1 , iCo-1:iCo+1)))                                                       &
                  ),    &
               sp)
          ! count - 1 to exclude referendce cell (iRo, iCo)
          validcount(iRo-1_i4,iCo-1_i4) = count(maske(iRo-1_i4:iRo+1_i4 , iCo-1_i4:iCo+1_i4)) - 1_i4
        end do
    end do
    !
    ! normalize every pixel to number of valid neighbouring cells (defined by maske) [0..8] --> [0..1]
    validmaske = (maske(2:size(maske, dim=1) - 1_i4, 2:size(maske, dim=2) - 1_i4) .and. (validcount > 0_i4))
    noValidPixels = count(validmaske)
    if (noValidPixels .GT. 0_i4) then
       NNDVMatrix = merge(NNDVMatrix/real(validcount, sp),NNDVMatrix, validmaske)
       ! average over all pixels
       NNDV_sp = 1.0_sp - sum(NNDVMatrix, mask=validmaske) / noValidPixels
       if (present(valid)) valid = .TRUE.
    ! case if maske is full of .false. or no valid neighbouring cells available for every pixel (validcount(:,:)=0)
    else
       NNDV_sp = 0.0_sp
       if (present(valid)) valid = .FALSE.
    end if
    !
  END FUNCTION NNDV_sp

 FUNCTION NNDV_dp(mat1, mat2, mask, valid)

    IMPLICIT NONE

    REAL(dp),    DIMENSION(:,:),                                            INTENT(IN)  :: mat1, mat2
    LOGICAL,     DIMENSION(:,:),                                  OPTIONAL, INTENT(IN)  :: mask
    LOGICAL,                                                      OPTIONAL, INTENT(OUT) :: valid
    REAL(dp)                                                                            :: NNDV_dp

    INTEGER(i4)                                                                         :: iCo, iRo
    INTEGER(i4)                                                                         :: noValidPixels
    INTEGER(i4), DIMENSION(size(shape(mat1)) )                                          :: shapemask
    INTEGER(i4), DIMENSION(size(mat1, dim=1), size(mat1, dim=2))                        :: validcount
    REAL(dp),    DIMENSION(size(mat1, dim=1)+2_i4, size(mat1, dim=2)+2_i4)              :: bufferedMat1, bufferedMat2
    REAL(dp)   , DIMENSION(size(mat1, dim=1), size(mat1, dim=2))                        :: NNDVMatrix              
    LOGICAL,     DIMENSION(size(mat1, dim=1)+2_i4, size(mat1, dim=2)+2_i4)              :: maske
    LOGICAL,     DIMENSION(size(mat1, dim=1), size(mat1, dim=2))                        :: validmaske

    ! check if input has all the same dimensions
    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(mat1)
    end if
    !
    if (any(shape(mat1) .NE. shape(mat2)))  &
         stop 'NNDV_dp: shapes of input matrix 1 and input matrix 2 are not matching'
    if (any(shape(mat1) .NE. shapemask))    &
         stop 'NNDV_dp: shapes of input matrices and mask are not matching'
    !
    ! additional 2 rows and 2 cols added for checking the border cells without crashing the search agorithm
    ! so the search windows can always cover 9 cells even if it is checking a border cell
    ! buffer rows are initialized as false values within mask, i.e. they are not considered for the calculation
    ! of the criteria
    !
    ! initialize mask with default=.false.
    maske = .false.
    if (present(mask)) then
       maske(2:(size(maske,dim=1)-1_i4), 2:(size(maske,dim=2)-1_i4)) = mask
    else
       maske(2:(size(maske,dim=1)-1_i4), 2:(size(maske,dim=2)-1_i4)) = .true.
    end if
    !
    ! initialize bufferedMat1 & bufferedMat2 
    bufferedMat1 = 0.0_dp
    bufferedMat2 = 0.0_dp
    bufferedMat1(2:(size(maske,dim=1)-1), 2:(size(maske,dim=2)-1)) = mat1
    bufferedMat2(2:(size(maske,dim=1)-1), 2:(size(maske,dim=2)-1)) = mat2
    !
    ! initialize NNDV
    NNDVMatrix = 0.0_dp
    !
    NNDV_dp = 0.0_dp
    do iCo = 2_i4, size(bufferedMat1, dim=2) - 1
       do iRo = 2_i4, size(bufferedMat1, dim=1) - 1
          if (.NOT. maske(iRo,iCo)) cycle
          NNDVMatrix(iRo-1_i4,iCo-1_i4) =   &
               real(    &
                  abs(  &
                     count((bufferedMat1(iRo-1:iRo+1 , iCo-1:iCo+1) - bufferedMat1(iRo,iCo) > epsilon(0.0_dp)) .AND. &
                           (maske(iRo-1:iRo+1 , iCo-1:iCo+1))) -                                                     &
                     count((bufferedMat2(iRo-1:iRo+1 , iCo-1:iCo+1) - bufferedMat2(iRo,iCo) > epsilon(0.0_dp)) .AND. &
                           (maske(iRo-1:iRo+1 , iCo-1:iCo+1)))                                                       &
                  ),    &
               dp)
          ! count - 1 to exclude referendce cell (iRo, iCo)
          validcount(iRo-1_i4,iCo-1_i4) = count(maske(iRo-1_i4:iRo+1_i4 , iCo-1_i4:iCo+1_i4)) - 1_i4
        end do
    end do
    !
    ! normalize every pixel to number of valid neighbouring cells (defined by maske) [0..8] --> [0..1]
    validmaske = (maske(2:size(maske, dim=1) - 1_i4, 2:size(maske, dim=2) - 1_i4) .and. (validcount > 0_i4))
    noValidPixels = count(validmaske)
    if (noValidPixels .GT. 0_i4) then
       NNDVMatrix = merge(NNDVMatrix/real(validcount, dp),NNDVMatrix, validmaske)
       ! average over all pixels
       NNDV_dp = 1.0_dp - sum(NNDVMatrix, mask=validmaske) / noValidPixels
       if (present(valid)) valid = .TRUE.
    ! case if maske is full of .false. or no valid neighbouring cells available for every pixel (validcount(:,:)=0)
    else
       NNDV_dp = 0.0_dp
       if (present(valid)) valid = .FALSE.
    end if
    !
  END FUNCTION NNDV_dp

  ! ----------------------------------------------------------------------------------------------------------------

  FUNCTION PD_sp(mat1, mat2, mask, valid)

    IMPLICIT NONE

    REAL(sp),    DIMENSION(:,:),                                            INTENT(IN)  :: mat1, mat2
    LOGICAL,     DIMENSION(:,:)                                 , OPTIONAL, INTENT(IN)  :: mask
    LOGICAL,                                                      OPTIONAL, INTENT(OUT) :: valid
    REAL(sp)                                                                            :: PD_sp

    INTEGER(i4)                                                                         :: iCo, iRo
    INTEGER(i4)                                                                         :: noValidPixels
    INTEGER(i4), DIMENSION(size(shape(mat1)) )                                          :: shapemask
    INTEGER(i4), DIMENSION(size(mat1, dim=1), size(mat1, dim=2))                        :: validcount
    REAL(sp),    DIMENSION(size(mat1, dim=1)+2_i4, size(mat1, dim=2)+2_i4)              :: bufferedMat1, bufferedMat2
    REAL(sp),    DIMENSION(size(mat1, dim=1), size(mat1, dim=2))                        :: PDMatrix
    LOGICAL,     DIMENSION(size(mat1, dim=1)+2_i4, size(mat1, dim=2)+2_i4)              :: maske
    LOGICAL,     DIMENSION(size(mat1, dim=1), size(mat1, dim=2))                        :: validmaske
    ! check if input has all the same dimensions
    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(mat1)
    end if
    !
    if (any(shape(mat1) .NE. shape(mat2)))  &
         stop 'PD_sp: shapes of input matrix 1 and input matrix 2 are not matching'
    if (any(shape(mat1) .NE. shapemask))    &
         stop 'PD_sp: shapes of input matrices and mask are not matching'
    !
    ! additional 2 rows and 2 cols added for checking the border cells without crashing the search agorithm
    ! so the search windows can always cover 9 cells even if it is checking a border cell
    ! buffer rows are initialized as false values within mask, i.e. they are not considered for the calculation
    ! of the criteria
    !
    ! initialize mask with default=.false.
    maske = .false.
    if (present(mask)) then
       maske(2:(size(maske,dim=1)-1_i4), 2:(size(maske,dim=2)-1_i4)) = mask
    else
       maske(2:(size(maske,dim=1)-1_i4), 2:(size(maske,dim=2)-1_i4)) = .true.
    end if
    !
    ! initialize bufferedMat1 & bufferedMat2 
    bufferedMat1 = 0.0_sp
    bufferedMat2 = 0.0_sp
    bufferedMat1(2:(size(maske,dim=1)-1), 2:(size(maske,dim=2)-1)) = mat1
    bufferedMat2(2:(size(maske,dim=1)-1), 2:(size(maske,dim=2)-1)) = mat2
    !
    ! initialize PD
    PDMatrix = 0.0_sp
    do iCo = 2_i4, size(bufferedMat1, dim=2)-1_i4
       do iRo = 2_i4, size(bufferedMat1, dim=1)-1_i4
          ! no calculation at unmasked values
          if (.NOT. maske(iRo,iCo)) cycle
          ! .NEQV. is the fortran standard for .XOR. 
          ! result is written to -1 column and row because of buffer
          PDMatrix(iRo-1_i4,iCo-1_i4) = &
            real( &   
              count( & 
                     ( &
                        ! determine higher neighbouring values in mat1
                        (bufferedMat1(iRo-1_i4:iRo+1_i4, iCo-1_i4:iCo+1_i4) -                         &
                         bufferedMat1(iRo,iCo)                              > epsilon(0.0_sp)) .NEQV. &
                        ! determine higher neighbouring values in mat2                                &
                        (bufferedMat2(iRo-1_i4:iRo+1_i4, iCo-1_i4:iCo+1_i4) -                         &
                         bufferedMat2(iRo,iCo)                              > epsilon(0.0_sp))        &
                     ) &
                  ! exclude unmasked values
                  .and. (maske(iRo-1_i4:iRo+1_i4 , iCo-1_i4:iCo+1_i4)) &
                  ), &
            sp )
          ! count - 1 to exclude reference cell / center pixel (iRo, iCo)
          validcount(iRo-1_i4,iCo-1_i4) = count(maske(iRo-1_i4:iRo+1_i4 , iCo-1_i4:iCo+1_i4)) - 1_i4
          !
       end do
    end do
    !
    ! normalize every pixel to number of valid neighbouring cells (defined by maske) [0..8] --> [0..1]
    validmaske = (maske(2:size(maske, dim=1) - 1_i4, 2:size(maske, dim=2) - 1_i4) .and. (validcount > 0_i4))
    noValidPixels = count(validmaske)
    if (noValidPixels .GT. 0_i4) then
       PDMatrix = merge(PDMatrix/real(validcount, sp),PDMatrix, validmaske)
       ! average over all pixels
       PD_sp = 1.0_sp - sum(PDMatrix, mask=validmaske) / noValidPixels
       if (present(valid)) valid = .TRUE.
    ! case if maske is full of .false. or no valid neighbouring cells available for every pixel (validcount(:,:)=0)
    else
       PD_sp = 0.0_sp
       if (present(valid)) valid = .FALSE.
    end if
    !
  END FUNCTION PD_sp

  FUNCTION PD_dp(mat1, mat2, mask, valid)

    IMPLICIT NONE

    REAL(dp),    DIMENSION(:,:),                                            INTENT(IN)  :: mat1, mat2
    LOGICAL,     DIMENSION(:,:)                                 , OPTIONAL, INTENT(IN)  :: mask
    LOGICAL,                                                      OPTIONAL, INTENT(OUT) :: valid
    REAL(dp)                                                                            :: PD_dp

    INTEGER(i4)                                                                         :: iCo, iRo
    INTEGER(i4)                                                                         :: noValidPixels
    INTEGER(i4), DIMENSION(size(shape(mat1)) )                                          :: shapemask
    INTEGER(i4), DIMENSION(size(mat1, dim=1), size(mat1, dim=2))                        :: validcount
    REAL(dp),    DIMENSION(size(mat1, dim=1)+2_i4, size(mat1, dim=2)+2_i4)              :: bufferedMat1, bufferedMat2
    REAL(dp),    DIMENSION(size(mat1, dim=1), size(mat1, dim=2))                        :: PDMatrix
    LOGICAL,     DIMENSION(size(mat1, dim=1)+2_i4, size(mat1, dim=2)+2_i4)              :: maske
    LOGICAL,     DIMENSION(size(mat1, dim=1), size(mat1, dim=2))                        :: validmaske
    ! check if input has all the same dimensions
    if (present(mask)) then
       shapemask = shape(mask)
    else
       shapemask = shape(mat1)
    end if
    !
    if (any(shape(mat1) .NE. shape(mat2)))  &
         stop 'PD_dp: shapes of input matrix 1 and input matrix 2 are not matching'
    if (any(shape(mat1) .NE. shapemask))    &
         stop 'PD_dp: shapes of input matrices and mask are not matching'
    !
    ! additional 2 rows and 2 cols added for checking the border cells without crashing the search agorithm
    ! so the search windows can always cover 9 cells even if it is checking a border cell
    ! buffer rows are initialized as false values within mask, i.e. they are not considered for the calculation
    ! of the criteria
    !
    ! initialize mask with default=.false.
    maske = .false.
    if (present(mask)) then
       maske(2:(size(maske,dim=1)-1_i4), 2:(size(maske,dim=2)-1_i4)) = mask
    else
       maske(2:(size(maske,dim=1)-1_i4), 2:(size(maske,dim=2)-1_i4)) = .true.
    end if
    !
    ! initialize bufferedMat1 & bufferedMat2 
    bufferedMat1 = 0.0_dp
    bufferedMat2 = 0.0_dp
    bufferedMat1(2:(size(maske,dim=1)-1), 2:(size(maske,dim=2)-1)) = mat1
    bufferedMat2(2:(size(maske,dim=1)-1), 2:(size(maske,dim=2)-1)) = mat2
    !
    ! initialize PD
    PDMatrix = 0.0_dp
    do iCo = 2_i4, size(bufferedMat1, dim=2)-1_i4
       do iRo = 2_i4, size(bufferedMat1, dim=1)-1_i4
          ! no calculation at unmasked values
          if (.NOT. maske(iRo,iCo)) cycle
          ! .NEQV. is the fortran standard for .XOR. 
          ! result is written to -1 column and row because of buffer
          PDMatrix(iRo-1_i4,iCo-1_i4) = &
            real( &   
              count( & 
                     ( &
                        ! determine higher neighbouring values in mat1
                        (bufferedMat1(iRo-1_i4:iRo+1_i4, iCo-1_i4:iCo+1_i4) -                         &
                         bufferedMat1(iRo,iCo)                              > epsilon(0.0_dp)) .NEQV. &
                        ! determine higher neighbouring values in mat2                                &
                        (bufferedMat2(iRo-1_i4:iRo+1_i4, iCo-1_i4:iCo+1_i4) -                         &
                         bufferedMat2(iRo,iCo)                              > epsilon(0.0_dp))        &
                     ) &
                  ! exclude unmasked values
                  .and. (maske(iRo-1_i4:iRo+1_i4 , iCo-1_i4:iCo+1_i4)) &
                  ), &
            dp )
          ! count - 1 to exclude referendce cell (iRo, iCo)
          validcount(iRo-1_i4,iCo-1_i4) = count(maske(iRo-1_i4:iRo+1_i4 , iCo-1_i4:iCo+1_i4)) - 1_i4
          !
       end do
    end do
    !
    ! normalize every pixel to number of valid neighbouring cells (defined by maske) [0..8] --> [0..1]
    validmaske = (maske(2:size(maske, dim=1) - 1_i4, 2:size(maske, dim=2) - 1_i4) .and. (validcount > 0_i4))
    noValidPixels = count(validmaske)
    if (noValidPixels .GT. 0_i4) then
       PDMatrix = merge(PDMatrix/real(validcount, dp),PDMatrix, validmaske)
       ! average over all pixels
       PD_dp = 1.0_dp - sum(PDMatrix, mask=validmaske) / noValidPixels
       if (present(valid)) valid = .TRUE.
    ! case if maske is full of .false. or no valid neighbouring cells available for every pixel (validcount(:,:)=0)
    else
       PD_dp = 0.0_dp
       if (present(valid)) valid = .FALSE.
    end if
    !
  END FUNCTION PD_dp
  !
END MODULE mo_spatialsimilarity
