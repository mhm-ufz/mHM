!> \file mo_signatures.f90

!> \brief Module with calculations for several hydrological signatures.

!> \details This module contains calculations for hydrological signatures.
!> It contains:\n
!> * Autocorrelation
!> * Rising and declining limb densities
!> * Flow duration curves
!> * Peak distribution

!> \authors Remko Nijzink,
!> \date March 2014
!  modified, Stephan Thober, Jan 2015 - uncommented unused subroutines that are incompatible
!                                       with multiple time periods
!            Juliane Mai,    Jun 2015 - revising all routines
!                                     - removing routines which are specially low flow or high flow
!                                     - adding optional mask argument to all routines
!                                     - using CHS routines like percentile and autocorr

MODULE mo_signatures

  USE mo_kind, ONLY: i4, sp, dp

  IMPLICIT NONE

  PUBLIC :: Autocorrelation         ! Autocorrelation function
  PUBLIC :: FlowDurationCurve       ! Flow duration curve (i.e. CDF of runoff)
  PUBLIC :: Limb_densities          ! Rising and declining limb densities
  PUBLIC :: Moments                 ! Moments of data and log-transformed data, e.g. mean and standard deviation.
  PUBLIC :: PeakDistribution        ! Peak distribution parameter
  PUBLIC :: RunoffRatio             ! Runoff ratio (accumulated daily discharge [mm/d] / accumulated daily precipitation [mm/d])
  PUBLIC :: ZeroFlowRatio           ! Ratio of zero flow days to total observation days                             

  ! ------------------------------------------------------------------

CONTAINS

  !-------------------------------------------------------------------------------
  !     NAME
  !         Autocorrelation

  !     PURPOSE
  !>        \brief   Autocorrelation of a given data series.

  !>        \details Calculates the autocorrelation of a data series at given lags.
  !>                 An optional  argument for masking data points can be given.
  !>                 The function is basically a wrapper of the function autocorr
  !>                 from the module mo_corr.\n
  !>                 An optional mask of data points can be specified.

  !     CALLING SEQUENCE
  !         auto_correlation = Autocorrelation(data, lags, mask=mask)

  !     INTENT(IN)
  !>        \param[in] "real(dp), dimension(:)     :: data"   Array of data
  !>        \param[in] "integer(i4), dimension(:)  :: lags"   Array of lags where autocorrelation is requested 

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "logical, dimension(size(data,1))  :: mask"   Mask for data points given

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         \return     real(dp), dimension(size(lags,1)) :: Autocorrelation &mdash; autocorrelation of data at given lags

  !     RESTRICTIONS
  !>        Works only with 1d double precision input data.

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         Used as hydrologic signature with lag 1 in
  !            Euser, T., Winsemius, H. C., Hrachowitz, M., Fenicia, F., Uhlenbrook, S., & Savenije, H. H. G. (2013).
  !            A framework to assess the realism of model structures using hydrological signatures.
  !            Hydrology and Earth System Sciences, 17(5), 1893–1912. doi:10.5194/hess-17-1893-2013

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Jun 2015

  FUNCTION Autocorrelation(data, lags, mask)

    use mo_corr, only: autocorr

    IMPLICIT NONE

    real(dp),    dimension(:),                      intent(in)  :: data            ! Data series
    integer(i4), dimension(:),                      intent(in)  :: lags            ! Lags for autocorrelation
    logical,     dimension(size(data,1)), optional, intent(in)  :: mask            ! mask for data points
    real(dp),    dimension(size(lags,1))                        :: autocorrelation ! Autocorrelation of data at given lags

    if (present(mask)) then
       autocorrelation = autocorr(data, lags, mask=mask)
    else
       autocorrelation = autocorr(data, lags)
    end if

  END FUNCTION Autocorrelation

  !-------------------------------------------------------------------------------
  !     NAME
  !         FlowDurationCurves

  !     PURPOSE
  !>        \brief   Flow duration curves.

  !>        \details Calculates the flow duration curves for a given data vector. The Flow duration curve at a
  !>                 certain quantile x is the data point p where x% of the data points are above the value p.\n
  !>                 Hence the function percentile of the module mo_percentile is used. But percentile is
  !>                 determining the point p where x% of the data points are below that value. Therfore, the
  !>                 given quantiles are transformed by (1.0-quantile) to get the percentiles of exceedance probabilities.\n
  !>
  !>                 Optionally, the concavity index CI can be calculated [Zhang2014]. CI is defined by
  !>                         \f[ CI = \frac{q_{10\%}-q_{99\%}}{q_{1\%}-q_{99\%}} \f]
  !>                 where \f$ q_{x} \f$ is the data point where x% of the data points are above that value.
  !>                 Hence, exceedance probabilities are used.\n
  !>
  !>                 Optionally, the FDC mid-segment slope \f$FDC_{MSS}\f$ as used by Shafii et. al (2014) can be returned. 
  !>                 The \f$FDC_{MSS}\f$ is defined as
  !>                         \f[ FDC_{MSS} = \log(q_{m_1})-\log(q_{m_2}) \f]
  !>                 where \f$ m_1 \f$ and \f$ m_2 \f$ are the lowest and highest flow exceedance probabilities within the 
  !>                 midsegment of FDC. The settings \f$m_1=0.2\f$ and \f$0.7\f$ are used by Shafii et. al (2014) and are 
  !>                 implemented like that.\n
    !>
  !>                 Optionally, the FDC medium high-segment volume \f$FDC_{MHSV}\f$ as used by Shafii et. al (2014) can be 
  !>                 returned. The \f$FDC_{MHSV}\f$ is defined as
  !>                         \f[ FDC_{MHSV} = \sum_{h=1}^{H} q_h \f]
  !>                 where \f$ h=1,2,...,H \f$ are flow indeces located within the high-flow segment (exceedance probabilities  
  !>                 lower than \f$m_1\f$). \f$H\f$ is the index of the maximum flow. The settings \f$m_1=0.2\f$ is used here  
  !>                 to be consistent with the definitions of the low-segment (0.7-1.0) and the mid-segment (0.2-0.7).\n
  !>
  !>                 Optionally, the FDC high-segment volume \f$FDC_{HSV}\f$ as used by Shafii et. al (2014) can be returned. 
  !>                 The \f$FDC_{HSV}\f$ is defined as
  !>                         \f[ FDC_{HSV} = \sum_{h=1}^{H} q_h \f]
  !>                 where \f$ h=1,2,...,H \f$ are flow indeces located within the high-flow segment (exceedance probabilities  
  !>                 lower than \f$m_1\f$). \f$H\f$ is the index of the maximum flow. The settings \f$m_1=0.02\f$ is used by  
  !>                 Shafii et. al (2014) and is implemented like that.\n
  !>
  !>                 Optionally, the FDC low-segment volume \f$FDC_{LSV}\f$ as used by Shafii et. al (2014) can be returned. 
  !>                 The \f$FDC_{LSV}\f$ is defined as
  !>                         \f[ FDC_{LSV} = -\sum_{l=1}^{L} (\log(q_l) - \log(q_L)) \f]
  !>                 where \f$ l=1,2,...,L \f$ are flow indeces located within the low-flow segment (exceedance probabilities  
  !>                 larger than \f$m_1\f$). \f$L\f$ is the index of the minimum flow. The settings \f$m_1=0.7\f$ is used by  
  !>                 Shafii et. al (2014) and is implemented like that.\n
  !>
  !>                 An optional mask of data points can be specified.

  !     CALLING SEQUENCE
  !         flow_duration_curve = FlowDurationCurves(data, quantiles, mask=mask)

  !     INTENT(IN)
  !>        \param[in] "real(dp), dimension(:) :: data"              data series
  !>        \param[in] "real(dp), dimension(:) :: Quantiles"         Percentages of exceedance

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)  
  !         None
  
  !     INTENT(IN), OPTIONAL
  !>        \param[in] "logical, dimension(size(data,1)) :: mask"    mask of data array

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !>        \param[out] "real(dp), optional :: concavity_index"      concavity index as defined by Sauquet et al. (2011)
  !>        \param[out] "real(dp), optional :: mid_segment_slope"    mid-segment slope as defined by Shafii et al. (2014)
  !>        \param[out] "real(dp), optional :: mhigh_segment_volume" medium high-segment volume 
  !>        \param[out] "real(dp), optional :: high_segment_volume"  high-segment volume as defined by Shafii et al. (2014)
  !>        \param[out] "real(dp), optional :: low_segment_volume"   low-segment volume as defined by Shafii et al. (2014)

  !     RETURN
  !>        \return real(dp), dimension(size(quantiles,1)) :: FlowDurationCurve &mdash; Flow Duration Curve value at resp. quantile

  !     RESTRICTIONS
  !         Thresholds in mid_segment_slope, mhigh_segment_volume, high_segment_volume, low_segment_volume are hard coded.

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         FDC is used as hydrologic signature (quantiles not specified) in
  !            Euser, T., Winsemius, H. C., Hrachowitz, M., Fenicia, F., Uhlenbrook, S., & Savenije, H. H. G. (2013).
  !            A framework to assess the realism of model structures using hydrological signatures.
  !            Hydrology and Earth System Sciences, 17(5), 1893–1912. doi:10.5194/hess-17-1893-2013
  !         Concavity Index used as hydrologic signature in
  !            Zhang, Y., Vaze, J., Chiew, F. H. S., Teng, J., & Li, M. (2014).
  !            Predicting hydrological signatures in ungauged catchments using spatial interpolation, index model, and
  !            rainfall-runoff modelling.
  !            Journal of Hydrology, 517(C), 936–948. doi:10.1016/j.jhydrol.2014.06.032
  !         Concavity index is defined using exceedance probabilities by
  !            Sauquet, E., & Catalogne, C. (2011).
  !            Comparison of catchment grouping methods for flow duration curve estimation at ungauged sites in France.
  !            Hydrology and Earth System Sciences, 15(8), 2421–2435. doi:10.5194/hess-15-2421-2011
  !         mid_segment_slope, high_segment_volume, low_segment_volume used as hydrologic signature in
  !            Shafii, M., & Tolson, B. A. (2015).
  !            Optimizing hydrological consistency by incorporating hydrological signatures into model calibration objectives.
  !            Water Resources Research, 51(5), 3796–3814. doi:10.1002/2014WR016520
  
  !     HISTORY
  !>        \author Remko Nijzink, Juliane Mai
  !>        \date March 2014
  !         Modified   Juliane Mai, Jun 2015  - mask added
  !                                           - function instead of subroutine
  !                                           - use of percentile
  !                                           - add concavity_index
  !                    Juliane Mai, Jun 2015  - add mid_segment_slope, mhigh_segment_volume, high_segment_volume,
  !                                             low_segment_volume

  FUNCTION FlowDurationCurve(data, quantiles, mask, concavity_index, &
       mid_segment_slope, mhigh_segment_volume, high_segment_volume, low_segment_volume)

    use mo_percentile, only: percentile
    use mo_utils,      only: ge, le

    IMPLICIT NONE

    real(dp), dimension(:),                           intent(in) :: data                 ! data series
    real(dp), dimension(:),                           intent(in) :: quantiles            ! Percentages of exceedance (x-axis of FDC)
    logical,  dimension(:),                optional, intent(in)  :: mask                 ! mask for data
    real(dp),                              optional, intent(out) :: concavity_index      ! Concavity index as defined by
    !                                                                                    ! Sauquet et al. (2011)
    real(dp),                              optional, intent(out) :: mid_segment_slope    ! mid-segment slope as defined by
    !                                                                                    ! Shafii et al. (2014)
    real(dp),                              optional, intent(out) :: mhigh_segment_volume ! medium-high-segment volume
    !                                                                                    ! 0.0 <= exceed. prob. <= 0.2
    real(dp),                              optional, intent(out) :: high_segment_volume  ! high-segment volume as defined by
    !                                                                                    ! Shafii et al. (2014)
    !                                                                                    ! 0.0 <= exceed. prob. <= 0.02
    real(dp),                              optional, intent(out) :: low_segment_volume   ! low-segment volume as defined by
    !                                                                                    ! Shafii et al. (2014)
    !                                                                                    ! 0.7 <= exceed. prob. <= 1.0
    real(dp), dimension(size(quantiles,1))                       :: FlowDurationCurve    ! data point where x% of the data points
    !                                                                                    ! are above that value

    ! local variables
    logical, dimension(size(data,1)) :: maske            ! mask for data points
    real(dp)                         :: min_flow_value   ! minimal flow value
    real(dp)                         :: flow_value_thres ! flow value at a threshold quantile
    
    ! checking optionals
    if (present(mask)) then
       maske = mask
    else
       maske = .true.
    end if
    
    FlowDurationCurve = percentile(data, (1._dp-quantiles)*100._dp, mask=maske, mode_in=5)

    if (present(concavity_index)) then
       concavity_index = &
            ( percentile(data, (1._dp-0.10_dp)*100._dp, mask=maske, mode_in=5) - &
            percentile(data, (1._dp-0.99_dp)*100._dp, mask=maske, mode_in=5) ) / &
            ( percentile(data, (1._dp-0.01_dp)*100._dp, mask=maske, mode_in=5) - &
            percentile(data, (1._dp-0.99_dp)*100._dp, mask=maske, mode_in=5) )
    end if

    if (present(mid_segment_slope)) then
       ! mid-flows are defined to be between 0.2 and 0.7 by Shafii et. al (2014)
       mid_segment_slope = &
            log( percentile(data, (1._dp-0.2_dp)*100._dp, mask=maske, mode_in=5) ) - &
            log( percentile(data, (1._dp-0.7_dp)*100._dp, mask=maske, mode_in=5) ) 
    end if

    if (present(mhigh_segment_volume)) then
       ! medium high-flows are defined to be between 0.0 and 0.2 as to be constistent
       ! with the mid-segment (0.2-0.7) and low-segment (0.7-1.0) definitions
       flow_value_thres = percentile(data, (1._dp-0.2_dp)*100._dp, mask=maske, mode_in=5)
       mhigh_segment_volume = sum(data,mask=(maske .and. ge(data,flow_value_thres)))
       ! print*, 'flow_value_thres     = ',flow_value_thres
       ! print*, 'mhigh_segment_volume = ',mhigh_segment_volume
    end if

    if (present(high_segment_volume)) then
       ! high-flows are defined to be between 0.0 and 0.02 by Shafii et. al (2014)
       flow_value_thres = percentile(data, (1._dp-0.02_dp)*100._dp, mask=maske, mode_in=5)
       high_segment_volume = sum(data,mask=(maske .and. ge(data,flow_value_thres)))
    end if

    if (present(low_segment_volume)) then
       ! low-flows are defined to be between 0.7 and 1.0 by Shafii et. al (2014)
       min_flow_value   = minval(data, mask=maske)
       flow_value_thres = percentile(data, (1._dp-0.7)*100._dp, mask=maske, mode_in=5)
       low_segment_volume = -1.0_dp * &
            sum( log(data) - log(min_flow_value), mask=(maske .and. le(data,flow_value_thres)))
    end if

  END FUNCTION FlowDurationCurve

  !-------------------------------------------------------------------------------
  !     NAME
  !         Limb_densities

  !     PURPOSE
  !>        \brief   Calculates limb densities

  !>        \details Calculates rising and declinging limb densities. The peaks of the given series are
  !>                 first determined by looking for points where preceding and subsequent datapoint are lower.
  !>                 Second, the number of datapoints with rising values (nrise) and declining values (ndecline)
  !>                 are counted basically by comparing neighbors.\n

  !>                 The duration the data increase (nrise) divided by the number of peaks (npeaks)
  !>                 gives the rising limb density RLD
  !>                     \f[ RLD=t_{rise}/n_{peak} \f]
  !>                 whereas the duration the data decrease (ndecline) divided by the number of peaks (npeaks)
  !>                 gives the declining limb density DLD
  !>                     \f[ DLD=t_{fall}/n_{peak}. \f]
  
  !>                 An optional mask of data points can be specified.

  !     CALLING SEQUENCE
  !         call Limb_densities(data, mask=mask, RLD=rising_limb_density, DLD=declining_limb_density)

  !     INTENT(IN)
  !>        \param[in] "real(dp), dimension(:)           :: data"    data series

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "logical, dimension(size(data,1)) :: mask"    mask for data series

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !>        \param[out] "real(dp), optional              :: RLD"     rising    limb density
  !>        \param[out] "real(dp), optional              :: DLD"     declining limb density

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         Rising limb density used as hydrologic signature in
  !            Euser, T., Winsemius, H. C., Hrachowitz, M., Fenicia, F., Uhlenbrook, S., & Savenije, H. H. G. (2013).
  !            A framework to assess the realism of model structures using hydrological signatures.
  !            Hydrology and Earth System Sciences, 17(5), 1893–1912. doi:10.5194/hess-17-1893-2013

  !     HISTORY
  !>        \author Remko Nijzink
  !>        \date March 2014
  !         Modified   Juliane Mai, Jun 2015  - RLD and DLD as optional
  !                                           - optional mask for data can be given

  SUBROUTINE Limb_densities(data, mask, RLD, DLD)

    use mo_message,          only : message

    IMPLICIT NONE

    real(dp), dimension(:),                      intent(in)  :: data    ! data series (with rising and declining limbs)
    logical,  dimension(size(data,1)), optional, intent(in)  :: mask    ! mask for data series    
    real(dp),                          optional, intent(out) :: RLD     ! rising limb density
    real(dp),                          optional, intent(out) :: DLD     ! declining limb density

    ! local variables
    logical, dimension(size(data,1))    :: maske      ! mask for data points
    integer(i4)                         :: jj         ! Counter
    integer(i4)                         :: n_peak     ! Number of peaks
    integer(i4)                         :: n_decline  ! Number of declining data points
    integer(i4)                         :: n_rise     ! Number of rising data points
    logical, dimension(size(data,1))    :: goes_up    ! True if rising, False if declining or contains masked values
    real(dp)                            :: thres_rise ! Threshold that is has to rise at least to be detected as rising value

    ! checking optionals
    if (present(mask)) then
       maske = mask
    else
       maske = .true.
    end if

    if ( (.not. present(RLD)) .and. (.not. present(DLD)) ) then
       call message('mo_signatures: limb_densities: Neither RLD or DLD is specified in calling sequence.')
       stop
    end if
    
    ! initialize
    n_rise    = 0_i4
    n_decline = 0_i4
    n_peak    = 0_i4

    goes_up    = .False.
    thres_rise = 1.0_dp
    do jj=1, size(data,1)-1
       if (maske(jj) .and. maske(jj+1)) then
          if (data(jj) < data(jj+1) - thres_rise) then
             goes_up(jj) = .true.
             ! print*, jj, '  ', data(jj), '  ', data(jj+1)
          end if
       end if
    end do
    n_rise    = count(goes_up)
    n_decline = count(maske) - count(goes_up)

    ! write(*,*) 'goes_up = ', goes_up(1:178)

    ! peak is where goes_up switches from true to false
    n_peak = 0_i4
    do jj=1, size(data,1)-1
       if (maske(jj) .and. maske(jj+1)) then
          if (goes_up(jj) .and. .not.(goes_up(jj+1)) ) then
             n_peak = n_peak+1_i4
             ! print*, jj
          end if
       end if
    end do

    ! do jj=2, size(data,1)-1

    !    ! check for peak
    !    if (maske(jj-1) .and. maske(jj) .and. maske(jj+1)) then
    !       if ( (data(jj-1) .lt. data(jj)-1.0_dp) .and. (data(jj+1) .lt. data(jj)-1.0_dp) ) then
    !          n_peak = n_peak+1_i4
    !          write(*,*) jj-1
    !       end if
    !    end if

    !    ! check if data has rised
    !    if (maske(jj-1) .and. maske(jj)) then
    !       if (data(jj-1) .lt. data(jj)-1.0_dp) then
    !          n_rise = n_rise+1_i4
    !       else
    !          n_decline = n_decline+1_i4
    !       end if
    !    end if

    !    ! ! check if data has declined
    !    ! if (maske(jj-1) .and. maske(jj)) then
    !    !    if (data(jj-1) .gt. data(jj)) then
    !    !       n_decline = n_decline+1_i4
    !    !    end if
    !    ! end if
       
    ! end do

    ! write(*,*) 'n_peak = ', n_peak

    if (present(RLD)) then
       if (n_peak>0_i4) then
          RLD = real(n_rise,dp)/   real(n_peak,dp)
       else
          RLD = 0.0_dp
       end if
    end if
    
    if (present(DLD)) then
       if (n_peak>0_i4) then
          DLD = real(n_decline,dp)/real(n_peak,dp)
       else
          DLD = 0.0_dp
       end if
    end if
    
  END SUBROUTINE Limb_densities

    !-------------------------------------------------------------------------------
  !     NAME
  !         MaximumMonthlyFlow

  !     PURPOSE
  !>        \brief   Maximum of average flows per months.

  !>        \details Maximum of average flow per month is defined as
  !>                    \f[ max_{monthly flow} = Max( F(i), i=1,..12 ) \f]
  !>                 where \$f F(i) $\f is the average flow of month i.

  !     CALLING SEQUENCE
  !         max_monthly_flow = MaximumMonthlyFlow(data, yr_start=yr_start, mo_start=1, dy_start=1, mask=mask)

  !     INTENT(IN)
  !>        \param[in] "real(dp), dimension(:)     :: data"          array of data

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "logical, dimension(size(data,1))  :: mask"   mask for data points given
  !>        \param[in] "integer(i4)  :: yr_start"   year  of date of first data point given 
  !>        \param[in] "integer(i4)  :: mo_start"   month of date of first data point given (default: 1)
  !>        \param[in] "integer(i4)  :: dy_start"   month of date of first data point given (default: 1)

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        \return real(dp) :: MaximumMonthlyFlow &mdash; Maximum of average flow per month

  !     RESTRICTIONS
  !>        Works only with 1d double precision input data.\n
  !>        Assumes data are daily values.

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         used as hydrologic signature in
  !            Shafii, M., & Tolson, B. A. (2015).
  !            Optimizing hydrological consistency by incorporating hydrological signatures into model calibration objectives.
  !            Water Resources Research, 51(5), 3796–3814. doi:10.1002/2014WR016520

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Jun 2015

  FUNCTION MaximumMonthlyFlow(data, mask, yr_start, mo_start, dy_start)

    use mo_julian,  only: date2dec, dec2date
    use mo_message, only: message

    IMPLICIT NONE

    real(dp),    dimension(:),                      intent(in)  :: data               ! Data series
    logical,     dimension(size(data,1)), optional, intent(in)  :: mask               ! mask for data points
    integer(i4),                          optional, intent(in)  :: yr_start           ! year  of date of first data point given
    integer(i4),                          optional, intent(in)  :: mo_start           ! month of date of first data point given
    !                                                                                 ! DEFAULT: 1
    integer(i4),                          optional, intent(in)  :: dy_start           ! day   of date of first data point given
    !                                                                                 ! DEFAULT: 1
    real(dp)                                                    :: MaximumMonthlyFlow ! return: maximum of average monthly flow
    
    ! local variables
    logical,     dimension(size(data,1)) :: maske
    integer(i4)                          :: ii               ! counter
    integer(i4)                          :: yr, mo, dy, imo  ! date variables
    integer(i4), dimension(12)           :: counter          ! number of data points per month
    real(dp),    dimension(12)           :: flow_month       ! summed data points per months
    real(dp)                             :: ref_jul_day      ! julian day of one day before start day

    if (present(mask)) then
       maske = mask
    else
       maske = .true.
    end if

    if (.not. present(yr_start)) then
       call message('mo_signatures: MaximumMonthlyFlow: Year of of data point has to be given!')
       stop
    else
       yr = yr_start
    end if

    if (present(mo_start)) then
       mo = mo_start
    else
       mo = 1
    end if

    if (present(dy_start)) then
       dy = dy_start
    else
       dy = 1
    end if

    flow_month = 0.0_dp
    counter    = 0_i4
    ref_jul_day = date2dec(yy=yr, mm=mo, dd=dy) - 1.0_dp 
    
    do ii=1, size(data,1)
       if (maske(ii)) then
          ! determine current month
          call dec2date(ref_jul_day+real(ii,dp), mm=imo)
          ! add value
          counter(imo)    = counter(imo) + 1
          flow_month(imo) = flow_month(imo) + data(ii)
       end if
    end do

    if (any(counter == 0_i4)) then
       call message('mo_signatures: MaximumMonthlyFlow: There are months with no data points!')
       call message('                                   Aborted!')
       stop
    end if

    ! average
    MaximumMonthlyFlow = maxval(flow_month / real(counter,dp))

  END FUNCTION MaximumMonthlyFlow
  
  !-------------------------------------------------------------------------------
  !     NAME
  !         Moments

  !     PURPOSE
  !>        \brief   Moments of data and log-transformed data, e.g. mean and standard deviation.

  !>        \details Returns several moments of data series given, i.e.
  !>                     * mean               of data 
  !>                     * standard deviation of data 
  !>                     * median             of data 
  !>                     * maximum/ peak      of data 
  !>                     * mean               of log-transformed data 
  !>                     * standard deviation of log-transformed data 
  !>                     * median             of log-transformed data 
  !>                     * maximum/ peak      of log-transformed data 
  !>                 An optional mask of data points can be specified.

  !     CALLING SEQUENCE
  !         call moments(data, mask=mask, mean_data=mean_data, stddev_data=stddev_data, median_data=median_data, max_data=max_data,&
  !                                       mean_log=mean_log,   stddev_log=stddev_log,   median_log=median_log,   max_log=max_log)

  !     INTENT(IN)
  !>        \param[in] "real(dp), dimension(:)     :: data"          array of data

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "logical, dimension(size(data,1))  :: mask"   mask for data points given

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !>        \param[in]  "real(dp)  :: mean_data"    mean               of data
  !>        \param[in]  "real(dp)  :: stddev_data"  standard deviation of data
  !>        \param[out] "real(dp)  :: median_data"  median             of data
  !>        \param[out] "real(dp)  :: max_data"     maximum/ peak      of data
  !>        \param[out] "real(dp)  :: mean_log"     mean               of log-transformed data
  !>        \param[out] "real(dp)  :: stddev_log"   standard deviation of log-transformed data
  !>        \param[out] "real(dp)  :: median_log"   median             of log-transformed data
  !>        \param[out] "real(dp)  :: max_log"      maximum/ peak      of log-transformed data

  !     RETURN
  !         None

  !     RESTRICTIONS
  !>        Works only with 1d double precision input data.

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         mean_log and stddev_log used as hydrologic signature in
  !            Zhang, Y., Vaze, J., Chiew, F. H. S., Teng, J., & Li, M. (2014).
  !            Predicting hydrological signatures in ungauged catchments using spatial interpolation, index model, and
  !            rainfall-runoff modelling.
  !            Journal of Hydrology, 517(C), 936–948. doi:10.1016/j.jhydrol.2014.06.032
  !         mean_data, stddev_data, median_data, max_data, mean_log, and stddev_log used as hydrologic signature in
  !            Shafii, M., & Tolson, B. A. (2015).
  !            Optimizing hydrological consistency by incorporating hydrological signatures into model calibration objectives.
  !            Water Resources Research, 51(5), 3796–3814. doi:10.1002/2014WR016520

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Jun 2015

  SUBROUTINE Moments(data, mask, mean_data, stddev_data, median_data, max_data, mean_log, stddev_log, median_log, max_log)

    use mo_message,    only: message
    use mo_moment,     only: mean, stddev
    use mo_percentile, only: median

    IMPLICIT NONE

    real(dp),    dimension(:),                      intent(in)  :: data            ! Data series
    logical,     dimension(size(data,1)), optional, intent(in)  :: mask            ! mask for data points
    real(dp),                             optional, intent(out) :: mean_data       ! mean               of data                
    real(dp),                             optional, intent(out) :: stddev_data     ! standard deviation of data
    real(dp),                             optional, intent(out) :: median_data     ! median             of data                
    real(dp),                             optional, intent(out) :: max_data        ! max/ peak          of data                
    real(dp),                             optional, intent(out) :: mean_log        ! mean               of log-transformed data
    real(dp),                             optional, intent(out) :: stddev_log      ! standard deviation of log-transformed data
    real(dp),                             optional, intent(out) :: median_log      ! median             of log-transformed data
    real(dp),                             optional, intent(out) :: max_log         ! max/ peak          of log-transformed data

    ! local variables
    logical,     dimension(size(data,1)) :: maske
    real(dp),    dimension(size(data,1)) :: logdata

    if (present(mask)) then
       maske = mask
    else
       maske = .true.
    end if

    if ( .not.(present(mean_data)) .and. .not.(present(stddev_data)) .and. &
         .not.(present(median_data)) .and. .not.(present(max_data)) .and. &
         .not.(present(mean_log))  .and. .not.(present(stddev_log)) .and. &
         .not.(present(median_log))  .and. .not.(present(max_log)) ) then
       call message('mo_signatures: Moments: None of the optional output arguments is specified')
       stop
    end if

    if (present(mean_data))   mean_data   = mean(data, mask=maske)
    if (present(stddev_data)) stddev_data = stddev(data, mask=maske)
    if (present(median_data)) median_data = median(data, mask=maske)
    if (present(max_data))    max_data    = maxval(data, mask=maske)

    if (present(mean_log) .or. present(stddev_log)) then 
       where (data > 0.0_dp)
          logdata = log(data)
       elsewhere
          logdata = -9999.0_dp  ! will not be used, since mask is set to .false.
          maske    = .false.
       end where
    
       if (present(mean_log))   mean_log   = mean(logdata, mask=maske)
       if (present(stddev_log)) stddev_log = stddev(logdata, mask=maske)
       if (present(median_log)) median_log = median(logdata, mask=maske)
       if (present(max_log))    max_log    = maxval(logdata, mask=maske)
    end if

  END SUBROUTINE Moments

  !-------------------------------------------------------------------------------
  !     NAME
  !         PeakDistribution

  !     PURPOSE
  !>        \brief   Calculates the peak distribution.

  !>        \details First, the peaks of the time series given are identified. For the peak distribution
  !>                 only this subset of data points are considered. Second, the peak distribution at the
  !>                 quantiles given is calculated. Calculates the peak distribution at the quantiles given
  !>                 using mo_percentile. Since the exceedance probabilities are usually used in 
  !>                 hydrology the function percentile is used with (1.0-quantiles). \n
  !>                 
  !>                 Optionally, the slope of the peak distribution between 10th and 50th percentile, i.e.
  !>                    \f[ slope = \frac{\mathrm{peak\_{data}}_{0.1}-\mathrm{peak\_{data}}_{0.5}}{0.9-0.5} \f]
  !>                 can be returned.\n
  !>                 An optional mask for the data points can be given.

  !     CALLING SEQUENCE
  !         peak_distr_at_quantiles = PeakDistribution(data, quantiles, peaks, Q_quantiles)

  !     INTENT(IN)
  !>        \param[in] "real(dp), dimension(:)  :: data"                       data array
  !>        \param[in] "real(dp), dimension(:)  :: quantiles"                  requested quantiles for distribution

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "logical, dimension(size(data,1)), optional  :: mask"   mask of data array

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !>        \param[in] "real(dp), optional :: slope_peak_distribution"         slope of the Peak distribution between \n
  !>                                                                           10th and 50th percentile

  !     RETURN
  !>        \return real(dp), dimension(size(quantiles,1)) :: PeakDistribution &mdash; Distribution of peak values \n
  !>                                                                                   at resp. quantiles

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         slope_peak_distribution used as hydrologic signature in
  !            Euser, T., Winsemius, H. C., Hrachowitz, M., Fenicia, F., Uhlenbrook, S., & Savenije, H. H. G. (2013).
  !            A framework to assess the realism of model structures using hydrological signatures.
  !            Hydrology and Earth System Sciences, 17(5), 1893–1912. doi:10.5194/hess-17-1893-2013

  !     HISTORY
  !>        \author Remko Nijzink
  !>        \date March 2014
  !         Modified   Juliane Mai, Jun 2015  - mask added
  !                                           - function instead of subroutine
  !                                           - use of percentile

  FUNCTION PeakDistribution(data, quantiles, mask, slope_peak_distribution)

    USE mo_percentile, only: percentile

    IMPLICIT NONE

    real(dp), dimension(:),                          intent(in)  :: data                      ! data points      
    real(dp), dimension(:),                          intent(in)  :: quantiles                 ! percentages of occurence
    logical,  dimension(size(data,1)),     optional, intent(in)  :: mask                      ! mask of data points
    real(dp),                              optional, intent(out) :: slope_peak_distribution   ! slope of peak flow duration curve
    real(dp), dimension(size(quantiles,1))                       :: PeakDistribution          ! distribution of peaks in data
    !                                                                                         ! returns values of distribution at
    !                                                                                         ! given quantiles

    ! local variables
    integer(i4)                              :: ii, jj      ! counters
    logical, dimension(size(data,1))         :: maske       ! mask of data 
    integer(i4)                              :: n_peak      ! Number of peaks
    real(dp),                 dimension(2)   :: pp          ! array containing some quantiles
    real(dp),                 dimension(2)   :: data_pp     ! data points of quantiles pp 
    real(dp) ,   allocatable, dimension(:)   :: data_peak   ! series containing only peak values of original series data

    ! checking optionals
    if (present(mask)) then
       maske = mask
    else
       maske = .true.
    end if
    
    ! count peaks
    n_peak=0_i4
    do jj=2, size(data,1)-1
       if (maske(jj-1) .and. maske(jj) .and. maske(jj+1)) then
          if ( (data(jj-1) .le. data(jj)) .and. (data(jj+1) .le. data(jj)) ) then
             n_peak = n_peak+1_i4
          end if
       end if
    end do

    allocate(data_peak(n_peak))

    ! find peaks
    jj=0
    do ii=2, size(data,1)-1
       if (maske(ii-1) .and. maske(ii) .and. maske(ii+1)) then
          if( (data(ii-1) .le. data(ii)) .and. (data(ii+1) .le. data(ii)) ) then
             jj=jj+1_i4       
             data_peak(jj)=data(ii)
          end if
       end if
    end do

    if (present(slope_peak_distribution)) then
       ! calculate slope between 10% and 50% quantiles, per definition
       pp      = (/ 0.1_dp, 0.5_dp /)
       data_pp = percentile(data_peak,(1.0_dp-pp)*100._dp, mode_in=5)   ! (1-p) because exceedence probability is required
       slope_peak_distribution = (data_pp(1)-data_pp(2))/(0.9_dp-0.5_dp)
    end if

    PeakDistribution = percentile(data_peak,(1.0_dp-Quantiles)*100._dp, mode_in=5) 
    deallocate( data_peak)

  END FUNCTION PeakDistribution

    !-------------------------------------------------------------------------------
  !     NAME
  !         RunoffRatio

  !     PURPOSE
  !>        \brief   Runoff ratio (accumulated daily discharge [mm/d] / accumulated daily precipitation [mm/d]).

  !>        \details The runoff ratio is defined as
  !>                       \f[ runoff_ratio = \frac{\sum_{t=1}^{N} q_t}{\sum_{t=1}^{N} p_t}\f]
  !>                 where \f$p_t\f$ and \f$q_t\f$ are precipitation and discharge, respectively. \n 
  !>                 Therefore, precipitation over the entire basin is required and both discharge and precipitation
  !>                 have to be converted to the same units [mm/d].\n
  !>  
  !>                 Input discharge is given in [m**3/s] as this is mHM default while precipitation has to be given
  !>                 in [mm/km**2 / day].\n
  !>                 
  !>                 Either "precip_sum" or "precip_series" has to be specified.
  !>                 If "precip_series" is used the optional mask is also applied to precipitation values.
  !>                 The "precip_sum" is the accumulated "precip_series".\n
  !>
  !>                 Optionally, a mask for the data (=discharge) can be given. If optional "log_data" is set to .true.
  !>                 the runoff ratio will be calculated as
  !>                       \f[ runoff\_ratio = \frac{\sum_{t=1}^{N} \log(q_t)}{\sum_{t=1}^{N} p_t}\f]
  !>                 where \f$p_t\f$ and \f$q_t\f$ are precipitation and discharge, respectively. \n 

  !     CALLING SEQUENCE
  !         Runoff_Ratio = RunoffRatio(data, basin_area, mask=mask, precip_sum=precip_sum,       log_data=.False.)
  !            OR
  !         Runoff_Ratio = RunoffRatio(data, basin_area, mask=mask, precip_series=precip_series, log_data=.False.)

  !     INTENT(IN)
  !>        \param[in] "real(dp), dimension(:)     :: data"          array of data   [m**3/s]
  !>        \param[in] "real(dp)                   :: basin_area"    area of basin   [km**2]

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "logical, dimension(size(data,1))  :: mask"            mask for data points given
  !>        \param[in] "real(dp)                          :: precip_sum"      sum of daily precip. values of whole period 
  !>                                                                          [mm/km**2 / day]
  !>        \param[in] "real(dp), dimension(size(data,1)) :: precip_series"   daily precipitation values [mm/km**2 / day]
  !>        \param[in] "logical                           :: log_data"        ratio using logarithmic data

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         \return     real(dp), dimension(size(lags,1)) :: RunoffRation &mdash; Ratio of discharge and precipitation

  !     RESTRICTIONS
  !>        Works only with 1d double precision input data.

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         Used as hydrologic signature in
  !            Shafii, M., & Tolson, B. A. (2015).
  !            Optimizing hydrological consistency by incorporating hydrological signatures into model calibration objectives.
  !            Water Resources Research, 51(5), 3796–3814. doi:10.1002/2014WR016520

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Jun 2015

  FUNCTION RunoffRatio(data, basin_area, mask, precip_series, precip_sum, log_data)

    use mo_message, only: message

    IMPLICIT NONE

    real(dp),    dimension(:),                      intent(in)  :: data            ! Discharge series [m**3/s]
    real(dp),                                       intent(in)  :: basin_area      ! Basin area       [km**2]
    logical,     dimension(size(data,1)), optional, intent(in)  :: mask            ! mask for data points
    real(dp),    dimension(size(data,1)), optional, intent(in)  :: precip_series   ! series of daily precip. [mm/km**2 / d]
    !                                                                              ! optional mask will be applied
    real(dp),                             optional, intent(in)  :: precip_sum      ! sum    of daily precip. [mm/km**2 / d]
    !                                                                              ! optional mask will not be applied
    logical,                              optional, intent(in)  :: log_data        ! if log(data) is used in ratio
    !                                                                              ! DEFAULT: .false.
    real(dp)                                                    :: RunoffRatio     ! sum(data) / sum(precip)

    ! local variables
    logical, dimension(size(data,1)) :: maske    ! mask of data 
    logical                          :: log_dat  ! if logarithmic data are used --> sum(log(data)) / sum(precip)
    real(dp)                         :: sum_discharge
    real(dp)                         :: sum_precip

    ! checking optionals
    if (present(mask)) then
       maske = mask
    else
       maske = .true.
    end if

    if (present(log_data)) then
       log_dat = log_data
    else
       log_dat = .false.
    end if

    if ( (present(precip_series) .and. present(precip_sum)) .or. &
         (.not. present(precip_series) .and. .not. present(precip_sum)) ) then
       call message('mo_signatures: RunoffRatio: Exactly one precipitation information')
       call message('                            (precipitation series or sum of precipitation) ')
       call message('                            has to be specified!')
       stop
    end if

    if (present(mask) .and. present(precip_sum)) then
       call message('mo_signatures: RunoffRatio: Already aggregated precipitation (precip_sum) and' )
       call message('                            mask can not be used together.')
       call message('                            Precip_series should be used instead!')
       stop
    end if

    ! mhm output [m**3/s]  --> required [mm/d]
    !    [m**3/s] / [km**2] = [m**3/(s km**2)]
    ! => [m**3/(s km**2) * 60*60*24/1000**2] = [m/d]
    ! => [m**3/(s km**2) * 60*60*24*1000/1000**2] = [mm/d]
    ! => [m**3/(s km**2) * 86.4 ] = [mm/d]
    ! => discharge value [m**3/s] / catchment area [km**2] * 86.4 [km**2 s/m**3 * mm/d]
    if (log_dat) then
       sum_discharge = sum(log(data)*86.4_dp/basin_area,mask=maske)
    else
       sum_discharge = sum(data*86.4_dp/basin_area,mask=maske)
    end if

    if (present(precip_sum)) then
       sum_precip = precip_sum
    else
       sum_precip = sum(precip_series,mask=maske)
    end if

    RunoffRatio = sum_discharge / sum_precip

  END FUNCTION RunoffRatio

  !-------------------------------------------------------------------------------
  !     NAME
  !         ZeroFlowRatio

  !     PURPOSE
  !>        \brief   Ratio of zero values to total number of data points.

  !>        \details An optional mask of data points can be specified.

  !     CALLING SEQUENCE
  !         Zero_Flow_Ratio = ZeroFlowRatio(data, mask=mask)

  !     INTENT(IN)
  !>        \param[in] "real(dp), dimension(:)     :: data"          array of data

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "logical, dimension(size(data,1))  :: mask"   mask for data points given

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         \return     real(dp), dimension(size(lags,1)) :: ZeroFlowRatio &mdash; Ratio of zero values to total number of data points

  !     RESTRICTIONS
  !>        Works only with 1d double precision input data.

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         Used as hydrologic signature in
  !            Zhang, Y., Vaze, J., Chiew, F. H. S., Teng, J., & Li, M. (2014).
  !            Predicting hydrological signatures in ungauged catchments using spatial interpolation, index model, and
  !            rainfall-runoff modelling.
  !            Journal of Hydrology, 517(C), 936–948. doi:10.1016/j.jhydrol.2014.06.032

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Jun 2015

  FUNCTION ZeroFlowRatio(data, mask)

    use mo_message, only: message
    use mo_utils,   only: eq

    IMPLICIT NONE

    real(dp),    dimension(:),                      intent(in)  :: data            ! Data series
    logical,     dimension(size(data,1)), optional, intent(in)  :: mask            ! mask for data points
    real(dp)                                                    :: ZeroFlowRatio   ! Autocorrelation of data at given lags

    ! local variables
    logical, dimension(size(data,1)) :: maske   ! mask of data 
    integer(i4)                      :: nall    ! total number of data points
    integer(i4)                      :: nzero   ! number of zero data points

    ! checking optionals
    if (present(mask)) then
       maske = mask
    else
       maske = .true.
    end if

    nall  = count(maske)
    nzero = count(maske .and. (eq(data, 0.0_dp)))

    if (nall > 0) then
       ZeroFlowRatio = real(nzero, dp) / real(nall, dp)
    else
       call message('mo_signatures: ZeroFlowRatio: all data points are masked')
       stop
    end if

  END FUNCTION ZeroFlowRatio

END MODULE mo_signatures
