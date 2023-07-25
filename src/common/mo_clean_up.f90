!> \file    mo_clean_up.f90
!> \brief   \copybrief mo_clean_up
!> \details \copydetails mo_clean_up

!> \brief   Module to clean up after a mHM run.
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    May 2022
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_common
module mo_clean_up

  implicit none
  private
  public :: deallocate_global_variables

  contains

  !> \brief Deallocate all global variables.
  subroutine deallocate_global_variables()
    use mo_mhm_cli, only : set_verbosity_level
    use mo_common_run_variables, only : run_cfg
    use mo_global_variables, only : &
      couple_cfg, &
      meteo_handler, &
      L1_sm, &
      L1_sm_mask, &
      L1_neutronsdata, &
      L1_neutronsdata_mask, &
      L1_smObs, &
      L1_neutronsObs, &
      L1_etObs, &
      L1_twsaObs, &
      BFI_obs, &
      BFI_qBF_sum, &
      BFI_qT_sum, &
      L1_inter, &
      L1_snowPack, &
      L1_sealSTW, &
      L1_soilMoist, &
      L1_unsatSTW, &
      L1_satSTW, &
      L1_neutrons, &
      L1_pet_calc, &
      L1_temp_calc, &
      L1_prec_calc, &
      L1_aETSoil, &
      L1_aETCanopy, &
      L1_aETSealed, &
      L1_baseflow, &
      L1_infilSoil, &
      L1_fastRunoff, &
      L1_melt, &
      L1_percol, &
      L1_preEffect, &
      L1_rain, &
      L1_runoffSeal, &
      L1_slowRunoff, &
      L1_snow, &
      L1_Throughfall, &
      L1_total_runoff, &
      neutron_integral_AFast

    use mo_common_variables, only : &
      resolutionHydrology, &
      L0_Domain, &
      mhmFileRestartOut, &
      mrmFileRestartOut, &
      dirMorpho, &
      dirLCover, &
      dirOut, &
      fileLatLon, &
      level0, &
      level1, &
      l0_l1_remap, &
      L0_elev, &
      L0_LCover, &
      LCfilename, &
      LC_year_start, &
      LC_year_end, &
      global_parameters, &
      global_parameters_name, &
      domainMeta

    use mo_mpr_global_variables, only : &
      HorizonDepth_mHM, &
      GeoUnitList, &
      GeoUnitKar, &
      LAIBoundaries, &
      LAIUnitList, &
      LAILUT, &
      LAIPer, &
      L0_slope_emp, &
      L0_gridded_LAI, &
      L0_slope, &
      L0_asp, &
      L0_soilId, &
      L0_geoUnit, &
      dirgridded_LAI, &
      L1_fSealed, &
      L1_alpha, &
      L1_degDayInc, &
      L1_degDayMax, &
      L1_degDayNoPre, &
      L1_degDay, &
      L1_karstLoss, &
      L1_fAsp, &
      L1_petLAIcorFactor, &
      L1_HarSamCoeff, &
      L1_PrieTayAlpha, &
      L1_aeroResist, &
      L1_surfResist, &
      L1_fRoots, &
      L1_maxInter, &
      L1_kfastFlow, &
      L1_kSlowFlow, &
      L1_kBaseFlow, &
      L1_kPerco, &
      L1_soilMoistFC, &
      L1_soilMoistSat, &
      L1_soilMoistExp, &
      L1_jarvis_thresh_c1, &
      L1_tempThresh, &
      L1_unsatThresh, &
      L1_sealedThresh, &
      L1_wiltingPoint, &
      soilDB, &
      L1_No_Count, &
      L1_bulkDens, &
      L1_latticeWater, &
      L1_COSMICL3
    use mo_mrm_global_variables, only : &
      dirGauges, &
      dirTotalRunoff, &
      dirBankfullRunoff, &
      level11, &
      l0_l11_remap, &
      mRM_runoff, &
      domain_mrm, &
      L0_gaugeLoc, &
      L0_InflowGaugeLoc, &
      L0_fAcc, &
      L0_fDir, &
      L0_draSC, &
      L0_draCell, &
      L0_streamNet, &
      L0_floodPlain, &
      L0_noutlet, &
      L0_celerity, &
      L11_L1_ID, &
      L1_total_runoff_in, &
      L11_cellCoor, &
      L1_L11_ID, &
      L11_areaCell, &
      L11_fAcc, &
      L11_fDir, &
      L11_nOutlets, &
      L11_celerity, &
      L11_meandering, &
      L11_LinkIn_fAcc, &
      L11_rowOut, &
      L11_colOut, &
      L11_Qmod, &
      L11_qOUT, &
      L11_qTIN, &
      L11_qTR, &
      L11_fromN, &
      L11_toN, &
      L11_netPerm, &
      L11_fRow, &
      L11_fCol, &
      L11_tRow, &
      L11_tCol, &
      L11_rOrder, &
      L11_label, &
      L11_sink, &
      L11_length, &
      L11_aFloodPlain, &
      L11_slope, &
      L11_nLinkFracFPimp, &
      L11_K, &
      L11_xi, &
      L11_tsRout, &
      L11_C1, &
      L11_C2, &
      L11_bankfull_runoff_in, &
      L0_channel_depth, &
      L0_channel_elevation, &
      L0_river_head_mon_sum, &
      gauge, &
      InflowGauge, &
      riv_temp_pcs, &
      mrm_L0_slope => L0_slope

    use mo_common_mHM_mRM_variables, only : &
      resolutionRouting, &
      warmPer, &
      evalPer, &
      simPer, &
      warmingDays, &
      LCyearId, &
      mhmFileRestartIn, &
      mrmFileRestartIn

    ! mo_global_variables
    if ( allocated(L1_sm) ) deallocate(L1_sm)
    if ( allocated(L1_sm_mask) ) deallocate(L1_sm_mask)
    if ( allocated(L1_neutronsdata) ) deallocate(L1_neutronsdata)
    if ( allocated(L1_neutronsdata_mask) ) deallocate(L1_neutronsdata_mask)
    if ( allocated(L1_smObs) ) deallocate(L1_smObs)
    if ( allocated(L1_neutronsObs) ) deallocate(L1_neutronsObs)
    if ( allocated(L1_etObs) ) deallocate(L1_etObs)
    if ( allocated(L1_twsaObs) ) deallocate(L1_twsaObs)
    if ( allocated(BFI_obs) ) deallocate(BFI_obs)
    if ( allocated(BFI_qBF_sum) ) deallocate(BFI_qBF_sum)
    if ( allocated(BFI_qT_sum) ) deallocate(BFI_qT_sum)
    if ( allocated(L1_inter) ) deallocate(L1_inter)
    if ( allocated(L1_snowPack) ) deallocate(L1_snowPack)
    if ( allocated(L1_sealSTW) ) deallocate(L1_sealSTW)
    if ( allocated(L1_soilMoist) ) deallocate(L1_soilMoist)
    if ( allocated(L1_unsatSTW) ) deallocate(L1_unsatSTW)
    if ( allocated(L1_satSTW) ) deallocate(L1_satSTW)
    if ( allocated(L1_neutrons) ) deallocate(L1_neutrons)
    if ( allocated(L1_pet_calc) ) deallocate(L1_pet_calc)
    if ( allocated(L1_temp_calc) ) deallocate(L1_temp_calc)
    if ( allocated(L1_prec_calc) ) deallocate(L1_prec_calc)
    if ( allocated(L1_aETSoil) ) deallocate(L1_aETSoil)
    if ( allocated(L1_aETCanopy) ) deallocate(L1_aETCanopy)
    if ( allocated(L1_aETSealed) ) deallocate(L1_aETSealed)
    if ( allocated(L1_baseflow) ) deallocate(L1_baseflow)
    if ( allocated(L1_infilSoil) ) deallocate(L1_infilSoil)
    if ( allocated(L1_fastRunoff) ) deallocate(L1_fastRunoff)
    if ( allocated(L1_melt) ) deallocate(L1_melt)
    if ( allocated(L1_percol) ) deallocate(L1_percol)
    if ( allocated(L1_preEffect) ) deallocate(L1_preEffect)
    if ( allocated(L1_rain) ) deallocate(L1_rain)
    if ( allocated(L1_runoffSeal) ) deallocate(L1_runoffSeal)
    if ( allocated(L1_slowRunoff) ) deallocate(L1_slowRunoff)
    if ( allocated(L1_snow) ) deallocate(L1_snow)
    if ( allocated(L1_Throughfall) ) deallocate(L1_Throughfall)
    if ( allocated(L1_total_runoff) ) deallocate(L1_total_runoff)
    if ( allocated(neutron_integral_AFast) ) deallocate(neutron_integral_AFast)

    ! mo_common_variables
    if ( allocated(resolutionHydrology) ) deallocate(resolutionHydrology)
    if ( allocated(L0_Domain) ) deallocate(L0_Domain)
    if ( allocated(mhmFileRestartOut) ) deallocate(mhmFileRestartOut)
    if ( allocated(mrmFileRestartOut) ) deallocate(mrmFileRestartOut)
    if ( allocated(dirMorpho) ) deallocate(dirMorpho)
    if ( allocated(dirLCover) ) deallocate(dirLCover)
    if ( allocated(dirOut) ) deallocate(dirOut)
    if ( allocated(fileLatLon) ) deallocate(fileLatLon)
    if ( allocated(level0) ) deallocate(level0)
    if ( allocated(level1) ) deallocate(level1)
    if ( allocated(l0_l1_remap) ) deallocate(l0_l1_remap)
    if ( allocated(L0_elev) ) deallocate(L0_elev)
    if ( allocated(L0_LCover) ) deallocate(L0_LCover)
    if ( allocated(LCfilename) ) deallocate(LCfilename)
    if ( allocated(LC_year_start) ) deallocate(LC_year_start)
    if ( allocated(LC_year_end) ) deallocate(LC_year_end)
    if ( allocated(global_parameters) ) deallocate(global_parameters)
    if ( allocated(global_parameters_name) ) deallocate(global_parameters_name)
    if ( allocated(domainMeta%indices) ) deallocate(domainMeta%indices)
    if ( allocated(domainMeta%L0DataFrom) ) deallocate(domainMeta%L0DataFrom)
    if ( allocated(domainMeta%optidata) ) deallocate(domainMeta%optidata)
    if ( allocated(domainMeta%doRouting) ) deallocate(domainMeta%doRouting)

    ! mo_mpr_global_variables
    if ( allocated(HorizonDepth_mHM) ) deallocate(HorizonDepth_mHM)
    if ( allocated(GeoUnitList) ) deallocate(GeoUnitList)
    if ( allocated(GeoUnitKar) ) deallocate(GeoUnitKar)
    if ( allocated(LAIBoundaries) ) deallocate(LAIBoundaries)
    if ( allocated(LAIUnitList) ) deallocate(LAIUnitList)
    if ( allocated(LAILUT) ) deallocate(LAILUT)
    if ( allocated(LAIPer) ) deallocate(LAIPer)
    if ( allocated(L0_slope_emp) ) deallocate(L0_slope_emp)
    if ( allocated(L0_gridded_LAI) ) deallocate(L0_gridded_LAI)
    if ( allocated(L0_slope) ) deallocate(L0_slope)
    if ( allocated(L0_asp) ) deallocate(L0_asp)
    if ( allocated(L0_soilId) ) deallocate(L0_soilId)
    if ( allocated(L0_geoUnit) ) deallocate(L0_geoUnit)
    if ( allocated(dirgridded_LAI) ) deallocate(dirgridded_LAI)
    if ( allocated(L1_fSealed) ) deallocate(L1_fSealed)
    if ( allocated(L1_alpha) ) deallocate(L1_alpha)
    if ( allocated(L1_degDayInc) ) deallocate(L1_degDayInc)
    if ( allocated(L1_degDayMax) ) deallocate(L1_degDayMax)
    if ( allocated(L1_degDayNoPre) ) deallocate(L1_degDayNoPre)
    if ( allocated(L1_degDay) ) deallocate(L1_degDay)
    if ( allocated(L1_karstLoss) ) deallocate(L1_karstLoss)
    if ( allocated(L1_fAsp) ) deallocate(L1_fAsp)
    if ( allocated(L1_petLAIcorFactor) ) deallocate(L1_petLAIcorFactor)
    if ( allocated(L1_HarSamCoeff) ) deallocate(L1_HarSamCoeff)
    if ( allocated(L1_PrieTayAlpha) ) deallocate(L1_PrieTayAlpha)
    if ( allocated(L1_aeroResist) ) deallocate(L1_aeroResist)
    if ( allocated(L1_surfResist) ) deallocate(L1_surfResist)
    if ( allocated(L1_fRoots) ) deallocate(L1_fRoots)
    if ( allocated(L1_maxInter) ) deallocate(L1_maxInter)
    if ( allocated(L1_kfastFlow) ) deallocate(L1_kfastFlow)
    if ( allocated(L1_kSlowFlow) ) deallocate(L1_kSlowFlow)
    if ( allocated(L1_kBaseFlow) ) deallocate(L1_kBaseFlow)
    if ( allocated(L1_kPerco) ) deallocate(L1_kPerco)
    if ( allocated(L1_soilMoistFC) ) deallocate(L1_soilMoistFC)
    if ( allocated(L1_soilMoistSat) ) deallocate(L1_soilMoistSat)
    if ( allocated(L1_soilMoistExp) ) deallocate(L1_soilMoistExp)
    if ( allocated(L1_jarvis_thresh_c1) ) deallocate(L1_jarvis_thresh_c1)
    if ( allocated(L1_tempThresh) ) deallocate(L1_tempThresh)
    if ( allocated(L1_unsatThresh) ) deallocate(L1_unsatThresh)
    if ( allocated(L1_sealedThresh) ) deallocate(L1_sealedThresh)
    if ( allocated(L1_wiltingPoint) ) deallocate(L1_wiltingPoint)
    if ( allocated(soilDB%id) ) deallocate(soilDB%id)
    if ( allocated(soilDB%nHorizons) ) deallocate(soilDB%nHorizons)
    if ( allocated(soilDB%is_present) ) deallocate(soilDB%is_present)
    if ( allocated(soilDB%UD) ) deallocate(soilDB%UD)
    if ( allocated(soilDB%LD) ) deallocate(soilDB%LD)
    if ( allocated(soilDB%clay) ) deallocate(soilDB%clay)
    if ( allocated(soilDB%sand) ) deallocate(soilDB%sand)
    if ( allocated(soilDB%DbM) ) deallocate(soilDB%DbM)
    if ( allocated(soilDB%depth) ) deallocate(soilDB%depth)
    if ( allocated(soilDB%RZdepth) ) deallocate(soilDB%RZdepth)
    if ( allocated(soilDB%Wd) ) deallocate(soilDB%Wd)
    if ( allocated(soilDB%nTillHorizons) ) deallocate(soilDB%nTillHorizons)
    if ( allocated(soilDB%thetaS_Till) ) deallocate(soilDB%thetaS_Till)
    if ( allocated(soilDB%thetaS) ) deallocate(soilDB%thetaS)
    if ( allocated(soilDB%Db) ) deallocate(soilDB%Db)
    if ( allocated(soilDB%thetaFC_Till) ) deallocate(soilDB%thetaFC_Till)
    if ( allocated(soilDB%thetaFC) ) deallocate(soilDB%thetaFC)
    if ( allocated(soilDB%thetaPW_Till) ) deallocate(soilDB%thetaPW_Till)
    if ( allocated(soilDB%thetaPW) ) deallocate(soilDB%thetaPW)
    if ( allocated(soilDB%Ks) ) deallocate(soilDB%Ks)
    if ( allocated(L1_No_Count) ) deallocate(L1_No_Count)
    if ( allocated(L1_bulkDens) ) deallocate(L1_bulkDens)
    if ( allocated(L1_latticeWater) ) deallocate(L1_latticeWater)
    if ( allocated(L1_COSMICL3) ) deallocate(L1_COSMICL3)

    ! mo_mrm_global_variables
    if ( allocated(dirGauges) ) deallocate(dirGauges)
    if ( allocated(dirTotalRunoff) ) deallocate(dirTotalRunoff)
    if ( allocated(dirBankfullRunoff) ) deallocate(dirBankfullRunoff)
    if ( allocated(level11) ) deallocate(level11)
    if ( allocated(l0_l11_remap) ) deallocate(l0_l11_remap)
    if ( allocated(mRM_runoff) ) deallocate(mRM_runoff)
    if ( allocated(domain_mrm) ) deallocate(domain_mrm)
    if ( allocated(L0_gaugeLoc) ) deallocate(L0_gaugeLoc)
    if ( allocated(L0_InflowGaugeLoc) ) deallocate(L0_InflowGaugeLoc)
    if ( allocated(L0_fAcc) ) deallocate(L0_fAcc)
    if ( allocated(L0_fDir) ) deallocate(L0_fDir)
    if ( allocated(L0_draSC) ) deallocate(L0_draSC)
    if ( allocated(L0_draCell) ) deallocate(L0_draCell)
    if ( allocated(L0_streamNet) ) deallocate(L0_streamNet)
    if ( allocated(L0_floodPlain) ) deallocate(L0_floodPlain)
    if ( allocated(L0_noutlet) ) deallocate(L0_noutlet)
    if ( allocated(L0_celerity) ) deallocate(L0_celerity)
    if ( allocated(L11_L1_ID) ) deallocate(L11_L1_ID)
    if ( allocated(L1_total_runoff_in) ) deallocate(L1_total_runoff_in)
    if ( allocated(L11_cellCoor) ) deallocate(L11_cellCoor)
    if ( allocated(L1_L11_ID) ) deallocate(L1_L11_ID)
    if ( allocated(L11_areaCell) ) deallocate(L11_areaCell)
    if ( allocated(L11_fAcc) ) deallocate(L11_fAcc)
    if ( allocated(L11_fDir) ) deallocate(L11_fDir)
    if ( allocated(L11_nOutlets) ) deallocate(L11_nOutlets)
    if ( allocated(L11_celerity) ) deallocate(L11_celerity)
    if ( allocated(L11_meandering) ) deallocate(L11_meandering)
    if ( allocated(L11_LinkIn_fAcc) ) deallocate(L11_LinkIn_fAcc)
    if ( allocated(L11_rowOut) ) deallocate(L11_rowOut)
    if ( allocated(L11_colOut) ) deallocate(L11_colOut)
    if ( allocated(L11_Qmod) ) deallocate(L11_Qmod)
    if ( allocated(L11_qOUT) ) deallocate(L11_qOUT)
    if ( allocated(L11_qTIN) ) deallocate(L11_qTIN)
    if ( allocated(L11_qTR) ) deallocate(L11_qTR)
    if ( allocated(L11_fromN) ) deallocate(L11_fromN)
    if ( allocated(L11_toN) ) deallocate(L11_toN)
    if ( allocated(L11_netPerm) ) deallocate(L11_netPerm)
    if ( allocated(L11_fRow) ) deallocate(L11_fRow)
    if ( allocated(L11_fCol) ) deallocate(L11_fCol)
    if ( allocated(L11_tRow) ) deallocate(L11_tRow)
    if ( allocated(L11_tCol) ) deallocate(L11_tCol)
    if ( allocated(L11_rOrder) ) deallocate(L11_rOrder)
    if ( allocated(L11_label) ) deallocate(L11_label)
    if ( allocated(L11_sink) ) deallocate(L11_sink)
    if ( allocated(L11_length) ) deallocate(L11_length)
    if ( allocated(L11_aFloodPlain) ) deallocate(L11_aFloodPlain)
    if ( allocated(L11_slope) ) deallocate(L11_slope)
    if ( allocated(L11_nLinkFracFPimp) ) deallocate(L11_nLinkFracFPimp)
    if ( allocated(L11_K) ) deallocate(L11_K)
    if ( allocated(L11_xi) ) deallocate(L11_xi)
    if ( allocated(L11_tsRout) ) deallocate(L11_tsRout)
    if ( allocated(L11_C1) ) deallocate(L11_C1)
    if ( allocated(L11_C2) ) deallocate(L11_C2)
    if ( allocated(L11_bankfull_runoff_in) ) deallocate(L11_bankfull_runoff_in)
    if ( allocated(L0_channel_depth) ) deallocate(L0_channel_depth)
    if ( allocated(L0_channel_elevation) ) deallocate(L0_channel_elevation)
    if ( allocated(L0_river_head_mon_sum) ) deallocate(L0_river_head_mon_sum)
    if ( allocated(mrm_L0_slope) ) deallocate(mrm_L0_slope)
    if ( allocated(gauge%domainId) ) deallocate(gauge%domainId)
    if ( allocated(gauge%fname) ) deallocate(gauge%fname)
    if ( allocated(gauge%gaugeId) ) deallocate(gauge%gaugeId)
    if ( allocated(gauge%Q) ) deallocate(gauge%Q)
    if ( allocated(gauge%T) ) deallocate(gauge%T)
    if ( allocated(InflowGauge%domainId) ) deallocate(InflowGauge%domainId)
    if ( allocated(InflowGauge%fname) ) deallocate(InflowGauge%fname)
    if ( allocated(InflowGauge%gaugeId) ) deallocate(InflowGauge%gaugeId)
    if ( allocated(InflowGauge%Q) ) deallocate(InflowGauge%Q)
    if ( allocated(InflowGauge%T) ) deallocate(InflowGauge%T)
    call riv_temp_pcs%clean_up()

    ! mo_common_mHM_mRM_variables
    if ( allocated(resolutionRouting) ) deallocate(resolutionRouting)
    if ( allocated(warmPer) ) deallocate(warmPer)
    if ( allocated(evalPer) ) deallocate(evalPer)
    if ( allocated(simPer) ) deallocate(simPer)
    if ( allocated(warmingDays) ) deallocate(warmingDays)
    if ( allocated(LCyearId) ) deallocate(LCyearId)
    if ( allocated(mhmFileRestartIn) ) deallocate(mhmFileRestartIn)
    if ( allocated(mrmFileRestartIn) ) deallocate(mrmFileRestartIn)

    ! mo_common_run_variables
    call run_cfg%clean_up()

    ! meteo handler clean up
    call meteo_handler%clean_up()

    ! coupling config clean up
    call couple_cfg%clean_up()

    ! reset verbosity
    call set_verbosity_level()

  end subroutine deallocate_global_variables

end module mo_clean_up
