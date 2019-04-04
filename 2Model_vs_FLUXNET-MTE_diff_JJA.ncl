;Function: 6 panel plots of sensible heat flux, latent heat flux, and Bowen Ratio CAM3.5, SPCAM3.5, and FLUXNET-MTE, 1996-2011 JJA mean.
;Date: 20170615
;Note: please_find_me_figure14

load "/export/home/hongcheq/Programming/NCL/funs_and_procs/Calendar.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"    
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

f_CAM = addfile("/lustre/DATA/pritchard/hongcheq/AMIP_surface_fluxes/SHFLX.LHFLX.spcam3.5sld_NoSP.cam2.h0.1996-2015.nc","r")

f_SPCAM = addfile( "/lustre/DATA/pritchard/hongcheq/AMIP_surface_fluxes/SHFLX.LHFLX.spcam3.5sld_YesSP.cam2.h0.1996-2015.nc","r")

f_MTE = addfile("/lustre/DATA/pritchard/hongcheq/AMIP_surface_fluxes/FLUXNET-MTE_latent_sensible_heat_flux_1982-2011.nc","r")

lat_min = -30.0
lat_max = 60.0

LH_crit = 10.0   ; when calculating Bowen ratio, only consider LHFLX > latent_critical

SHFLX_CAM = f_CAM->SHFLX(0:191,{lat_min:lat_max},:)  ; 1996-2011
LHFLX_CAM = f_CAM->LHFLX(0:191,{lat_min:lat_max},:)  ; 1996-2011

BR_CAM = SHFLX_CAM ; metadata
BR_CAM@units = "unitless"
BR_CAM@long_name = "Bowen Ratio (Sensible heat flux / latent heat flux)"
BR_CAM = where(abs(LHFLX_CAM) .gt. LH_crit, SHFLX_CAM / LHFLX_CAM,BR_CAM@_FillValue) 

SHFLX_SPCAM = f_SPCAM->SHFLX(0:191,{lat_min:lat_max},:)
LHFLX_SPCAM = f_SPCAM->LHFLX(0:191,{lat_min:lat_max},:)

BR_SPCAM = SHFLX_SPCAM ; metadata
BR_SPCAM@units = "unitless"
BR_SPCAM@long_name = "Bowen Ratio (Sensible heat flux / latent heat flux)"
BR_SPCAM = where(abs(LHFLX_SPCAM) .gt. LH_crit, SHFLX_SPCAM / LHFLX_SPCAM,BR_SPCAM@_FillValue) 

SHFLX_MTE = f_MTE->SH_MTE_regrid(168:359,{lat_min:lat_max},:) ; 1982-2011
LHFLX_MTE = f_MTE->LH_MTE_regrid(168:359,{lat_min:lat_max},:) ; 1982-2011

BR_MTE = SHFLX_MTE ; metadata
BR_MTE@units = "unitless"
BR_MTE@long_name = "Bowen ratio from FLUXNET-MTE"
LHFLX_MTE_replace = LHFLX_MTE ; metadata
LHFLX_MTE_replace = where(LHFLX_MTE .eq. 0.0, LHFLX_MTE@_FillValue, LHFLX_MTE)   
BR_MTE = where(abs(LHFLX_MTE) .gt. LH_crit, SHFLX_MTE / LHFLX_MTE_replace, BR_MTE@_FillValue) 

;========= Calculation of JJA mean for every single year during 1996-2011, 16 years ===

;================ CAM3.5 ========
SHFLX_CAM_16yr = new(dimsizes(SHFLX_CAM),double)
SHFLX_CAM_16yr(0,:,:) = SHFLX_CAM(0,:,:) ; metadata
printVarSummary(SHFLX_CAM_16yr)

LHFLX_CAM_16yr = new(dimsizes(LHFLX_CAM),double)
LHFLX_CAM_16yr(0,:,:) = LHFLX_CAM(0,:,:) ; metadata

BR_CAM_16yr = new(dimsizes(BR_CAM),double)
BR_CAM_16yr(0,:,:) = BR_CAM(0,:,:)  ; metadata

do iyear = 0,15
SHFLX_CAM_16yr(iyear,:,:) = (SHFLX_CAM(5+12*iyear,:,:) + SHFLX_CAM(6+12*iyear,:,:) + SHFLX_CAM(7+12*iyear,:,:) ) / 3.0
LHFLX_CAM_16yr(iyear,:,:) = (LHFLX_CAM(5+12*iyear,:,:) + LHFLX_CAM(6+12*iyear,:,:) + LHFLX_CAM(7+12*iyear,:,:) ) / 3.0
BR_CAM_16yr(iyear,:,:) = (BR_CAM(5+12*iyear,:,:) + BR_CAM(6+12*iyear,:,:) + BR_CAM(7+12*iyear,:,:) ) / 3.0
end do ; iyear

printVarSummary(SHFLX_CAM_16yr)
printVarSummary(LHFLX_CAM_16yr)
printVarSummary(BR_CAM_16yr)

;============== SPCAM3.5 ===========

SHFLX_SPCAM_16yr = new(dimsizes(SHFLX_SPCAM),double)
SHFLX_SPCAM_16yr(0,:,:) = SHFLX_SPCAM(0,:,:) ; metadata

LHFLX_SPCAM_16yr = new(dimsizes(LHFLX_SPCAM),double)
LHFLX_SPCAM_16yr(0,:,:) = LHFLX_SPCAM(0,:,:) ; metadata

BR_SPCAM_16yr = new(dimsizes(BR_SPCAM),double)
BR_SPCAM_16yr(0,:,:) = BR_SPCAM(0,:,:)  ; metadata

do iyear = 0,15
SHFLX_SPCAM_16yr(iyear,:,:) = (SHFLX_SPCAM(5+12*iyear,:,:) + SHFLX_SPCAM(6+12*iyear,:,:) + SHFLX_SPCAM(7+12*iyear,:,:) ) / 3.0
LHFLX_SPCAM_16yr(iyear,:,:) = (LHFLX_SPCAM(5+12*iyear,:,:) + LHFLX_SPCAM(6+12*iyear,:,:) + LHFLX_SPCAM(7+12*iyear,:,:) ) / 3.0
BR_SPCAM_16yr(iyear,:,:) = (BR_SPCAM(5+12*iyear,:,:) + BR_SPCAM(6+12*iyear,:,:) + BR_SPCAM(7+12*iyear,:,:) ) / 3.0
end do ; iyear

printVarSummary(SHFLX_SPCAM_16yr)
printVarSummary(LHFLX_SPCAM_16yr)
printVarSummary(BR_SPCAM_16yr)

;============== FLUXNET-MTE JJA mean for every single year in the 16 years ========================

SHFLX_MTE_16yr = new(dimsizes(SHFLX_MTE),double)
SHFLX_MTE_16yr(0,:,:) = SHFLX_MTE(0,:,:) ; metadata

LHFLX_MTE_16yr = new(dimsizes(LHFLX_MTE),double)
LHFLX_MTE_16yr(0,:,:) = LHFLX_MTE(0,:,:) ; metadata

BR_MTE_16yr = new(dimsizes(BR_MTE),double)
BR_MTE_16yr(0,:,:) = BR_MTE(0,:,:)  ; metadata

do iyear = 0,15
SHFLX_MTE_16yr(iyear,:,:) = (SHFLX_MTE(5+12*iyear,:,:) + SHFLX_MTE(6+12*iyear,:,:) + SHFLX_MTE(7+12*iyear,:,:) ) / 3.0
LHFLX_MTE_16yr(iyear,:,:) = (LHFLX_MTE(5+12*iyear,:,:) + LHFLX_MTE(6+12*iyear,:,:) + LHFLX_MTE(7+12*iyear,:,:) ) / 3.0
BR_MTE_16yr(iyear,:,:) = (BR_MTE(5+12*iyear,:,:) + BR_MTE(6+12*iyear,:,:) + BR_MTE(7+12*iyear,:,:) ) / 3.0
end do ; iyear

printVarSummary(SHFLX_MTE_16yr)
printVarSummary(LHFLX_MTE_16yr)
printVarSummary(BR_MTE_16yr)

;============= ttest of CAM3.5 vs MTE ====================

prob_CAM_SHFLX = ttest(dim_avg_n(SHFLX_CAM_16yr,0),dim_variance_n(SHFLX_CAM_16yr,0),16,dim_avg_n(SHFLX_MTE_16yr,0),dim_variance_n(SHFLX_MTE_16yr,0),16,False,False)
prob_CAM_SHFLX!0 = "lat"
prob_CAM_SHFLX&lat = SHFLX_CAM_16yr&lat
prob_CAM_SHFLX!1 = "lon"
prob_CAM_SHFLX&lon = SHFLX_CAM_16yr&lon

prob_CAM_LHFLX = ttest(dim_avg_n(LHFLX_CAM_16yr,0),dim_variance_n(LHFLX_CAM_16yr,0),16,dim_avg_n(LHFLX_MTE_16yr,0),dim_variance_n(LHFLX_MTE_16yr,0),16,False,False)
prob_CAM_LHFLX!0 = "lat"
prob_CAM_LHFLX&lat = LHFLX_CAM_16yr&lat
prob_CAM_LHFLX!1 = "lon"
prob_CAM_LHFLX&lon = LHFLX_CAM_16yr&lon

prob_CAM_BR = ttest(dim_avg_n(BR_CAM_16yr,0),dim_variance_n(BR_CAM_16yr,0),16,dim_avg_n(BR_MTE_16yr,0),dim_variance_n(BR_MTE_16yr,0),16,False,False)
prob_CAM_BR!0 = "lat"
prob_CAM_BR&lat = BR_CAM_16yr&lat
prob_CAM_BR!1 = "lon"
prob_CAM_BR&lon = BR_CAM_16yr&lon

;============ ttest of SPCAM3.5 vs MTE ===================
prob_SPCAM_SHFLX = ttest(dim_avg_n(SHFLX_SPCAM_16yr,0),dim_variance_n(SHFLX_SPCAM_16yr,0),16,dim_avg_n(SHFLX_MTE_16yr,0),dim_variance_n(SHFLX_MTE_16yr,0),16,False,False)
prob_SPCAM_SHFLX!0 = "lat"
prob_SPCAM_SHFLX&lat = SHFLX_SPCAM_16yr&lat
prob_SPCAM_SHFLX!1 = "lon"
prob_SPCAM_SHFLX&lon = SHFLX_SPCAM_16yr&lon

prob_SPCAM_LHFLX = ttest(dim_avg_n(LHFLX_SPCAM_16yr,0),dim_variance_n(LHFLX_SPCAM_16yr,0),16,dim_avg_n(LHFLX_MTE_16yr,0),dim_variance_n(LHFLX_MTE_16yr,0),16,False,False)
prob_SPCAM_LHFLX!0 = "lat"
prob_SPCAM_LHFLX&lat = LHFLX_SPCAM_16yr&lat
prob_SPCAM_LHFLX!1 = "lon"
prob_SPCAM_LHFLX&lon = LHFLX_SPCAM_16yr&lon

prob_SPCAM_BR = ttest(dim_avg_n(BR_SPCAM_16yr,0),dim_variance_n(BR_SPCAM_16yr,0),16,dim_avg_n(BR_MTE_16yr,0),dim_variance_n(BR_MTE_16yr,0),16,False,False)
prob_SPCAM_BR!0 = "lat"
prob_SPCAM_BR&lat = BR_SPCAM_16yr&lat
prob_SPCAM_BR!1 = "lon"
prob_SPCAM_BR&lon = BR_SPCAM_16yr&lon

;============== Land sea mask of prob_CAM_SHFLX, prob_CAM_LHFLX, prob_CAM_BR, prob_SPCAM_SHFLX, prob_SPCAM_LHFLX, prob_SPCAM_BR =========

a    = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
lsdata = a->LSMASK

lsm  = landsea_mask(lsdata,prob_CAM_SHFLX&lat,prob_CAM_SHFLX&lon)
prob_CAM_SHFLX = mask(prob_CAM_SHFLX,lsm.eq.0,False)

lsm  = landsea_mask(lsdata,prob_CAM_LHFLX&lat,prob_CAM_LHFLX&lon)
prob_CAM_LHFLX = mask(prob_CAM_LHFLX,lsm.eq.0,False)

lsm  = landsea_mask(lsdata,prob_CAM_BR&lat,prob_CAM_BR&lon)
prob_CAM_BR = mask(prob_CAM_BR,lsm.eq.0,False)

lsm  = landsea_mask(lsdata,prob_SPCAM_SHFLX&lat,prob_SPCAM_SHFLX&lon)
prob_SPCAM_SHFLX = mask(prob_SPCAM_SHFLX,lsm.eq.0,False)

lsm  = landsea_mask(lsdata,prob_SPCAM_LHFLX&lat,prob_SPCAM_LHFLX&lon)
prob_SPCAM_LHFLX = mask(prob_SPCAM_LHFLX,lsm.eq.0,False)

lsm  = landsea_mask(lsdata,prob_SPCAM_BR&lat,prob_SPCAM_BR&lon)
prob_SPCAM_BR = mask(prob_SPCAM_BR,lsm.eq.0,False)

;========== Calculation of JJA mean in CAM3.5 ============
SHFLX_CAM_mon = clmMonTLL(SHFLX_CAM)
SHFLX_CAM_JJA_mean = SHFLX_CAM_mon(0,:,:) ;metadata
SHFLX_CAM_JJA_mean = (SHFLX_CAM_mon(5,:,:) + SHFLX_CAM_mon(6,:,:) + SHFLX_CAM_mon(7,:,:)) / 3.0

LHFLX_CAM_mon = clmMonTLL(LHFLX_CAM)
LHFLX_CAM_JJA_mean = LHFLX_CAM_mon(0,:,:) ;metadata
LHFLX_CAM_JJA_mean = (LHFLX_CAM_mon(5,:,:) + LHFLX_CAM_mon(6,:,:) + LHFLX_CAM_mon(7,:,:)) / 3.0

BR_CAM_mon = clmMonTLL(BR_CAM)
BR_CAM_JJA_mean = BR_CAM_mon(0,:,:) ;metadata
BR_CAM_JJA_mean = (BR_CAM_mon(5,:,:) + BR_CAM_mon(6,:,:) + BR_CAM_mon(7,:,:)) / 3.0
printVarSummary(BR_CAM_JJA_mean)

;========== Calculation of JJA mean in SPCAM3.5 ===========
SHFLX_SPCAM_mon = clmMonTLL(SHFLX_SPCAM)
SHFLX_SPCAM_JJA_mean = SHFLX_SPCAM_mon(0,:,:) ;metadata
SHFLX_SPCAM_JJA_mean = (SHFLX_SPCAM_mon(5,:,:) + SHFLX_SPCAM_mon(6,:,:) + SHFLX_SPCAM_mon(7,:,:)) / 3.0

LHFLX_SPCAM_mon = clmMonTLL(LHFLX_SPCAM)
LHFLX_SPCAM_JJA_mean = LHFLX_SPCAM_mon(0,:,:) ;metadata
LHFLX_SPCAM_JJA_mean = (LHFLX_SPCAM_mon(5,:,:) + LHFLX_SPCAM_mon(6,:,:) + LHFLX_SPCAM_mon(7,:,:)) / 3.0

BR_SPCAM_mon = clmMonTLL(BR_SPCAM)
BR_SPCAM_JJA_mean = BR_SPCAM_mon(0,:,:) ;metadata
BR_SPCAM_JJA_mean = (BR_SPCAM_mon(5,:,:) + BR_SPCAM_mon(6,:,:) + BR_SPCAM_mon(7,:,:)) / 3.0
printVarSummary(BR_SPCAM_JJA_mean)

;========== Calculation of JJA mean in FLUXNET-MTE ===========
SHFLX_MTE_mon = clmMonTLL(SHFLX_MTE)
SHFLX_MTE_JJA_mean = SHFLX_MTE_mon(0,:,:) ;metadata
SHFLX_MTE_JJA_mean = (SHFLX_MTE_mon(5,:,:) + SHFLX_MTE_mon(6,:,:) + SHFLX_MTE_mon(7,:,:)) / 3.0

LHFLX_MTE_mon = clmMonTLL(LHFLX_MTE)
LHFLX_MTE_JJA_mean = LHFLX_MTE_mon(0,:,:) ;metadata
LHFLX_MTE_JJA_mean = (LHFLX_MTE_mon(5,:,:) + LHFLX_MTE_mon(6,:,:) + LHFLX_MTE_mon(7,:,:)) / 3.0

BR_MTE_mon = clmMonTLL(BR_MTE)
BR_MTE_JJA_mean = BR_MTE_mon(0,:,:) ;metadata
BR_MTE_JJA_mean = (BR_MTE_mon(5,:,:) + BR_MTE_mon(6,:,:) + BR_MTE_mon(7,:,:)) / 3.0
printVarSummary(BR_MTE_JJA_mean)

;================== Land sea mask of CAM3.5, and SPCAM3.5 variables =====
a    = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
lsdata = a->LSMASK

lsm  = landsea_mask(lsdata,SHFLX_CAM_JJA_mean&lat,SHFLX_CAM_JJA_mean&lon)
SHFLX_CAM_JJA_mean = mask(SHFLX_CAM_JJA_mean,lsm.eq.0,False)

lsm  = landsea_mask(lsdata,LHFLX_CAM_JJA_mean&lat,LHFLX_CAM_JJA_mean&lon)
LHFLX_CAM_JJA_mean = mask(LHFLX_CAM_JJA_mean,lsm.eq.0,False)

lsm  = landsea_mask(lsdata,BR_CAM_JJA_mean&lat,BR_CAM_JJA_mean&lon)
BR_CAM_JJA_mean = mask(BR_CAM_JJA_mean,lsm.eq.0,False)

lsm  = landsea_mask(lsdata,SHFLX_SPCAM_JJA_mean&lat,SHFLX_SPCAM_JJA_mean&lon)
SHFLX_SPCAM_JJA_mean = mask(SHFLX_SPCAM_JJA_mean,lsm.eq.0,False)

lsm  = landsea_mask(lsdata,LHFLX_SPCAM_JJA_mean&lat,LHFLX_SPCAM_JJA_mean&lon)
LHFLX_SPCAM_JJA_mean = mask(LHFLX_SPCAM_JJA_mean,lsm.eq.0,False)

lsm  = landsea_mask(lsdata,BR_SPCAM_JJA_mean&lat,BR_SPCAM_JJA_mean&lon)
BR_SPCAM_JJA_mean = mask(BR_SPCAM_JJA_mean,lsm.eq.0,False)

;========================== Calculation of Diff variables ===================
diff_SHFLX_CAM = SHFLX_CAM_JJA_mean ; meta data
diff_SHFLX_CAM = SHFLX_CAM_JJA_mean - SHFLX_MTE_JJA_mean

diff_LHFLX_CAM = LHFLX_CAM_JJA_mean ; meta data
diff_LHFLX_CAM = LHFLX_CAM_JJA_mean - LHFLX_MTE_JJA_mean

diff_BR_CAM = BR_CAM_JJA_mean
diff_BR_CAM = BR_CAM_JJA_mean - BR_MTE_JJA_mean

diff_SHFLX_SPCAM = SHFLX_SPCAM_JJA_mean ; meta data
diff_SHFLX_SPCAM = SHFLX_SPCAM_JJA_mean - SHFLX_MTE_JJA_mean

diff_LHFLX_SPCAM = LHFLX_SPCAM_JJA_mean ; meta data
diff_LHFLX_SPCAM = LHFLX_SPCAM_JJA_mean - LHFLX_MTE_JJA_mean

diff_BR_SPCAM = BR_SPCAM_JJA_mean
diff_BR_SPCAM = BR_SPCAM_JJA_mean - BR_MTE_JJA_mean

;=============================== Calculation of lat-weighted RMSE, global mean of SHFLX, LHFLX, Bowen ratio in CAM against FLUXNET-MTE =====
clat = cos(0.01745329*diff_SHFLX_CAM&lat)
rmsd_cam_SH = wgt_arearmse(SHFLX_MTE_JJA_mean,SHFLX_CAM_JJA_mean,clat,1.0,0)
rmsd_cam_LH = wgt_arearmse(LHFLX_MTE_JJA_mean,LHFLX_CAM_JJA_mean,clat,1.0,0)
rmsd_cam_BR = wgt_arearmse(BR_MTE_JJA_mean,BR_CAM_JJA_mean,clat,1.0,0)

mean_cam_SH = wgt_areaave(diff_SHFLX_CAM,clat,1,0)
mean_cam_LH = wgt_areaave(diff_LHFLX_CAM,clat,1,0)
mean_cam_BR = wgt_areaave(diff_BR_CAM,clat,1,0)
;============================== Calculation of lat-weighted RMSE, global mean of SHFLX, LHFLX, Bowen Ratio in SPCAM against FLUXNET-MTE ======

rmsd_spcam_SH = wgt_arearmse(SHFLX_MTE_JJA_mean,SHFLX_SPCAM_JJA_mean,clat,1.0,0)
rmsd_spcam_LH = wgt_arearmse(LHFLX_MTE_JJA_mean,LHFLX_SPCAM_JJA_mean,clat,1.0,0)
rmsd_spcam_BR = wgt_arearmse(BR_MTE_JJA_mean,BR_SPCAM_JJA_mean,clat,1.0,0)

mean_spcam_SH = wgt_areaave(diff_SHFLX_SPCAM,clat,1,0)
mean_spcam_LH = wgt_areaave(diff_LHFLX_SPCAM,clat,1,0)
mean_spcam_BR = wgt_areaave(diff_BR_SPCAM,clat,1,0)
;======== Visualization of nine panels ==========

wks = gsn_open_wks("pdf","T_test_Diff_1996-2011JJA_mean_sensible_latent_heat_flux_Bowen_ratio_CAM3.5_SPCAM3.5_FLUXNET-MTE")

plot = new(6,"graphic")

res = True
res@cnFillOn = True
res@cnLinesOn = False
res@gsnDraw = False
res@gsnFrame = False

;res@mpMinLatF = -60.0
;res@mpMaxLatF = 90.0

res@mpMinLatF = -30.0
res@mpMaxLatF = 60.0

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -30.0
res@cnMaxLevelValF = 30.0
res@cnLevelSpacingF = 6.0

res@tiMainString = "(a) CAM3.5 - FLUXNET-MTE, Sensible Heat Flux"
res@gsnLeftString = "Mean="+sprintf("%5.2f",mean_cam_SH)
res@gsnCenterString = "RMSE="+sprintf("%5.2f",rmsd_cam_SH)
res@tiMainFontHeightF   = 0.023
res@gsnStringFontHeightF = 0.02
res@gsnRightString = "W/m~S~2~N~"
plot(0) = gsn_csm_contour_map(wks,diff_SHFLX_CAM,res)

res@tiMainString = "(b) SPCAM3.5 - FLUXNET-MTE, Sensible Heat Flux"
res@gsnLeftString = "Mean="+sprintf("%5.2f",mean_spcam_SH)
res@gsnCenterString = "RMSE="+sprintf("%5.2f",rmsd_spcam_SH)
res@gsnRightString = "W/m~S~2~N~"
plot(1) = gsn_csm_contour_map(wks,diff_SHFLX_SPCAM,res)

res@cnMinLevelValF = -30.0
res@cnMaxLevelValF = 30.0
res@cnLevelSpacingF = 6.0

res@tiMainString = "(c) CAM3.5 - FLUXNET-MTE, Latent Heat Flux"
res@gsnLeftString = "Mean="+sprintf("%5.2f",mean_cam_LH)
res@gsnCenterString = "RMSE="+sprintf("%5.2f",rmsd_cam_LH)
res@gsnRightString = "W/m~S~2~N~"
plot(2) = gsn_csm_contour_map(wks,diff_LHFLX_CAM,res)

res@tiMainString = "(d) SPCAM3.5 - FLUXNET-MTE, Latent Heat Flux"
res@gsnLeftString = "Mean="+sprintf("%5.2f",mean_spcam_LH)
res@gsnCenterString = "RMSE="+sprintf("%5.2f",rmsd_spcam_LH)
res@gsnRightString = "W/m~S~2~N~"
plot(3) = gsn_csm_contour_map(wks,diff_LHFLX_SPCAM,res)

res@cnMinLevelValF = -3.0
res@cnMaxLevelValF = 3.0
res@cnLevelSpacingF = 0.6

res@tiMainString  = "(e) CAM3.5 - FLUXNET-MTE, Bowen Ratio"
res@gsnLeftString = "Mean="+sprintf("%5.2f",mean_cam_BR)
res@gsnCenterString = "RMSE="+sprintf("%5.2f",rmsd_cam_BR)
res@gsnRightString = "unitless"
plot(4) = gsn_csm_contour_map(wks,diff_BR_CAM,res)

res@tiMainString = "(f) SPCAM3.5 - FLUXNET-MTE, Bowen Ratio"
res@gsnLeftString = "Mean="+sprintf("%5.2f",mean_spcam_BR)
res@gsnCenterString = "RMSE="+sprintf("%5.2f",rmsd_spcam_BR)
res@gsnRightString = "unitless"
plot(5) = gsn_csm_contour_map(wks,diff_BR_SPCAM,res)

;====================== 
;----Controling FDR to address the multiple testing problem (Wilks, 2016)---------
p_result_FDR = new(6,double)

do i_prob_FDR = 0,5
if (i_prob_FDR .eq. 0)
prob_1d = ndtooned(prob_CAM_SHFLX)
else if (i_prob_FDR .eq. 1)
prob_1d = ndtooned(prob_SPCAM_SHFLX)
else if (i_prob_FDR .eq. 2)
prob_1d = ndtooned(prob_CAM_LHFLX)
else if (i_prob_FDR .eq. 3)
prob_1d = ndtooned(prob_SPCAM_LHFLX)
else if (i_prob_FDR .eq. 4)
prob_1d = ndtooned(prob_CAM_BR)
else if (i_prob_FDR .eq. 5)
prob_1d = ndtooned(prob_SPCAM_BR)
end if
end if
end if
end if
end if
end if

num_total = num(.not.ismissing(prob_1d))
p_val = new(num_total,double)

i_num2 = 0
do i_num = 0, dimsizes(prob_1d) - 1
if (.not.ismissing(prob_1d(i_num)))
p_val(i_num2) = prob_1d(i_num)
i_num2 = i_num2+1
end if
end do

qsort(p_val)
printVarSummary(p_val)
;print(p_val)

p_result = new(1,double)

i_count = 0
alpha_FDR = 0.10   ; (See Wilks, 2016. For weak spacial correlation,alpha_FDR = alpha_global; alpha_FDR = Alpha_FDR = 2 alpha_global)

do while (p_val(i_count) .le. (i_count+1.0)*alpha_FDR/num_total)
p_result_FDR(i_prob_FDR) = p_val(i_count)
i_count = i_count + 1
end do
print(p_result_FDR)
delete(p_val)
end do  ; do i_prob_FDR = 0,5 

;========================

;============== PLOT 2 overlay, adding significant stippling ====
;=========================== PLOT 2 overlay, adding significant stippling =======
;siglvl = 0.05

plot2  = new(6,graphic)

res2 = True

res2@gsnDraw = False
res2@gsnFrame = False

res2@cnLevelSelectionMode = "ManualLevels"
res2@cnMinLevelValF = 0.00
res2@cnMaxLevelValF = 1.05
res2@cnLevelSpacingF = 0.05

res2@cnInfoLabelOn = False
res2@cnLinesOn = False
res2@cnLineLabelsOn = False

res2@cnFillScaleF = 0.6

plot2(0) = gsn_csm_contour(wks,prob_CAM_SHFLX,res2)
opt = True
opt@gsnShadeFillType = "pattern"
opt@gsnShadeLow = 17
plot2(0) = gsn_contour_shade(plot2(0),p_result_FDR(0),30,opt)
overlay(plot(0),plot2(0))

plot2(1) = gsn_csm_contour(wks,prob_SPCAM_SHFLX,res2)
plot2(1) = gsn_contour_shade(plot2(1),p_result_FDR(1),30,opt)
overlay(plot(1),plot2(1))

plot2(2) = gsn_csm_contour(wks,prob_CAM_LHFLX,res2)
plot2(2) = gsn_contour_shade(plot2(2),p_result_FDR(2),30,opt)
overlay(plot(2),plot2(2))

plot2(3) = gsn_csm_contour(wks,prob_SPCAM_LHFLX,res2)
plot2(3) = gsn_contour_shade(plot2(3),p_result_FDR(3),30,opt)
overlay(plot(3),plot2(3))

plot2(4) = gsn_csm_contour(wks,prob_CAM_BR,res2)
plot2(4) = gsn_contour_shade(plot2(4),p_result_FDR(4),30,opt)
overlay(plot(4),plot2(4))

plot2(5) = gsn_csm_contour(wks,prob_SPCAM_BR,res2)
plot2(5) = gsn_contour_shade(plot2(5),p_result_FDR(5),30,opt)
overlay(plot(5),plot2(5))

;=========================== panel plot ==========

pres = True
pres@txString = "Anomalies Against FLUXNET-MTE, 1996-2011 JJA"

pres@amJust = "BottomLeft"
;pres@gsnPanelFigureStrings = (/"(a)","(b)","(c)","(d)","(e)","(f)"/)


gsn_panel(wks,plot,(/3,2/),pres)




end
