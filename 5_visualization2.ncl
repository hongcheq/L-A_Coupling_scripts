;Function: Forcing and response metric, Model VS obs, Using Pearson's R
;Date: 20170622

;Note: 20170706, updated with GPCP v2.3 (previous analysis is based on GPCP v2.2)
;Note: please_find_me_figure2

begin

;tag1 = "CAM3.5"
;tag1 = "SPCAM3.5"

f_obs = addfile("/lustre/DATA/pritchard/hongcheq/Paul_metric/obs_forcing_response.nc","r")
f_CAM   = addfile("/lustre/DATA/pritchard/hongcheq/Paul_metric/CAM3.5model_forcing_response.nc","r")
f_SPCAM = addfile("/lustre/DATA/pritchard/hongcheq/Paul_metric/SPCAM3.5model_forcing_response.nc","r")

;sig_test = "yes"
;sig_test = "no"

lat_min = -60.0
lat_max = 60.0
;lat_max = 90.0

obs_forcing_sig = f_obs->forcing_sig({lat_min:lat_max},:)
obs_response_sig = f_obs->response_sig({lat_min:lat_max},:)

cam_forcing_sig = f_CAM->forcing_sig({lat_min:lat_max},:)
cam_response_sig = f_CAM->response_sig({lat_min:lat_max},:)

spcam_forcing_sig = f_SPCAM->forcing_sig({lat_min:lat_max},:)
spcam_response_sig = f_SPCAM->response_sig({lat_min:lat_max},:)

lon = f_CAM->lon
;lon@units = "degrees_east"
lat = f_CAM->lat({lat_min:lat_max})
;lat@units = "degrees_north"
;
;obs_forcing_tau_sig!0 = "lat"
;obs_forcing_tau_sig&lat = lat
;obs_forcing_tau_sig!1 = "lon"
;obs_forcing_tau_sig&lon = lon
;
;obs_response_tau_sig!0 = "lat"
;obs_response_tau_sig&lat = lat
;obs_response_tau_sig!1 = "lon"
;obs_response_tau_sig&lon = lon
;
;cam_forcing_tau_sig!0 = "lat"
;cam_forcing_tau_sig&lat = lat
;cam_forcing_tau_sig!1 = "lon"
;cam_forcing_tau_sig&lon = lon
;
;cam_response_tau_sig!0 = "lat"
;cam_response_tau_sig&lat = lat
;cam_response_tau_sig!1 = "lon"
;cam_response_tau_sig&lon = lon
;
;spcam_forcing_tau_sig!0 = "lat"
;spcam_forcing_tau_sig&lat = lat
;spcam_forcing_tau_sig!1 = "lon"
;spcam_forcing_tau_sig&lon = lon
;
;spcam_response_tau_sig!0 = "lat"
;spcam_response_tau_sig&lat = lat
;spcam_response_tau_sig!1 = "lon"
;spcam_response_tau_sig&lon = lon
;


;;======================= NaN values from matlab=========
;
;if (any(isnan_ieee(obs_forcing_tau_sig))) then 
;    obs_forcing_tau_sig@_FillValue = 1.e20 ; or whatever value you want to use 
;    replace_ieeenan (obs_forcing_tau_sig, obs_forcing_tau_sig@_FillValue, 0) 
;end if 
;
;if (any(isnan_ieee(obs_response_tau_sig))) then 
;    obs_response_tau_sig@_FillValue = 1.e20 ; or whatever value you want to use 
;    replace_ieeenan (obs_response_tau_sig, obs_response_tau_sig@_FillValue, 0) 
;end if 
;
;if (any(isnan_ieee(cam_forcing_tau_sig))) then 
;    cam_forcing_tau_sig@_FillValue = 1.e20 ; or whatever value you want to use 
;    replace_ieeenan (cam_forcing_tau_sig, cam_forcing_tau_sig@_FillValue, 0) 
;end if 
;
;if (any(isnan_ieee(cam_response_tau_sig))) then 
;    cam_response_tau_sig@_FillValue = 1.e20 ; or whatever value you want to use 
;    replace_ieeenan (cam_response_tau_sig, cam_response_tau_sig@_FillValue, 0) 
;end if 
;
;if (any(isnan_ieee(spcam_forcing_tau_sig))) then 
;    spcam_forcing_tau_sig@_FillValue = 1.e20 ; or whatever value you want to use 
;    replace_ieeenan (spcam_forcing_tau_sig, spcam_forcing_tau_sig@_FillValue, 0) 
;end if 
;
;if (any(isnan_ieee(spcam_response_tau_sig))) then 
;    spcam_response_tau_sig@_FillValue = 1.e20 ; or whatever value you want to use 
;    replace_ieeenan (spcam_response_tau_sig, spcam_response_tau_sig@_FillValue, 0) 
;end if 
;
;======================================================

;siglvl = 0.10 ; significant level
;
;rad    = 4.0*atan(1.0)/180.0
;wgty = cos(f_CAM->lat({-60:90}) * rad)
;
;obs_forcing_tau_sig_ave = wgt_areaave_Wrap(obs_forcing_tau_sig({-60:90},:),wgty,1.0,0)
;obs_response_tau_sig_ave = wgt_areaave_Wrap(obs_response_tau_sig({-60:90},:),wgty,1.0,0)
;;===============
;
;cam_forcing_tau_sig_ave = wgt_areaave_Wrap(cam_forcing_tau_sig({-60:90},:),wgty,1.0,0)
;cam_response_tau_sig_ave = wgt_areaave_Wrap(cam_response_tau_sig({-60:90},:),wgty,1.0,0)
;
;spcam_forcing_tau_sig_ave = wgt_areaave_Wrap(spcam_forcing_tau_sig({-60:90},:),wgty,1.0,0)
;spcam_response_tau_sig_ave = wgt_areaave_Wrap(spcam_response_tau_sig({-60:90},:),wgty,1.0,0)
;
;;;when you have both positive and negative value, you can not determine the amplitude buy doing the subtraction
;;diff_forcing = forcing_obs        ; metadata    
;;diff_response = response_obs          ; metadata
;;diff_forcing = forcing_model - forcing_obs
;;diff_response = response_model - response_obs
;
;
;;diff_forcing_ave = wgt_areaave_Wrap(diff_forcing({-60:90},:),wgty,1.0,0)
;;diff_response_ave = wgt_areaave_Wrap(diff_response({-60:90},:),wgty,1.0,0)
;
;============================= Pattern Correlation and RMSE ==========

clat = cos(0.01745329*lat)        ; cos(lat) as weight

p_cor_cam_forcing_sig = pattern_cor(obs_forcing_sig,cam_forcing_sig,clat,0) ; cos(lat) as weight, centered pattern correlation
rmsd_cam_forcing_sig = wgt_arearmse(obs_forcing_sig,cam_forcing_sig,clat,1.0,0)
p_cor_cam_response_sig = pattern_cor(obs_response_sig,cam_response_sig,clat,0) ; cos(lat) as weight, centered pattern correlation
rmsd_cam_response_sig = wgt_arearmse(obs_response_sig,cam_response_sig,clat,1.0,0)

p_cor_spcam_forcing_sig = pattern_cor(obs_forcing_sig,spcam_forcing_sig,clat,0) ; cos(lat) as weight, centered pattern correlation
rmsd_spcam_forcing_sig = wgt_arearmse(obs_forcing_sig,spcam_forcing_sig,clat,1.0,0)
p_cor_spcam_response_sig = pattern_cor(obs_response_sig,spcam_response_sig,clat,0) ; cos(lat) as weight, centered pattern correlation
rmsd_spcam_response_sig = wgt_arearmse(obs_response_sig,spcam_response_sig,clat,1.0,0)

;=============================
wks = gsn_open_wks("pdf","land-atmosphere_coupling_strength_forcing_response_OBS_CAM3.5_SPCAM3.5.pdf")

plot = new(6,graphic)

res = True
res@gsnDraw = False
res@gsnFrame = False
res@cnInfoLabelOn = False
res@cnFillOn = True
res@cnFillMode = "RasterFill"
res@lbLabelBarOn = False
res@cnLinesOn = False

;res@cnFillPalette = "MPL_RdBu"
;res@cnFillPalette = "temp_diff_18lev"

res@cnLevelSelectionMode =  "ExplicitLevels"   
res@cnLevels    = (/ -0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8/) 

res@mpMaxLatF = lat_max
res@mpMinLatF = lat_min

res@gsnLeftString = "GRACE,GPCP"
res@tiMainString =  "(a) Observations: Forcing metric"
plot(0) = gsn_csm_contour_map(wks,obs_forcing_sig,res)

res@gsnLeftString = "GRACE,GPCP"
res@tiMainString =  "(b) Observations: Response metric"
plot(1) = gsn_csm_contour_map(wks,obs_response_sig,res)

res@gsnLeftString = "Pattern Correlation = " + sprintf("%5.4f",p_cor_cam_forcing_sig)
res@gsnRightString = "RMSE = "+sprintf("%5.2f",rmsd_cam_forcing_sig)
res@tiMainString =  "(c) CAM3.5: Forcing metric"
plot(2) = gsn_csm_contour_map(wks,cam_forcing_sig,res)

res@gsnLeftString = "Pattern Correlation = " + sprintf("%5.4f",p_cor_cam_response_sig)
res@gsnRightString = "RMSE = "+sprintf("%5.2f",rmsd_cam_response_sig)
res@tiMainString =  "(d) CAM3.5: Response metric"
plot(3) = gsn_csm_contour_map(wks,cam_response_sig,res)

res@gsnLeftString = "Pattern Correlation = " + sprintf("%5.4f",p_cor_spcam_forcing_sig)
res@gsnRightString = "RMSE = "+sprintf("%5.2f",rmsd_spcam_forcing_sig)
res@tiMainString =  "(e) SPCAM3.5: Forcing metric"
plot(4) = gsn_csm_contour_map(wks,spcam_forcing_sig,res)

res@gsnLeftString = "Pattern Correlation = " + sprintf("%5.4f",p_cor_spcam_response_sig)
res@gsnRightString = "RMSE = "+sprintf("%5.2f",rmsd_spcam_response_sig)
res@tiMainString =  "(f) SPCAM3.5: Response metric"
plot(5) = gsn_csm_contour_map(wks,spcam_response_sig,res)
;====================================================


;================================
resP = True
resP@txString = "Land-Atmosphere Coupling on Seasonal Timescales, 2003--2013"
resP@gsnPanelLabelBar = True

resP@amJust = "TopLeft"
;resP@gsnPanelFigureStrings = (/"(a)","(b)","(c)","(d)","(e)","(f)"/)

gsn_panel(wks,plot,(/3,2/),resP)

;================== Another difference map for the reviewer's own interest =====

diff_1 = spcam_forcing_sig ; metadata
printVarSummary(obs_forcing_sig)

diff_1 = spcam_forcing_sig - cam_forcing_sig
diff_2 = spcam_response_sig 
diff_2 = spcam_response_sig - cam_response_sig

wks2 = gsn_open_wks("pdf","Fig2_Difference_map.pdf") 
res@gsnFrame = True
res@gsnDraw = True

res@lbLabelBarOn = True
res@cnLevelSelectionMode =  "ExplicitLevels"
delete(res@cnLevels)
res@cnLevels    = (/ -0.3,-0.2,-0.1,0,0.1,0.2,0.3/)

delete(res@gsnLeftString)
delete(res@gsnRightString)

res@tiMainString = "Forcing metric, SPCAM3.5 minus CAM3.5"
plot2_a = gsn_csm_contour_map(wks2,diff_1,res) 

res@tiMainString = "Response metric, SPCAM3.5 minus CAM3.5"
plot2_b = gsn_csm_contour_map(wks2,diff_2,res) 


end
