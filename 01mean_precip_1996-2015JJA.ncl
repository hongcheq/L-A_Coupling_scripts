;Function: Check the 1996-2015 JJA precipitation difference in CAM3.5 and SPCAM3.5 vs GPCP
;Date: 2017/01/14

;Date: 20170730
;Note: please_find_me_figure_s1

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
load "/export/home/hongcheq/Programming/NCL/funs_and_procs/Calendar.ncl"

begin

;f_gpcp = addfile("/lustre/DATA/pritchard/hongcheq/GPCP/mon/update/precip.mon.mean.nc","r")
;Note: the previous line indicates that this analysis already used newer version of GPCP v2.3, but let's hard wire to the latest version of GPCPv2.3 (you may see no changes since both of them are using GPCPv2.3)

f_gpcp = addfile("/lustre/DATA/pritchard/hongcheq/GPCP/mon/v2.3/precip.mon.mean.nc","r")

f_cam = addfile("/lustre/DATA/pritchard/hongcheq/AMIP_precipitation/PRECT.spcam3.5sld_NoSP.cam2.h0.1996-2015_JJA_monthly.nc","r")
f_spcam = addfile("/lustre/DATA/pritchard/hongcheq/AMIP_precipitation/PRECT.spcam3.5sld_YesSP.cam2.h0.1996-2015_JJA_monthly.nc","r")

precip_06_gpcp = f_gpcp->precip(209:443:12,:,:)   ; 199601--201512  only June
precip_07_gpcp = f_gpcp->precip(210:443:12,:,:)   ; 199601--201512  only July
precip_08_gpcp = f_gpcp->precip(211:443:12,:,:)   ; 199601--201512  only August
;print(Calendar(precip_06_gpcp&time))
;print(Calendar(precip_07_gpcp&time))
;print(Calendar(precip_08_gpcp&time))

precip_cam = f_cam->PRECT
precip_cam = precip_cam * 8.64e7 ; change m/sec to mm/day
precip_cam@units = "mm/day"

;print(Calendar(precip_cam&time))

precip_spcam = f_spcam->PRECT
precip_spcam = precip_spcam * 8.64e7 ; change m/sec to mm/day
precip_spcam@units = "mm/day"
;print(Calendar(precip_spcam&time))

precip_06_gpcp_mean = dim_avg_n_Wrap(precip_06_gpcp,0)
precip_07_gpcp_mean = dim_avg_n_Wrap(precip_07_gpcp,0)
precip_08_gpcp_mean = dim_avg_n_Wrap(precip_08_gpcp,0)

;printVarSummary(precip_06_gpcp_mean)
;printVarSummary(precip_07_gpcp_mean)
;printVarSummary(precip_08_gpcp_mean)

precip_JJA_gpcp = precip_06_gpcp_mean   ;  meta data
precip_JJA_gpcp = 1/3.0 * (precip_06_gpcp_mean +  precip_07_gpcp_mean + precip_08_gpcp_mean)

;printVarSummary(precip_JJA_gpcp)

precip_JJA_cam = dim_avg_n_Wrap(precip_cam,0)
;printVarSummary(precip_JJA_cam)

precip_JJA_spcam = dim_avg_n_Wrap(precip_spcam,0)
;printVarSummary(precip_JJA_spcam)

;==============================================
;============regrid GPCP to model grid=====
precip_JJA_gpcp_regrid = linint2(precip_JJA_gpcp&lon,precip_JJA_gpcp&lat,precip_JJA_gpcp,True,precip_JJA_cam&lon,precip_JJA_cam&lat,0)

precip_JJA_gpcp_regrid!0 = "lat"
precip_JJA_gpcp_regrid&lat = precip_JJA_cam&lat
precip_JJA_gpcp_regrid!1 = "lon"
precip_JJA_gpcp_regrid&lon = precip_JJA_cam&lon
precip_JJA_gpcp_regrid@units = precip_JJA_gpcp@units
precip_JJA_gpcp_regrid@long_name = "1996-2015 JJA mean precipitation"

;printVarSummary(precip_JJA_gpcp_regrid)

;=======================================
precip_cam_gpcp = precip_JJA_cam ; meta data
precip_cam_gpcp = precip_JJA_cam - precip_JJA_gpcp_regrid ; both units are mm/day


;printVarSummary(precip_cam_gpcp)

precip_spcam_gpcp = precip_JJA_spcam ; meta data
precip_spcam_gpcp = precip_JJA_spcam - precip_JJA_gpcp_regrid ; both units are mm/day
;printVarSummary(precip_spcam_gpcp)

;======================== Visualation =========

wks = gsn_open_wks("pdf","mean_precipitation_1996-2015JJA_SPCAM3.5_vs_CAM3.5_vs_GPCP")

plot = new(4,graphic)

res = True
res@gsnDraw = False
res@gsnFrame = False
res@cnLinesOn = False
res@cnInfoLabelOn = False

res@cnFillOn = True

res@mpMinLatF = 6.0
res@mpMaxLatF = 35.0
res@mpMinLonF = 33.0
res@mpMaxLonF = 64.0

res@gsnLeftString = "Precipitation"

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 0.0
res@cnMaxLevelValF = 10.0
res@cnLevelSpacingF = 1.0

res@tiMainString = "CAM3.5"
plot(0) = gsn_csm_contour_map(wks,precip_JJA_cam,res)

res@tiMainString = "SPCAM3.5"
plot(1) = gsn_csm_contour_map(wks,precip_JJA_spcam,res)

;================
;============== Add rectangles ========
lon1a = 246.0
lat1a = 27.0

lon1b = 257.0
lat1b = 35.0

lon2a = 43.0
lat2a = 16.0

lon2b = 54.0
lat2b = 25.0

lon3a = 257.0
lat3a = 24.0

lon3b = 259.0
lat3b = 33.0

lon4a = 259.0
lat4a = 23.0

lon4b = 261.0
lat4b = 31.0

xpts1 = (/lon1a,lon1b,lon1b,lon1a,lon1a/)
ypts1 = (/lat1a,lat1a,lat1b,lat1b,lat1a/)

xpts2 = (/lon2a,lon2b,lon2b,lon2a,lon2a/)
ypts2 = (/lat2a,lat2a,lat2b,lat2b,lat2a/)

xpts3 = (/lon3a,lon3b,lon3b,lon3a,lon3a/)
ypts3 = (/lat3a,lat3a,lat3b,lat3b,lat3a/)

xpts4 = (/lon4a,lon4b,lon4b,lon4a,lon4a/)
ypts4 = (/lat4a,lat4a,lat4b,lat4b,lat4a/)
;-------------------------

res@cnMinLevelValF = -5.0
res@cnMaxLevelValF = 5.0
res@cnLevelSpacingF = 1.0

res@tiMainString = "CAM3.5 - GPCP"

precip_JJA_cam_Arabia = ndtooned(precip_JJA_cam({lat2a:lat2b},{lon2a:lon2b}))
precip_JJA_gpcp_Arabia = ndtooned(precip_JJA_gpcp_regrid({lat2a:lat2b},{lon2a:lon2b}))
RMSE_CAM_GPCP = dim_rmsd(precip_JJA_cam_Arabia,precip_JJA_gpcp_Arabia)
res@gsnCenterString = "RMSE="+sprintf("%5.2f",RMSE_CAM_GPCP)

plot(2) = gsn_csm_contour_map(wks,precip_cam_gpcp,res)

;-------------------------
dum = new(16,graphic)  ; each line must be associated with a unique line

respts = True
;respts@gsLineColor = "red"
respts@gsLineColor = "black"
respts@gsLineThicknessF = 2.0

;do ipts = 0,3
;dum(ipts) = gsn_add_polyline(wks,plot(2),xpts1(ipts:ipts+1),ypts1(ipts:ipts+1),respts)
;end do ; ipts
;
do ipts = 0,3
dum(ipts+4) = gsn_add_polyline(wks,plot(2),xpts2(ipts:ipts+1),ypts2(ipts:ipts+1),respts)
end do ; ipts
;
;do ipts = 0,3
;dum(ipts+8) = gsn_add_polyline(wks,plot(2),xpts3(ipts:ipts+1),ypts3(ipts:ipts+1),respts)
;end do ; ipts
;
;do ipts = 0,3
;dum(ipts+12) = gsn_add_polyline(wks,plot(2),xpts4(ipts:ipts+1),ypts4(ipts:ipts+1),respts)
;end do ; ipts
;
;==============
precip_JJA_spcam_Arabia = ndtooned(precip_JJA_spcam({lat2a:lat2b},{lon2a:lon2b}))
RMSE_SPCAM_GPCP = dim_rmsd(precip_JJA_spcam_Arabia,precip_JJA_gpcp_Arabia)
res@gsnCenterString = "RMSE="+sprintf("%5.2f",RMSE_SPCAM_GPCP)

res@tiMainString = "SPCAM3.5 - GPCP"
plot(3) = gsn_csm_contour_map(wks,precip_spcam_gpcp,res)

;==============
dum2 = new(16,graphic)

;do ipts = 0,3
;dum2(ipts) = gsn_add_polyline(wks,plot(3),xpts1(ipts:ipts+1),ypts1(ipts:ipts+1),respts)
;end do ; ipts
; 
do ipts = 0,3
dum2(ipts+4) = gsn_add_polyline(wks,plot(3),xpts2(ipts:ipts+1),ypts2(ipts:ipts+1),respts)
end do ; ipts 
;
;do ipts = 0,3
;dum2(ipts+8) = gsn_add_polyline(wks,plot(3),xpts3(ipts:ipts+1),ypts3(ipts:ipts+1),respts)
;end do ; ipts
;
;do ipts = 0,3
;dum2(ipts+12) = gsn_add_polyline(wks,plot(3),xpts4(ipts:ipts+1),ypts4(ipts:ipts+1),respts)
;end do ; ipts
;
;;============ Panel plot -----------

resP = True
resP@txString = "Mean(precipitation), 1996-2015 JJA"

resP@gsnPanelFigureStrings = (/"(a)","(b)","(c)","(d)"/)
resP@amJust = "TopLeft"


gsn_panel(wks,plot,(/2,2/),resP)


end
