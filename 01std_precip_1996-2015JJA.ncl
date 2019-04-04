;Function: Check the 1996-2015 JJA precipitation difference in CAM3.5 and SPCAM3.5 vs GPCP
;Date: 2017/01/14

;Date: 20170706. Note: this GPCP 1-Degree Daily combinaitions v1.2 uses v2.2, not v2.3.
;Note: please_find_me_figure_s1

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
load "/export/home/hongcheq/Programming/NCL/funs_and_procs/Calendar.ncl"

begin

f_gpcp = addfile("/lustre/DATA/pritchard/hongcheq/GPCP/daily/bilinear_Regrid_64x128_GPCP_1DD_v1.2_199610-201510.nc","r")
f_cam = addfile("/lustre/DATA/pritchard/hongcheq/AMIP_precipitation/PRECT.spcam3.5sld_NoSP.cam2.h1.1996-2015_JJA_daily.nc","r")
f_spcam = addfile("/lustre/DATA/pritchard/hongcheq/AMIP_precipitation/PRECT.spcam3.5sld_YesSP.cam2.h1.1996-2015_JJA_daily.nc","r")

precip_gpcp = f_gpcp->precip
;print(Calendar(precip_gpcp&time))

p_obs = new((/19,92,64,128/),double)  ; 1997-2015 19 years of JJA

p_obs_daily = new((/19*92,64,128/),double) 

do iyear = 0,18
;the # corresponding to June 1st from year 1997-2015
interval06 = (/243,608,973,1339,1704,2069,2434,2800,3165,3530,3895,4261,4626,4991,5356,5722,6087,6452,6817/)
p_obs(iyear,:,:,:) = precip_gpcp(interval06(iyear):interval06(iyear)+91,:,:)

end do ; do iyear

do ilat = 0,63
   do ilon = 0,127
     do iyear = 0,18
          p_obs_daily(iyear*92:iyear*92+91,ilat,ilon) = p_obs(iyear,:,ilat,ilon)
     end do
   end do
end do

p_obs_daily = p_obs_daily * 8.64e7  ; change m/sec to mm/day
p_obs_daily@units = "mm/day"
;printVarSummary(p_obs_daily)


;======================
precip_cam = new((/1748,64,128/),double)
precip_spcam = new((/1748,64,128/),double)

precip_cam = f_cam->PRECT(92:1839,:,:)      ; do not use 1996 data to keep consistency with GPCP data
precip_cam = precip_cam * 8.64e7 ; change m/sec to mm/day
precip_cam@units = "mm/day"

precip_spcam = f_spcam->PRECT(92:1839,:,:) 
precip_spcam = precip_spcam * 8.64e7 ; change m/sec to mm/day
precip_spcam@units = "mm/day"
;print(Calendar(precip_spcam&time))

;==================================
precip_JJA_gpcp = dim_stddev_n_Wrap(p_obs_daily,0)

precip_JJA_cam = dim_stddev_n_Wrap(precip_cam,0)
;printVarSummary(precip_JJA_cam)

precip_JJA_spcam = dim_stddev_n_Wrap(precip_spcam,0)
;printVarSummary(precip_JJA_spcam)

;==============================================

;=======================================
precip_cam_gpcp = precip_JJA_cam ; meta data
;printVarSummary(precip_JJA_gpcp)
;printVarSummary(precip_JJA_cam)
precip_cam_gpcp = precip_JJA_cam - precip_JJA_gpcp ; both units are mm/day


;printVarSummary(precip_cam_gpcp)

precip_spcam_gpcp = precip_JJA_spcam ; meta data
precip_spcam_gpcp = precip_JJA_spcam - precip_JJA_gpcp ; both units are mm/day
;printVarSummary(precip_spcam_gpcp)

;======================== Visualation =========

wks = gsn_open_wks("pdf","standard_deviation_precipitation_1996-2015JJA_SPCAM3.5_vs_CAM3.5_vs_GPCP")

plot = new(4,graphic)

res = True
res@gsnDraw = False
res@gsnFrame = False
res@cnLinesOn = False
res@cnInfoLabelOn = False

res@cnFillOn = True

;res@mpMinLatF = -60.0
;res@mpMaxLatF = 90.0
res@mpMinLatF = 6.0
res@mpMaxLatF = 35.0
res@mpMinLonF = 33.0
res@mpMaxLonF = 64.0

res@gsnLeftString = "Precipitation"

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 0.0
res@cnMaxLevelValF = 16.0
res@cnLevelSpacingF = 2.0

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
;---------------------------------

res@cnMinLevelValF = -7.0
res@cnMaxLevelValF = 7.0
res@cnLevelSpacingF = 1.0

res@tiMainString = "CAM3.5 - GPCP"

precip_JJA_cam_Arabia = ndtooned(precip_JJA_cam({lat2a:lat2b},{lon2a:lon2b}))
precip_JJA_gpcp_Arabia = ndtooned(precip_JJA_gpcp({lat2a:lat2b},{lon2a:lon2b}))
RMSE_CAM_GPCP = dim_rmsd(precip_JJA_cam_Arabia,precip_JJA_gpcp_Arabia)
res@gsnCenterString = "RMSE="+sprintf("%5.2f",RMSE_CAM_GPCP)

plot(2) = gsn_csm_contour_map(wks,precip_cam_gpcp,res)

;----------------------------
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
;============ Panel plot -----------

resP = True
resP@txString = "STD(precip), 1997-2015 JJA"

resP = True
resP@gsnPanelFigureStrings = (/"(e)","(f)","(g)","(h)"/)
resP@amJust = "TopLeft"

gsn_panel(wks,plot,(/2,2/),resP)


end
