;Function: GLACE JJA ens mean precipitation, CAM3.5 vs SPCAM3.5
;Date: 2017/06/23
;Note: please_find_me_figure8

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
load "/export/home/hongcheq/Programming/NCL/funs_and_procs/Calendar.ncl"

begin

tag1 = "NoSP"
tag1 = "YesSP"

pathwrite_NoSP = "/gdata/pritchard2/hongcheq/Data/GLACE/Post-processing/SPCAM3.5_NoSP_Ensemble_write1994"
pathwrite_YesSP = "/gdata/pritchard2/hongcheq/Data/GLACE/Post-processing/SPCAM3.5_YesSP_Ensemble_write1994"
filewrite = "/ens_writePRECT.nc"

f_NoSP = addfile(pathwrite_NoSP+filewrite,"r")
f_YesSP = addfile(pathwrite_YesSP+filewrite,"r")

var_NoSP = new((/16*92,64,128/),double)
var_YesSP = new((/16*92,64,128/),double)

do iens = 0, 15
   do iday = 0, 91
     var_NoSP(iday+iens*92,:,:) = f_NoSP->var(iens,iday,:,:)
     var_YesSP(iday+iens*92,:,:) = f_YesSP->var(iens,iday,:,:)
   end do
end do

var_NoSP = var_NoSP * 8.64e7 ;
var_NoSP@units = "mm/day"
var_YesSP = var_YesSP * 8.64e7 ; 
var_YesSP@units = "mm/day"

precip_W_NoSP = dim_avg_n_Wrap(var_NoSP,0)
precip_W_YesSP = dim_avg_n_Wrap(var_YesSP,0)

printVarSummary(precip_W_NoSP)
printVarSummary(precip_W_YesSP)

diff_W = precip_W_YesSP  ; meta data  
diff_W = precip_W_YesSP - precip_W_NoSP         ; ENS(Writing), YesSP - NoSP

printVarSummary(diff_W)

;============ Student t-test ===============
;siglvl = 0.05

;prob = ttest(precip_W_NoSP,dim_variance_n_Wrap(var_NoSP,0),16*92,precip_W_YesSP,dim_variance_n_Wrap(var_YesSP,0),16*92,False,False)

; Your precipitation data is daily.
;Assuming synoptic temporal scale (weekly) precipitation are independent, reducing original sample size by a factor of 7 to account for the autocorrelatin, i.e.,  16*92/7 = 210
;prob = ttest(precip_W_NoSP,dim_variance_n_Wrap(var_NoSP,0),210,precip_W_YesSP,dim_variance_n_Wrap(var_YesSP,0),210,False,False)
prob = ttest(precip_W_NoSP,dim_variance_n_Wrap(var_NoSP,0),16*92/7,precip_W_YesSP,dim_variance_n_Wrap(var_YesSP,0),16*92/7,False,False)

prob!0 = "lat"
prob&lat = diff_W&lat
prob!1 = "lon"
prob&lon = diff_W&lon

printVarSummary(prob)
;--------------------
;----Controling FDR to address the multiple testing problem (Wilks, 2016)---------
;res@mpMinLatF = 17.0
;res@mpMaxLatF = 45.0
;res@mpMinLonF = 237.0
;res@mpMaxLonF = 270.0

SWUS_lat1 = 17.0
SWUS_lat2 = 45.0
SWUS_lon1 = 237.0
SWUS_lon2 = 270.0

prob_SWUS = prob({SWUS_lat1:SWUS_lat2},{SWUS_lon1:SWUS_lon2})

prob_1d = ndtooned(prob_SWUS)
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
alpha_FDR = 0.15   ; (See Wilks, 2016. For weak spacial correlation,alpha_FDR = alpha_global; alpha_FDR = Alpha_FDR = 2 alpha_global)

do while (p_val(i_count) .le. (i_count+1.0)*alpha_FDR/num_total)
p_result_FDR = p_val(i_count)
i_count = i_count + 1
end do
print(p_result_FDR)

;=============== Visualization
wks = gsn_open_wks("pdf","SWUS_GLACE_Diff_mean_precip_ENS_W_R_SPCAM3.5_CAM3.5")

plot = new(3,graphic)

res = True
res@gsnDraw = False
res@gsnFrame = False
res@cnLinesOn = False
res@cnInfoLabelOn = False

res@cnFillOn = True

;res@cnFillMode = "RasterFill"

res@mpMinLatF = -60.0
res@mpMaxLatF = 90.0

; US
;res@mpMinLatF = 20.0
;res@mpMaxLatF = 50.0
;res@mpMinLonF = 230.0
;res@mpMaxLonF = 270.0

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF =  0.0
res@cnMaxLevelValF =  10.0
res@cnLevelSpacingF = 1.0

res@mpOutlineBoundarySets       = "National"
res@gsnLeftString = ""
res@tiMainString = "SW US & N Mexico"

;res@mpMinLatF = 06.0
;res@mpMaxLatF = 35.0
;res@mpMinLonF = 33.0
;res@mpMaxLonF = 64.0

res@mpMinLatF = 17.0
res@mpMaxLatF = 45.0
res@mpMinLonF = 237.0
res@mpMaxLonF = 270.0

res@tiMainString = "(a) CAM3.5"
plot(0) = gsn_csm_contour_map_ce(wks,precip_W_NoSP,res)

res@tiMainString = "(b) SPCAM3.5"
plot(1) = gsn_csm_contour_map_ce(wks,precip_W_YesSP,res)

res@cnMinLevelValF = -10.0
res@cnMaxLevelValF = 10.0
res@cnLevelSpacingF = 1.0

res@tiMainString = "(c) SPCAM3.5 - CAM3.5"
plot(2) = gsn_csm_contour_map_ce(wks,diff_W,res)

;================== Add rectangles ============
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

dum = new(16,graphic)  ; each line must be associated with a unique line
dum1 = new(16,graphic)  ; each line must be associated with a unique line
dum2 = new(16,graphic)  ; each line must be associated with a unique line

respts = True
;respts@gsLineColor = "red"
respts@gsLineColor = "black"
respts@gsLineThicknessF = 2.0

;========== only certain segments of lines need to be drawn over SW US & N Mexico

SW_US_x = (/lon1a,lon1b,lon3a,lon3b,lon4a,lon4b,lon4b,lon4a,lon3b,lon3a,lon1b,lon1a,lon1a/)
SW_US_y = (/lat1a,lat1a,lat3a,lat3a,lat4a,lat4a,lat4b,lat4b,lat3b,lat3b,lat1b,lat1b,lat1a/)

do ipts = 0,3
dum(ipts) = gsn_add_polyline(wks,plot(0),xpts2(ipts:ipts+1),ypts2(ipts:ipts+1),respts)
dum1(ipts) = gsn_add_polyline(wks,plot(1),xpts2(ipts:ipts+1),ypts2(ipts:ipts+1),respts)
dum2(ipts) = gsn_add_polyline(wks,plot(2),xpts2(ipts:ipts+1),ypts2(ipts:ipts+1),respts)
end do ; ipts

do ipts = 0,11
dum(ipts+4) = gsn_add_polyline(wks,plot(0),SW_US_x(ipts:ipts+1),SW_US_y(ipts:ipts+1),respts)
dum1(ipts+4) = gsn_add_polyline(wks,plot(1),SW_US_x(ipts:ipts+1),SW_US_y(ipts:ipts+1),respts)
dum2(ipts+4) = gsn_add_polyline(wks,plot(2),SW_US_x(ipts:ipts+1),SW_US_y(ipts:ipts+1),respts)
end do ; ipts


;==============================
;
;do ipts = 0,3
;dum(ipts) = gsn_add_polyline(wks,plot(0),xpts1(ipts:ipts+1),ypts1(ipts:ipts+1),respts)
;dum1(ipts) = gsn_add_polyline(wks,plot(1),xpts1(ipts:ipts+1),ypts1(ipts:ipts+1),respts)
;dum2(ipts) = gsn_add_polyline(wks,plot(2),xpts1(ipts:ipts+1),ypts1(ipts:ipts+1),respts)
;end do ; ipts
;
;do ipts = 0,3
;dum(ipts+4) = gsn_add_polyline(wks,plot(0),xpts2(ipts:ipts+1),ypts2(ipts:ipts+1),respts)
;dum1(ipts+4) = gsn_add_polyline(wks,plot(1),xpts2(ipts:ipts+1),ypts2(ipts:ipts+1),respts)
;dum2(ipts+4) = gsn_add_polyline(wks,plot(2),xpts2(ipts:ipts+1),ypts2(ipts:ipts+1),respts)
;end do ; ipts
;
;do ipts = 0,3
;dum(ipts+8) = gsn_add_polyline(wks,plot(0),xpts3(ipts:ipts+1),ypts3(ipts:ipts+1),respts)
;dum1(ipts+8) = gsn_add_polyline(wks,plot(1),xpts3(ipts:ipts+1),ypts3(ipts:ipts+1),respts)
;dum2(ipts+8) = gsn_add_polyline(wks,plot(2),xpts3(ipts:ipts+1),ypts3(ipts:ipts+1),respts)
;end do ; ipts
;
;do ipts = 0,3
;dum(ipts+12) = gsn_add_polyline(wks,plot(0),xpts4(ipts:ipts+1),ypts4(ipts:ipts+1),respts)
;dum1(ipts+12) = gsn_add_polyline(wks,plot(1),xpts4(ipts:ipts+1),ypts4(ipts:ipts+1),respts)
;dum2(ipts+12) = gsn_add_polyline(wks,plot(2),xpts4(ipts:ipts+1),ypts4(ipts:ipts+1),respts)
;end do ; ipts
;;============== PLOT 2 overlay, adding significant stippling ====
;=========================== PLOT 2 overlay, adding significant stippling =======
;siglvl = 0.05
siglvl = p_result_FDR ; Controlling False Discovery Rate (Wilks, 2016)

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

plot2 = gsn_csm_contour(wks,prob,res2)
opt = True
opt@gsnShadeFillType = "pattern"
opt@gsnShadeLow = 17
plot2 = gsn_contour_shade(plot2,siglvl,30,opt)
overlay(plot(2),plot2)

;;============================

resP = True
resP@txString = "ENS(W), JJA mean precipitation"

resP@amJust = "TopLeft"
;resP@gsnPanelFigureStrings = (/"(a)","(b)","(c)"/)

gsn_panel(wks,plot,(/1,3/),resP)



end
