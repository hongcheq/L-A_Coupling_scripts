;Function: precipitation PDF differences between SP and NoSP in the ENS(writing)
;Date: 03/22/2017
;Update: 08/04/2017, ENS(W), SW US Bowen ratio as a function of H2OSOI at levsoi=0,1,...,9. in two models (one blue, one red). This is in response to Gabe's comment.
;update: SM binned from  0 to 1 at increment of 0.1, box-and-whisker plot of CAM (blue) and SPCAM(red), only levsoi = 0
;Note: please_find_me_figure12, modifications according to Mike's suggestions

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
load "/export/home/hongcheq/Programming/NCL/funs_and_procs/Calendar.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"


begin

total=16
year = "1994"
Label1 = "write"
;Label2 = "read"


;================
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

;================

location = "SW_US"

Label = Label1

do levsoi_index = 0,0
levsoi_x = levsoi_index

CAM_yval = new((/10,5/),double)   ; number of bins X 5 points in constructing a box-plot
SPCAM_yval = new((/10,5/),double)

CAM_SM_std_yval = new(10,double)
SPCAM_SM_std_yval = new(10,double)

CAM_mean_10bins = new(10,double)
SPCAM_mean_10bins = new(10,double)

do itag = 0,1
if (itag .eq. 0) then
tag1 = "NoSP"
tag2 = "nocrm"
end if
if (itag .eq. 1) then
tag1 = "YesSP"
tag2 = "yescrm"
end if

path1 = "/gdata/pritchard2/hongcheq/Data/GLACE/SPCAM3.5_"+tag1+year
path2 = (/"/spcam3.5_"+tag2+"_SM_"+Label,"/spcam3.5_"+tag2+"_SM_Ensemble_"+Label/)
path3 = (/"","/ens_1","/ens_2","/ens_3","/ens_4","/ens_5",\
             "/ens_6","/ens_7","/ens_8","/ens_9","/ens_10",\
              "/ens_11","/ens_12","/ens_13","/ens_14","/ens_15"/)
path4 = "/spcam3.5sld_"+tag1+".cam2.h1."+year+"-06-01-00000.nc"
path4b = "/spcam3.5sld_"+tag1+".clm2.h1."+year+"-06-01-00000.nc"

;============
;================== Reading soil moisture ===========

var_SM = new((/16*92,64,128/),double)       ; ens*time,lat,lon

do ifile = 0,total-1
if (ifile .eq. 0) then
f = addfile(path1+path2(0)+path3(ifile)+path4,"r")
fl = addfile(path1+path2(0)+path3(ifile)+path4b,"r")
else
f = addfile(path1+path2(1)+path3(ifile)+path4,"r")
fl = addfile(path1+path2(1)+path3(ifile)+path4b,"r")
end if

land = "True"
;land = "False"

if (land .eq. "True") then
var_SM(ifile*92:ifile*92+91,:,:) = fl->H2OSOI(1:92,levsoi_x,:,:)   ; for soil moisture variables, there're four dimensions, including(time, levsoil, lat, lon)
end if

end do ; do ifile = 0,total-1

var_SM!1 = "lat"
var_SM&lat = f->lat
var_SM!2 = "lon"
var_SM&lon = f->lon

printVarSummary(var_SM)

;=====================================

;========================= ET calculated directly from CLM outputs. (LHFLX from CAM outputs contain land + ocean points, your location 1 contains oceanic points causing spurious high LHFLX values)

f_FCEV = addfile("/gdata/pritchard2/hongcheq/Data/GLACE/Post-processing/SPCAM3.5_"+tag1+"_Ensemble_"+Label+"1994/ens_"+Label+"FCEV_regression.nc","r")
f_FCTR = addfile("/gdata/pritchard2/hongcheq/Data/GLACE/Post-processing/SPCAM3.5_"+tag1+"_Ensemble_"+Label+"1994/ens_"+Label+"FCTR_regression.nc","r")
f_FGEV = addfile("/gdata/pritchard2/hongcheq/Data/GLACE/Post-processing/SPCAM3.5_"+tag1+"_Ensemble_"+Label+"1994/ens_"+Label+"FGEV_regression.nc","r")

f_lonlat = addfile("/lustre/DATA/pritchard/hongcheq/GLACE/Subsurface/SPCAM3.5_NoSP1994/spcam3.5_nocrm_SM_write/spcam3.5sld_NoSP.cam2.h1.1994-06-01-00000.nc","r")

var_ET = f_FCEV->var + f_FCTR->var + f_FGEV->var

var_ET!1 = "lat"
var_ET&lat = f_lonlat->lat
var_ET!2 = "lon"
var_ET&lon = f_lonlat->lon

var_ET@long_name = "total ET"
var_ET@units = f_FCEV->var@units

printVarSummary(var_ET)

;================================

f_FSH = addfile("/gdata/pritchard2/hongcheq/Data/GLACE/Post-processing/SPCAM3.5_"+tag1+"_Ensemble_"+Label+"1994/ens_"+Label+"FSH_regression.nc","r")

var_SH = f_FSH->var

var_SH!1 = "lat"
var_SH&lat = f_lonlat->lat
var_SH!2 = "lon"
var_SH&lon = f_lonlat->lon

printVarSummary(var_SH)


;var = var_SH / where(var_ET .ne. 0. , var_ET, var_ET@_FillValue)  ; Bowen ratio = sensible heat flux / latent heat flux

var = var_SH / where(abs(var_ET) .gt. 10.0 , var_ET, var_ET@_FillValue)  ; Bowen ratio = sensible heat flux / latent heat flux, only |ET|>10 watts/m2 are calculated

var!1 = "lat"
var&lat = f_lonlat->lat
var!2 = "lon"
var&lon = f_lonlat->lon

var@units = "unitless"
var@long_name = "Bowen Ratio"

VAR = "BowenRatio"

printVarSummary(var)


;=========================
opt = True

bin_min = -1.0
bin_max = 6.0
bin_num = 30
bin_width = ( bin_max - bin_min  ) / bin_num

opt@bin_min = bin_min
opt@bin_max = bin_max

box_num = bin_num

;opt@bin_min = 0.0
;opt@bin_max = 2.0

var_box_1 = ndtooned(var(:,{lat1a:lat1b},{lon1a:lon1b}))  ; 1st box over SW US
var_box_3 = ndtooned(var(:,{lat3a:lat3b},{lon3a:lon3b}))  ; 1st box over SW US
var_box_4 = ndtooned(var(:,{lat4a:lat4b},{lon4a:lon4b}))  ; 1st box over SW US

var_box_13 = array_append_record(var_box_1,var_box_3,0)
var_box_134 = array_append_record(var_box_13,var_box_4,0)

var_box = var_box_134

var_SM_box_1 = ndtooned(var_SM(:,{lat1a:lat1b},{lon1a:lon1b}))
var_SM_box_3 = ndtooned(var_SM(:,{lat3a:lat3b},{lon3a:lon3b}))
var_SM_box_4 = ndtooned(var_SM(:,{lat4a:lat4b},{lon4a:lon4b}))

var_SM_box_13 = array_append_record(var_SM_box_1,var_SM_box_3,0)
var_SM_box_134 = array_append_record(var_SM_box_13,var_SM_box_4,0)

var_SM_box = var_SM_box_134

;var_box = where((var_SM_box .ge. 0.15 .and. var_SM_box .le. 0.25), var_box, var_box@_FillValue)

;var_SM_box = where((var_SM_box .ge. 0.15 .and. var_SM_box .le. 0.25), var_SM_box, var_SM_box@_FillValue)

if (itag .eq. 0) then
var_pdf1 = pdfx(var_box,box_num,opt)   ; this pdf function actually calculates "relative frequency", you need to further divide it by the bin width of X-axis to obtain PDF so that the integral can be 1. And the Y-axis units becomes [% / (units of X-axis)], in this case, its only [%]
var_pdf1 = var_pdf1 / bin_width

end if
if (itag .eq. 1) then
var_pdf2 = pdfx(var_box,box_num,opt)
var_pdf2 = var_pdf2 / bin_width
end if

;============= seperate H2OSOI and Bowen ratio into [p0,p10), [p10,p20), ... , [p90,p100) bins =======

var_SM_box_temp = var_SM_box(ind(.not.ismissing(var_SM_box)))   ; only retain non-missing values

printVarSummary(var_SM_box)
printVarSummary(var_SM_box_temp)     ; actually you do not have missing values here. var_SM_box_temp and var_SM_box are the same size.

n_sample = dimsizes(var_SM_box_temp)

print(n_sample)
; 30912

qsort(var_SM_box_temp)

p0 = var_SM_box_temp(0)
p10 = var_SM_box_temp(n_sample * 10 / 100 - 1)
p20 = var_SM_box_temp(n_sample * 20 / 100 - 1)
p30 = var_SM_box_temp(n_sample * 30 / 100 - 1)
p40 = var_SM_box_temp(n_sample * 40 / 100 - 1)
p50 = var_SM_box_temp(n_sample * 50 / 100 - 1)
p60 = var_SM_box_temp(n_sample * 60 / 100 - 1)
p70 = var_SM_box_temp(n_sample * 70 / 100 - 1)
p80 = var_SM_box_temp(n_sample * 80 / 100 - 1)
p90 = var_SM_box_temp(n_sample * 90 / 100 - 1)
p100 = var_SM_box_temp(n_sample - 1)

temp2 = stat_dispersion(var_SM_box_temp,False)

print(p0)
print(p10)
print(p20)
print(p30)
print(p40)
print(p50)
print(p60)
print(p70)
print(p80)
print(p90)
print(p100)

print(temp2(3))   ; lower dectile 
print(temp2(8))  ; median
  ; After examination, your percentiles with 10 increment works!

if (itag .eq. 0) then
bin_boundary_cam = (/p0,p10,p20,p30,p40,p50,p60,p70,p80,p90,p100/)
CAM_BR_134_mean = avg(var_box)
CAM_BR_134_mean_gt3 = avg(var_box(ind(var_box .gt. 3.0)))
end if  ; if (itag .eq. 0)

if (itag .eq. 1) then
bin_boundary_spcam = (/p0,p10,p20,p30,p40,p50,p60,p70,p80,p90,p100/)
SPCAM_BR_134_mean = avg(var_box)
SPCAM_BR_134_mean_gt3 = avg(var_box(ind(var_box .gt. 3.0)))
end if  ; if (itag .eq. 0)

if (itag .eq. 0) then

print(dimsizes(var_SM_box))
SM_cam = new((/10,dimsizes(var_SM_box)/),double)
BR_cam = new((/10,dimsizes(var_SM_box)/),double)

do bin_index = 0,9
SM_cam(bin_index,:) = where((var_SM_box .ge. bin_boundary_cam(0+bin_index) .and. var_SM_box .lt. bin_boundary_cam(1+bin_index)) , var_SM_box  , SM_cam@_FillValue)
BR_cam(bin_index,:) = where((var_SM_box .ge. bin_boundary_cam(0+bin_index) .and. var_SM_box .lt. bin_boundary_cam(1+bin_index)) , var_box  , BR_cam@_FillValue)
end do

do bin_index = 0,9
stats_CAM = stat_dispersion(BR_cam(bin_index,:),False)
IQR_CAM = stats_CAM(10) - stats_CAM(6)

;CAM_yval(bin_index,0) = stats_CAM(2)
CAM_yval(bin_index,0) = stats_CAM(6) - 1.5 * IQR_CAM
CAM_yval(bin_index,1) = stats_CAM(6)
CAM_yval(bin_index,2) = stats_CAM(8)
;;CAM_yval(bin_index,2) = stats_CAM(0)
CAM_yval(bin_index,3) = stats_CAM(10)
CAM_yval(bin_index,4) = stats_CAM(10) + 1.5 * IQR_CAM
;CAM_yval(bin_index,4) = stats_CAM(14)

CAM_mean_10bins(bin_index) = stats_CAM(0)    ; mean of Bowen ratio for this bin in CAM

;--------- Calculate STD of soil moisture within each bin ------
stats_CAM_SM = stat_dispersion(SM_cam(bin_index,:),False)
CAM_SM_std_yval(bin_index) = stats_CAM_SM(1)

end do

end if 

if (itag .eq. 1) then
SM_spcam = new((/10,dimsizes(var_SM_box)/),double)
BR_spcam = new((/10,dimsizes(var_SM_box)/),double)

do bin_index = 0,9
SM_spcam(bin_index,:) = where((var_SM_box .ge. bin_boundary_cam(0+bin_index) .and. var_SM_box .lt. bin_boundary_cam(1+bin_index)) , var_SM_box  , SM_spcam@_FillValue)
BR_spcam(bin_index,:) = where((var_SM_box .ge. bin_boundary_cam(0+bin_index) .and. var_SM_box .lt. bin_boundary_cam(1+bin_index)) , var_box  , BR_spcam@_FillValue)
end do

do bin_index = 0,9
stats_SPCAM = stat_dispersion(BR_spcam(bin_index,:),False)
IQR_SPCAM = stats_SPCAM(10) - stats_SPCAM(6)

;SPCAM_yval(bin_index,0) = stats_SPCAM(2)
SPCAM_yval(bin_index,0) = stats_SPCAM(6) - 1.5 * IQR_SPCAM
SPCAM_yval(bin_index,1) = stats_SPCAM(6)
SPCAM_yval(bin_index,2) = stats_SPCAM(8)
;SPCAM_yval(bin_index,2) = stats_SPCAM(0)
SPCAM_yval(bin_index,3) = stats_SPCAM(10)
SPCAM_yval(bin_index,4) = stats_SPCAM(10) + 1.5 * IQR_SPCAM
;SPCAM_yval(bin_index,4) = stats_SPCAM(14)

SPCAM_mean_10bins(bin_index) = stats_SPCAM(0)    ; mean of Bowen ratio for this bin in SPCAM

;------------ Calculate STD of soil moisture within each bin -----
stats_SPCAM_SM = stat_dispersion(SM_spcam(bin_index,:),False)
SPCAM_SM_std_yval(bin_index) = stats_SPCAM_SM(1)

end do

end if


;=============== calculate mean and std of SM and Bowen ratio ========
if (itag .eq. 0 ) then
stats_SM1 = stat_dispersion(var_SM_box,False)
SM_mean1 = stats_SM1(0)
SM_std1 = stats_SM1(1)

stats_BR1 = stat_dispersion(var_box,False)
BR_mean1 = stats_BR1(0)
BR_std1 = stats_BR1(1)

;----In response to Gabe's comments Bowen ratio as a function of SM------
BR_CAM = var_box
SM_CAM = var_SM_box
;----------
end if

if (itag .eq. 1 ) then
stats_SM2 = stat_dispersion(var_SM_box,False)
SM_mean2 = stats_SM2(0)
SM_std2 = stats_SM2(1)

stats_BR2 = stat_dispersion(var_box,False)
BR_mean2 = stats_BR2(0)
BR_std2 = stats_BR2(1)

;------------- save Bowen ratio and SM for SPCAM
BR_SPCAM = var_box
SM_SPCAM = var_SM_box

end if

end do ; itag


;================= Visualization =================

;========================

wks1 = gsn_open_wks("x11","TEST_Binning_by_cam_percentiles_Box_whisker_plot_"+levsoi_index+"_In_response_to_Mike_Bowen_ratio_vs_SM_CAM3.5_vs_SPCAM3.5_"+location+"_"+VAR+"_1")

;plot = new(3,graphic)

res1 = True
res1@gsnFrame = False
res1@gsnDraw = False
res1@tiMainString = "(a)"

res1@trYMinF = 0.0
res1@trYMaxF = 1.1

res1@tiYAxisString = "Volumetric soil moisture (0.007m), mm~S~3~N~/mm~S~3~N~"
res1@tiXAxisString = "Xth percentile"

res1@xyMarkLineMode = "MarkLines"
res1@xyMarkers = (/4,6/)
res1@xyLineColors = (/"blue","red"/)
res1@xyMarkerColors = (/"blue","red"/)

x_axis_plot0 = ispan(0,100,10)

plot = gsn_csm_xy(wks1,x_axis_plot0,(/bin_boundary_cam,bin_boundary_spcam/),res1)
draw(plot)
frame(wks1)
;================ box plot =========

wks2 = gsn_open_wks("x11","TEST_Binning_by_cam_percentiles_Box_whisker_plot_"+levsoi_index+"_In_response_to_Mike_Bowen_ratio_vs_SM_CAM3.5_vs_SPCAM3.5_"+location+"_"+VAR+"_2")
;-------- CAM3.5
res = True
res@gsnFrame = False
res@gsnDraw = False

res@tmXBLabels = (/"["+sprintf("%5.3f",p0)+","+sprintf("%5.3f",p10)+")","["+sprintf("%5.3f",p10)+","+sprintf("%5.3f",p20)+")","["+sprintf("%5.3f",p20)+","+sprintf("%5.3f",p30)+")","["+sprintf("%5.3f",p30)+","+sprintf("%5.3f",p40)+")","["+sprintf("%5.3f",p40)+","+sprintf("%5.3f",p50)+")","["+sprintf("%5.3f",p50)+","+sprintf("%5.3f",p60)+")","["+sprintf("%5.3f",p60)+","+sprintf("%5.3f",p70)+")","["+sprintf("%5.3f",p70)+","+sprintf("%5.3f",p80)+")","["+sprintf("%5.3f",p80)+","+sprintf("%5.3f",p90)+")","["+sprintf("%5.3f",p90)+","+sprintf("%5.3f",p100)+")"/)
;res@tmXBLabels = (/"[0,0.1)","[0.1,0.2)","[0.2,0.3)","[0.3,0.4)","[0.4,0.5)","[0.5,0.6)","[0.6,0.7)","[0.7,0.8)","[0.8,0.9)","[0.9,1.0)"/)
res@tmXBLabelFontHeightF = 0.017
res@tmXBLabelAngleF = 45

res@tiMainString = "(b)"
;res@tiMainString = "EN SW US"

res@tiXAxisString = "Volumetric soil moisture (0.007m)"
;res@tiXAxisString = "H2OSOI (levsoi="+levsoi_index+")"
res@tiYAxisString = "Bowen ratio"

res@trYMaxF = 25.0
res@trYMinF = -5.0

color_str = new(10,string)
color_str = "blue"

opti = True
opti@boxColors = color_str

x_axis = (/1,6,11,16,21,26,31,36,41,46/)

plot2 = boxplot(wks2,x_axis,CAM_yval,opti,res,False)

;--------- SPCAM3.5

x_axis2 = x_axis + 1.5

color_str = "red"
opti@boxColors = color_str

res@tmXBLabels = (/"","","","","","","","","",""/)
plot2_temp = boxplot(wks2,x_axis2,SPCAM_yval,opti,res,False)

overlay(plot2,plot2_temp)

draw(plot2)
frame(wks2)

;------------ panel (c)

diff_mean_10bins = SPCAM_mean_10bins - CAM_mean_10bins
perc_mean_10bins = diff_mean_10bins/CAM_mean_10bins * 100.0

bar_labels = sprintf("%4.1f",perc_mean_10bins)+"%"
print(bar_labels)

print(diff_mean_10bins)

wks3 = gsn_open_wks("x11","TEST_Binning_by_cam_percentiles_Box_whisker_plot_"+levsoi_index+"_In_response_to_Mike_Bowen_ratio_vs_SM_CAM3.5_vs_SPCAM3.5_"+location+"_"+VAR+"_3")
res3 = True

res3@gsnFrame = False
res3@gsnDraw = False

res3@tiMainString = "(c)"

res3@trYMinF = -0.8
res3@trYMaxF = 0.5

res3@gsnXYBarChart = True
res3@gsnYRefLine = 0.

res3@gsnAboveYRefLineColor = "Grey"
res3@gsnBelowYRefLineColor = "Black"

res3@tmXBMode = "Explicit"
res3@tmXBLabelFontHeightF = 0.017
res3@tmXBLabelAngleF = 45

res3@tmXBValues = x_axis
res3@tmXBLabels = (/"["+sprintf("%5.3f",p0)+","+sprintf("%5.3f",p10)+")","["+sprintf("%5.3f",p10)+","+sprintf("%5.3f",p20)+")","["+sprintf("%5.3f",p20)+","+sprintf("%5.3f",p30)+")","["+sprintf("%5.3f",p30)+","+sprintf("%5.3f",p40)+")","["+sprintf("%5.3f",p40)+","+sprintf("%5.3f",p50)+")","["+sprintf("%5.3f",p50)+","+sprintf("%5.3f",p60)+")","["+sprintf("%5.3f",p60)+","+sprintf("%5.3f",p70)+")","["+sprintf("%5.3f",p70)+","+sprintf("%5.3f",p80)+")","["+sprintf("%5.3f",p80)+","+sprintf("%5.3f",p90)+")","["+sprintf("%5.3f",p90)+","+sprintf("%5.3f",p100)+")"/)

res3@tiYAxisString = "Bowen ratio difference, SPCAM-CAM"
res3@tiXAxisString = "Volumetric soil moisture (0.007m)"

plot_new = gsn_csm_xy(wks3,x_axis,diff_mean_10bins,res3)

wks4 = gsn_open_wks("x11","TEST_Binning_by_cam_percentiles_Box_whisker_plot_"+levsoi_index+"_In_response_to_Mike_Bowen_ratio_vs_SM_CAM3.5_vs_SPCAM3.5_"+location+"_"+VAR+"_4")
res3@tiMainString = "(d)"
res3@trYMinF = -12.0
res3@trYMaxF = 20.0
res3@tiYAxisString = "Bowen ratio difference in %, (SPCAM-CAM)/CAM"
plot_wks4_new = gsn_csm_xy(wks4,x_axis,perc_mean_10bins,res3)

;**********************************************************
; add text labels
;**********************************************************  
  txres                     = True                
  txres@txFontHeightF       = 0.014               
  txres@txJust              = "CenterLeft"        
  
;  do i_x = 0,9
;    gsn_text(wks3,plot_new,bar_labels(i_x),x_axis(i_x)-1.8,0.45,txres) ; add labels
;  end do
draw(plot_new)
frame(wks3)

;  do i_x = 0,9
;    gsn_text(wks4,plot_wks4_new,bar_labels(i_x),x_axis(i_x)-1.8,0.45,txres) ; add labels
;  end do
draw(plot_wks4_new)
frame(wks4)


;==============================

;
;
;  tres               = True                    ; text mods desired
;  tres@txFontHeightF = 0.011                   ; change font height
;  
;  gsn_text_ndc(wks,"Difference of Bowen ratio mean (SPCAM - CAM)",0.80,0.72, tres)
;  gsn_text_ndc(wks,"Grouped by tenth percentiles",0.80,0.70, tres)
;  gsn_text_ndc(wks,sprintf("%5.2f",diff_mean_10bins(0)),0.63,0.68, tres)
;  gsn_text_ndc(wks,sprintf("%5.2f",diff_mean_10bins(1)),0.665,0.68, tres)
;  gsn_text_ndc(wks,sprintf("%5.2f",diff_mean_10bins(2)),0.70,0.68, tres)
;  gsn_text_ndc(wks,sprintf("%5.2f",diff_mean_10bins(3)),0.735,0.68, tres)
;  gsn_text_ndc(wks,sprintf("%5.2f",diff_mean_10bins(4)),0.77,0.68, tres)
;  gsn_text_ndc(wks,sprintf("%5.2f",diff_mean_10bins(5)),0.805,0.68, tres)
;  gsn_text_ndc(wks,sprintf("%5.2f",diff_mean_10bins(6)),0.84,0.68, tres)
;  gsn_text_ndc(wks,sprintf("%5.2f",diff_mean_10bins(7)),0.875,0.68, tres)
;  gsn_text_ndc(wks,sprintf("%5.2f",diff_mean_10bins(8)),0.91,0.68, tres)
;  gsn_text_ndc(wks,sprintf("%5.2f",diff_mean_10bins(9)),0.95,0.68, tres)
;  
;  gsn_text_ndc(wks,"Over all groups",0.80,0.66, tres)
;  gsn_text_ndc(wks,sprintf("%5.2f",SPCAM_BR_134_mean - CAM_BR_134_mean),0.80,0.64, tres)
;
;======================

;============== another figure as the response to Reviewer 2 ==========
;wks5 = gsn_open_wks("x11","TEST_Binning_by_cam_percentiles_Box_whisker_plot_"+levsoi_index+"_In_response_to_Mike_Bowen_ratio_vs_SM_CAM3.5_vs_SPCAM3.5_"+location+"_"+VAR+"_5")
;res5 = True
;
;res5@gsnFrame = False
;res5@gsnDraw = False
;
;res5@tiMainString = "(f)"
;
;;res5@trYMinF = -0.8
;;res5@trYMaxF = 0.5
;
;;res5@gsnXYBarChart = True
;;res5@gsnYRefLine = 0.
;
;;res5@gsnAboveYRefLineColor = "Grey"
;;res5@gsnBelowYRefLineColor = "Black"
;
;res5@tmXBMode = "Explicit"
;res5@tmXBLabelFontHeightF = 0.017
;res5@tmXBLabelAngleF = 45
;
;res5@tmXBValues = x_axis
;res5@tmXBLabels = (/"["+sprintf("%5.3f",p0)+","+sprintf("%5.3f",p10)+")","["+sprintf("%5.3f",p10)+","+sprintf("%5.3f",p20)+")","["+sprintf("%5.3f",p20)+","+sprintf("%5.3f",p30)+")","["+sprintf("%5.3f",p30)+","+sprintf("%5.3f",p40)+")","["+sprintf("%5.3f",p40)+","+sprintf("%5.3f",p50)+")","["+sprintf("%5.3f",p50)+","+sprintf("%5.3f",p60)+")","["+sprintf("%5.3f",p60)+","+sprintf("%5.3f",p70)+")","["+sprintf("%5.3f",p70)+","+sprintf("%5.3f",p80)+")","["+sprintf("%5.3f",p80)+","+sprintf("%5.3f",p90)+")","["+sprintf("%5.3f",p90)+","+sprintf("%5.3f",p100)+")"/)
;
;res5@tiYAxisString = "Standard deviation, mm~S~3~N~/mm~S~3~N~"
;res5@tiXAxisString = "Volumetric soil moisture (0.007m)"
;
;plot_new_5 = gsn_csm_xy(wks5,x_axis,(/CAM_SM_std_yval,SPCAM_SM_std_yval/),res5)
;
;draw(plot_new_5)
;frame(wks5)

;=========== 
asciiwrite("CAM_SoilMoistre_STD_within_each_bin.txt",sprintf("%4.3E",CAM_SM_std_yval))
asciiwrite("SPCAM_SoilMoistre_STD_within_each_bin.txt",sprintf("%4.3E",SPCAM_SM_std_yval))

print(CAM_SM_std_yval)
print(SPCAM_SM_std_yval)


;============================
resP = True
;resP@gsnPanelFigureStrings = (/"(b)"/)
resP@amJust = "TopRight"

;resP@gsnPanelFigureStrings = (/"(c)"/)
;gsn_panel(wks,plot,(/1,3/),resP)
;gsn_panel(wks,plot(0:1),(/2,1/),resP)


end do ; do levsoi_index



end
