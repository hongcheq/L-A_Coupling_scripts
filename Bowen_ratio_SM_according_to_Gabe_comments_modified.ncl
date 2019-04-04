;Function: precipitation PDF differences between SP and NoSP in the ENS(writing)
;Date: 03/22/2017
;Update: 08/04/2017, ENS(W), SW US Bowen ratio as a function of H2OSOI at levsoi=0,1,...,9. in two models (one blue, one red). This is in response to Gabe's comment.
;update: SM binned from  0 to 1 at increment of 0.1, box-and-whisker plot of CAM (blue) and SPCAM(red), only levsoi = 0
;Note: please_find_me_figure12

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

;============= seperate H2OSOI and Bowen ratio into [0,0.1), [0.1,0.2), ... , [0.9,1) bins =======
if (itag .eq. 0) then

print(dimsizes(var_SM_box))
SM_cam = new((/10,dimsizes(var_SM_box)/),double)
BR_cam = new((/10,dimsizes(var_SM_box)/),double)

do bin_index = 0,9
SM_cam(bin_index,:) = where((var_SM_box .ge. (0.0+0.1*bin_index) .and. var_SM_box .lt. (0.1+0.1*bin_index)) , var_SM_box  , SM_cam@_FillValue)
BR_cam(bin_index,:) = where((var_SM_box .ge. (0.0+0.1*bin_index) .and. var_SM_box .lt. (0.1+0.1*bin_index)) , var_box  , BR_cam@_FillValue)
end do

do bin_index = 0,9
stats_CAM = stat_dispersion(BR_cam(bin_index,:),False)
IQR_CAM = stats_CAM(10) - stats_CAM(6)

CAM_yval(bin_index,0) = stats_CAM(6) - 1.5 * IQR_CAM
CAM_yval(bin_index,1) = stats_CAM(6)
CAM_yval(bin_index,2) = stats_CAM(8)
CAM_yval(bin_index,3) = stats_CAM(10)
CAM_yval(bin_index,4) = stats_CAM(10) + 1.5 * IQR_CAM
end do

end if 

if (itag .eq. 1) then
SM_spcam = new((/10,dimsizes(var_SM_box)/),double)
BR_spcam = new((/10,dimsizes(var_SM_box)/),double)

do bin_index = 0,9
SM_spcam(bin_index,:) = where((var_SM_box .ge. (0.0+0.1*bin_index) .and. var_SM_box .lt. (0.1+0.1*bin_index)) , var_SM_box  , SM_spcam@_FillValue)
BR_spcam(bin_index,:) = where((var_SM_box .ge. (0.0+0.1*bin_index) .and. var_SM_box .lt. (0.1+0.1*bin_index)) , var_box  , BR_spcam@_FillValue)
end do

do bin_index = 0,9
stats_SPCAM = stat_dispersion(BR_spcam(bin_index,:),False)
IQR_SPCAM = stats_SPCAM(10) - stats_SPCAM(6)

SPCAM_yval(bin_index,0) = stats_SPCAM(6) - 1.5 * IQR_SPCAM
SPCAM_yval(bin_index,1) = stats_SPCAM(6)
SPCAM_yval(bin_index,2) = stats_SPCAM(8)
SPCAM_yval(bin_index,3) = stats_SPCAM(10)
SPCAM_yval(bin_index,4) = stats_SPCAM(10) + 1.5 * IQR_SPCAM
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

wks = gsn_open_wks("png","Box_whisker_plot_"+levsoi_index+"_In_response_to_Gabe_Bowen_ratio_vs_SM_CAM3.5_vs_SPCAM3.5_"+location+"_"+VAR+".png")

;res = True
;res@gsnDraw = False
;res@gsnFrame = False
;
;res@tiMainString =  "ENS(W), "+location
;res@tiYAxisString = "Bowen ratio"
;res@tiXAxisString = "H2OSOI (levsoi="+levsoi_x +")"
;
;;res@trYMinF = -3.0
;;res@trYMaxF = 18.0
;;res@trXMinF = 0.0
;;res@trXMaxF = 0.5
;
;;======= add a legend
; res@pmLegendDisplayMode    = "Always"              ; turn on legend
;
; res@pmLegendSide           = "Top"                 ; Change location of 
; res@pmLegendParallelPosF   = .80                   ; move units right
; res@pmLegendOrthogonalPosF = -0.35                  ; more neg = down
;
; res@pmLegendWidthF         = 0.08                  ; Change width and
; res@pmLegendHeightF        = 0.15                  ; height of legend.
; res@lgLabelFontHeightF     = .015                   ; change font height
; res@lgPerimOn              = False                 ; no box around
;; labels for the legend
; res@xyExplicitLegendLabels = (/"CAM3.5","SPCAM3.5"/)
;
;res@xyMarkLineModes = (/"Markers","Markers"/)
;res@xyMarkers = (/16, 16/)
;res@xyMarkerSizeF = 0.003
;res@xyMarkerColors = (/"blue","red"/)
;==========================


;plot = gsn_csm_xy(wks,(/SM_CAM,SM_SPCAM/),(/BR_CAM,BR_SPCAM/),res)

;================ box plot =========
;-------- CAM3.5
res = True
res@tmXBLabels = (/"[0,0.1)","[0.1,0.2)","[0.2,0.3)","[0.3,0.4)","[0.4,0.5)","[0.5,0.6)","[0.6,0.7)","[0.7,0.8)","[0.8,0.9)","[0.9,1.0)"/)
res@tmXBLabelFontHeightF = 0.01

res@tiMainString = "(c) Bowen ratio grouped by volumetric soil moisture"
;res@tiMainString = "EN SW US"

res@tiXAxisString = "Volumetric soil moisture (0.007m)"
;res@tiXAxisString = "H2OSOI (levsoi="+levsoi_index+")"
res@tiYAxisString = "Bowen ratio"

color_str = new(10,string)
color_str = "blue"

opti = True
opti@boxColors = color_str

x_axis = (/1,6,11,16,21,26,31,36,41,46/)

plot = boxplot(wks,x_axis,CAM_yval,opti,res,False)


;--------- SPCAM3.5

x_axis2 = x_axis + 1.5

color_str = "red"
opti@boxColors = color_str

res@tmXBLabels = (/"","","","","","","","","",""/)
plot2 = boxplot(wks,x_axis2,SPCAM_yval,opti,res,False)

overlay(plot,plot2)

;txres = True
;txres@txFontHeightF = 0.012
;dum_temp3 = gsn_add_text(wks,plot,"  CAM: avg(H2OSOI_3)="+sprintf("%5.3f",SM_mean1),4.5,90.,txres)
;dum_temp4 = gsn_add_text(wks,plot,"  CAM: std(H2OSOI_3)="+sprintf("%5.3f",SM_std1),4.5,87.,txres)
;dum_temp5 = gsn_add_text(wks,plot,"  CAM: avg(Bowen_R)="+sprintf("%5.3f",BR_mean1),4.5,84.,txres)
;dum_temp6 = gsn_add_text(wks,plot,"  CAM: std(Bowen_R)="+sprintf("%5.3f",BR_std1),4.5,81.,txres)


;dum_temp7 = gsn_add_text(wks,plot,"SPCAM: avg(H2OSOI_3)="+sprintf("%5.3f",SM_mean2),4.5,75.,txres)
;dum_temp8 = gsn_add_text(wks,plot,"SPCAM: std(H2OSOI_3)="+sprintf("%5.3f",SM_std2),4.5,72.,txres)
;dum_temp9 = gsn_add_text(wks,plot,"SPCAM: avg(Bowen_R)="+sprintf("%5.3f",BR_mean2),4.5,69.,txres)
;dum_temp10 = gsn_add_text(wks,plot,"SPCAM: std(Bowen_R)="+sprintf("%5.3f",BR_std2),4.5,66.,txres)


resP = True
;resP@gsnPanelFigureStrings = (/"(b)"/)
resP@amJust = "TopRight"

;resP@gsnPanelFigureStrings = (/"(c)"/)
gsn_panel(wks,plot,(/1,1/),resP)


end do ; do levsoi_index



end
