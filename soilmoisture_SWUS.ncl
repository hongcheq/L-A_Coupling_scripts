;Function: precipitation PDF differences between SP and NoSP in the ENS(writing)
;Date: 03/22/2017
;Note: please_find_me_figure9

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

subfigure = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)"/)

var_cam = new((/10,5/),double)
var_spcam = new((/10,5/),double)

do itag = 0,1
if (itag .eq. 0) then
tag1 = "NoSP"
tag2 = "nocrm"
end if
if (itag .eq. 1) then
tag1 = "YesSP"
tag2 = "yescrm"
end if


do ilevsoi = 0,9

path1 = "/gdata/pritchard2/hongcheq/Data/GLACE/SPCAM3.5_"+tag1+year
path2 = (/"/spcam3.5_"+tag2+"_SM_"+Label,"/spcam3.5_"+tag2+"_SM_Ensemble_"+Label/)
path3 = (/"","/ens_1","/ens_2","/ens_3","/ens_4","/ens_5",\
              "/ens_6","/ens_7","/ens_8","/ens_9","/ens_10",\
              "/ens_11","/ens_12","/ens_13","/ens_14","/ens_15"/)
path4 = "/spcam3.5sld_"+tag1+".cam2.h1."+year+"-06-01-00000.nc"
path4b = "/spcam3.5sld_"+tag1+".clm2.h1."+year+"-06-01-00000.nc"

var = new((/16*92,64,128/),double)       ; ens*time,lat,lon

do ifile = 0,total-1
if (ifile .eq. 0) then
f = addfile(path1+path2(0)+path3(ifile)+path4,"r")
fl = addfile(path1+path2(0)+path3(ifile)+path4b,"r")
else
f = addfile(path1+path2(1)+path3(ifile)+path4,"r")
fl = addfile(path1+path2(1)+path3(ifile)+path4b,"r")
end if

VAR = (/"H2OSOI"/)   
VAR2 = (/"H2OSOI"/) 

land = "True"
;land = "False"

if (land .eq. "True") then
var(ifile*92:ifile*92+91,:,:) = fl->$VAR$(1:92,ilevsoi,:,:)   ; for soil moisture variables, there're four dimensions, including(time, levsoil, lat, lon)
end if

if (land .eq. "False") then
var(ifile*92:ifile*92+91,:,:) = f->$VAR$(1:92,:,:)
end if

end do ; do ifile = 0,total-1

var!1 = "lat"
var&lat = f->lat
var!2 = "lon"
var&lon = f->lon

printVarSummary(var)

var_box_1 = ndtooned(var(:,{lat1a:lat1b},{lon1a:lon1b}))  ; 1st box over SW US
var_box_3 = ndtooned(var(:,{lat3a:lat3b},{lon3a:lon3b}))  ;  box in the middle SW US
var_box_4 = ndtooned(var(:,{lat4a:lat4b},{lon4a:lon4b}))  ; rightmost box over  SW US

var_box_13 = array_append_record(var_box_1,var_box_3,0)
var_box_134 = array_append_record(var_box_13,var_box_4,0)

var_box = var_box_134 ; make it consistent with the following code
printVarSummary(var_box)

if (itag .eq. 0) then
stats_temp = stat_dispersion(var_box,False)
IQR =  stats_temp(10) - stats_temp(6)
var_cam(ilevsoi,0) = stats_temp(6) - 1.5 * IQR ;  1.5 interquartile range between the end of the whiskers and the lower quartile
var_cam(ilevsoi,1) = stats_temp(6) ; lower quartile
var_cam(ilevsoi,2) = stats_temp(8) ; median
var_cam(ilevsoi,3) = stats_temp(10) ; upper quartile
var_cam(ilevsoi,4) = stats_temp(10) + 1.5 * IQR ; maximum

end if

if (itag .eq. 1) then

stats_temp = stat_dispersion(var_box,False)
IQR =  stats_temp(10) - stats_temp(6)

var_spcam(ilevsoi,0) = stats_temp(6) - 1.5 * IQR 
var_spcam(ilevsoi,1) = stats_temp(6) ; lower quartile
var_spcam(ilevsoi,2) = stats_temp(8) ; median
var_spcam(ilevsoi,3) = stats_temp(10) ; upper quartile
var_spcam(ilevsoi,4) = stats_temp(10) + 1.5 * IQR 

end if

end do ; do ilevsoi

end do ; do itag

;print(var_cam(0,4))
;print(var_spcam(0,4))
;exit

;================= Visualization =================

;========================

wks = gsn_open_wks("pdf","CAM3.5_vs_SPCAM3.5_"+location+"_"+VAR+".pdf")

res = True
res@gsnDraw = False
res@gsnFrame = False

res@tiMainString = "Volumetric soil moisture, SW US"
;res@tiMainString = "Volumetric soil moisture, "+location
;res@tiMainString = VAR+ ", "+location
;res@tiYAxisString = var@units
res@tiYAxisString = "mm~S~3~N~/mm~S~3~N~"
;res@tiYAxisString = var@long_name + ", "+var@units
res@tiXAxisString = fl->levsoi@long_name + ", m"

res@trYMinF = -0.25
res@trYMaxF = 0.7

opti = True
color_str_cam = new(6,string)
color_str_cam = "blue"

opti@boxColors = color_str_cam

;res@xyMarkLineMode = "MarkLines"
;res@xyMarkers = (/16,6/)

;;======= add a legend
; res@pmLegendDisplayMode    = "Always"              ; turn on legend

; res@pmLegendSide           = "Top"                 ; Change location of 
; res@pmLegendParallelPosF   = .20                   ; move units right
; res@pmLegendOrthogonalPosF = -0.35                  ; more neg = down

; res@pmLegendWidthF         = 0.08                  ; Change width and
; res@pmLegendHeightF        = 0.15                  ; height of legend.
; res@lgLabelFontHeightF     = .015                   ; change font height
;res@lgPerimOn              = False                 ; no box around
;; labels for the legend
; res@xyExplicitLegendLabels = (/"CAM3.5","SPCAM3.5"/)
;;==========================
;;res@xyLineColors = (/"blue","red"/)

print(var_cam)
print(var_spcam)

x_axis = (/0,1,2,3,4,5,6,7,8,9/)

res@tmXBLabels = (/"0.007","0.028","0.062","0.119","0.212","0.366","0.620","1.038","1.728","2.865"/)
res@tmXBLabelAngleF = 90.0

plot = boxplot(wks,x_axis, var_cam,opti,res,False)

res2 = True
res2@tmXBLabels = (/"","","","","","","","","",""/)  ; so that the second would not have x axis tickmark

color_str_spcam = new(6,string)
color_str_spcam = "red"

opti@boxColors = color_str_spcam

plot2 = boxplot(wks,x_axis+0.3, var_spcam,opti,res2,False)

overlay(plot,plot2)


;plot = gsn_csm_xy(wks,fl->levsoi,(/var_cam,var_spcam/),res)

resP = True
;resP@gsnPanelFigureStrings = subfigure(ilevsoi)
;resP@gsnPanelFigureStrings = "(a)"
resP@amJust = "TopRight"

gsn_panel(wks,plot,(/1,1/),resP)





end
