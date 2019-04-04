;Function: precipitation PDF differences between SP and NoSP in the ENS(writing)
;Date: 03/22/2017
;Note: please_find_me_figure6

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

location = "location2"

Label = Label1

do itag = 0,1
if (itag .eq. 0) then
tag1 = "NoSP"
tag2 = "nocrm"
end if
if (itag .eq. 1) then
tag1 = "YesSP"
tag2 = "yescrm"
end if


var = new((/16*92,64,128/),double)       ; ens*time,lat,lon

do ifile = 0,total-1

f = addfile("/lustre/DATA/pritchard/hongcheq/GLACE_revision_PBL/"+tag1+".ENS"+ifile+".PBLH_Calculated.nc","r")

VAR = (/"PBLH"/)   
VAR2 = (/"PBLH"/) 

var(ifile*92:ifile*92+91,:,:) = f->$VAR$(:,:,:)

end do ; do ifile = 0,total-1

var!1 = "lat"
var&lat = f->lat
var!2 = "lon"
var&lon = f->lon

printVarSummary(var)

opt = True

bin_min = 100.0
bin_max = 1000.0   
;bin_n = 10
bin_n = 10
bin_width = (bin_max-bin_min)/bin_n

opt@bin_min = bin_min
opt@bin_max = bin_max

var_box = ndtooned(var(:,{lat2a:lat2b},{lon2a:lon2b}))  ; 1st box over SW US

if (itag .eq. 0) then
var_pdf1 = pdfx(var_box,bin_n,opt)
var_pdf1 = var_pdf1 / bin_width
end if
if (itag .eq. 1) then
var_pdf2 = pdfx(var_box,bin_n,opt)
var_pdf2 = var_pdf2 / bin_width
end if

end do ; itag
;================= Visualization =================

;========================

wks = gsn_open_wks("pdf","CAM3.5_vs_SPCAM3.5_"+location+"_"+VAR+"_PDF")

res = True
res@gsnDraw = False
res@gsnFrame = False

res@tiMainString = "PBL Height"
;res@tiMainString = VAR + " PDF, "+"Arabian Peninsula"
;res@tiYAxisString = "PDF, [% / ("+var@units+")]"
;res@tiXAxisString = var@units
res@tiYAxisString = "PDF, [% / m]"
res@tiXAxisString = "m"

res@trYMinF = 0.0
res@trYMaxF = 0.6

;======= add a legend
 res@pmLegendDisplayMode    = "Always"              ; turn on legend

 res@pmLegendSide           = "Top"                 ; Change location of 
 res@pmLegendParallelPosF   = .20                   ; move units right
 res@pmLegendOrthogonalPosF = -0.35                  ; more neg = down

 res@pmLegendWidthF         = 0.08                  ; Change width and
 res@pmLegendHeightF        = 0.15                  ; height of legend.
 res@lgLabelFontHeightF     = .015                   ; change font height
 res@lgPerimOn              = False                 ; no box around
; labels for the legend
 res@xyExplicitLegendLabels = (/"CAM3.5","SPCAM3.5"/)
;==========================
res@xyLineThicknesses = (/3.0,3.0/)

plot = gsn_csm_xy(wks,var_pdf1@bin_center,(/var_pdf1,var_pdf2/),res)

resP = True
;resP@gsnPanelFigureStrings = (/"(d)"/)
resP@amJust = "TopRight"

gsn_panel(wks,plot,(/1,1/),resP)






end
