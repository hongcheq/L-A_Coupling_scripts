;Date: 20160825 modified
;Function: Compute Omega_reading minus Omega_writing for the first link SM-precip
;Note: please_find_me_figure1
;Note: please_find_me_figure3

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
load "/export/home/hongcheq/Programming/NCL/funs_and_procs/Calendar.ncl"

begin

;VAR = (/"PRECT"/)
;VAR = (/"PRECC"/)
;VAR = (/"PRECL"/)
VAR = (/"TSA"/)
;VAR = (/"TS"/)
;;VAR = (/"ET"/) ; please do not use this line for calculating ET anymore, use 02***.ncl instead.
year = "1994"

tag1 = "YesSP"

low_fre = "no"
;low_fre = "yes"

if (low_fre .eq. "no") then
f1 = addfile("/gdata/pritchard2/hongcheq/Data/GLACE/Post-processing/SPCAM3.5_"+tag1+"_Ensemble_read"+year+"/Omega_"+VAR+"_read.nc","r")
f2 = addfile("/gdata/pritchard2/hongcheq/Data/GLACE/Post-processing/SPCAM3.5_"+tag1+"_Ensemble_write"+year+"/Omega_"+VAR+"_write.nc","r")
else
f1 = addfile("/gdata/pritchard2/hongcheq/Data/GLACE/Post-processing/SPCAM3.5_"+tag1+"_Ensemble_read"+year+"/Omega_"+VAR+"_read_low_fre.nc","r")
f2 = addfile("/gdata/pritchard2/hongcheq/Data/GLACE/Post-processing/SPCAM3.5_"+tag1+"_Ensemble_write"+year+"/Omega_"+VAR+"_write_low_fre.nc","r")
end if


;if (VAR .eq. "ET") then
;;Omega_RW = f1->Omega ; metadata first
;;Omega_RW  = (f1->Omega - f2->Omega ) * f2->sigma   ; SM -- ET link, see Guo paper GLACE II analysis Fig. 1
;;Omega_RW@long_name = "(Omega_ET(R)-Omega_ET(W) ) * sigma_ET(W), the first link of SM -- Precip, ie SM-ET"
;Omega_RW@units = f2->sigma@units ; Omega is dimensionless
;Omega_RW = Omega_RW  *10.0e4
;unit_label = f2->sigma@units+"*1/(10^4)"
;Omega_RW@units = f2->sigma@units + unit_label
;;unit_label = f2->sigma@units
;;Omega_RW@units =  unit_label
;else
Omega_RW = f1->Omega ; copy meta data
Omega_RW = f1->Omega - f2->Omega 
; set _FillValue to Antarctic land points
;Omega_RW({-90:-60},:) = Omega_RW@_FillValue
Omega_RW({-90:-30},:) = Omega_RW@_FillValue
Omega_RW({60:90},:) = Omega_RW@_FillValue


;end if

printVarSummary(Omega_RW)

;=========== land/sea mask ===============

a = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
lsdata = a->LSMASK
lsm  = landsea_mask(lsdata,Omega_RW&lat,Omega_RW&lon)
Omega_RW2 = mask(Omega_RW,lsm.eq.0,False) ; set ocean to _FillValue
Omega_RW2 = mask(Omega_RW2,lsm.eq.2,False) ; set lake to _FillValue
Omega_RW2 = mask(Omega_RW2,lsm.eq.4,False) ; set ice shelf to _FillValue  ; (NCL only fraction over Antarctica)

Omega_RW2!0 = "lat"
Omega_RW2&lat = Omega_RW&lat
Omega_RW2!1 = "lon"
Omega_RW2&lon = Omega_RW&lon

printVarSummary(Omega_RW2)
if (VAR .eq. "ET") then
Omega_RW2@long_name = "SPCAM3.5_"+tag1+", " + VAR + Omega_RW@long_name
else
Omega_RW2@long_name = "Omega_R-Omega_W_"+VAR+"_SPCAM3.5_"+tag1
end if
;===================================
;======================================
if (low_fre .eq. "no") then
newfile = "./Delta_Omega_R--W_"+VAR+".nc"
else
newfile = "./Delta_Omega_R--W_"+VAR+"_low_fre.nc"
end if

system("rm -f "+newfile)
ncdf = addfile(newfile,"c")

ncdf->Delta_Omega = Omega_RW2

;==================== write the global mean into a txt file ============
;===================================
if (VAR .eq. "PRECT") then
rad = 4.0*atan(1.0)/180.0
clat = cos(Omega_RW2&lat*rad)

Omega_RW2_ave = wgt_areaave(Omega_RW2,clat,1.0,0)
print(Omega_RW2_ave)

if (low_fre .eq. "no") then
asciiwrite("LOG_SM-precip_strength_01SM-Precip_Diff_Omega_Reading_Writing"+low_fre+".txt",Omega_RW2_ave)
else
asciiwrite("LOG_SM-precip_strength_01SM-Precip_Diff_Omega_Reading_Writing"+low_fre+"_low_fre.txt",Omega_RW2_ave)
end if

end if
;========================================
wks = gsn_open_wks("pdf","Omega_R-Omega_W_"+VAR+"_SPCAM3.5_"+tag1+low_fre+"_low_fre.pdf")

res = True
res@cnFillOn = True
res@cnLinesOn = False

;res@mpMinLatF = -60.0
;res@mpMaxLatF = 90.0

res@mpMinLatF = -30.0
res@mpMaxLatF = 60.0

;res@mpMinLatF = 0.0
;res@mpMaxLatF = 50.0
;res@mpMinLonF = 230.0
;res@mpMaxLonF = 270.0

res@cnFillPalette = "BlueWhiteOrangeRed"

res@cnLevelSelectionMode = "ExplicitLevels"

if (VAR .eq. "PRECT") then

if (low_fre .eq. "no") then
res@cnLevels = (/-0.1,-0.05,0.03,0.06,0.09,0.12,0.15,0.18/); for precipitation
else
res@cnLevels = (/0.03,0.06,0.09,0.12,0.15,0.18/); for precipitation
end if

end if
if (VAR .eq. "TS" .or. VAR .eq. "TSA") then
res@cnLevels = (/-0.1, -0.05, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3/) ; for temperature
end if
if (low_fre .eq. "no") then
res@cnFillColors = (/17, 40, 122, 153, 169, 185, 201, 217, 233/) ; for precipitation
else
res@cnFillColors = (/122, 153, 169, 185, 201, 217, 233/) ; for precipitation
end if
;;;res@cnLevels = (/0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45/) ; for SM-ET link, Delta(Omega) * sigma_ET(W)
;res@cnLevels = (/0.1,0.8,1.5,2.2,2.9,3.6,4.3,5.0/)
;res@cnFillColors = (/122,153,158,163,173,183,193,213,233/) ; for SM-ET link

;res@tiMainString = "Omega(Reading) - Omega(Writing), based on "+VAR
if (VAR .eq. "ET") then
res@gsnLeftString = unit_label
else
res@gsnLeftString = ""
end if

;;;==============Greek letters title=======================
;txres = True
;txres@txFontHeightF = 0.025
;txres@txJust = "CenterLeft"
;
;txres@txFuncCode  = ":"
;
;;change the subscripts according to the variable name
;if (VAR .eq. "TS") then
;gsn_text_ndc(wks,"[ :F8:W:B::F10:Ts:E:"+"  (Read) - "+":F8:W:B::F10:Ts:E:"+"  (Write) ],   "+tag1,.20,.75,txres)
;end if
;if (VAR .eq. "TSA") then
;gsn_text_ndc(wks,"[ :F8:W:B::F10:Tsa:E:"+"  (Read) - "+":F8:W:B::F10:Tsa:E:"+"  (Write) ],   "+tag1,.20,.75,txres)
;end if
;if (VAR .eq. "PRECT") then
;gsn_text_ndc(wks,"[ :F8:W:B::F10:pre:E:"+"  (Read) - "+":F8:W:B::F10:pre:E:"+"  (Write) ],   "+tag1,.20,.75,txres)
;end if
;
;======================================
res@pmLabelBarOrthogonalPosF = 0.2   ; label bar position

res@tmXBLabelFontHeightF = 0.015
res@tmYLLabelFontHeightF = 0.015

res@gsnFrame = False
res@gsnDraw = False

res@cnLineLabelsOn = False ; no contour labels
;res@lbLabelBarOn = False ; no color label bar

res@mpOutlineBoundarySets       = "National"

res@gsnLeftString = "(b) SPCAM3.5"
;res@gsnLeftString = "(d) SPCAM3.5"

plot = gsn_csm_contour_map(wks,Omega_RW2,res)

;;;========== Add marker
;lat_loc1 = 31.0
;lon_loc1 = 250.0
;
;lat_loc2 = 22.0
;lon_loc2 = 48.0
;
;;lat_loc1 = 34.0
;;lon_loc1 = 257.0
;
;;lat_loc2 = 24.0
;;lon_loc2 = 46.0
;
;pmres = True
;pmres@gsMarkerColor = "black"
;pmres@gsMarkerSizeF = 0.01
;pmres@gsLineThicknessF = 6.0
;pmres@gsMarkerIndex = 4
;dum0 = gsn_add_polymarker(wks,plot,lon_loc1,lat_loc1,pmres)
;dum00 = gsn_add_polymarker(wks,plot,lon_loc2,lat_loc2,pmres)
;;=======================

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

dum = new(16,graphic)  ; each line must be associated with a unique line


respts = True
;respts@gsLineColor = "red"
respts@gsLineColor = "black"
respts@gsLineThicknessF = 2.0

;========== only certain segments of lines need to be drawn over SW US & N Mexico

SW_US_x = (/lon1a,lon1b,lon3a,lon3b,lon4a,lon4b,lon4b,lon4a,lon3b,lon3a,lon1b,lon1a,lon1a/)
SW_US_y = (/lat1a,lat1a,lat3a,lat3a,lat4a,lat4a,lat4b,lat4b,lat3b,lat3b,lat1b,lat1b,lat1a/)

if (VAR .eq. "TSA") then

do ipts = 0,3
dum(ipts) = gsn_add_polyline(wks,plot,xpts2(ipts:ipts+1),ypts2(ipts:ipts+1),respts)
end do ; ipts

do ipts = 0,11
dum(ipts+4) = gsn_add_polyline(wks,plot,SW_US_x(ipts:ipts+1),SW_US_y(ipts:ipts+1),respts)
end do ; ipts

;do ipts = 0,1
;dum(ipts+8) = gsn_add_polyline(wks,plot,xpts3(ipts:ipts+1),ypts3(ipts:ipts+1),respts)
;end do ; ipts

;do ipts = 0,2
;dum(ipts+12) = gsn_add_polyline(wks,plot,xpts4(ipts:ipts+1),ypts4(ipts:ipts+1),respts)
;end do ; ipts

end if ; if VAR
;======================================

resP = True
;resP@gsnFrame = False
;resP@gsnDraw = False
if (VAR .eq. "PRECT") then
;resP@gsnPanelFigureStrings = (/"(b)"/)
resP@txString = ""
end if
if (VAR .eq. "TSA") then
;resP@gsnPanelFigureStrings = (/"(b)"/)
resP@txString = ""
end if

resP@amJust = "TopLeft"

gsn_panel(wks,plot,(/1,1/),resP)
;======================================
;draw(wks)
;frame(wks)

end
