;Function: regression of Precip - SM, and SM - TSA at two chosen locations
;Date:  2017/06/19
;Note: please_find_me_figure11

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
load "/export/home/hongcheq/Programming/NCL/funs_and_procs/Calendar.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

;======================
;----------------------------------------------------------------------
; Function to attach a labelbar outside of an XY plot
;----------------------------------------------------------------------
function attach_labelbar(wks,plot,labels,colors)
local lbres, vph, vpw, nboxes
begin
  nboxes = dimsizes(labels)
  getvalues plot                      ; Get plot size for use in
    "vpHeightF" : vph                 ; creating labelbar.
    "vpWidthF"  : vpw
  end getvalues

  lbres                    = True          ; labelbar only resources
  lbres@lbAutoManage       = False          ; Necessary to control sizes
  lbres@lbPerimOn          = False
  lbres@vpWidthF           = 0.2 * vpw     ; labelbar width
  lbres@vpHeightF          = vph           ; labelbar height
  lbres@lbFillColors       = colors
  lbres@lbMonoFillPattern  = True          ; Solid fill pattern
  lbres@lbLabelFontHeightF = 0.02          ; font height. default is small
  lbres@lbLabelJust        = "CenterLeft"  ; left justify labels
  lbres@lbBoxLinesOn       = False

  lbid = gsn_create_labelbar(wks,nboxes,labels,lbres)
;
; Now, create some annotation resources indicating how we want to
; attach the labelbar to the plot.
;                                                                 
  amres                = True
  amres@amParallelPosF = 0.61     ; Move away from plot
  annoid = gsn_add_annotation(plot,lbid,amres)

  return(annoid)
end


;============ Main code
begin

tag1 = "NoSP"
tag2 = "SW_US"

;lat_loc = 31.0
;lon_loc = 250.0

Label = "write"
;Label = "read"

f_PRECT = addfile("/gdata/pritchard2/hongcheq/Data/GLACE/Post-processing/SPCAM3.5_"+tag1+"_Ensemble_"+Label+"1994/ens_"+Label+"PRECT_regression.nc","r")
f_H2OSOI = addfile("/gdata/pritchard2/hongcheq/Data/GLACE/Post-processing/SPCAM3.5_"+tag1+"_Ensemble_"+Label+"1994/ens_"+Label+"H2OSOI_regression_sfc_layer.nc","r")
;f_SOILLIQ = addfile("/gdata/pritchard2/hongcheq/Data/GLACE/Post-processing/SPCAM3.5_"+tag1+"_Ensemble_"+Label+"1994/ens_"+Label+"BTRAN_regression.nc","r")
;f_LHFLX = addfile("/gdata/pritchard2/hongcheq/Data/GLACE/Post-processing/SPCAM3.5_"+tag1+"_Ensemble_"+Label+"1994/ens_"+Label+"LHFLX_regression.nc","r")
;f_LHFLX = addfile("/gdata/pritchard2/hongcheq/Data/GLACE/Post-processing/SPCAM3.5_"+tag1+"_Ensemble_"+Label+"1994/ens_"+Label+"NLH_regression.nc","r")

f_TSA = addfile("/gdata/pritchard2/hongcheq/Data/GLACE/Post-processing/SPCAM3.5_"+tag1+"_Ensemble_"+Label+"1994/ens_"+Label+"TSA_regression.nc","r")

f_Bowen = addfile(tag1+"_SW_US_BowenRatio.nc","r")

BR = f_Bowen->BR

printVarSummary(BR)

;=================
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

H2OSOI_3_loc1_1 = f_H2OSOI->var(:,{lat1a:lat1b},{lon1a:lon1b})    ;levsoi = 3
TSA_loc1_1 = f_TSA->var(:,{lat1a:lat1b},{lon1a:lon1b})

H2OSOI_3_loc1_3 = f_H2OSOI->var(:,{lat3a:lat3b},{lon3a:lon3b})    ;levsoi = 3
TSA_loc1_3 = f_TSA->var(:,{lat3a:lat3b},{lon3a:lon3b})

H2OSOI_3_loc1_4 = f_H2OSOI->var(:,{lat4a:lat4b},{lon4a:lon4b})    ;levsoi = 3
TSA_loc1_4 = f_TSA->var(:,{lat4a:lat4b},{lon4a:lon4b})

H2OSOI_3_loc1_13 = array_append_record(ndtooned(H2OSOI_3_loc1_1),ndtooned(H2OSOI_3_loc1_3),0)
H2OSOI_3_loc1_134 = array_append_record(H2OSOI_3_loc1_13,ndtooned(H2OSOI_3_loc1_4),0)

TSA_loc1_13 = array_append_record(ndtooned(TSA_loc1_1),ndtooned(TSA_loc1_3),0)
TSA_loc1_134 = array_append_record(TSA_loc1_13,ndtooned(TSA_loc1_4),0)

H2OSOI_3_loc1_1d = H2OSOI_3_loc1_134
TSA_loc1_1d = TSA_loc1_134

printVarSummary(TSA_loc1_1d)
;============== another 2D variable capsuling TSA_loc1_1d so that Bowen -1 --- 6, can be divided ======
nlevels = 43
BR_min = 0
BR_max = 7
levels = fspan(BR_min,BR_max,nlevels)
print(levels)

data2d = new((/nlevels-1,dimsizes(TSA_loc1_1d)/),typeof(TSA_loc1_1d))
printVarSummary(data2d)

do icolor = 0,nlevels-2
   do iloc = 0,dimsizes(TSA_loc1_1d)-1
      if (.not.ismissing(BR(iloc))) then
      if (BR(iloc) .ge. levels(icolor) .and. BR(iloc) .le.levels(icolor+1)) then
          data2d(icolor,iloc) = TSA_loc1_1d(iloc)
      end if
      end if

   end do ;  do iloc
end do ; do icolors

printVarSummary(data2d)

;==============statistics of the data
stat_SM = stat_dispersion(H2OSOI_3_loc1_1d,False)
SM_mean = stat_SM(0)
SM_lower_quartile = stat_SM(6)
SM_median = stat_SM(8)
SM_upper_quartile = stat_SM(10)
SM_std = stat_SM(1)

;================= mean and std of latent heat flux =---
stat_TSA = stat_dispersion(TSA_loc1_1d,False)
TSA_mean = stat_TSA(0)
TSA_std = stat_TSA(1)

;========= Set soil moisture segment min and max to hone in

;SM_seg_min = 0.15
SM_seg_min = 0.15
SM_seg_max = 0.25

;=========== Hone in SM regimes where SP vs NoSP are hugely different =======
H2OSOI_3_loc1_1d_seg = where(H2OSOI_3_loc1_1d .ge. SM_seg_min .and. H2OSOI_3_loc1_1d .le. SM_seg_max, H2OSOI_3_loc1_1d,H2OSOI_3_loc1_1d@_FillValue)
;=====================================

;=========== Hone in SM regimes where SP vs NoSP are hugely different =======
TSA_loc1_1d_seg = where(H2OSOI_3_loc1_1d .ge. SM_seg_min .and. H2OSOI_3_loc1_1d .le. SM_seg_max, TSA_loc1_1d,TSA_loc1_1d@_FillValue)
;=====================================

rc_SM_T = regline(H2OSOI_3_loc1_1d,TSA_loc1_1d)

T_loc1 = new((/3,dimsizes(TSA_loc1_1d)/),double)

T_loc1(0,:) = TSA_loc1_1d
T_loc1(1,:) = rc_SM_T * ( H2OSOI_3_loc1_1d - rc_SM_T@xave )  + rc_SM_T@yave
;====================================

rc_SM_T_seg = regline(H2OSOI_3_loc1_1d_seg,TSA_loc1_1d_seg)

T_loc1(2,:) = rc_SM_T_seg * ( H2OSOI_3_loc1_1d_seg - rc_SM_T_seg@xave )  + rc_SM_T_seg@yave

;===============
wks = gsn_open_wks("pdf",tag1+tag2+"_Regression_precip_SM_TSA.pdf")

res = True
res@gsnFrame = False
res@gsnDraw = False

res@xyMarkLineModes = (/"Markers"/)
res@xyMarkers = 16
res@xyMarkerSizeF = 0.002
res@xyLineThicknesses = (/1/)

;=======================================
res@tiMainString = "(a) CAM3.5, SW US"
;res@tiMainString = "H2OSOI(levsoi=0) vs TSA, "+tag1+", "+tag2
;res@gsnLeftString = "R=" + sprintf("%5.3f", rc_SM_T) + ", ("+sprintf("%5.3f", rc_SM_T_seg)+ "@["+SM_seg_min+","+SM_seg_max+"])"
;res@gsnLeftString = "R=" + sprintf("%5.3f", rc_SM_T)
;res@tiXAxisString = "soil liquid moisture "+"(levsoi = "+sprintf("%5.3f",SOILLIQ_3_loc1@levsoi)+"), "+ SOILLIQ_3_loc1@units
;res@gsnRightString = "mm3/mm3"
res@tiXAxisString = "volumetric soil moisture (0.007m), mm~S~3~N~/mm~S~3~N~"
res@tiYAxisString = TSA_loc1_1@long_name+", "+TSA_loc1_1@units
res@lbLabelBarOn = False

;============== axis ranges ====
res@trXMinF = 0.0
res@trXMaxF = 0.5
res@trYMinF = 285.0
res@trYMaxF = 314.0
;=========================
gsn_define_colormap(wks,"NCV_bright")
;gsn_define_colormap(wks,"NCV_rainbow2")
;gsn_define_colormap(wks,"NCV_jaisnd")

;======================add some stats of precip =========
;res@gsnCenterString = sprintf("%5.3f",SM_mean) + "," + sprintf("%5.3f",SM_std)
;res@gsnRightString = sprintf("%5.3f",SM_lower_quartile) + "," + sprintf("%5.3f",SM_median) + "," + sprintf("%5.3f",SM_upper_quartile)

;============= Scatter plots, but colored by Bowen ratio of every sample ======
nstep = 256/nlevels
colors = ispan(nlevels,256,nstep)
printVarSummary(colors)

;colors = (/"darkgoldenrod","darkgreen","coral4", \ 
;             "darkslateblue","limegreen",     \
;             "lightgoldenrod","darkseagreen1","lightsteelblue1"/)

overlays = new(nlevels-2,graphic)

do ioverlay = 0,nlevels-2
res@xyMarkerColors = colors(ioverlay)
;res@xyMarkerColor = (/"red"/)
if (ioverlay .eq. 0) then
    plot = gsn_csm_xy(wks,H2OSOI_3_loc1_1d,data2d(ioverlay,:),res)
else
overlays(ioverlay-1) = gsn_csm_xy(wks,H2OSOI_3_loc1_1d,data2d(ioverlay,:),res)
overlay(plot,overlays(ioverlay-1))
end if

end do ; do ioverlay


delete(res@xyMarkLineModes)
delete(res@xyLineThicknesses)

res@xyDashPatterns = 1
res@xyMarkLineModes = (/"Lines","Lines"/)
res@xyLineColors = (/"blue","red"/)
res@xyLineThicknesses = (/3,3/)

;plot2 = gsn_csm_xy(wks,H2OSOI_3_loc1_1d,T_loc1(1,:),res)

; overlay(plot,plot2)


;================= color bar

labels  = new(dimsizes(colors),string)  ; labels in # nlevels 
printVarSummary(labels)

;labels = "" + ispan(BR_min,BR_max,1)   
labels = sprintf("%4.2f",fspan(BR_min,BR_max,nlevels))

;labelbar = attach_labelbar(wks,plot,labels,colors)

;=================
txres = True
txres@txFontHeightF = 0.018

;dum_temp3 = gsn_add_text(wks,plot,"avg(SM)="+sprintf("%5.3f",SM_mean),0.4,313.0,txres)
;dum_temp4 = gsn_add_text(wks,plot,"std(SM)="+sprintf("%5.3f",SM_std),0.4,312.0,txres)

;dum_temp5 = gsn_add_text(wks,plot,"avg(TSA)="+sprintf("%5.3f",TSA_mean),0.4,311.0,txres)
;dum_temp6 = gsn_add_text(wks,plot,"std(TSA)="+sprintf("%5.3f",TSA_std),0.4,310.0,txres)

;============ Panel plot =======

resP = True

;resP@txString = tag2

resP@amJust = "TopLeft"
;resP@gsnPanelFigureStrings = (/"(a)"/)

gsn_panel(wks,plot,(/1,1/),resP)


end
