;Function: Differences of JJA mean Bowen Ratio between CAM3.5, SPCAM3.5 in GLACE experiments. Mean of 16 cases in SPCAM3.5 minus mean of 16 cases in CAM3.5
;Date: 2017/07/21
;Update: only show land, and add statistical significance test.
;Note: please_find_me_figure13_response_to_reviewer's_comments

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
load "/export/home/hongcheq/Programming/NCL/funs_and_procs/Calendar.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

;system("ncrcat -D 2 *NoSPPRECT.nc -o 16member_mean_NoSPPRECT.nc")
;system("ncrcat -D 2 *YesSPPRECT.nc -o 16member_mean_YesSPPRECT.nc")

;=========================
;BR_CAM_mean = new((/64,128/),double)
;BR_SPCAM_mean = new((/64,128/),double)

min_lat = -30.0
max_lat = 60.0

   nameCAM = "16member_mean_NoSPPRECT.nc"
   fCAM = addfile(nameCAM,"r")
  
   PRECT_CAM = fCAM->PRECT(:,{min_lat:max_lat},:)
   PRECT_CAM_mean = dim_avg_n_Wrap(PRECT_CAM,0)


;   BR_CAM_mean = BR_CAM_mean / 16.0
;   BR_CAM_mean!0 = "lat"
;   BR_CAM_mean&lat = fCAM->lat
;   BR_CAM_mean!1 = "lon"
;   BR_CAM_mean&lon = fCAM->lon

printVarSummary(PRECT_CAM_mean) 

   nameSPCAM = "16member_mean_YesSPPRECT.nc"
   fSPCAM=addfile(nameSPCAM,"r")

   PRECT_SPCAM=fSPCAM->PRECT(:,{min_lat:max_lat},:)
   PRECT_SPCAM_mean = dim_avg_n_Wrap(PRECT_SPCAM,0)


;   BR_SPCAM_mean = BR_SPCAM_mean / 16.0
;   BR_SPCAM_mean!0 = "lat"
;   BR_SPCAM_mean&lat = fCAM->lat
;   BR_SPCAM_mean!1 = "lon"
;   BR_SPCAM_mean&lon = fCAM->lon
;
   Diff = PRECT_SPCAM_mean ; meta data
   Diff = PRECT_SPCAM_mean - PRECT_CAM_mean
   ;================= Student t test ========

;prob = ttest(BR_CAM_mean,dim_variance_n_Wrap(BR_CAM,0),1344,BR_SPCAM_mean,dim_variance_n_Wrap(BR_SPCAM,0),1344,False,False)

prob = ttest(PRECT_CAM_mean,dim_variance_n_Wrap(PRECT_CAM,0),1344/7,PRECT_SPCAM_mean,dim_variance_n_Wrap(PRECT_SPCAM,0),1344/7,False,False)

prob!0 = "lat"
prob&lat = PRECT_CAM&lat
prob!1 = "lon"
prob&lon = PRECT_CAM&lon

;============= Land sea mask for Diff and prob variables =======

a    = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
lsdata = a->LSMASK
lsm  = landsea_mask(lsdata,Diff&lat,Diff&lon)
Diff = mask(Diff,lsm.eq.0,False)

lsm_prob = landsea_mask(lsdata,prob&lat,prob&lon)
prob = mask(prob,lsm.eq.0,False)

;----Controling FDR to address the multiple testing problem (Wilks, 2016)---------
prob_1d = ndtooned(prob)
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
alpha_FDR = 0.10   ; (See Wilks, 2016. For weak spacial correlation,alpha_FDR = alpha_global; alpha_FDR = Alpha_FDR = 2 alpha_global)

do while (p_val(i_count) .le. (i_count+1.0)*alpha_FDR/num_total)
p_result_FDR = p_val(i_count)
i_count = i_count + 1
end do
print(p_result_FDR)


   ;============== plot Mean of 16 ensemble members in SPCAM3.5 minus Mean of 16 ensemble members in CAM3.5  ======
   wks  = gsn_open_wks("pdf","Mean_of_16members_in_SPCAM--CAM_PRECT_Mean_Diff_refined")
 
   res = True
   
   res@gsnDraw = False
   res@gsnFrame = False
   
   res@cnInfoLabelOn = False
   res@cnFillOn = True
   
   res@mpMaxLatF = max_lat
   res@mpMinLatF = min_lat

   res@tiMainString = "Difference of Ensemble Mean, PRECT, 1994 JJA, SPCAM3.5 - CAM3.5"
   res@cnLevelSelectionMode = "ManualLevels"
   res@cnMinLevelValF = -10.0
   res@cnMaxLevelValF = 10.0
   res@cnLevelSpacingF = 2.0

   plot = gsn_csm_contour_map(wks, Diff, res)

;===================== Overlay with studen t test ===========

;siglvl = 0.05
siglvl = p_result_FDR

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
plot2 = gsn_contour_shade(plot2,siglvl,1000,opt)
overlay(plot,plot2)

;================= Print some stats of the Diff ========
opt = True
opt@PrintStat = True
stats = stat_dispersion(Diff,opt)

;=============

resP = True
resP@txString = ""
resP@amJust = "TopLeft"
resP@gsnPanelFigureStrings = (/""/)

gsn_panel(wks,plot,(/1,1/),resP)










end
