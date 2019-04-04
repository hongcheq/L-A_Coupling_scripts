;Function: Box plot of  JJA mean Bowen Ratio between CAM3.5, SPCAM3.5 in GLACE experiments, 16 x 16 cases
;Date: Aug 03rd, 2017
;Note: please_find_me_figure12

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
load "/export/home/hongcheq/Programming/NCL/funs_and_procs/Calendar.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

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


Diff = new((/16,16/),double)

path0 = "~/Programming/NCL/scripts/GLACE_Bowen_ratio_difference_between_SPCAM_CAM_sensitivity_to_initial_perturbations/"

; Note: in the path0, in this file Diff_SPCAM_CAM_GLACE_perturbations1.ncl, Bowen ratio is already calculated only when abs (LHFLX) > 10 watt/m2

do iCAM = 0,15
   nameCAM = path0+iCAM+"NoSPBowenRatio.nc"
   fCAM = addfile(nameCAM,"r")
   do iSPCAM = 0,15
   nameSPCAM = path0+iSPCAM+"YesSPBowenRatio.nc"
   fSPCAM=addfile(nameSPCAM,"r")

   BR_CAM = fCAM->BR
   BR_SPCAM=fSPCAM->BR

   BR_CAM_mean = dim_avg_n_Wrap(BR_CAM,0)
   BR_SPCAM_mean = dim_avg_n_Wrap(BR_SPCAM,0)

   
CAM_var_box_1 = avg(BR_CAM_mean({lat1a:lat1b},{lon1a:lon1b}))  ; 1st box over SW US
CAM_var_box_3 = avg(BR_CAM_mean({lat3a:lat3b},{lon3a:lon3b}))  ; 1st box over SW US
CAM_var_box_4 = avg(BR_CAM_mean({lat4a:lat4b},{lon4a:lon4b}))  ; 1st box over SW US

CAM_SWUS_avg = (CAM_var_box_1 + CAM_var_box_3 + CAM_var_box_4)  / 3.0


SPCAM_var_box_1 = avg(BR_SPCAM_mean({lat1a:lat1b},{lon1a:lon1b}))  ; 1st box over SW US
SPCAM_var_box_3 = avg(BR_SPCAM_mean({lat3a:lat3b},{lon3a:lon3b}))  ; 1st box over SW US
SPCAM_var_box_4 = avg(BR_SPCAM_mean({lat4a:lat4b},{lon4a:lon4b}))  ; 1st box over SW US
   
SPCAM_SWUS_avg = (SPCAM_var_box_1 + SPCAM_var_box_3 + SPCAM_var_box_4) / 3.0

   
   Diff(iCAM,iSPCAM) = SPCAM_SWUS_avg - CAM_SWUS_avg

   end do

end do   

printVarSummary(Diff)


;======= Calculation of stats for a box plot
var_diff = new((/1,5/),double)

Diff_1d = ndtooned(Diff)

stats_temp = stat_dispersion(Diff_1d,False)
IQR =  stats_temp(10) - stats_temp(6)
var_diff(0,0) = stats_temp(6) - 1.5 * IQR ; 1.5 interquartile range between the end of the whiskers and the lower quartile
var_diff(0,1) = stats_temp(6) ; lower quartile
var_diff(0,2) = stats_temp(8) ; median
var_diff(0,3) = stats_temp(10) ; upper quartile
var_diff(0,4) = stats_temp(10) + 1.5 * IQR



;============== box plot of the Diff(SPCAM - CAM) of areal & temperal mean Bowen ratio over SW US ======
wks  = gsn_open_wks("png","Diff_SPCAM_minus_CAM_Box_plot_Bowen_ratio.png")

res = True
res@tiMainString = "(b) Bowen ratio difference, SPCAM3.5 - CAM3.5"
;res@tiMainString = "Bowen ratio difference, SPCAM3.5 - CAM3.5, SW US"   
res@tmXBLabels = (/"Difference"/)

res@trYMinF = -1.5
res@trYMaxF = 1.5

x_axis = (/1/)

yval = var_diff

opti = True
opti@boxWidth = 0.80

plot = boxplot(wks,x_axis,yval,opti,res,False)

;===========
resP = True

;resP@gsnPanelFigureStrings = (/"(b)"/)
resP@amJust = "TopRight"
gsn_panel(wks,plot,(/1,1/),resP)


end
