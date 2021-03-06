;Function:compute terrestrial segment of land-atmosphere coupling with drizzling removal
;         
;Date:20170622
;Note: please_find_me_figure10

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

drizzle_threshold = (/1.0e-1,1.0e-2,1.0e-3,1.0e-4,1.0e-5/)

sub_label = (/"(a)","(c)","(e)","(g)","(i)"/)

do idt = 0,dimsizes(drizzle_threshold) -1
;do idt = 0,0                  ; 0.1mm/day
do ilevsoi = 3,3
;do ilevsoi = 0,9

levsoi = ilevsoi    ; soil level

tag = "CAM3.5"
tag2 = "NoSP"

path1 = "/lustre/DATA/pritchard/hongcheq/AMIP_terrestrial_coupling_strength/"
file_LH_JJA = "LHFLX.spcam3.5sld_"+tag2+".cam2.h1.1996-2015_JJA_daily.nc"
file_SM_JJA = "SOILLIQ.spcam3.5sld_"+tag2+".clm2.h1.1996-2015_JJA_daily.nc"
file_P_JJA = "PRECT.spcam3.5sld_"+tag2+".cam2.h1.1996-2015_JJA_daily.nc"

;;;file_LH_DJF = "LHFLX.cam3.5.cam2.h0.0000_0019_20years_DJF_daily.nc"
;;;file_SM_DJF = "SOILLIQ.cam3.5.clm2.h0.0000_0019_20years_DJF_daily.nc"

f_LH_JJA = addfile(path1+file_LH_JJA,"r")
f_SM_JJA = addfile(path1+file_SM_JJA,"r")
f_P_JJA = addfile(path1+file_P_JJA,"r")

LHFLX_JJA = f_LH_JJA->LHFLX
SOILLIQ_JJA = f_SM_JJA->SOILLIQ(:,levsoi,:,:)   ; levsoi = 0.1188651 m depth
PRECT_JJA = f_P_JJA->PRECT   ; its unit is meter/sec

PRECT_JJA = PRECT_JJA * 8.64*10^7 ;
PRECT_JJA@units = "mm/day"


LHFLX_JJA@_FillValue = 1.0e36
SOILLIQ_JJA@_FillValue = 1.0e36

printVarSummary(PRECT_JJA)
printVarSummary(LHFLX_JJA)
printVarSummary(SOILLIQ_JJA)


do itime = 0,1839 
   do ilat = 0,63
      do ilon = 0,127
           if (PRECT_JJA(itime,ilat,ilon) .lt. drizzle_threshold(idt) ) then         
              LHFLX_JJA(itime,ilat,ilon) = LHFLX_JJA@_FillValue
              SOILLIQ_JJA(itime,ilat,ilon) = SOILLIQ_JJA@_FillValue
           end if
      end do
   end do
end do

;;;f_LH_DJF = addfile(path1+file_LH_DJF,"r")
;;;f_SM_DJF = addfile(path1+file_SM_DJF,"r")
;;;LHFLX_DJF = f_LH_DJF->LHFLX
;;;SOILLIQ_DJF = f_SM_DJF->SOILLIQ(:,3,:,:) 

;============= calculate linear regression, standard deviation of SOILLIQ, correlation of SOILLIQ and LHFLX ========
rc_JJA = new((/dimsizes(LHFLX_JJA&lat),dimsizes(LHFLX_JJA&lon)/),double)
;;;rc_DJF = new((/dimsizes(LHFLX_DJF&lat),dimsizes(LHFLX_DJF&lon)/),double)

std_SM_JJA = new((/dimsizes(LHFLX_JJA&lat),dimsizes(LHFLX_JJA&lon)/),double)
;;;std_SM_DJF = new((/dimsizes(LHFLX_JJA&lat),dimsizes(LHFLX_JJA&lon)/),double)

corr_JJA = new((/dimsizes(LHFLX_JJA&lat),dimsizes(LHFLX_JJA&lon)/),double)
;;;corr_DJF = new((/dimsizes(LHFLX_JJA&lat),dimsizes(LHFLX_JJA&lon)/),double)

tc_JJA = new((/dimsizes(LHFLX_JJA&lat),dimsizes(LHFLX_JJA&lon)/),double)
;;;tc_DJF = new((/dimsizes(LHFLX_JJA&lat),dimsizes(LHFLX_JJA&lon)/),double)

do ilat = 0,dimsizes(LHFLX_JJA&lat) - 1
	do ilon = 0, dimsizes(LHFLX_JJA&lon) - 1
		rc_JJA(ilat,ilon)  = regline(SOILLIQ_JJA(:,ilat,ilon),LHFLX_JJA(:,ilat,ilon))
;;;             rc_DJF(ilat,ilon)  = regline(SOILLIQ_DJF(:,ilat,ilon),LHFLX_DJF(:,ilat,ilon))
		std_SM_JJA(ilat,ilon) = dim_stddev_n_Wrap(SOILLIQ_JJA(:,ilat,ilon),0) 
;;;		std_SM_DJF(ilat,ilon) = dim_stddev_n_Wrap(SOILLIQ_DJF(:,ilat,ilon),0)   
		corr_JJA(ilat,ilon) = escorc(LHFLX_JJA(:,ilat,ilon), SOILLIQ_JJA(:,ilat,ilon))
;;;		corr_DJF(ilat,ilon) = escorc(LHFLX_DJF(:,ilat,ilon), SOILLIQ_DJF(:,ilat,ilon))

		pr_JJA = rtest(corr_JJA(ilat,ilon),dimsizes(LHFLX_JJA&time),0)     ; 99% confidence level
;;;		pr_DJF = rtest(corr_DJF(ilat,ilon),dimsizes(LHFLX_DJF&time),0)
		
 		;; significance level set to 99%
		if (pr_JJA .lt. 0.01 ) then        
		tc_JJA(ilat,ilon) = std_SM_JJA(ilat,ilon) * rc_JJA(ilat,ilon) 
		end if
;;;		if (pr_DJF .lt. 0.01) then
;;;		tc_DJF(ilat,ilon) = std_SM_DJF(ilat,ilon) * rc_DJF(ilat,ilon)	
;;;		end if
	end do
end do

tc_JJA@units = LHFLX_JJA@units
tc_JJA!0 = "lat"
tc_JJA&lat = LHFLX_JJA&lat
tc_JJA!1 = "lon"
tc_JJA&lon = LHFLX_JJA&lon

;;;tc_DJF@units = LHFLX_DJF@units
;;;tc_DJF!0 = "lat"
;;;tc_DJF&lat = LHFLX_DJF&lat
;;;tc_DJF!1 = "lon"
;;;tc_DJF&lon = LHFLX_DJF&lon

;=============================
system("rm -f "+tag+"_levsoi"+levsoi+"_terretrial_coupling_strength.nc")
ncdf = addfile(tag+"_levsoi"+levsoi+"_terretrial_coupling_strength.nc","c")
filedimdef(ncdf,"time",-1,True)
ncdf->tc_JJA = tc_JJA
;;;ncdf->tc_DJF = tc_DJF

;========= Plot ===============

wks = gsn_open_wks("png",idt+"_"+tag+"_levsoi"+levsoi+"terrestrial_coupling_strength.png")
;;;plot = new(2,graphic)

res = True
res@gsnDraw = False
res@gsnFrame = False
res@cnInfoLabelOn = False
res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@lbLabelBarOn = False

res@mpMaxLatF = 60.0
res@mpMinLatF = -30.0

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -15.0
res@cnMaxLevelValF = 25.0
res@cnLevelSpacingF = 5.0
res@cnFillColors = (/ "slateblue4","slateblue1","steelblue3","steelblue1","darkseagreen2","yellow","gold","tan3","sienna2","firebrick" /) 

;res@tiMainString = "JJA"
res@gsnLeftString = sub_label(idt)+" "+tag
;plot(0) = gsn_csm_contour_map(wks,tc_JJA,res)
plot = gsn_csm_contour_map(wks,tc_JJA,res)
;;;res@tiMainString = "DJF"
;;;plot(1) = gsn_csm_contour_map(wks,tc_DJF,res)
;=================
resP = True
;resP@txString = "Terrestrial Coupling (drizzle_removal_threshold="+drizzle_threshold(idt)+"), AMIP 1996-2015, "+tag
resP@txString = "Terrestrial Coupling (with drizzle removal), AMIP 1996-2015 JJA"

;resP@gsnPanelLabelBar = True

resP@amJust = "TopLeft"
;resP@gsnPanelFigureStrings = (/"(a)"/)
; resP@gsnPanelFigureStrings = sub_label(idt)

gsn_panel(wks,plot,(/1,1/),resP)
;gsn_panel(wks,plot,(/2,1/),resP)

end do ; do ilevsoi
end do ; do idt
end
