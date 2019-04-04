;Date: 08/30/2016
;Function: Check one point, at one soil level to see its H2OSOI values for 16 ensemble from and 16 ensemble from reading.
;Note: please_find_me_figure_s1  ; (in the revised version, this is figure s1)

begin


wks = gsn_open_wks("pdf","SOILLIQ_ens") 
plot = new(20, graphic)

maxarray = (/10.0,20.0,30.0,60.0,100.0,200.0,300.0,400.0,500.0,700.0/)

do isoil = 0,9
;==============  Writing ensemble ============
total=16
year = "1994"
Label = "write"

tag1 = "NoSP"
tag2 = "nocrm"

path1 = "/gdata/pritchard2/hongcheq/Data/GLACE/SPCAM3.5_"+tag1+year
path2 = (/"/spcam3.5_"+tag2+"_SM_"+Label,"/spcam3.5_"+tag2+"_SM_Ensemble_"+Label/)
path3 = (/"","/ens_1","/ens_2","/ens_3","/ens_4","/ens_5",\
              "/ens_6","/ens_7","/ens_8","/ens_9","/ens_10",\
              "/ens_11","/ens_12","/ens_13","/ens_14","/ens_15"/)
path4 = "/spcam3.5sld_"+tag1+".cam2.h1."+year+"-06-01-00000.nc"
path4b = "/spcam3.5sld_"+tag1+".clm2.h1."+year+"-06-01-00000.nc"

;point_lon = 240.0
;point_lat = 40.0
;point_lon = 15.0
;point_lat = 50.0
point_lon = 260.0
point_lat = 50.0
;point_lon = 5.0
;point_lat = 50.0
;point_lon = 280.0
;point_lat = 40.0

Ens_write = new((/16,100,10/),double)
printVarSummary(Ens_write)

do ifile = 0,total-1
if (ifile .eq. 0) then
f = addfile(path1+path2(0)+path3(ifile)+path4,"r")
fl = addfile(path1+path2(0)+path3(ifile)+path4b,"r")
else
f = addfile(path1+path2(1)+path3(ifile)+path4,"r")
fl = addfile(path1+path2(1)+path3(ifile)+path4b,"r")
end if

SOILLIQ = fl->SOILLIQ(:,:,:,:)    ; to add lon and lat dimension to it
SOILLIQ!2 = "lat"
SOILLIQ&lat = f->lat
SOILLIQ!3 = "lon"
SOILLIQ&lon = f->lon
printVarSummary(SOILLIQ)

Ens_write(ifile,:,:) = SOILLIQ(:,:,{point_lat},{point_lon})

printVarSummary(Ens_write)

Ens_write@long_name = "soil liquid water at lon="+point_lon+",lat="+point_lat

end do ; ifile = 0,total-1

res = True
res@gsnDraw = False
res@gsnFrame = False
res@tiMainString = "ENS(W), levsoi = "+isoil
res@trYMinF = 0.0
res@trYMaxF = maxarray(isoil)

plot(0+isoil*2) = gsn_csm_xy(wks,Ens_write&time,Ens_write(:,:,isoil),res)

;====================== Reading Ensemble ================

Label = "read"

path1 = "/gdata/pritchard2/hongcheq/Data/GLACE/SPCAM3.5_"+tag1+year
path2 = (/"/spcam3.5_"+tag2+"_SM_"+Label,"/spcam3.5_"+tag2+"_SM_Ensemble_"+Label/)
path3 = (/"","/ens_1","/ens_2","/ens_3","/ens_4","/ens_5",\
              "/ens_6","/ens_7","/ens_8","/ens_9","/ens_10",\
              "/ens_11","/ens_12","/ens_13","/ens_14","/ens_15"/)
path4 = "/spcam3.5sld_"+tag1+".cam2.h1."+year+"-06-01-00000.nc"
path4b = "/spcam3.5sld_"+tag1+".clm2.h1."+year+"-06-01-00000.nc"

Ens_read = new((/16,100,10/),double)

do ifile = 0,total-1
if (ifile .eq. 0) then
f = addfile(path1+path2(0)+path3(ifile)+path4,"r")
fl = addfile(path1+path2(0)+path3(ifile)+path4b,"r")
else
f = addfile(path1+path2(1)+path3(ifile)+path4,"r")
fl = addfile(path1+path2(1)+path3(ifile)+path4b,"r")
end if

SOILLIQ2 = fl->SOILLIQ(:,:,:,:)    ; to add lon and lat dimension to it
SOILLIQ2!2 = "lat"
SOILLIQ2&lat = f->lat
SOILLIQ2!3 = "lon"
SOILLIQ2&lon = f->lon
printVarSummary(SOILLIQ2)

Ens_read(ifile,:,:) = SOILLIQ2(:,:,{point_lat},{point_lon})

printVarSummary(Ens_read)

Ens_read@long_name =  "soil liquid water at lon="+point_lon+",lat="+point_lat

end do ; ifile = 0,total-1

res@tiMainString = "ENS(S), levsoi = "+isoil

plot(1+isoil*2) = gsn_csm_xy(wks,Ens_read&time,Ens_read(:,:,isoil),res)

end do  ; do isoil = 0,9

;=============== Panel ===============
resP = True

resP@txString = "Sanity check @ (lon="+point_lon+", lat="+point_lat+")"
gsn_panel(wks,plot,(/5,4/),resP)


end 
