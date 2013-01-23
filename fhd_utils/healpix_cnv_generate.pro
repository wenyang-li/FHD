FUNCTION healpix_cnv_generate,obs,file_path_fhd,nside=nside,mask=mask,radius=radius,$
    restore_last=restore_last,silent=silent,_Extra=extra

;vis_path_default,data_directory,filename,file_path,obs=obs
IF Keyword_Set(restore_last) AND (file_test(file_path_fhd+'_hpxcnv'+'.sav') EQ 0) THEN BEGIN 
    IF ~Keyword_Set(silent) THEN print,file_path_fhd+'_hpxcnv'+'.sav' +' Not found. Recalculating.' 
    restore_last=0
ENDIF
IF Keyword_Set(restore_last) THEN BEGIN
    IF ~Keyword_Set(silent) THEN print,'Saved Healpix grid map restored'
    restore,file_path_fhd+'_hpxcnv'+'.sav'
    nside=hpx_cnv.nside
    RETURN,hpx_cnv
ENDIF

t00=Systime(1)
astr=(*obs.bin).astr
dimension=obs.dimension
elements=obs.elements
IF N_Elements(radius) EQ 0 THEN radius=obs.degpix*(dimension>elements)/4.
;all angles in DEGREES
;uses RING index scheme
IF not Keyword_Set(nside) THEN BEGIN
    pix_sky=4.*!Pi*!RaDeg^2./Product(astr.cdelt)
    Nside=2.^(Ceil(ALOG(Sqrt(pix_sky/12.))/ALOG(2))) ;=1024. for 0.1119 degrees/pixel
ENDIF
npix=nside2npix(nside)

ang2vec,obs.obsdec,obs.obsra,cen_coords,/astro

Query_disc,nside,cen_coords,radius,hpx_inds0,ninds,/deg

pix2vec_ring,nside,hpx_inds0,pix_coords
vec2ang,pix_coords,pix_dec,pix_ra,/astro
ad2xy,pix_ra,pix_dec,astr,xv_hpx,yv_hpx
pix_coords=0
pix_ra=0
pix_dec=0

;NOTE: slightly more restrictive boundary here ('LT' and 'GT' instead of 'LE' and 'GE') 
pix_i_use=where((xv_hpx GT 0) AND (xv_hpx LT dimension-1) AND (yv_hpx GT 0) AND (yv_hpx LT elements-1),n_hpx_use)
xv_hpx=xv_hpx[pix_i_use]
yv_hpx=yv_hpx[pix_i_use]
IF Keyword_Set(mask) THEN BEGIN
    hpx_mask00=mask[Floor(xv_hpx),Floor(yv_hpx)]
    hpx_mask01=mask[Floor(xv_hpx),Ceil(yv_hpx)]
    hpx_mask10=mask[Ceil(xv_hpx),Floor(yv_hpx)]
    hpx_mask11=mask[Ceil(xv_hpx),Ceil(yv_hpx)]
    hpx_mask=Temporary(hpx_mask00)*Temporary(hpx_mask01)*Temporary(hpx_mask10)*Temporary(hpx_mask11)
    pix_i_use2=where(hpx_mask,n_hpx_use)
    xv_hpx=xv_hpx[pix_i_use2]
    yv_hpx=yv_hpx[pix_i_use2]
    pix_i_use=pix_i_use[pix_i_use2]
ENDIF ;ELSE BEGIN
;    mask=fltarr(dimension,elements)
;    mask[Floor(xv_hpx),Floor(yv_hpx)]=1
;    mask[Floor(xv_hpx),Ceil(yv_hpx)]=1
;    mask[Ceil(xv_hpx),Floor(yv_hpx)]=1
;    mask[Ceil(xv_hpx),Ceil(yv_hpx)]=1
;ENDELSE
hpx_inds=hpx_inds0[pix_i_use]

x_frac=xv_hpx-Floor(xv_hpx)
y_frac=yv_hpx-Floor(yv_hpx)
image_inds=Long64(Floor(xv_hpx)+dimension*Floor(yv_hpx))
corner_inds=Long64([0,1,dimension,dimension+1])

min_bin=Min(Floor(xv_hpx)+dimension*Floor(yv_hpx))>0L
max_bin=Max(Ceil(xv_hpx)+dimension*Ceil(yv_hpx))<(dimension*elements-1L)
h00=histogram(Floor(xv_hpx)+dimension*Floor(yv_hpx),min=min_bin,max=max_bin,/binsize,reverse_ind=ri00)
h01=histogram(Floor(xv_hpx)+dimension*Ceil(yv_hpx),min=min_bin,max=max_bin,/binsize,reverse_ind=ri01)
h10=histogram(Ceil(xv_hpx)+dimension*Floor(yv_hpx),min=min_bin,max=max_bin,/binsize,reverse_ind=ri10)
h11=histogram(Ceil(xv_hpx)+dimension*Ceil(yv_hpx),min=min_bin,max=max_bin,/binsize,reverse_ind=ri11)
htot=h00+h01+h10+h11
inds=where(htot,n_img_use)

n_arr=htot[inds]
;hpx_inds=

i_use=inds+min_bin
sa=Ptrarr(n_img_use,/allocate)
ija=Ptrarr(n_img_use,/allocate)

FOR i=0L,n_img_use-1L DO BEGIN
    ind0=inds[i]
    sa0=fltarr(n_arr[i])
    ija0=Lon64arr(n_arr[i])
    bin_i=Total([0,h00[ind0],h01[ind0],h10[ind0],h11[ind0]],/cumulative)-1
    IF h00[ind0] GT 0 THEN BEGIN
        bi=0
        inds1=ri00[ri00[ind0]:ri00[ind0+1]-1]
        sa0[bin_i[bi]+1:bin_i[bi+1]]=x_frac[inds1]*y_frac[inds1]
        ija0[bin_i[bi]+1:bin_i[bi+1]]=inds1
    ENDIF
    IF h01[ind0] GT 0 THEN BEGIN
        bi=1
        inds1=ri01[ri01[ind0]:ri01[ind0+1]-1]
        sa0[bin_i[bi]+1:bin_i[bi+1]]=x_frac[inds1]*(1.-y_frac[inds1])
        ija0[bin_i[bi]+1:bin_i[bi+1]]=inds1
    ENDIF
    IF h10[ind0] GT 0 THEN BEGIN
        bi=2
        inds1=ri10[ri10[ind0]:ri10[ind0+1]-1]
        sa0[bin_i[bi]+1:bin_i[bi+1]]=(1.-x_frac[inds1])*y_frac[inds1]
        ija0[bin_i[bi]+1:bin_i[bi+1]]=inds1
    ENDIF
    IF h11[ind0] GT 0 THEN BEGIN
        bi=3
        inds1=ri11[ri11[ind0]:ri11[ind0+1]-1]
        sa0[bin_i[bi]+1:bin_i[bi+1]]=(1.-x_frac[inds1])*(1.-y_frac[inds1])
        ija0[bin_i[bi]+1:bin_i[bi+1]]=inds1
    ENDIF
    *sa[i]=sa0
    *ija[i]=ija0
        
ENDFOR

hpx_cnv={nside:nside,ija:ija,sa:sa,i_use:i_use,inds:hpx_inds}

save,hpx_cnv,filename=file_path_fhd+'_hpxcnv'+'.sav'
;print,Systime(1)-t00 ; currently takes 4 - 8 seconds (43MB file). Takes ~ 0.17s to apply
RETURN,hpx_cnv
END