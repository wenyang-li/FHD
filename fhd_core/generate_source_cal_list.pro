FUNCTION generate_source_cal_list,obs,psf,catalog_path=catalog_path,calibration_spectral_index=calibration_spectral_index,$
    max_calibration_sources=max_calibration_sources,calibration_flux_threshold=calibration_flux_threshold,$
    no_restrict_cal_sources=no_restrict_cal_sources,no_extend=no_extend,mask=mask,$
    allow_sidelobe_cal_sources=allow_sidelobe_cal_sources,_Extra=extra
;catalog=getvar_savefile(catalog_path,'catalog')
RESTORE,catalog_path,/relaxed ;catalog
IF N_Elements(calibration_flux_threshold) EQ 0 THEN calibration_flux_threshold=0.
IF Keyword_Set(allow_sidelobe_cal_sources) THEN IF N_Elements(no_restrict_cal_sources) EQ 0 THEN no_restrict_cal_sources=1
IF N_Elements(no_restrict_cal_sources) EQ 0 THEN no_restrict_cal_sources=1
astr=obs.astr
dimension=obs.dimension
elements=obs.elements
degpix=obs.degpix

FoV=!RaDeg/obs.kpix
freq_arr=psf.freq
freq_use=Mean(freq_arr)/1E6
n_pol=obs.n_pol

ra0=astr.crval[0]
dec0=astr.crval[1]
angs=angle_difference(dec0,ra0,catalog.dec,catalog.ra,/degree)
i_use=where(Abs(angs) LE FoV/2.,n_use)

IF Keyword_Set(no_restrict_cal_sources) THEN BEGIN
    fft_alias_range=0.
    IF Keyword_Set(allow_sidelobe_cal_sources) THEN cal_beam_threshold=0.01 ELSE cal_beam_threshold=0.05
ENDIF ELSE BEGIN
    fft_alias_range=dimension/4.
    cal_beam_threshold=0.2
ENDELSE

IF n_use GT 0 THEN BEGIN
    catalog=catalog[i_use]
    source_list=catalog
;    source_list.ra=catalog.ra
;    source_list.dec=catalog.dec
    ad2xy,source_list.ra,source_list.dec,astr,x_arr,y_arr
    source_list.x=x_arr
    source_list.y=y_arr
    IF N_Elements(calibration_spectral_index) EQ 0 THEN calibration_spectral_index=catalog.alpha
    FOR i=0,7 DO source_list.flux.(i)=catalog.flux.(i)*(freq_use/catalog.freq)^calibration_spectral_index
    source_list.alpha=calibration_spectral_index
;    source_list.StoN=catalog.StoN
    
    beam=fltarr(dimension,elements)
    beam_arr=Ptrarr(n_pol<2)
    FOR pol_i=0,(n_pol<2)-1 DO BEGIN
        beam_arr[pol_i]=Ptr_new(beam_image(psf,obs,pol_i=pol_i,/fast)>0.)
        beam+=*beam_arr[pol_i]^2.
    ENDFOR
    beam=Sqrt(beam/(n_pol<2))
    IF Keyword_Set(allow_sidelobe_cal_sources) THEN beam_i=where(beam GT cal_beam_threshold) $
        ELSE beam_i=region_grow(beam,dimension/2.+dimension*elements/2.,threshold=[Max(beam)/2.<cal_beam_threshold,Max(beam)>1.])
    beam_mask=fltarr(dimension,elements) & beam_mask[beam_i]=1.
    IF N_Elements(mask) EQ N_Elements(beam_mask) THEN beam_mask*=mask

    src_use=where((x_arr GE fft_alias_range) AND (x_arr LE dimension-1-fft_alias_range) AND (y_arr GE fft_alias_range) $
        AND (y_arr LE elements-1-fft_alias_range) AND (source_list.flux.I GT calibration_flux_threshold),n_src_use)
    
    IF n_src_use EQ 0 THEN RETURN,source_comp_init(n_sources=0,freq=obs.freq_center);
    src_use2=where(beam_mask[Round(x_arr[src_use]),Round(y_arr[src_use])],n_src_use)
    IF n_src_use GT 0 THEN src_use=src_use[src_use2]
    source_list=source_list[src_use]
    beam_list=Ptrarr(n_pol<2)
;    FOR pol_i=0,(n_pol<2)-1 DO BEGIN
;        beam_list[pol_i]=Ptr_new((*beam_arr[pol_i])[source_list.x,source_list.y])
;        source_list.flux.(pol_i)=source_list.flux.I*(*beam_list[pol_i])
;    ENDFOR
    
    influence=source_list.flux.I*beam[source_list.x,source_list.y]
    
    order=Reverse(sort(influence))
    source_list=source_list[order]
    source_list.id=Lindgen(n_src_use)
    IF Keyword_Set(no_extend) THEN source_list.extend=Ptrarr(n_src_use) ELSE BEGIN
        extend_i=where(Ptr_valid(source_list.extend),n_extend)
        FOR ext_i=0L,n_extend-1 DO BEGIN
            extend_list=*source_list[extend_i[ext_i]].extend
            ad2xy,extend_list.ra,extend_list.dec,astr,x_arr,y_arr
            extend_list.x=x_arr
            extend_list.y=y_arr
;            FOR pol_i=0,(n_pol<2)-1 DO extend_list.flux.(pol_i)=extend_list.flux.I*(*beam_list[pol_i])[extend_i[ext_i]]
            *source_list[extend_i[ext_i]].extend=extend_list
        ENDFOR
    ENDELSE
ENDIF ELSE RETURN,source_comp_init(n_sources=0,freq=obs.freq_center)

IF Keyword_Set(max_calibration_sources) THEN IF N_Elements(source_list) GT max_calibration_sources $
    THEN source_list=source_list[0:max_calibration_sources-1]
RETURN,source_list
END