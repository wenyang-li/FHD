FUNCTION beam_power,antenna1,antenna2,ant_pol1=ant_pol1,ant_pol2=ant_pol2,freq_i=freq_i,psf_dim=psf_dim,$
    psf_image_dim=psf_image_dim,psf_intermediate_res=psf_intermediate_res,psf_resolution=psf_resolution,$
    beam_mask_electric_field=beam_mask_electric_field,beam_mask_threshold=beam_mask_threshold,$
    xvals_uv_superres=xvals_uv_superres,yvals_uv_superres=yvals_uv_superres,zen_int_x=zen_int_x,zen_int_y=zen_int_y,$
    interpolate_beam_threshold=interpolate_beam_threshold,debug_beam_clip_grow=debug_beam_clip_grow,$
    debug_beam_conjugate=debug_beam_conjugate, debug_beam_clip_floor=debug_beam_clip_floor,$
    debug_clip_beam_mask=debug_clip_beam_mask,debug_gauss=debug_gauss, _Extra=extra
    
icomp = Complex(0, 1)
freq_center=antenna1.freq[freq_i]
dimension_super=(size(xvals_uv_superres,/dimension))[0]

Jones1=antenna1.Jones[*,*,freq_i]
Jones2=antenna2.Jones[*,*,freq_i]

beam_ant1=DComplex(*(antenna1.response[ant_pol1,freq_i]))
beam_ant2=DComplex(Conj(*(antenna2.response[ant_pol2,freq_i])))
beam_norm=1.

IF ant_pol1 NE ant_pol2 THEN BEGIN
    power_zenith_beam1=abs((*Jones1[0,ant_pol1])*(Conj(*Jones2[0,ant_pol1]))+(*Jones1[1,ant_pol1])*(Conj(*Jones2[1,ant_pol1])))
    power_zenith1=Interpolate(power_zenith_beam1,zen_int_x,zen_int_y,cubic=-0.5)
    power_zenith_beam2=abs((*Jones1[0,ant_pol2])*(Conj(*Jones2[0,ant_pol2]))+(*Jones1[1,ant_pol2])*(Conj(*Jones2[1,ant_pol2])))
    power_zenith2=Interpolate(power_zenith_beam2,zen_int_x,zen_int_y,cubic=-0.5)
    power_zenith=Sqrt(power_zenith1*power_zenith2)
    
    power_beam1=(*Jones1[0,ant_pol1]*beam_ant1)*(Conj(*Jones2[0,ant_pol1])*beam_ant2)+$
           (*Jones1[1,ant_pol1]*beam_ant1)*(Conj(*Jones2[1,ant_pol1])*beam_ant2) 
    power_beam2=(*Jones1[0,ant_pol2]*beam_ant1)*(Conj(*Jones2[0,ant_pol2])*beam_ant2)+$
           (*Jones1[1,ant_pol2]*beam_ant1)*(Conj(*Jones2[1,ant_pol2])*beam_ant2) 
    power_beam=Sqrt(power_beam1*power_beam2)
    debug_point=1
ENDIF ELSE BEGIN
    power_zenith_beam=abs((*Jones1[0,ant_pol1])*(Conj(*Jones2[0,ant_pol2]))+(*Jones1[1,ant_pol1])*(Conj(*Jones2[1,ant_pol2])))
    power_zenith=Interpolate(power_zenith_beam,zen_int_x,zen_int_y,cubic=-0.5)
    ;power_zenith=abs((*Jones1[0,ant_pol1])[zen_int_x,zen_int_y]*(Conj(*Jones2[0,ant_pol2]))[zen_int_x,zen_int_y]+$
    ;             (*Jones1[1,ant_pol1])[zen_int_x,zen_int_y]*(Conj(*Jones2[1,ant_pol2]))[zen_int_x,zen_int_y])
    power_beam=(*Jones1[0,ant_pol1]*beam_ant1)*(Conj(*Jones2[0,ant_pol2])*beam_ant2)+$
               (*Jones1[1,ant_pol1]*beam_ant1)*(Conj(*Jones2[1,ant_pol2])*beam_ant2)           
ENDELSE
;power_beam_zenith=Abs(power_beam[zen_int_x,zen_int_y])
;power_zenith_norm=power_zenith/power_beam_zenith/Max(Abs(power_beam))
psf_base_single=dirty_image_generate(power_beam/power_zenith,/no_real)

psf_base_superres=Interpolate(psf_base_single,xvals_uv_superres,yvals_uv_superres,cubic=-0.5)

uv_mask_superres=Fltarr(dimension_super,dimension_super)
psf_mask_threshold_use = Max(Abs(psf_base_superres))/beam_mask_threshold
IF ant_pol1 NE ant_pol2 THEN BEGIN
    seed_i=where(Abs(psf_base_superres) GE Max(Abs(psf_base_superres))/2.,n_seed)
    beam_i=region_grow(Abs(psf_base_superres),seed_i,$
        thresh=[psf_mask_threshold_use,Max(Abs(psf_base_superres))])
ENDIF ELSE BEGIN
    beam_i=region_grow(Abs(psf_base_superres),dimension_super*(1.+dimension_super)/2.,$
        thresh=[psf_mask_threshold_use,Max(Abs(psf_base_superres))])
ENDELSE
uv_mask_superres[beam_i]=1

IF Keyword_Set(debug_beam_clip_grow) THEN BEGIN
    mask_dist_test = morph_distance(uv_mask_superres, neighbor=3, /background)
    uv_mask_superres_int=Fltarr(dimension_super,dimension_super)
    uv_mask_superres_int[where(mask_dist_test LT psf_resolution)] = 1
    inds_ring = where(uv_mask_superres_int - uv_mask_superres, n_ring)
    
    psf_base_real = real_part(psf_base_superres) > 0
    psf_base_imag = imaginary(psf_base_superres)
    IF n_ring GT 0 THEN BEGIN
        inds_cut = where(psf_base_real LE 0, n_cut)
        IF n_cut GT 0 THEN psf_base_imag[inds_cut] = 0
    ENDIF
    psf_base_superres_int = psf_base_real + icomp*psf_base_imag
    
    psf_base_superres = psf_base_superres_int
    uv_mask_superres = uv_mask_superres_int
ENDIF

psf_base_superres*=psf_intermediate_res^2. ;FFT normalization correction in case this changes the total number of pixels
psf_base_superres/=beam_norm
psf_val_ref=Total(psf_base_superres)

IF Keyword_Set(debug_clip_beam_mask) THEN BEGIN
    xvals_i=Reform(meshgrid(psf_dim,psf_dim,1)*psf_resolution,psf_dim^2.)
    yvals_i=Reform(meshgrid(psf_dim,psf_dim,2)*psf_resolution,psf_dim^2.)
    beam_mask = Fltarr(size(psf_base_superres[xvals_i,yvals_i], /dimension)) + 1
    FOR i=0,psf_resolution-1 DO FOR j=0,psf_resolution-1 DO beam_mask *= $
        uv_mask_superres[xvals_i+i,yvals_i+j]
    FOR i=0,psf_resolution-1 DO FOR j=0,psf_resolution-1 DO $
        uv_mask_superres[xvals_i+i,yvals_i+j] *= beam_mask
ENDIF

IF Keyword_Set(interpolate_beam_threshold) THEN BEGIN
    psf_amp = interpol_2d(abs(psf_base_superres*uv_mask_superres), uv_mask_superres) > 0
    psf_phase = Fltarr(size(psf_base_superres, /dimension))
    psf_phase[beam_i] = Atan(psf_base_superres[beam_i], /phase)
    psf_phase = interpol_2d(psf_phase, uv_mask_superres)
    psf_base_superres = psf_amp*Cos(psf_phase) + icomp*psf_amp*Sin(psf_phase)
ENDIF ELSE psf_base_superres*=uv_mask_superres

IF Keyword_Set(debug_beam_clip_floor) THEN BEGIN
    i_use = where(abs(psf_base_superres))
    psf_amp = abs(psf_base_superres)
    psf_phase = Atan(psf_base_superres, /phase)
    
    psf_floor = psf_mask_threshold_use*(psf_intermediate_res^2.)/beam_norm
    psf_amp[i_use] -= psf_floor
    psf_base_superres = psf_amp*Cos(psf_phase) + icomp*psf_amp*Sin(psf_phase)
ENDIF

if keyword_set(debug_gauss) then begin
  max_ind = where(abs(psf_base_superres) EQ max(abs(psf_base_superres)))
  max_ind_col = max_ind mod N_elements(psf_base_superres[0,*])
  max_val = abs(psf_base_superres[max_ind])
  uvbeam_input_1D = FLTARR(N_elements(psf_base_superres[0,*]))
  max_val_arr = FLTARR(N_elements(psf_base_superres[0,*]),N_elements(psf_base_superres[0,*]))
  max_val_arr[*,*] = max_val

  gauss = exp(-(INDGEN(1400,/float)-699.)^2./(2*200^2.))
  uvbeam_input = gauss # transpose(gauss) * max_val_arr
  singauss = FLTARR(1400)
  singauss[100:1299]=sin(2.*!Pi*(1/1200.)*INDGEN(1200,/float))*.0009
  singauss = singauss*gauss
  uvbeam_imag = (Complex(0,1)*(singauss)) # transpose(gauss)
  uvbeam_input = uvbeam_input + uvbeam_imag
  
  sub_level = min(abs(uvbeam_input[*,699]))
  uvbeam_input2 = uvbeam_input - sub_level
  zeroed = where(uvbeam_input2 LT 0,n_count)
  if n_count GT 0 then uvbeam_input2[zeroed] = 0.
  uvbeam_input=uvbeam_input2
  
  psf_base_superres = uvbeam_input
endif else $
psf_base_superres*=psf_val_ref/Total(psf_base_superres)
IF Keyword_Set(debug_beam_conjugate) THEN psf_base_superres=Conj(psf_base_superres)
RETURN,psf_base_superres
END