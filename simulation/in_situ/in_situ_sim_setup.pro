PRO in_situ_sim_setup, in_situ_sim_input, vis_arr, vis_weights, flag_calibration, n_pol=n_pol, enhance_eor=enhance_eor, $
    eor_savefile=eor_savefile, file_path_vis=file_path_vis, file_path_fhd=file_path_fhd,sim_noise=sim_noise, $
    hdr=hdr, params=params, calibration_catalog_file_path=calibration_catalog_file_path, $
    diffuse_calibrate=diffuse_calibrate,transfer_calibration=transfer_calibration,freq_start=freq_start, $
    freq_end=freq_end,tile_flag_list=tile_flag_list,deproject_w_term=deproject_w_term,dft_threshold=dft_threshold, $
    remove_sim_flags=remove_sim_flags,_Extra=extra
    
  print, "Performing in-situ simulation"
  
  if ~keyword_set(n_pol) then n_pol=2
  if n_pol EQ 2 then pol_name=['XX','YY'] else pol_name=['XX','YY','XY','YX']
  
  vis_model_arr=PTRARR(n_pol,/allocate)
  obs_id = file_basename(file_path_vis, '.uvfits')
  
  ;restore model visibilities given the in_situ_sim_input to act as the input data visibilities
  if isa(in_situ_sim_input,'String') then begin
    for pol_i=0, n_pol-1 do begin
      vis_model_arr[pol_i] = GETVAR_SAVEFILE(in_situ_sim_input+'/vis_data/'+obs_id+'_vis_model_'+pol_name[pol_i]+'.sav', 'vis_model_ptr')
      print, "Using " + in_situ_sim_input+'/vis_data/'+obs_id+'_vis_model_'+pol_name[pol_i]+'.sav as input model'
    endfor
  endif
  
  ;***Begin in-situ model making to act as input data visibilities if read-in is not available
  IF *vis_model_arr[0] EQ !NULL then begin
    print, "Read-in file not found/provided in in_situ_sim_input. Creating model"
    
    ;Note: explicitly reference dft_threshold here to remove it from EXTRA, which would be passed on to lower-level routines
    obs=fhd_struct_init_obs(file_path_vis,hdr,params,n_pol=n_pol,dft_threshold=dft_threshold,_Extra=extra)
    n_pol=obs.n_pol
    n_freq=obs.n_freq
    
    IF Keyword_Set(deproject_w_term) THEN vis_arr=simple_deproject_w_term(obs,params,vis_arr,direction=deproject_w_term)
    
    ;Read in or construct a new beam model. Also sets up the structure PSF
    psf=beam_setup(obs,status_str,antenna,file_path_fhd=file_path_fhd,restore_last=0,silent=silent,timing=t_beam,no_save=no_save,_Extra=extra)
    jones=fhd_struct_init_jones(obs,status_str,file_path_fhd=file_path_fhd,restore=0,mask=beam_mask,_Extra=extra)
    
    vis_weights=vis_flag_basic(vis_weights,obs,params,n_pol=n_pol,n_freq=n_freq,freq_start=freq_start,$
      freq_end=freq_end,tile_flag_list=tile_flag_list,vis_ptr=vis_arr,unflag_all=1,no_frequency_flagging=1,_Extra=extra)
    vis_weights_update,vis_weights,obs,psf,params,_Extra=extra
    
    IF Keyword_Set(calibration_catalog_file_path) THEN catalog_use=calibration_catalog_file_path
    IF ~Keyword_Set(calibration_source_list) THEN $
      calibration_source_list=generate_source_cal_list(obs,psf,catalog_path=catalog_use,_Extra=extra)
    skymodel_cal=fhd_struct_init_skymodel(obs,source_list=calibration_source_list,catalog_path=catalog_use,return_cal=1,diffuse_model=diffuse_calibrate,_Extra=extra)
    cal=fhd_struct_init_cal(obs,params,skymodel_cal,source_list=calibration_source_list,$
      catalog_path=catalog_use,transfer_calibration=transfer_calibration,_Extra=extra)
      
    vis_model_arr=vis_source_model(cal.skymodel,obs,status_str,psf,params,vis_weights,cal,jones,model_uv_arr=model_uv_arr,fill_model_vis=1,$
      timing=model_timing,silent=silent,error=error,/calibration_flag,spectral_model_uv_arr=spectral_model_uv_arr,_Extra=extra)
    for pol_i=0, n_pol-1 do begin
      vis_weights_use=0>*vis_weights[pol_i]<1
      *vis_model_arr[pol_i]=*vis_model_arr[pol_i]*vis_weights_use
    endfor
    
    print, "Saving input model visibilities to " + file_dirname(file_path_fhd) +'/sim_outputs/'+obs_id+'_input_model.sav'
    file_mkdir, file_dirname(file_path_fhd) +'/sim_outputs'
    save, vis_model_arr, filename=file_dirname(file_path_fhd) +'/sim_outputs/'+obs_id+'_input_model.sav'
    
    fhd_log_settings,file_path_fhd,obs=obs,psf=psf,cal=cal,antenna=antenna,skymodel=skymodel,cmd_args=cmd_args,/overwrite,sub_dir='sim_outputs'
    
    undefine, psf, jones, skymodel_cal, cal, calibration_source_list
    
  endif
  ;***End in-situ model making to act as input data visibilities if read-in is not available
  
  vis_arr=temporary(vis_model_arr)
  
  ;Remove all weighting to remove pfb effects and flagged channels. Remove calibration flagging.
  if keyword_set(remove_sim_flags) then for pol_i=0, n_pol-1 do (*vis_weights[pol_i])[*,*]=1.
  flag_calibration=0
  
  ;restore EoR visibilities
  If keyword_set(eor_savefile) then begin
    vis_eor=PTRARR(n_pol,/allocate)
    
    ;*Search for the specified eor savefile
    if (total(file_test(eor_savefile)) GT 0) AND (file_test(eor_savefile,/directory) EQ 0) then begin
      size_savefile=(size(eor_savefile))[1]
      
      void = GETVAR_SAVEFILE(eor_savefile[0], names=names)
      vis_varname = names[where(strmatch(names, '*vis*') EQ 1,n_count)]  ;assumption: visibilities in sav file have "vis" in the name
      
      ;Restore visibilities that are in different polarization save files, or restore the all-pol save file
      if size_savefile EQ n_pol then $
        for pol_i=0,n_pol-1 do vis_eor[pol_i] = GETVAR_SAVEFILE(eor_savefile[pol_i], vis_varname) $
      else vis_eor = GETVAR_SAVEFILE(eor_savefile, vis_varname)
      
    endif else begin
      for pol_i=0,n_pol-1 do begin
        if total(file_test(eor_savefile + '/' + obs_id + '_vis_' + pol_name[pol_i] + '.sav')) GT 0 then begin
          void = GETVAR_SAVEFILE(eor_savefile + '/' + obs_id + '_vis_' + pol_name[pol_i] + '.sav', names=names)
          vis_varname = names[where(strmatch(names, '*vis*') EQ 1,n_count)]  ;assumption: visibilities in sav file have "vis" in the name
          
          ;Restore visibilities that are in the <obsid>_vis_XX/vis_YY format
          vis_eor[pol_i] = GETVAR_SAVEFILE(eor_savefile + '/' + obs_id + '_vis_' + pol_name[pol_i] + '.sav', vis_varname)
        endif else message, "eor_savefile not found! Tried " + eor_savefile + " and " + eor_savefile + '/' + obs_id + "_vis_" + pol_name[pol_i] + ".sav (for all pol)"
      endfor
    endelse
    ;*End of search for the specified eor savefile
    
    ;Boost the eor signal by a specified amount
    If keyword_set(enhance_eor) then begin
      print, "Enhancing input EoR by "+enhance_eor+"x"
      for pol_i=0,n_pol-1 do *vis_eor[pol_i]=*vis_eor[pol_i]*enhance_eor
    endif
    
    ;Combine the calibrated visibilities in the correct format for the script
    for pol_i=0,n_pol-1 do *vis_arr[pol_i] = *vis_arr[pol_i]+*vis_eor[pol_i]
    
  endif
  
  ;Optionally add noise to the visibilities (from a save file)
  If keyword_set(sim_noise) then begin
  
    if total(file_test(sim_noise)) GT 0 then begin
      ;Restore or create the noise visibilities
      void = getvar_savefile(sim_noise,names=names)
      if strmatch(names, 'obs', /FOLD_CASE) then begin
        obs_restore = getvar_savefile(sim_noise,strjoin((names[where(strmatch(names, 'obs', /FOLD_CASE) EQ 1)])[0]))
        vis_noise = vis_noise_simulation(cal_sim_input, vis_arr, obs_id, obs_restore, n_pol=n_pol, file_path_fhd=file_path_fhd)
      endif else vis_noise = getvar_savefile(sim_noise,strjoin(names[0]))
      
    endif else message, 'Need to specify either a filepath to an obs with vis_noise calculation or a noise simulation.'
    
    ;Add the noise to the visibilities, but keeping zeroed visibilities fully zero
    for pol_i=0, n_pol-1 do begin
      nonzero_i = where(real_part(*vis_arr[pol_i]) NE 0)
      (*vis_arr[pol_i])[nonzero_i] = (*vis_arr[pol_i])[nonzero_i]+(*vis_noise[pol_i])[nonzero_i]
    endfor
    
  endif
  
  return
END
