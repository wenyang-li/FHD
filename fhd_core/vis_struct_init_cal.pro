FUNCTION vis_struct_init_cal,obs,params,gain_arr_ptr=gain_arr_ptr,n_pol=n_pol,n_freq=n_freq,n_tile=n_tile,n_time=n_time,source_list=source_list,$
    tile_A=tile_A,tile_B=tile_B,freq=freq,bin_offset=bin_offset,tile_names=tile_names,u_loc=u_loc,v_loc=v_loc,n_cal_src=n_cal_src,$
    galaxy_cal=galaxy_cal,min_cal_baseline=min_cal_baseline,max_cal_baseline=max_cal_baseline,n_vis_cal=n_vis_cal,$
    cal_time_average=cal_time_average,ref_antenna=ref_antenna,cal_convergence_threshold=cal_convergence_threshold,$
    calibration_origin=calibration_origin,_Extra=extra
IF N_Elements(tile_A) EQ 0 THEN tile_A=(*obs.baseline_info).tile_A
IF N_Elements(tile_B) EQ 0 THEN tile_B=(*obs.baseline_info).tile_B
IF N_Elements(freq) EQ 0 THEN freq=(*obs.baseline_info).freq
IF N_Elements(bin_offset) EQ 0 THEN bin_offset=(*obs.baseline_info).bin_offset
IF N_Elements(tile_names) EQ 0 THEN tile_names=(*obs.baseline_info).tile_names
IF N_Elements(u_loc) EQ 0 THEN u_loc=params.uu
IF N_Elements(v_loc) EQ 0 THEN v_loc=params.vv

IF N_Elements(n_vis_cal) EQ 0 THEN n_vis_cal=obs.n_vis
IF N_Elements(n_pol) EQ 0 THEN n_pol=obs.n_pol
IF N_Elements(n_freq) EQ 0 THEN n_freq=obs.n_freq
IF N_Elements(n_tile) EQ 0 THEN n_tile=obs.n_tile
IF N_Elements(n_time) EQ 0 THEN n_time=N_Elements(bin_offset)
IF N_Elements(n_cal_src) EQ 0 THEN n_cal_src=-1
IF N_Elements(source_list) EQ 0 THEN source_comp_init,source_list,n_sources=n_cal_src,freq=obs.freq_center
IF N_Elements(galaxy_cal) EQ 0 THEN galaxy_cal=0
IF N_Elements(min_cal_baseline) EQ 0 THEN min_cal_baseline=obs.min_baseline
IF N_Elements(max_cal_baseline) EQ 0 THEN max_cal_baseline=obs.max_baseline
IF N_Elements(cal_time_average) EQ 0 THEN cal_time_average=1 ;time average visibilities before calculating calibration solutions by default
IF N_Elements(min_cal_solutions) EQ 0 THEN min_cal_solutions=5
IF N_Elements(max_cal_iter) EQ 0 THEN max_cal_iter=10L
IF N_Elements(ref_antenna) EQ 0 THEN ref_antenna=1L
IF N_Elements(cal_convergence_threshold) EQ 0 THEN cal_convergence_threshold=1E-3
IF N_Elements(calibration_origin) EQ 0 THEN calibration_origin='Calibration not transferred'

IF N_Elements(gain_arr_ptr) EQ 0 THEN BEGIN
    gain_arr=Complexarr(n_freq,n_tile)+1 
    gain_arr_ptr=Ptrarr(n_pol,/allocate)
    FOR pol_i=0,n_pol-1 DO *gain_arr_ptr[pol_i]=gain_arr
ENDIF

cal_struct={n_pol:n_pol,n_freq:n_freq,n_tile:n_tile,n_time:n_time,uu:u_loc,vv:v_loc,source_list:source_list,max_iter:max_cal_iter,$
    tile_A:tile_A,tile_B:tile_B,tile_names:tile_names,bin_offset:bin_offset,freq:freq,gain:gain_arr_ptr,n_cal_src:n_cal_src,$
    galaxy_cal:galaxy_cal,min_cal_baseline:min_cal_baseline,max_cal_baseline:max_cal_baseline,n_vis_cal:n_vis_cal,$
    time_avg:cal_time_average,min_solns:min_cal_solutions,ref_antenna:ref_antenna,conv_thresh:cal_convergence_threshold,$
    cal_origin:calibration_origin}
RETURN,cal_struct
END