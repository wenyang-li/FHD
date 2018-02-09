FUNCTION uvfits_params_simulate,hdr,params_in,sim_baseline_uu=sim_baseline_uu,sim_baseline_vv=sim_baseline_vv,$
    sim_baseline_ww=sim_baseline_ww,sim_baseline_i=sim_baseline_i,sim_baseline_time=sim_baseline_time,$
    sim_tile_locations_x=sim_tile_locations_x,sim_tile_locations_y=sim_tile_locations_y,sim_tile_i=sim_tile_i,n_time=n_time,$
    sim_antenna1=sim_antenna1, sim_antenna2=sim_antenna2, _Extra=extra

;    params={uu:uu_arr,vv:vv_arr,ww:ww_arr,baseline_arr:baseline_arr,time:time}
n_tile=hdr.n_tile
n_baseline=hdr.nbaselines ;excludes time axis!
IF N_Elements(n_time) EQ 0 THEN n_time=1.
freq_arr=hdr.freq_arr 
freq_use=Median(freq_arr)

tile_sim=(N_Elements(sim_tile_locations_x)+N_Elements(sim_tile_locations_y)) GT 0
baseline_sim=(N_Elements(sim_baseline_uu)+N_Elements(sim_baseline_vv)+$
              N_Elements(sim_baseline_ww)+N_Elements(sim_baseline_time)) GT 0
IF tile_sim THEN baseline_sim=0
IF baseline_sim THEN BEGIN
  n_baseline=Max([N_Elements(sim_baseline_uu),N_Elements(sim_baseline_vv),$
      N_Elements(sim_baseline_ww),N_Elements(sim_baseline_time)])
;  IF n_baseline GT hdr.nbaselines THEN BEGIN
    n_baseline/=n_time
    params_in_flag=0 
    n_tile_check=Ceil((-1.+Sqrt(1+8.*n_baseline))/2.) ;solution to n_baselines=n_tile*(n_tile-1)/2 + n_tile
    n_tile=n_tile_check
    hdr.n_tile=n_tile_check
    hdr.nbaselines=n_baseline
;  ENDIF ELSE params_in_flag=Keyword_Set(params_in) 
ENDIF ELSE params_in_flag=Keyword_Set(params_in)

IF params_in_flag THEN BEGIN
    default_uu=params_in.uu
    default_vv=params_in.vv
    default_ww=params_in.ww
    default_i=params_in.baseline_arr
    default_time=params_in.time
    
    n_time=N_Elements(Uniq(default_time))
    n_bt=N_Elements(default_uu)
    n_baseline=Ceil(n_bt/n_time)
ENDIF ELSE BEGIN
    n_default=N_Elements(sim_baseline_uu)>N_Elements(sim_baseline_vv)>N_Elements(sim_baseline_ww)>N_Elements(sim_baseline_i)>N_Elements(sim_baseline_time)
    IF n_default EQ 0 THEN n_default=n_baseline*n_time 
    n_mod=Long(2.^(Ceil(Alog(n_tile+1)/Alog(2.))))
    IF N_Elements(default_uu) EQ 0 THEN default_uu=((Findgen(n_default)+1) mod n_mod)/freq_use 
    IF N_Elements(default_vv) EQ 0 THEN default_vv=(Float(Floor((Findgen(n_default)+1) / n_mod)))/freq_use
    IF N_Elements(default_ww) EQ 0 THEN default_ww=Fltarr(n_default)
    IF N_Elements(default_i) EQ 0 THEN BEGIN
        default_i_single=Lonarr(n_default/n_time)
        ii=0L
        FOR a_i=1,n_tile DO FOR b_i=a_i,n_tile DO BEGIN
            IF ii GE n_default/n_time THEN CONTINUE
            default_i_single[ii]=b_i+a_i*n_mod
            ii+=1L
        ENDFOR
        default_i=default_i_single
        FOR t_i=1,n_time-1 DO default_i=[default_i,default_i_single]
    ENDIF
    IF N_Elements(default_time) EQ 0 THEN default_time=Floor(Lindgen(n_default)/n_baseline)
ENDELSE

IF baseline_sim THEN BEGIN
    IF N_Elements(sim_baseline_uu) EQ 0 THEN sim_baseline_uu=default_uu
    IF N_Elements(sim_baseline_vv) EQ 0 THEN sim_baseline_vv=default_vv
    IF N_Elements(sim_baseline_ww) EQ 0 THEN sim_baseline_ww=default_ww
    IF N_Elements(sim_baseline_i) EQ 0 THEN sim_baseline_i=default_i
    IF N_Elements(sim_baseline_time) EQ 0 THEN sim_baseline_time=default_time
    n_use=Minmax([N_Elements(sim_baseline_uu),N_Elements(sim_baseline_vv),N_Elements(sim_baseline_ww),$
               N_Elements(sim_baseline_i),N_Elements(sim_baseline_time)])
    IF n_use[0] NE n_use[1] THEN BEGIN
        n_use=n_use[0]
        sim_baseline_uu=sim_baseline_uu[0:n_use-1]
        sim_baseline_vv=sim_baseline_vv[0:n_use-1]
        sim_baseline_ww=sim_baseline_ww[0:n_use-1]
        sim_baseline_i=sim_baseline_i[0:n_use-1]
        sim_baseline_time=sim_baseline_time[0:n_use-1]
    ENDIF
ENDIF 

IF tile_sim THEN BEGIN

ENDIF

IF N_Elements(sim_baseline_uu) EQ 0 THEN sim_baseline_uu=default_uu
IF N_Elements(sim_baseline_vv) EQ 0 THEN sim_baseline_vv=default_vv
IF N_Elements(sim_baseline_ww) EQ 0 THEN sim_baseline_ww=default_ww
IF N_Elements(sim_baseline_i) EQ 0 THEN sim_baseline_i=default_i
IF N_Elements(sim_baseline_time) EQ 0 THEN sim_baseline_time=default_time

; Define default antenna numbers
name_mod=2.^((Ceil(Alog(Sqrt(hdr.nbaselines*2.-hdr.n_tile))/Alog(2.)))>Floor(Alog(Min(sim_baseline_i))/Alog(2.)))        
IF N_Elements(sim_antenna1) EQ 0 THEN sim_antenna1=Long(Floor(sim_baseline_i/name_mod)) ;tile numbers start from 1      
IF N_Elements(sim_antenna2) EQ 0 THEN sim_antenna2=Long(Fix(sim_baseline_i mod name_mod))

params={uu:sim_baseline_uu,vv:sim_baseline_vv,ww:sim_baseline_ww,$
        baseline_arr:sim_baseline_i,time:sim_baseline_time,$
        antenna1:sim_antenna1, antenna2:sim_antenna2}
RETURN,params
END