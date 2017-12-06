FUNCTION vis_param_extract,params,header


uu_arr=Double(reform(params[header.uu_i,*]))
vv_arr=Double(reform(params[header.vv_i,*]))
ww_arr=Double(reform(params[header.ww_i,*]))
time=reform(params[header.date_i,*])
IF header.baseline_i GE 0 THEN $
    baseline_arr=reform(params[header.baseline_i,*]) $
    ELSE baseline_arr = Lonarr(N_Elements(time))
IF Tag_exist(header, "ant1_i") THEN $
    IF header.ant1_i GE 0 THEN ant1_arr = reform(params[header.ant1_i,*])
IF N_Elements(ant1_arr) EQ 0 THEN ant1_arr = Lonarr(N_Elements(time))
IF Tag_exist(header, "ant2_i") THEN $
    IF header.ant2_i GE 0 THEN ant2_arr = reform(params[header.ant2_i,*])
IF N_Elements(ant2_arr) EQ 0 THEN ant2_arr = Lonarr(N_Elements(time))

struct={uu:uu_arr,vv:vv_arr,ww:ww_arr,baseline_arr:baseline_arr,time:time, antenna1:ant1_arr, antenna2:ant2_arr}
RETURN,struct
END