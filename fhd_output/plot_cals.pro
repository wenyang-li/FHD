pro plot_cals,cal,obs,cal_res=cal_res,vis_baseline_hist=vis_baseline_hist,file_path_base=file_path_base,no_ps=no_ps
; Make plot of the cal solutions, save to png
; cgPS_Open/cgPS_Close write .ps first, then converts to png. Supply .png
; filename to automatically overwrite .ps.

no_ps=1
IF Keyword_Set(no_ps) THEN ext_name='.png' ELSE ext_name='.ps'
phase_filename=file_path_base+'_cal_phase'+ext_name
amp_filename=file_path_base+'_cal_amp'+ext_name
res_filename=file_path_base+'_cal_residual'+ext_name
vis_hist_filename=file_path_base+'_cal_hist'+ext_name
IF file_test(file_dirname(file_path_base),/directory) EQ 0 THEN file_mkdir,file_dirname(file_path_base)

IF N_Elements(cal_res) EQ 0 THEN IF tag_exist(cal,'gain_residual') THEN cal_res=cal.gain_residual
tile_names = cal.tile_names
n_tiles=obs.n_tile
n_pol=cal.n_pol
obs2=*obs.baseline_info
tile_use=obs2.tile_use
freq_use=obs2.freq_use
freq_i_use=where(freq_use,nf_use) & IF nf_use EQ 0 THEN freq_i_use=lindgen(obs.n_freq)

gains0 = *cal.gain[0] ; save on typing
IF n_pol GT 1 THEN gains1 = *cal.gain[1]

gains0=gains0[freq_i_use,*]
IF n_pol GT 1 THEN gains1=gains1[freq_i_use,*]
freq=cal.freq[freq_i_use]
freq=freq/10^6. ; in MHz

plot_pos=calculate_plot_positions(.8, nplots=128, /no_colorbar, ncol=17, nrow=9, plot_margin = [.05, .05, .02, .25]); use 1 extra row/column to leave room for title and axes
plot_pos=reform(plot_pos.plot_pos,17,9,4)
; Now remove unwanted row/column
plot_pos = plot_pos[1:*,*,*]
plot_pos = plot_pos[*,1:*,*]
plot_pos = reform(plot_pos,128,4)
; Shift a half a width over and up
width = plot_pos[1,0]-plot_pos[0,0]
height = abs(plot_pos[16,1]-plot_pos[0,1])
plot_pos[*,0] -= width/2
plot_pos[*,1] += height/2
plot_pos[*,2] -= width/2
plot_pos[*,3] += height/2

cgPS_Open,phase_filename,scale_factor=2,/quiet,/nomatch

n_baselines=obs.nbaselines
tile_A=obs2.tile_A[0:n_baselines-1]
tile_B=obs2.tile_B[0:n_baselines-1]
tile_exist=(histogram(tile_A,min=1,/bin,max=(max(tile_A)>max(tile_B)))+histogram(tile_B,min=1,/bin,max=(max(tile_A)>max(tile_B))))<1

ytickv=[-!pi,0,!pi]
ytickname=['-!9p!X','0','!9p!X']
yrange=[-1.5*!pi,1.5*!pi]
xtickv=[ceil(min(freq)/10)*10,floor(max(freq)/10)*10]
;xtickv=[min(freq),max(freq)]
xtickname=strtrim(round(xtickv),2)
xrange=[min(freq)-(max(freq)-min(freq))/8,max(freq)+(max(freq)-min(freq))/8]

FOR tile_i=0L,n_tiles-1 DO BEGIN
    tile_name=tile_names[tile_i]
;    rec=Floor(tile_name/10)
;    tile=tile_name mod 10

    IF tile_exist[tile_i] EQ 0 THEN BEGIN
      ; no tile found... must have been flagged in pre-processing
      axiscolor='grey'
      cgplot,1,title=strtrim(tile_name,2),XTICKFORMAT="(A1)",YTICKFORMAT="(A1)",position=plot_pos[tile_i,*],$
        /noerase,charsize=.5,axiscolor=axiscolor
    ENDIF ELSE BEGIN
      IF tile_use[tile_i] EQ 0 THEN axiscolor='red' ELSE axiscolor='black'
      IF tile_i EQ cal.ref_antenna THEN axiscolor='blue'
      IF ~(tile_i mod 16) THEN BEGIN
        IF (tile_i gt (n_tiles-17)) THEN BEGIN
          ; both axes
          cgplot,freq,phunwrap(atan(gains0[*,tile_i],/phase)),color='blue',title=strtrim(tile_name,2),$
            xticks=1,xtickv=xtickv,xtickname=xtickname,yticks=2,ytickv=ytickv,ytickname=ytickname,position=plot_pos[tile_i,*],$
            yticklen=0.04,yrange=yrange,xrange=xrange,charsize=.5,/noerase,axiscolor=axiscolor,psym=3
        ENDIF ELSE BEGIN
          ; just the y-axis
          cgplot,freq,phunwrap(atan(gains0[*,tile_i],/phase)),color='blue',title=strtrim(tile_name,2),$
            xticks=1,xtickv=xtickv,XTICKFORMAT="(A1)",yticks=2,ytickv=ytickv,ytickname=ytickname,position=plot_pos[tile_i,*],$
            yticklen=0.04,yrange=yrange,xrange=xrange,charsize=.5,/noerase,axiscolor=axiscolor,psym=3
        ENDELSE
      ENDIF ELSE BEGIN
        IF (tile_i gt (n_tiles-17)) THEN BEGIN
          ; just x-axis
          cgplot,freq,phunwrap(atan(gains0[*,tile_i],/phase)),color='blue',title=strtrim(tile_name,2),$
            xticks=1,xtickv=xtickv,yticks=2,ytickv=ytickv,YTICKFORMAT="(A1)",position=plot_pos[tile_i,*],$
            yticklen=0.04,yrange=yrange,xrange=xrange,charsize=.5,/noerase,axiscolor=axiscolor,psym=3
        ENDIF ELSE BEGIN
          ; No axes
          cgplot,freq,phunwrap(atan(gains0[*,tile_i],/phase)),color='blue',title=strtrim(tile_name,2),$
            xticks=1,xtickv=xtickv,XTICKFORMAT="(A1)",yticks=2,ytickv=ytickv,YTICKFORMAT="(A1)",position=plot_pos[tile_i,*],yrange=yrange,xrange=xrange,$
            yticklen=0.04,charsize=.5,/noerase,axiscolor=axiscolor,psym=3
        ENDELSE
      ENDELSE
      IF n_pol GT 1 THEN cgoplot,freq,phunwrap(atan(gains1[*,tile_i],/phase)),color='red',psym=3
    ENDELSE
ENDFOR
cgtext,.4,max(plot_pos[*,3]+height/4),obs.obsname,/normal
cgPS_Close,/png,Density=75,Resize=100.,/allow_transparent,/nomessage
    
cgPS_Open,amp_filename,scale_factor=2,/quiet,/nomatch

IF n_pol GT 1 THEN max_amp = mean(abs([gains0,gains1])) + 2*stddev(abs([gains0,gains1])) $
    ELSE max_amp = Mean(abs(gains0)) + 2*stddev(abs(gains0))
yrange=[0,max_amp]
ytickv=[0,max_amp/2,max_amp]

FOR tile_i=0L,n_tiles-1 DO BEGIN
    tile_name=tile_names[tile_i]
;    rec=Floor(tile_name/10)
;    tile=tile_name mod 10
    
    IF tile_exist[tile_i] EQ 0  THEN BEGIN
      ; no tile found... must have been flagged in pre-processing
      axiscolor='grey'
      cgplot,1,title=strtrim(tile_name,2),XTICKFORMAT="(A1)",YTICKFORMAT="(A1)",position=plot_pos[tile_i,*],$
        /noerase,charsize=.5,axiscolor=axiscolor
    ENDIF ELSE BEGIN
      IF tile_use[tile_i] EQ 0 THEN axiscolor='red' ELSE axiscolor='black'
      IF ~(tile_i mod 16) THEN BEGIN
        IF (tile_i gt (n_tiles-17)) THEN BEGIN
          ; both axes
          cgplot,freq,abs(gains0[*,tile_i]),color='blue',title=strtrim(tile_name,2),$
            xticks=1,xtickv=xtickv,xtickname=xtickname,yticks=2,ytickv=ytickv,position=plot_pos[tile_i,*],$
            yticklen=0.04,yrange=yrange,xrange=xrange,charsize=.5,/noerase,axiscolor=axiscolor,psym=3
        ENDIF ELSE BEGIN
          ; just the y-axis
          cgplot,freq,abs(gains0[*,tile_i]),color='blue',title=strtrim(tile_name,2),$
            xticks=1,xtickv=xtickv,XTICKFORMAT="(A1)",yticks=2,ytickv=ytickv,position=plot_pos[tile_i,*],$
            yticklen=0.04,yrange=yrange,xrange=xrange,charsize=.5,/noerase,axiscolor=axiscolor,psym=3
        ENDELSE
      ENDIF ELSE BEGIN
        IF (tile_i gt (n_tiles-17)) THEN BEGIN
          ; just x-axis
          cgplot,freq,abs(gains0[*,tile_i]),color='blue',title=strtrim(tile_name,2),$
            xticks=1,xtickv=xtickv,yticks=2,ytickv=ytickv,YTICKFORMAT="(A1)",position=plot_pos[tile_i,*],$
            yticklen=0.04,yrange=yrange,xrange=xrange,charsize=.5,/noerase,axiscolor=axiscolor,psym=3
        ENDIF ELSE BEGIN
          ; No axes
          cgplot,freq,abs(gains0[*,tile_i]),color='blue',title=strtrim(tile_name,2),$
            xticks=1,xtickv=xtickv,XTICKFORMAT="(A1)",yticks=2,ytickv=ytickv,YTICKFORMAT="(A1)",position=plot_pos[tile_i,*],yrange=yrange,xrange=xrange,$
            yticklen=0.04,charsize=.5,/noerase,axiscolor=axiscolor,psym=3
        ENDELSE
      ENDELSE
      IF n_pol GT 1 THEN cgoplot,freq,abs(gains1[*,tile_i]),color='red',psym=3
    ENDELSE
ENDFOR

cgtext,.4,max(plot_pos[*,3]+height/4),obs.obsname,/normal
cgPS_Close,/png,Density=75,Resize=100.,/allow_transparent,/nomessage

IF Keyword_Set(cal_res) THEN BEGIN
    res_type=Size(*cal_res[0],/type)
    IF res_type EQ 8 THEN gains0r=*cal_res.gain[0] ELSE gains0r=*cal_res[0]
    gains0r=*cal_res.gain[0]
    gains0r=gains0r[freq_i_use,*]
    
    gains0_orig=gains0r+gains0
    gains0r=Abs(gains0_orig)-Abs(gains0)
    IF n_pol GT 1 THEN BEGIN    
        IF res_type EQ 8 THEN gains1r=*cal_res.gain[1] ELSE gains1r=*cal_res[1]
        gains1r=gains1r[freq_i_use,*]
        gains1_orig=gains1r+gains1
        gains1r=Abs(gains1_orig)-Abs(gains1)
        max_amp = mean(abs([gains0r,gains1r])) + 2*stddev(abs([gains0r,gains1r]))
    ENDIF ELSE max_amp = mean(abs(gains0r)) + 2*stddev(abs(gains0r))
    
;    gains0r=*cal_res.gain[0]
;    gains0r=gains0r[freq_i_use,*]
;    sign0r=Real_part(gains0r)*weight_invert(Abs(real_part(gains0r)))
;    gains0r=Abs(gains0r)*sign0r
;    IF n_pol GT 1 THEN BEGIN    
;        gains1r=*cal_res.gain[1]
;        gains1r=gains1r[freq_i_use,*]
;        sign1r=Real_part(gains1r)*weight_invert(Abs(real_part(gains1r)))
;        gains1r=Abs(gains1r)*sign1r
;        max_amp = mean(abs([gains0r,gains1r])) + 2*stddev(abs([gains0r,gains1r]))
;    ENDIF ELSE max_amp = mean(abs(gains0r)) + 2*stddev(abs(gains0r))
    
    IF max_amp GT 0 THEN BEGIN
        cgPS_Open,res_filename,scale_factor=2,/quiet,/nomatch
        yrange=[-max_amp,max_amp]
        ytickv=[-max_amp,0,max_amp]
        
        FOR tile_i=0L,n_tiles-1 DO BEGIN
            tile_name=tile_names[tile_i]
;            rec=Floor(tile_name/10)
;            tile=tile_name mod 10
    
            IF tile_exist[tile_i] EQ 0  THEN BEGIN
              ; no tile found... must have been flagged in pre-processing
              axiscolor='grey'
              cgplot,1,title=strtrim(tile_name,2),XTICKFORMAT="(A1)",YTICKFORMAT="(A1)",position=plot_pos[tile_i,*],$
                /noerase,charsize=.5,axiscolor=axiscolor
            ENDIF ELSE BEGIN
              IF tile_use[tile_i] EQ 0 THEN axiscolor='red' ELSE axiscolor='black'
              IF ~(tile_i mod 16) THEN BEGIN
                IF (tile_i gt (n_tiles-17)) THEN BEGIN
                  ; both axes
                  cgplot,freq,gains0r[*,tile_i],color='blue',title=strtrim(tile_name,2),$
                    xticks=1,xtickv=xtickv,xtickname=xtickname,yticks=2,ytickv=ytickv,position=plot_pos[tile_i,*],$
                    yticklen=0.04,yrange=yrange,xrange=xrange,charsize=.5,/noerase,axiscolor=axiscolor,psym=3
                ENDIF ELSE BEGIN
                  ; just the y-axis
                  cgplot,freq,gains0r[*,tile_i],color='blue',title=strtrim(tile_name,2),$
                    xticks=1,xtickv=xtickv,XTICKFORMAT="(A1)",yticks=2,ytickv=ytickv,position=plot_pos[tile_i,*],$
                    yticklen=0.04,yrange=yrange,xrange=xrange,charsize=.5,/noerase,axiscolor=axiscolor,psym=3
                ENDELSE
              ENDIF ELSE BEGIN
                IF (tile_i gt (n_tiles-17)) THEN BEGIN
                  ; just x-axis
                  cgplot,freq,gains0r[*,tile_i],color='blue',title=strtrim(tile_name,2),$
                    xticks=1,xtickv=xtickv,yticks=2,ytickv=ytickv,YTICKFORMAT="(A1)",position=plot_pos[tile_i,*],$
                    yticklen=0.04,yrange=yrange,xrange=xrange,charsize=.5,/noerase,axiscolor=axiscolor,psym=3
                ENDIF ELSE BEGIN
                  ; No axes
                  cgplot,freq,gains0r[*,tile_i],color='blue',title=strtrim(tile_name,2),$
                    xticks=1,xtickv=xtickv,XTICKFORMAT="(A1)",yticks=2,ytickv=ytickv,YTICKFORMAT="(A1)",position=plot_pos[tile_i,*],yrange=yrange,xrange=xrange,$
                    yticklen=0.04,charsize=.5,/noerase,axiscolor=axiscolor,psym=3
                ENDELSE
              ENDELSE
              IF n_pol GT 1 THEN cgoplot,freq,gains1r[*,tile_i],color='red',psym=3
            ENDELSE
        ENDFOR
        
        cgtext,.4,max(plot_pos[*,3]+height/4),obs.obsname,/normal
        cgPS_Close,/png,Density=75,Resize=100.,/allow_transparent,/nomessage
    ENDIF
ENDIF

IF size(vis_baseline_hist,/type) EQ 8 THEN BEGIN
   ratio=vis_baseline_hist.vis_res_ratio_mean ; just save some typing
   sigma=vis_baseline_hist.vis_res_sigma
   base_len=vis_baseline_hist.baseline_length

   cgPS_Open,filename=vis_hist_filename,/quiet,/nomatch
   !p.multi=[0,2,1]
   FOR pol=0,1 DO BEGIN
      cgplot,base_len,ratio[pol,*],color='red',/xlog,yrange=[0,max(ratio+sigma)]
      cgoplot,base_len,ratio[pol,*]+sigma[pol,*],linestyle=2
      cgoplot,base_len,ratio[pol,*]-sigma[pol,*],linestyle=2
   ENDFOR
   cgPS_Close,/png,Density=75,Resize=100.,/allow_transparent,/nomessage
ENDIF

END    
