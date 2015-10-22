;+
; :Copyright: (c) 2014, Sullivan, I., Morales, M., Hazelton, B.
;All rights reserved.
;Please acknowledge use of this software by citing:
;Sullivan I. S., Morales M. F., Hazelton B. J. et al
;	"Fast Holographic Deconvolution: a new technique for precision radio interferometry"
;	Astrophysical Journal 759 17 (2012)
;
;Redistribution and use in source and binary forms, with or without
;modification, are permitted provided that the following conditions are met:
;
;* Redistributions of source code must retain the above copyright notice, this
;  list of conditions and the following disclaimer.
;
;* Redistributions in binary form must reproduce the above copyright notice,
;  this list of conditions and the following disclaimer in the documentation
;  and/or other materials provided with the distribution.
;
; THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
; AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
; IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
; DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
; FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
; DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
; SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
; CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
; OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
; OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
;-

FUNCTION fhd_struct_init_antenna, obs, beam_model_version=beam_model_version,$
    psf_resolution=psf_resolution, psf_intermediate_res=psf_intermediate_res,$
    psf_image_resolution=psf_image_resolution, timing=timing,$
    psf_dim=psf_dim, psf_max_dim=psf_max_dim, beam_offset_time=beam_offset_time, _Extra=extra
;+
; :Description:
;    Builds the beam model for individual antenna elements. Generates the full tile response for a phased array.
;
; :Returns:
;   FHD **antenna** structure
; :Params:
;    obs : in, required, type=structure
;       FHD **obs** structure
;
; :Keywords:
;    psf_resolution : in, type=Float, default=16.0
;    beam_model_version : in, type=Int, default=1 
;    psf_image_resolution : in/out, type=Float, default=10.0
;    beam_offset_time : in, optional, type=Float, default=0.0
;    psf_max_dim : in, optional, type=Int
;    psf_dim : out, type=Int
;    psf_intermediate_res : out, type=Int
;    timing : out, optional, type=Double
;    _Extra
;
; :History:
;-
compile_opt idl2,strictarrsubs   
t0=Systime(1)

IF N_Elements(beam_model_version) EQ 0 THEN beam_model_version=1
instrument=obs.instrument
tile_gain_fn=instrument+'_beam_setup_gain' ;mwa_beam_setup_gain paper_beam_setup_gain hera_beam_setup_gain
tile_init_fn=instrument+'_beam_setup_init' ;mwa_beam_setup_init paper_beam_setup_init hera_beam_setup_init
n_tiles=obs.n_tile
n_freq=obs.n_freq
n_pol=obs.n_pol
n_ant_pol=2 ;use as default, since almost all instruments have two instrumental polarizations (either linear or circular)

obsra=obs.obsra
obsdec=obs.obsdec
zenra=obs.zenra
zendec=obs.zendec
obsx=obs.obsx
obsy=obs.obsy
;phasera=obs.phasera
;phasedec=obs.phasedec
Jdate=obs.Jd0
IF Keyword_Set(beam_offset_time) THEN Jdate_use=Jdate+beam_offset_time/24./3600. ELSE Jdate_use=Jdate
frequency_array=(*obs.baseline_info).freq
freq_bin_i=(*obs.baseline_info).fbin_i
nfreq_bin=Max(freq_bin_i)+1

tile_A=(*obs.baseline_info).tile_A
tile_B=(*obs.baseline_info).tile_B
nbaselines=obs.nbaselines
ant_names=tile_A[uniq(tile_A[0:nbaselines-1],Sort(tile_A[0:nbaselines-1]))]


dimension=obs.dimension
elements=obs.elements
kbinsize=obs.kpix
;kx_span=kbinsize*dimension ;Units are # of wavelengths
;ky_span=kx_span
degpix=obs.degpix
astr=obs.astr

speed_light=299792458. ;speed of light, in meters/second
IF N_Elements(psf_resolution) EQ 0 THEN psf_resolution=16.0 ;=32? ;super-resolution factor
IF N_Elements(psf_image_resolution) EQ 0 THEN psf_image_resolution=10.0
Eq2Hor,obsra,obsdec,Jdate_use,obsalt,obsaz,lat=obs.lat,lon=obs.lon,alt=obs.alt ; this may or may not include refraction
obsalt=Float(obsalt)
obsaz=Float(obsaz)
obsza=90.-obsalt

freq_center=fltarr(nfreq_bin)
FOR fi=0L,nfreq_bin-1 DO BEGIN
    fi_i=where(freq_bin_i EQ fi,n_fi)
    IF n_fi EQ 0 THEN freq_center[fi]=Interpol(frequency_array,freq_bin_i,fi) ELSE freq_center[fi]=Median(frequency_array[fi_i])
ENDFOR

;initialize antenna structure
antenna_str={n_pol:n_ant_pol,antenna_type:instrument,names:ant_names,model_version:beam_model_version,freq:freq_center,nfreq_bin:nfreq_bin,$
    n_ant_elements:0,Jones:Ptrarr(n_ant_pol,n_ant_pol,nfreq_bin),coupling:Ptrarr(n_ant_pol,nfreq_bin),gain:Ptrarr(n_ant_pol),coords:Ptrarr(3),$
    delays:Ptr_new(),size_meters:0.,height:0.,response:Ptrarr(n_ant_pol,nfreq_bin),group_id:Lonarr(n_ant_pol)-1}
    
;update structure with instrument-specific values, and return as a structure array, with an entry for each tile/antenna
;first, update to include basic configuration data
antenna=Call_function(tile_init_fn,obs,antenna_str,_Extra=extra) ;mwa_beam_setup_init

psf_dim=Ceil((Max(antenna.size_meters)*2.*Max(frequency_array)/speed_light)/kbinsize/Cos(obsza*!DtoR))  
psf_dim=Ceil(psf_dim/2.)*2. ;dimension MUST be even

IF Keyword_Set(psf_max_dim) THEN BEGIN
    psf_max_dim=Ceil(psf_max_dim/2.)*2 ;dimension MUST be even
    IF psf_max_dim LT psf_dim THEN print,'Warning! PSF dim cut to '+Strn(psf_max_dim)+', fit dim was '+Strn(psf_dim)
    psf_dim=psf_dim<psf_max_dim
ENDIF

psf_intermediate_res=(Ceil(Sqrt(psf_resolution)/2)*2.)<psf_resolution
psf_image_dim=psf_dim*psf_image_resolution*psf_intermediate_res ;use a larger box to build the model than will ultimately be used, to allow higher resolution in the initial image-space beam model
psf_superres_dim=psf_dim*psf_resolution
psf_scale=dimension*psf_intermediate_res/psf_image_dim

;xvals_uv_superres=meshgrid(psf_superres_dim,psf_superres_dim,1)/(Float(psf_resolution)/psf_intermediate_res)-Floor(psf_dim/2)*psf_intermediate_res+Floor(psf_image_dim/2)
;yvals_uv_superres=meshgrid(psf_superres_dim,psf_superres_dim,2)/(Float(psf_resolution)/psf_intermediate_res)-Floor(psf_dim/2)*psf_intermediate_res+Floor(psf_image_dim/2)

xvals_celestial=meshgrid(psf_image_dim,psf_image_dim,1)*psf_scale-psf_image_dim*psf_scale/2.+obsx
yvals_celestial=meshgrid(psf_image_dim,psf_image_dim,2)*psf_scale-psf_image_dim*psf_scale/2.+obsy
;turn off refraction for speed, then make sure it is also turned off in Eq2Hor below
apply_astrometry, obs, x_arr=xvals_celestial, y_arr=yvals_celestial, ra_arr=ra_arr, dec_arr=dec_arr, /xy2ad, /ignore_refraction
;xy2ad,xvals_celestial,yvals_celestial,astr,ra_arr,dec_arr
valid_i=where(Finite(ra_arr),n_valid)
ra_use=ra_arr[valid_i]
dec_use=dec_arr[valid_i]

;NOTE: Eq2Hor REQUIRES Jdate_use to have the same number of elements as RA and Dec for precession!!
;;NOTE: The NEW Eq2Hor REQUIRES Jdate_use to be a scalar! They created a new bug when they fixed the old one
Eq2Hor,ra_use,dec_use,Jdate_use,alt_arr1,az_arr1,lat=obs.lat,lon=obs.lon,alt=obs.alt,precess=1,/nutate, refract=0
za_arr=fltarr(psf_image_dim,psf_image_dim)+90. & za_arr[valid_i]=90.-alt_arr1
az_arr=fltarr(psf_image_dim,psf_image_dim) & az_arr[valid_i]=az_arr1

;now, update antenna structure to include gains
antenna=Call_function(tile_gain_fn,obs,antenna,za_arr=za_arr,az_arr=az_arr,psf_image_dim=psf_image_dim,Jdate_use=Jdate_use,_Extra=extra) ;mwa_beam_setup_gain

;Finally, update antenna structure to include the response of each antenna
antenna=general_antenna_response(obs,antenna,za_arr=za_arr,az_arr=az_arr)

timing=Systime(1)-t0
RETURN,antenna
END