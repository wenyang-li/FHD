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

;+
; :Returns:
;   antenna structure with response tag (pointer array) updated with single antenna response or phased array response
;
; :Params:
;    obs : in, required, type=structure
;       FHD **obs** structure
;    antenna : in, required, type=structure array
;       FHD **antenna** structure
;
; :Keywords:
;    za_arr : in, required, type=fltarr
;       Array of zenith angles (degrees)
;    az_arr : in, required, type=fltarr
;       Array of azimuth angles (degrees)       
;
; :History:
;-
FUNCTION general_antenna_response, obs, antenna, za_arr=za_arr, az_arr=az_arr

n_ant=obs.n_tile
n_ant_pol=antenna[0].n_pol ;this needs to be the same for all antennas!
double_flag=Tag_exist(obs, 'double_precision') ? obs.double_precision : 0 
icomp=Complex(0, 1, double=double_flag)
Pi=double_flag ? !DPi : !Pi
c_light_vacuum=299792458.

response=Ptrarr(n_ant_pol,n_ant)
nfreq_bin=antenna[0].nfreq_bin
freq_center=antenna[0].freq
IF N_Elements(psf_image_dim) EQ 0 THEN psf_image_dim = (size(za_arr,/dimension))[0]

proj_east=Sin(za_arr*!DtoR)*Sin(az_arr*!DtoR) & proj_east_use=Reform(proj_east,(psf_image_dim)^2.)
proj_north=Sin(za_arr*!DtoR)*Cos(az_arr*!DtoR) & proj_north_use=Reform(proj_north,(psf_image_dim)^2.)
proj_z=Cos(za_arr*!DtoR) & proj_z_use=Reform(proj_z,(psf_image_dim)^2.)

;;initialize antenna structure
;antenna_str={n_pol:n_ant_pol,antenna_type:instrument,model_version:beam_model_version,freq:freq_center,nfreq_bin:nfreq_bin,$
;    n_ant_elements:0,Jones:Ptrarr(n_ant_pol,n_ant_pol,nfreq_bin),coupling:Ptrarr(n_ant_pol,nfreq_bin),gain:Ptrarr(n_ant_pol),coords:Ptrarr(3),$
;    delays:Ptr_new(),size_meters:0.,height:0.,group_id:Lonarr(n_ant_pol)-1}

group_arr=antenna.group_id
FOR pol_i=0,n_ant_pol-1 DO BEGIN
    g_hist=histogram(group_arr[pol_i,*],min=0,/binsize,reverse_ind=g_ri)
    n_group=N_Elements(g_hist)
    FOR grp_i=0L,n_group-1 DO BEGIN
        ng=g_hist[grp_i]
        IF ng EQ 0 THEN CONTINUE
        g_inds=g_ri[g_ri[grp_i]:g_ri[grp_i+1]-1]
        ref_i=g_inds[0]
        n_ant_elements=antenna[ref_i].n_ant_elements
        coupling=antenna[ref_i].coupling
        gain=antenna[ref_i].gain
        xc_arr=*(antenna[ref_i].coords[0])
        yc_arr=*(antenna[ref_i].coords[1])
        zc_arr=*(antenna[ref_i].coords[2])
        delays=*(antenna[ref_i].delays)
        
        ;phase of each dipole for the source (relative to the beamformer settings)
        D_d=(proj_east_use#xc_arr+proj_north_use#yc_arr+proj_z_use#zc_arr)
        D_d=Reform(D_d,psf_image_dim,psf_image_dim,n_ant_elements)
        
        response_grp=Ptrarr(nfreq_bin)
        
        FOR freq_i=0L,nfreq_bin-1 DO BEGIN
            response=Complexarr(psf_image_dim, psf_image_dim, double=double_flag)
            Kconv=(2.*Pi)*(freq_center[freq_i]/c_light_vacuum) 
            antenna_gain_arr=Exp(-icomp*Kconv*D_d)
            voltage_delay=Exp(icomp*2.*Pi*delays*(freq_center[freq_i])*Reform((*gain[pol_i])[freq_i,*])) 
            meas_current=(*coupling[pol_i,freq_i])#voltage_delay
            zenith_norm=Mean((*coupling[pol_i,freq_i])#Replicate(1.,n_ant_elements))
            meas_current/=zenith_norm
            
            FOR ii=0L,n_ant_elements-1 DO BEGIN
                response+=antenna_gain_arr[*,*,ii]*meas_current[ii]/n_ant_elements
            ENDFOR
            response_grp[freq_i]=Ptr_new(response)
        ENDFOR
        FOR g_i=0L,ng-1 DO antenna[g_inds[g_i]].response[pol_i,*]=response_grp
    ENDFOR
ENDFOR
RETURN,antenna
END