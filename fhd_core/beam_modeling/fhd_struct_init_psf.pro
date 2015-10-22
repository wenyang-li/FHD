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
; :Description:
;    Initializes the structure containing gridded PSF data for a visibility. 
;
; :Returns:
;   FHD **psf** structure
;
;
; :Keywords:
;    beam_ptr : in, optional, type=Pointer
;    xvals
;    yvals
;    fbin_i
;    psf_resolution
;    psf_dim
;    n_pol
;    n_freq
;    freq_cen
;    pol_norm
;    freq_norm
;    group_arr
;    complex_flag : in, optional
;
; :History:
;-
FUNCTION fhd_struct_init_psf,beam_ptr=beam_ptr,complex_flag=complex_flag,$
    xvals=xvals,yvals=yvals,fbin_i=fbin_i,psf_resolution=psf_resolution,psf_dim=psf_dim,$
    n_pol=n_pol,n_freq=n_freq,freq_cen=freq_cen,pol_norm=pol_norm,freq_norm=freq_norm,group_arr=group_arr
    
IF N_Elements(n_pol) EQ 0 THEN n_pol=1 ELSE n_pol=Fix(n_pol)
IF N_Elements(n_freq) EQ 0 THEN n_freq=1 ELSE n_freq=Fix(n_freq)
IF N_Elements(freq_cen) EQ 0 THEN freq_cen=Fltarr(n_freq) ELSE freq_cen=Float(freq_cen)
IF N_Elements(beam_ptr) EQ 0 THEN beam_ptr=Ptr_new() ;target will have dimensions (npol,nfreq,nbaselines)
IF N_Elements(psf_resolution) EQ 0 THEN psf_resolution=1. ELSE psf_resolution=Float(psf_resolution);over-resolution
IF N_Elements(psf_dim) EQ 0 THEN psf_dim=1. ELSE psf_dim=Float(psf_dim)
IF N_Elements(xvals) EQ 0 THEN xvals=Ptrarr(1) ;will have dimensions of (resolution,resolution)
IF N_Elements(pol_norm) EQ 0 THEN pol_norm=replicate(1.,n_pol) ELSE pol_norm=Float(pol_norm)
IF N_Elements(freq_norm) EQ 0 THEN freq_norm=replicate(1.,n_freq) ELSE freq_norm=Float(freq_norm)
IF N_Elements(fbin_i) EQ 0 THEN fbin_i=Lonarr(1) ELSE fbin_i=Long(fbin_i)
IF N_Elements(complex_flag) EQ 0 THEN complex_flag=1
IF N_Elements(group_arr) EQ 0 THEN group_arr=Lonarr(size(beam_arr,/dimension))

struct={beam_ptr:beam_ptr,xvals:xvals,yvals:yvals,pnorm:pol_norm,fnorm:freq_norm,id:group_arr,$
    fbin_i:fbin_i,resolution:psf_resolution,dim:psf_dim,complex_flag:complex_flag,n_pol:n_pol,n_freq:n_freq,freq:freq_cen}
RETURN,struct
END