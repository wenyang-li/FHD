pro mwa_apply_pfb_aliasing_matrix_then_inverse_bandpass_correction,vis_arr

	pfb = read_table('/nfs/mwa-09/r1/abrahamn/FHD/instrument_config/mwa_pfbmat_highband_80kHz.txt') ; 1st dim is output chan, 2nd dim is input chan

	bandpasscor = 1./total(pfb,2) ; sum over input channels
	numvis = n_elements(*vis_arr[0])/384L
	bandpasscormat = replicate(1,numvis)##bandpasscor

	for i=0,n_elements(vis_arr)-1 do *vis_arr[i]=(pfb#(*vis_arr[i]))*bandpasscormat
end


