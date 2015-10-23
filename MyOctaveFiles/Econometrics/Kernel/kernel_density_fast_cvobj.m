function cvscore = kernel_density_fast_cvobj(logbandwidth, data_to_fit, data)
		dens = kernel_density(data_to_fit, data, exp(logbandwidth));
		dens = dens + eps; # some kernels can assign zero density
		cvscore = -mean(log(dens));
endfunction
