function cvscore = kernel_density_cvscore(bandwidth, data, kernel)
		dens = kernel_density(data, data, exp(bandwidth), kernel, true, 0, 0, chol(cov(data)), kernel);
		dens = dens + eps; # some kernels can assign zero density
		cvscore = -mean(log(dens));
endfunction

