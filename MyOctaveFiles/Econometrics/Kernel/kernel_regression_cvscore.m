function cvscore = kernel_regression_cvscore(bandwidth, data, depvar, kernel, prewhiten)
	if (nargin < 4) kernel = "__kernel_normal"; endif
	if (nargin < 5) prewhiten = false; endif
	
	fit = kernel_regression(data, depvar, data, exp(bandwidth), kernel, prewhiten, true);
	cvscore = depvar - fit;
	cvscore = cvscore'*cvscore;
endfunction

