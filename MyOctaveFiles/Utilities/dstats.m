# ordinary decriptive statistics
# can call with or without row labels

function [info, clabels] = dstats(x, rlabels=0,silent=false);

info = [mean(x);std(x); skew(x); kurtosis(x); min(x); max(x)]; 
info = info';

clabels = char('mean','st. dev.','skew','kurtosis','min','max');
if !silent
	if !ischar(rlabels) prettyprint_c(info, clabels); endif
	if ischar(rlabels) prettyprint(info, rlabels, clabels); endif
endif

endfunction
