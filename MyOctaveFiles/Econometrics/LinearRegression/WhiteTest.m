# Perform White's test for heteroscedasticity
# Inputs:
#	e: the regression residuals
#	z: the variables suspected to influence the error variance
function [value, pvalue] = WhiteTest(e, z)
	
	n = rows(e);
	
	# If z doesn't include a constant, insert it
	test = min(std(z)) == 0;
	if test z = [ones(n,1) z]; endif
	k = columns(z) - 1;
	
	esq = e .^ 2;
	rsq = rsquare(esq, z); # the r-square of the artificial regression
	# test statistic and p-value
	value = n*rsq;
	pvalue = 1 - chi2cdf(value,k); 
	result = [value pvalue];

	# print it out
	rlabels = char("White's test");
	clabels = char("Value","p-value");
	prettyprint(result, rlabels, clabels); 

endfunction
