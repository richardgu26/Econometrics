# Perform Breuch-Godfrey test for autocorrelation
# Inputs:
#	e: the regression residuals
#	p: the number of lags to use
#	x: the conditioning variables in the model
function [value, pvalue] = BreuschGodfreyTest(e, p, x)
	n = rows(e);
	# get the lagged residuals
	elag = lags(e,p);
	
	# drop appropriate rows
	e = e(p+1:n,:);
	x = [x elag];
	x = x(p+1:n,:);
	n = rows(x);
	# the artifical regression
	rsq = rsquare(e, x); # the r-square of the artificial regression

	# test statistic and p-value
	value = n*rsq;
	pvalue = 1 - chi2cdf(value,p); 
	result = [value pvalue];

	# print it out
	rlabels = char("Breusch-Godfrey test");
	clabels = char("Value","p-value");
	prettyprint(result, rlabels, clabels); 

endfunction
