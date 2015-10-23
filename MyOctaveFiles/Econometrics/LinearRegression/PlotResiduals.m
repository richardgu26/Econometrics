# Plots residuals. If "outfile" is provided as an argument
# the plot is written to and eps file with the provided name".

function residuals = PlotResiduals(y, x, outfile)

	[b, sigsq, e] = ols(y,x);
	plot(e, "or;Residuals;")
	title("Regression residuals");
	if nargin > 2 print(outfile, "-depsc2"); endif
endfunction
