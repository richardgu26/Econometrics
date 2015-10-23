% octave-matlab compatibility function
function [b s r] = regress(y,x)
	[b s r] = ols(y,x);
end
