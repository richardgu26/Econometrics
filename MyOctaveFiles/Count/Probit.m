function [logdensity, score, prediction] = Probit(theta, data)
	y = data(:,1);
	x = data(:,2:columns(data));
	prob_1 = normcdf(x*theta);
	logdensity = y .* log(prob_1) + (1 - y) .* log(1 - prob_1);
	score = "na";
	prediction = prob_1;
endfunction
