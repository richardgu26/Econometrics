function [logdensity, score, prediction] = Logit(theta, data)
	y = data(:,1);
	x = data(:,2:columns(data));
	prob_1 = 1 ./ (1 + exp(-x*theta));
	logdensity = y .* log(prob_1) + (1 - y) .* log(1 - prob_1);
	score = y - prob_1;
	score = diag(score)*x;
#	score = "na"; # uncomment to use numeric score vector
	prediction = prob_1;
endfunction
