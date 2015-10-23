function [f] = mc_exp(z)
lim = 709; # exp overflows at x = 710
test = z > lim;
percentage = mean(test);
if sum(test) > 0
	warning("mc_exp: %f of the arguments were bounded", percentage);
endif	
z = lim*test + (1-test).*z;

f = exp(z);

endfunction
