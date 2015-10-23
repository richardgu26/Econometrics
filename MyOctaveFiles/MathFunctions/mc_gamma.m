function [f] = mc_gamma(z)

test = z > 170;
percentage = mean(test);
if sum(test) > 0
	warning("mc_gamma: %f of the arguments were bounded", percentage);
endif	
z = 170*test + (1-test).*z;

f = gamma(z);

endfunction
