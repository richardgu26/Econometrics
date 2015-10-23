function [f] = mc_lgamma(z)
lim = 1e305; # lgamma overflows past this
test = z > lim;
percentage = mean(test);
if sum(test) > 0
	warning("mc_lgamma: %f of the arguments were top bounded", percentage);
endif	
z = lim*test + (1-test).*z;

# lim = 1e-5; # must be positive
# test = z < lim;
# percentage = mean(test);
# if sum(test) > 0
# 	warning("mc_lgamma: %f of the arguments were bottom bounded", percentage);
# endif	
# z = lim*test + (1-test).*z;


f = lgamma(z);

endfunction
