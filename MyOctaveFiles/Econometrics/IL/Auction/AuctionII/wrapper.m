# Does Monte Carlo for Auction model estimated by CUE-II
function contrib = wrapper(args)
	n = args{1}; # sample size
	S = args{2}; # number of simulations for II

	theta0 = [0.5; 0.5];
	# get  Z
	randdraws1 = rand(n,1);
	randdraws2 = rand(n,6);
	randdraws2 = min(randdraws2')'; % min of 6 uniforms, used to get bid
	Z = aux_stat(theta0, randdraws1, randdraws2);

	# get CUE-II
	randdraws1 = rand(n,S);
	randdraws2 = zeros(n,S);
	for i = 1:S
		r = rand(n,6);
		randdraws2(:,i) = min(r')'; % min of 6 uniforms, used to get bid
	endfor

	theta = theta0 + 0.1*randn(2,1);
	% first exact id
	[ii1, obj_value1, convergence1, iters] = samin("ii_obj", {theta, randdraws1, randdraws2, Z, 0}, {[-1; 0],[3;2],3,1,0.5,20000,3,1e-10,1e-4,2,1});
	% now over id
	[ii2, obj_value2, convergence2, iters] = samin("ii_obj", {ii1, randdraws1, randdraws2, Z, 1}, {[-1; 0],[3;2],3,1,0.5,20000,3,1e-10,1e-4,2,1});
	contrib = [convergence1 convergence2 obj_value1 obj_value2  ii1' ii2'];
endfunction
