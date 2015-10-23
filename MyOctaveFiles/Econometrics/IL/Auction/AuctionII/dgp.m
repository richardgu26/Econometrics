% randdraws1: nX1 uniforms, to generate x
% randdraws2: nX1, each entry is the minimum of 6 uniforms, to generate bid

function data = dgp(theta, randdraws1, randdraws2)
	% the model
	theta1 = theta(1,:);
	theta2 = theta(2,:);
	n = rows(randdraws1);
	N = 6;
	# quality of good
	x = randdraws1;
	# valuations drawn from exponetial mean phi
	phi = exp(theta1 + theta2*x);
	# highest valuation
	v = -log(randdraws2).*phi;
	# get winning bid
	z = v./phi;
	D = exp(-5*z).*(60*exp(5*z) + 300*phi .* exp(4*z) - 300*phi .* exp(3*z) ...
	+ 200*phi .* exp(2*z) - 75*phi .* exp(z) + 12*phi)/60 - 137*phi/60;
	b = v - D ./ ((1 - exp(-v./phi)).^(N-1));
	b = b.*(b>0);
	data = [b x];
endfunction

