% This example code shows how CU-II can be used to estimate the auction model in the paper
% "Indirect Likelihood Inference", by Michael Creel and Dennis Kristensen.
% This is a simple version that is for illustrative purposes. The actual code
% used to get the results for the paper implements parallelization, to allow
% work with large numbers of simulations. The intention of this code is to clearly
% show the ideas, without worrying about computational performance.

% usage: execute Auction_CUII_Simple from the Matlab or Octave prompt.

function Auction_CUII_Simple()

	mc_reps = 100;	% number of Monte Carlo replication
	n = 80;		% sample size
	S = 50;		% number of simulations for CU-II moments and cov
			% set low here for speed, increase for accuracy

	% true parameter values
	theta0 = [0.5; 0.5];

	results = zeros(mc_reps, 3);

	for j = 1:mc_reps
		# get  Z
		randdraws1 = rand(n,1);
		randdraws2 = rand(n,6);
		randdraws2 = min(randdraws2')'; % min of 6 uniforms, used to get bid
		Z = aux_stat(theta0, randdraws1, randdraws2);

		# get CUE-II
		randdraws1 = rand(n,S);
		randdraws2 = zeros(n,S);
		for i = 1:S
			r = rand(n,3);
			randdraws2(:,i) = min(r')'; % min of 6 uniforms, used to get bid
		end

		theta = theta0 + 0.1*randn(2,1);
		% over identified estimator
		[ii, obj_value, convergence] = samin("ii_obj", {theta, randdraws1, randdraws2, Z}, {[-1; 0],[3;2],3,1,0.5,20000,3,1e-10,1e-4,2,1});
		results(j,:) = [convergence ii'];
		printf('%d of %d done so far\n', j, mc_reps);
	end

	% drop and non-converged: there are unlikely to be any
	test = (results(:,1) == 1);
	results = results(test,:);

	% display basic results
	thetahat1 = results(:,2);
	thetahat2 = results(:,3);
	fprintf('CU-II estimation results\n');
	fprintf('Sample size: %d\n', n);
	fprintf('Number of simulations for moments and cov: %d\n', S);

	fprintf('posterior mean theta1: %f\n', mean(thetahat1));
	fprintf('posterior mean theta2: %f\n', mean(thetahat2));
	fprintf('RMSE theta1: %f\n', sqrt(mean((thetahat1-theta0(1,:)).^2)));
	fprintf('RMSE theta2: %f\n', sqrt(mean((thetahat2-theta0(2,:)).^2)));
	hist(thetahat1, 50);
	title('thetahat1: true value is 0.5');
	figure;
	hist(thetahat2, 50);
	title('thetahat2: true value is 0.5');
end

% randdraws1: nX1 uniforms, to generate x
% randdraws2: nX1, each entry is the minimum of 6 uniforms, to generate bid

function data = dgp(theta, randdraws1, randdraws2)
	% the model
	theta1 = theta(1,:);
	theta2 = theta(2,:);
	n = size(randdraws1,1);
	N = 6;
	% quality of good
	x = randdraws1;
	% valuations drawn from exponetial mean phi
	phi = exp(theta1 + theta2*x);
	% highest valuation
	v = -log(randdraws2).*phi;
	% get winning bid
	z = v./phi;
	D = exp(-5*z).*(60*exp(5*z) + 300*phi .* exp(4*z) - 300*phi .* exp(3*z) ...
	+ 200*phi .* exp(2*z) - 75*phi .* exp(z) + 12*phi)/60 - 137*phi/60;
	b = v - D ./ ((1 - exp(-v./phi)).^(N-1));
	b = b.*(b>0);
	data = [b x];
end

% returns a stack of aux stats, each one in a row
% used for CUE-II (S cols. in randdraws)
function Z = aux_stat(theta, randdraws1, randdraws2)
	S = size(randdraws1,2);
	Z = zeros(S, 6);
	for s = 1:S
		data = dgp(theta, randdraws1(:,s), randdraws2(:,s));
		b = data(:,1);
		% bound bid for numeric stability
		b = (b>0.01).*b + (b<0.01)*0.01;
		y = log(b); % use log of bid
		x = [ones(size(y,1),1) data(:,2)];
		bhat = x\y;
		e = y - x*bhat;
		sig = log(e'*e/(size(y,1)-2));
		m1 = mean(log(b));
		m2 = std(log(b));
		m3 = mean((log(b)-m1).^3);
		Z(s,:) = [bhat; sig; m1; m2; m3]';
	end
end


function Z = auxstat(data)
	x = data(:,1);
	b = data(:,2);
	% top bound bid to control large outliers
	b = b.*(b>0);
	% the aux stat
	b = 0.01.*(b<0.01) + b.*(b>0.01);
	y = log(b);
	z = [ones(size(x,1),1) x];
	bhat = z\y;
	e = y-z*bhat;
	s = log(e'*e/(size(y,1)-2));
	m1 = mean(log(b));
	m2 = std(log(b));
	m3 = mean((log(b)-m1).^3);
	Z = [bhat' s m1 m2 m3];
end

% CU-GMM objective function for indirect inference
function [obj_value score] = ii_obj(theta, randdraws1, randdraws2, Z)

	% simulated aux stats
	Zs = aux_stat(theta, randdraws1, randdraws2);

	% moment conditions
	Zbar = mean(Zs);
	m = Z - Zbar; % row vector

	% weight matrix
	e = Zs - Z;
	V = e'*e/size(e,1);
	W = inv(V);

	obj_value = m*W*m';
	score = '';
end


