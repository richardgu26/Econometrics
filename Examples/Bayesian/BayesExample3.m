% Bayesian estimation of Nerlove model
% uniform prior
% uniform random walk proposal
% allows imposition of homogeneity of degree 1
% this is the simple model which is obviously 
% misspecified. Should add quadratic term in lnQ

function BayesExample3
	close all;
	load nerlove.data;
	data = data(:,2:end);
	data = log(data);
	n = rows(data);
	y = data(:,1);
	x = [ones(n,1) data(:,2:5)];

	S = 5000;
	theta = [-5; 0.7; 0.5; .4; 0.4];  % initial parameter values (use OLS with HOD1 to get these)
	tuning = [0.3; 0.05; 0.05; 0.05; 0.1];  % tunes the acceptance rate

	thetas = zeros(S, 5);
	accepts = zeros(S, 1);

	for s = 1:S
		thetastar =  proposal(theta, tuning);
		% in following 2 lines, beta_k is 1 - beta_l - beta_f
		theta2 = [theta(1:4,:); 1 - theta(3,:) - theta(4,:); theta(5,:)];
		thetastar2 = [thetastar(1:4,:); 1 - thetastar(3,:) - thetastar(4,:); thetastar(5,:)];
		% log likelihoods
		logn = NormalLF(thetastar2, y, x);
		logd = NormalLF(theta2, y, x);
		num = prior(thetastar)*proposal_density(theta, thetastar, tuning);
		den = prior(theta)*proposal_density(thetastar, theta, tuning);
		% note how ratio of likelihoods is computed as exp( ) of difference of logs
		% this avoids a divide by zero error, as exp(logd) is very small
		crit = exp(logn-logd)*num/den;
		accept = crit > rand(1,1);
		if accept
			theta = thetastar;
		end
		thetas(s,:) = theta;
		accepts(s,:) = accept;
	end

    thetas = thetas(1001:end,:); # drop burnin
	thetas(:,1) = thetas(:,1)/10; % scale constant to make a nice plot
	plot(thetas);
	legend('const10', 'bq', 'bL','bF','sig', 'location', 'southwest');
	thetas = thetas(1001:end,:);
	pm = mean(thetas);
	fprintf('posterior mean %f\n', pm);
	fprintf('acceptance rate %f\n', mean(accepts));
	print -dsvg BayesNerlove1.svg;
	
	% posterior for bQ
	figure;
	bq = thetas(:,2);
	x = linspace(min(bq), max(bq));
	x = x';
	y = kernel_density(x, bq);
	plot(x,y);
	print -dsvg BayesNerlove2.svg;
end


% the prior is a product density over the 5 free parameters
function p = prior(theta)
	pconst = unifpdf(theta(1,:),-10,0); % uniform for constant
	pq = unifpdf(theta(2,:),0,1); % uniform(0,1) for betaq (rules out decreasing RTS)
	pl = unifpdf(theta(3,:), 0, 0.5);
	pf = unifpdf(theta(4,:),0, 0.5);
	psig = unifpdf(theta(5,:), 0, 1); %sig
	p = pconst*pq*pl*pf*psig;
end

% the proposal density: random walk uniform
function f = proposal_density(thetastar, theta, tuning)
	f = 1; % because the proposal is uniform, we don't need to compute this, it's constant
end
% the proposal: random walk uniform
function thetastar = proposal(theta, tuning)
	thetastar = theta;
	whichp = unidrnd(5);
	%thetastar(whichp,:) = theta(whichp,:) + tuning(whichp,:)*rand - tuning(whichp,:)/2;
	thetastar = theta + tuning.*rand(size(tuning)) - tuning/2;
end


