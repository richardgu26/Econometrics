% Example of Bayesian estimation
% sampling from exponential(theta) distribution
% lognormal prior

% Shows how MCMC can be used to get posterior mean

% explore different sample sizes, different true thetas


function BayesExample2
	close all;
	n = 100;   % sample size
	truetheta = 3; % true theta
	y = exprnd(ones(n,1)*truetheta); % sample from exponential(theta)

	S = 1000;
	theta = 2; % start value for theta
	tuning = 0.25; % tunes the acceptance rate. Lowering increases acceptance 
                % try to get acceptance rate to be about 0.4 or so
	thetas = zeros(S,1);
	accepts = zeros(S,1);
	for s = 1:S
		thetastar =  proposal(theta, tuning);
		num = likelihood(y, thetastar)*prior(thetastar)*proposal_density(theta, thetastar, tuning);
		den = likelihood(y, theta)*prior(theta)*proposal_density(thetastar, theta, tuning);
		crit = num / den;
		accept = crit > rand(1,1);
		if accept
			theta = thetastar;
		end
		thetas(s,:) = theta;
		accepts(s,:) = accept;
	end
	plot(thetas);
	hold on;
	pm = mean(thetas(500:S,:));
	plot([0; S], [pm; pm], 'g');
	plot([0; S], [truetheta; truetheta], 'r');
	legend('chain', 'posterior mean', 'true theta');
	fprintf('posterior mean %f\n', pm);
	fprintf('acceptance rate %f\n', mean(accepts));
	%print -dpng BayesExample2.png;
end


% the prior is lognormal
function p = prior(theta)
	p = lognpdf(theta, 1, 1);
end

% the likelihood function
function dens = likelihood(y, theta)
	dens = (1./theta).*exp(-y./theta);
	dens = prod(dens); % independent obsn
end

% the proposal density: random walk lognormal
function f = proposal_density(thetastar, theta, tuning)
	f = lognpdf(thetastar, log(theta), tuning);
end
% the proposal: random walk lognormal
function thetastar = proposal(theta, tuning)
	thetastar = lognrnd(log(theta), tuning);
end


