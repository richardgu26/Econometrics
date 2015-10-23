% Example of Bayesian estimation
% sampling from exponential(theta) distribution
% lognormal prior

% Shows how likelihood, joint, marginal likelihood,
% and posterior are formed, and how posterior can be
% integrated to get posterior mean

% explore different sample sizes, different true thetas


function BayesExample1

    close all;

    n = 30;   % sample size
    theta = 3; % true theta
    y = exprnd(theta*ones(n,1)); % sample from exponential(theta)

    % make plots
    thetas =linspace(0.01,10,1000);
    delta = thetas(:,2)-thetas(:,1);
    thetas = thetas;
    p = prior(thetas);
    post = posterior(y, thetas);
    plot(thetas', [p' post']);

    hold on;

    % get posterior mean by another crude numeric integration
    posteriormean = sum(thetas.*post*delta,2);
    priormean = exp(1.5); % this is the mean of lognormal exp(mu+sig^2/2)
    fprintf('Posterior mean: %f\n', posteriormean);
    h = 1.1*max(post); % height of lines, a little more than height of posterior
    plot([theta; theta], [0; h], 'r');
    plot([priormean; priormean], [0; h], 'g');
    plot([posteriormean; posteriormean], [0; h], 'c');
    legend('prior', 'posterior', 'true theta', 'prior mean', 'posterior mean');
    %print -dpng BayesExampleN10.png;
end

% the prior is lognormal
function p = prior(theta)
    p = lognpdf(theta, 1, 1);
end

% the likelihood function
function dens = likelihood(y, theta)
    theta = repmat(theta, size(y,1), 1);
    y = repmat(y, 1, size(theta,2));
    dens = (1./theta).*exp(-y./theta);
    dens = prod(dens); % independent obsn
end

% joint is prior X likelihood
function dens = joint(y, theta)
    l = likelihood(y, theta);
    p = prior(theta);
    dens = l.*p;
end


% compute marginal likelihood of Y by integrating out theta (crude, only illustrative)
function dens = marginal(y);
    thetas = linspace(0.01,10,1000); % set up a grid
    delta = thetas(2)-thetas(1); % the step size of the grid
    R = size(thetas,2);
    temp = zeros(R,1);
    % evaluate joint on grid
    for r = 1:R
        theta = thetas(:,r);
        temp(r,:) = joint(y, theta);
    end
    % crude numeric integration, sum height X width
    dens = sum(delta*temp);
end

% the posterior, by Bayes' Law
function dens = posterior(y, theta)
    m = marginal(y);
    j = joint(y, theta);
    dens = j ./ m;
end

