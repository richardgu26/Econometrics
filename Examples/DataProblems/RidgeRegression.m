% ridge regression as example of Bayesian estimation

n = 100;
a = 0.95; % controls collinearity: close to 0 indep close to 1 highly collinear

for i = 1:1000

	% generate regressors that are collinear
	x1 = rand(n,1);
	x2 = a*x1 + (1-a)*randn(n,1);
	x = [ones(n,1) x1 x2];

	% true params and dep var
	b = ones(3,1);
	e = randn(n,1);
	y = x*b+e;

	% OLS
	bols = ols(y,x);

	% ridge regression
	k = 1;
	yy = [y; k*zeros(3,1)];
	xx = [x; k*eye(3)];
	bridge = ols(yy,xx);
	
	% store results
	bs(i,:) = [bols' bridge'];

end

% histograms to compare 
bs = bs - [b' b']; % take away true, so centered around 0
hist(bs(:,2),30);
%print('ridge_example_ols.png','-dpng');
figure
hist(bs(:,5),30);
%print('ridge_example_ridge.png','-dpng');

