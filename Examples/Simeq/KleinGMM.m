% Estimates the Klein consumption equation by
% GMM
1;
function m = KleinMoments(theta, data, momentargs)
	k = momentargs{1}; # use this so that data can hold dep, indeps, and instr
	y = data(:,1);
	x = data(:,2:k+1);
	z = data(:, k+2:columns(data));
	e = y -x*theta;
	m = diag(e)*z;
endfunction	

	
load klein.data;
data = klein;


# construct missing lags, and drop first row that has missing data
profits = data(:,3);
output = data(:,7);
data = [data lag(profits,1) lag(output,1)];
data = data(2:rows(data),:);

n = rows(data);

# define instruments
exogs = [1, 6, 8, 9, 10, 11, 12];
exogs = data(:,exogs);
exogs = [ones(n,1) exogs];

# CONSUMPTION
printf("CONSUMPTION EQUATION\n");
# define variables in consumption equation
y = data(:,2);
profits = data(:,3);
lagprofits = data(:,11);
wp = data(:,4);
wg = data(:,8);
wages = wp + wg;

# regressors in consumption equation
x = [profits lagprofits wages];
x = [ones(n,1) x];

# GMM estimation
k = columns(x);
momentargs = {k};
data = [y x exogs];
moments = 'KleinMoments';
theta = ols(y,x); % OLS start values, inconsistent, but ok (?)
weight = 1;
names = char("Constant", "Profits", "Lagged Profits", "Wages");
# initial consistent estimate: only used to get moment covariance (needed for t-stats) no screen output
[theta, obj_value, convergence] = gmm_estimate(theta, data, weight, moments, momentargs);

# moment covariance assuming no autocorrelation
m = feval(moments, theta, data, momentargs);
momentcov = cov(m);
weight = inv(momentcov);
# estimation results using efficient weight
gmmtitle = "Klein model 1 GMM example";
gmm_results(theta, data, weight, moments, momentargs);


# moment covariance assuming autocorrelation
# note: if there really is autocorrelation,
# then lagged endogs need to be dropped as 
# instruments. This is not done here, as this
# is just meant as an example of use of NW
# covariance estimator
m = feval(moments, theta, data, momentargs);
lags = 2;
momentcov = NeweyWest(m, lags);
weight = inv(momentcov);

# estimation results using efficient weight
gmmtitle = "Klein model 1 GMM example";
gmm_results(theta, data, weight, moments, momentargs);


