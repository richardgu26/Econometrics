# Estimates the Klein consumption equation by
# OLS, plots residuals, and does
# Breusch-Godfrey test

load klein.data;

# construct missing lags, and drop first row that has missing data
profits = data(:,3);
output = data(:,7);
data = [data lag(profits,1) lag(output,1)];
data = data(2:rows(data),:);

n = rows(data);

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


# OLS residuals to estimate rho
[junk, junk, e] = ols(y, x);
elag = lag(e,1);
e = e(2:rows(e),:);
elag = elag(2:rows(elag),:);

# estimated rho
[rho, junk, junk] = ols(e, elag);
rho
# THE DATA TRANSFORMATION
# special treatment 1st obsn.
y1 = y(1,:)*sqrt(1-rho^2);
x1 = x(1,:)*sqrt(1-rho^2);
# the other obns.
y = y - rho*lag(y,1);
x = x - rho*lag(x,1);
# now copy in the special treatment 1st obsn.
y(1,:) = y1;
x(1,:) = x1;

names = char("Constant", "Profits", "Lagged Profits", "Wages");
[junk, junk, e] = mc_ols(y, x, names);
PlotResiduals(y,x,"KleinResiduals2.eps");

# check for remaining autocorrelation
BreuschGodfreyTest(e, 1, x);

