# Estimates the Klein consumption equation by
# 2SLS, plots residuals, and does
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


# OLS estimation
names = char("Constant", "Profits", "Lagged Profits", "Wages");
[junk, junk, e] = mc_ols(y, x, names);
PlotResiduals(y,x,"KleinResiduals.eps");
pause(3);
BreuschGodfreyTest(e, 1, x);
 
