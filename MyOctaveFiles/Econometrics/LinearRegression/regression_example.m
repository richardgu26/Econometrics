# this is to test the LinRegLib translation

n = 100;
x = [ones(n,1) rand(n,2)];
theta = ones(3,1);
y = x*theta + 2*randn(n,1);
names = char("a","b","c");

R = [1, -1, 0];
r = 0;

[b, v, e] = mc_ols(y,x,names);
[b, v, e] = mc_olsr(y, x, R, r, names);

PlotResiduals(y,x);
pause(5);
