% Estimates the basic Nerlove Cobb-Douglas model

load nerlove.data;

data = data(:,2:6);
data = log(data);
n = rows(data);
y = data(:,1);
x = data(:,2:5);
x = [ones(n,1), x];


names = char('constant', 'output', 'labor', 'fuel', 'capital');

[b junk junk ess] = mc_ols(y,x,names, 0, 1);
