# Results for auction model
close all;
printf("\n\nDoing SBIL estimation using k-nearest neighbors nonparametric regression\n\n");
tic;

load simdata.320;
printf("%d rows read in\n", rows(simdata));
test = any(isnan(simdata'));
test = (test ==0);
simdata = simdata(test,:);
test = any(isinf(simdata'));
test = (test ==0);
simdata = simdata(test,:);
printf("%d valid rows \n", rows(simdata));

design = simdata(:,1); # 1 for design, 0 for from param space
insamp = !design;
design = !insamp; # make it boolean

theta = simdata(:,2:3);
theta_in = theta(insamp,:);
theta0 = theta(design,:);
theta0 = theta0(1,:);
Z = simdata(:,[4:9]);
[Z, m, s] = st_norm(Z);
Z = Z + m;
Z_in = Z(insamp,:);
Z_out = Z(design,:);

n = rows(Z_in);
k = floor(1.5*n^0.25);
n = rows(Z_out);
% all params
selection = 1:columns(Z);
Z_in1 = Z_in(:,selection);
Z_out1 = Z_out(:,selection);
contrib = knn_regression(Z_out1, theta_in, Z_in1, k, 1, 'false');
m = mean(contrib);
s = std(contrib);
e = contrib - repmat(theta0,rows(contrib),1);
b = mean(e);
e = e.^2;
mse = mean(e);
rmse = sqrt(mse);
lb_param_ub = [
-1     0.5   3	% beta0
0.0    0.5   2	% beta1
];
lb = lb_param_ub(:,1);
ub = lb_param_ub(:,3);
priormean = (ub+lb)'/2;
priorsdev = sqrt(((ub-lb).^2)/12);
priorsdev = priorsdev';
priorbias = priormean - theta0;
priorrmse = sqrt(priorbias.^2 + priorsdev.^2);
mae = mean(abs(e));
clabels = char("true", "mean", "pmean", "sdev.","psdev", "bias", "pbias","rmse", "prmse");
rlabels = char(
"beta1",
"beta2"
);
printf("\n\nSBIL estimation results\n");
prettyprint([theta0; m; priormean; s; priorsdev; b ; priorbias; rmse; priorrmse]', rlabels, clabels);
printf("number of Monte Carlo reps: %d\n", rows(Z_out));
printf("\n\n");
