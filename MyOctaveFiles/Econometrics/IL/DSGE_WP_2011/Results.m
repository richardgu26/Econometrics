# Results for DSGE model
# this one is nice because it does all dep vars at once, without loop
# no CIs in this version, could add them
printf("\n\nDoing SBIL estimation using k-nearest neighbors nonparametric regression\n\n");
tic;
load simdata.uniform;
test = sum(isnan(simdata),2);
test = (test ==0);
simdata = simdata(test,:);
printf("%d rows read in\n", rows(simdata));
design = simdata(:,2); # 1 for design, 0 for from param space
insamp = !design;
design = !insamp; # make it boolean
theta = simdata(:,3:9);
Z = simdata(:,10:columns(simdata));
Z = [Z(:,1:columns(Z)-4) Z(:,columns(Z))]; # we computed all st.devs., but it's better to drop some
Z = Z ./ repmat(std(Z),rows(Z),1);         # so as not to put too much weight on sigma, to detriment of
# other parameters
model_params0 = [0.33; 0.99; 0.023; 1.75; 0.95; 0.010448; 10];

# find which are in sample
Z_in = Z(insamp,:);
theta_in = theta(insamp,:);
n = rows(Z_in);
q = columns(Z_in);
#k = floor(0.5*n^0.25);
k = 10;   # this is for the small number of reps example, increase it for a real run
# find the out of sample (at Monte Carlo design point)
Z_out = Z(design,:);
theta_out = theta(design,:);
n = rows(Z_out);
contrib = zeros(n,columns(theta_in));
tolerance = 0;



idx = nearest_neighbors(Z_out, Z_in, k);
for i = 1:n
	keepers = theta_in(idx(i,:),:);
	thetahat = mean(keepers);
	contrib(i,:) = thetahat;
endfor
m = mean(contrib);
s = std(contrib);
e = contrib - repmat(model_params0',rows(contrib),1);
b = mean(e);
e = e.^2;
mse = mean(e);
rmse = sqrt(mse);
priormean = mean(theta_in);
priorsdev = std(theta_in);
priorbias = priormean - model_params0';
priorrmse = sqrt(priorbias.^2 + priorsdev.^2);
mae = mean(abs(e));
rlabels = char("true param", "SBIL mean", "prior mean", "SBIL st.dev.","prior sdev", "SBIL bias", "prior bias","SBIL rmse", "prior rmse");
clabels = char("alpha","beta","delta", "gam", "rho","sig", "eps");
printf("\n\nSBIL estimation results\n");
prettyprint([model_params0'; m; priormean; s; priorsdev; b ; priorbias; rmse; priorrmse], rlabels, clabels)
printf("number of Monte Carlo reps: %d\n", rows(Z_out));
printf("\n\n");

