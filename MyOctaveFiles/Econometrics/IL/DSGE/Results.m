# Results for DSGE 2 shock model
printf("\n\nDoing SBIL estimation using k-nearest neighbors nonparametric regression\n\n");
tic;

% read in simulations done at design point
FN=fopen('simdata_d', 'r');
ncols = fread(FN, 1, 'single');
simdata = fread(FN, Inf, 'single');
fclose(FN);
simdata = reshape(simdata, ncols, rows(simdata)/ncols);
simdata_d = simdata';
printf("at design: %d rows read in\n", rows(simdata_d));

% read in simulations done from parameter space
FN=fopen('simdata_p', 'r');
ncols = fread(FN, 1, 'single');
simdata = fread(FN, Inf, 'single');
fclose(FN);
simdata = reshape(simdata, ncols, rows(simdata)/ncols);
simdata_p = simdata';
printf("from parameter space: %d rows read in\n", rows(simdata_p));

% stack them together
simdata = [simdata_d; simdata_p];
% find in and out of sample
design = simdata(:,1); # 1 for design, 0 for from param space
insamp = !design;
design = !insamp; # make it boolean

% design
s = simdata(design,:);
printf("%d full size of design sample: \n", rows(s));
test = any(isnan(s'));
test = (test ==0);
s = s(test,:);
test = any(isinf(s'));
test = (test ==0);
s = s(test,:);
test = (s(:,12) != -1000);
s = s(test,:);
%test = max(abs(simdata'))' < 20;
%simdata = simdata(test,:);
printf("%d valid rows of design sample \n", rows(s));

% param space
simdata = simdata(insamp,:);
printf("%d full size of param space sample: \n", rows(simdata));
test = any(isnan(simdata'));
test = (test ==0);
simdata = simdata(test,:);
test = any(isinf(simdata'));
test = (test ==0);
simdata = simdata(test,:);
test = (simdata(:,12) != -1000);
simdata = simdata(test,:);
printf("%d valid rows of design sample \n", rows(simdata));
% drop off extremes to reduce NN computations: this is always clean, so don't report reduction is sample
test = max(abs(simdata'))' < 10;
simdata = simdata(test,:);
% put two together again
simdata = [s; simdata];
% update this
design = simdata(:,1); # 1 for design, 0 for from param space
insamp = !design;
design = !insamp; # make it boolean

% scale aux stats
Z = simdata(:,12:columns(simdata));
[Z, m] = st_norm(Z);
Z = Z + m;
Z_in = Z(insamp,:);
Z_out = Z(design,:);



% parameters
theta = simdata(:,2:11);
theta_in = theta(insamp,:);

% get information for psi
psi = theta_in(1:rows(theta_in),10);
psi_pmean = mean(psi);
psi_psdev = std(psi);

% true parameter values
theta0 = theta(design,:);
theta0 = theta0(1,:);

n = rows(Z_in);
k = floor(1.0*n^0.25);
n = rows(Z_out);

%%%%%%%%%%%%%% select statistics for each parameter and find neighbors %%%%%%%%%%%%%%%%%%%%%%%%%
% alpha beta delta nss
selection = [8:12]; % means of the variables
Z_in1 = Z_in(:,selection);
Z_out1 = Z_out(:,selection);
thetahat_abdn = knn_regression(Z_out1, theta_in, Z_in1, k, 3, 'false');

% gam
selection = [5 14 19:21]; % slope of MPL equation, plus correlations c with other variables
Z_in1 = Z_in(:,selection);
Z_out1 = Z_out(:,selection);
thetahat_gam = knn_regression(Z_out1, theta_in, Z_in1, k, 3, 'false');

% rho1
selection = [2 28:32]; % rho from production equatioon, plus autocovs. of output
Z_in1 = Z_in(:,selection);
Z_out1 = Z_out(:,selection);
thetahat_rho1 = knn_regression(Z_out1, theta_in, Z_in1, k, 3, 'false');

% sig1
selection = [3]; % sig from production equation
Z_in1 = Z_in(:,selection);
Z_out1 = Z_out(:,selection);
thetahat_sig1 = knn_regression(Z_out1, theta_in, Z_in1, k, 3, 'false');

% rho2
selection = [6 33:37]; % rho from MPL equation, plus autocovs. of cons.
Z_in1 = Z_in(:,selection);
Z_out1 = Z_out(:,selection);
thetahat_rho2 = knn_regression(Z_out1, theta_in, Z_in1, k, 3, 'false');

% sig2
selection = [7];  % sig from MPL equation
Z_in1 = Z_in(:,selection);
Z_out1 = Z_out(:,selection);
thetahat_sig2 = knn_regression(Z_out1, theta_in, Z_in1, k, 3, 'false');

thetahat = [thetahat_abdn(:,1) thetahat_abdn(:,2) thetahat_abdn(:,3) thetahat_gam(:,4) ...
		thetahat_rho1(:,5) thetahat_sig1(:,6) thetahat_rho2(:,7) thetahat_sig2(:,8)...
	       	thetahat_abdn(:,9)];

% recover estimate of psi
beta = thetahat(:,2);
delta = thetahat(:,3);
alpha = thetahat(:,1);
nss = thetahat(:,9);
gam = thetahat(:,4);

c1 = ((1./beta + delta - 1)./alpha).^(1./(1-alpha));
kss = nss./c1;
css = kss .* (c1.^(1-alpha) - delta);
c2 = (css).^(-gam./alpha);
psi = (1-alpha).*((c1./c2).^(-alpha));
thetahat = [thetahat psi];
contrib = thetahat;

%%%%%%%%%%%%%%%%% report results %%%%%%%%%%%%%%%%%%555
m = mean(contrib);
s = std(contrib);
e = contrib - repmat(theta0,rows(contrib),1);
b = mean(e);
e = e.^2;
mse = mean(e);
rmse = sqrt(mse);

% characteristics of prior
lb_ub = [
0.20   	0.4		% alpha
0.90    0.999		% beta
0.005    0.1		% delta
0     	3		% gam
0.0    	0.999		% rho1
0.001   0.1		% sigma1
0.0    	0.999		% rho2
0.001	0.1		% sigma2
6/24    12/24		% nss
];
lb = lb_ub(:,1);
ub = lb_ub(:,2);

% prior
priormean = (ub+lb)'/2;
priormean = [priormean psi_pmean];
priorsdev = sqrt(((ub-lb).^2)/12);
priorsdev = [priorsdev' psi_psdev];
priorbias = priormean - theta0;
priorrmse = sqrt(priorbias.^2 + priorsdev.^2);

% results
mae = mean(abs(e));
clabels = char("true", "mean", "pmean", "sdev.","psdev", "bias", "pbias","rmse", "prmse");
rlabels = char(
"alpha",
"beta",
"delta",
"gam",
"rho1",
"sig1",
"rho2",
"sig2",
"nss",
"psi"
);
printf("\n\nSBIL estimation results\n");
prettyprint([theta0; m; priormean; s; priorsdev; b ; priorbias; rmse; priorrmse]', rlabels, clabels);
printf("number of Monte Carlo reps: %d\n", rows(Z_out));
printf("\n\n");

% information for plots
% recover prior of psi, to plot it's prior
beta = theta_in(:,2);
delta = theta_in(:,3);
alpha = theta_in(:,1);
nss = theta_in(:,9);
gam= theta_in(4);
c1 = ((1./beta + delta - 1)./alpha).^(1./(1-alpha));
kss = nss./c1;
css = kss .* (c1.^(1-alpha) - delta);
c2 = (css).^(-gam./alpha);
priorpsi = (1-alpha).*((c1./c2).^(-alpha));
iss = delta.*kss;
yss = kss.^alpha .* nss.^(1-alpha);
css = yss - iss;
psi = css.^(-gam) .* (1-alpha) .* (kss.^alpha) .* (nss.^(-alpha));
% thetahat
thetahat = contrib;
#{
if expdesign==1
	save -ascii thetahat.d1 thetahat;
	save -ascii priorpsi.d1 priorpsi;
elseif expdesign==140
	save -ascii thetahat.d1_40 thetahat;
	save -ascii priorpsi.d1_40 priorpsi;
elseif expdesign==2
	save -ascii thetahat.d2 thetahat;
	save -ascii priorpsi.d2 priorpsi;
elseif expdesign==240
	save -ascii thetahat.d2_40 thetahat;
	save -ascii priorpsi.d2_40 priorpsi;
endif
#}
