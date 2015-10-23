% Simple example code to replicate results in
% "Indirect Likelihood Inference" (2013), by
% Michael Creel and Dennis Kristensen.
% 
design = 1;	% there are 2 designs: select one
		% you can also change parameters
		% as long as the values are inside
		% the support of the prior
S = 100;	% replications from param space to generate
		% this should be much larger to get 
		% reliable results
mc_reps = 20;  	% number of Monte Carlo replications

if design ==1
	% design 1
	lb_param_ub = [
	0.20   	0.33   	0.4		% alpha
	0.90    0.99   	0.999		% beta
	0.005    0.025   0.1		% delta
	0     	2      	3		% gam
	0.0    	0.9   	0.999		% rho1
	0.001   0.01 	0.1		% sigma1
	0.0    	0.7     0.999		% rho2
	0.001	0.005  	0.1		% sigma2
	6/24    8/24	12/24		% nss
	];
else design==2
	% design 2
	lb_param_ub = [
	0.20   	0.25   	0.4		% alpha
	0.90    0.97   	0.999		% beta
	0.005    0.04   0.1		% delta
	0     	1      	3		% gam
	0.0    	0.8   	0.999		% rho1
	0.001   0.02 	0.1		% sigma1
	0.0    	0.6     0.999		% rho2
	0.001	0.01  	0.1		% sigma2
	6/24    8/24	12/24		% nss
	];
end

model_params0 = lb_param_ub(:,2);
lb = lb_param_ub(:,1);
ub = lb_param_ub(:,3);

first_time = true;
for s = 1:mc_reps + S
	if ~first_time
		load state1;
		rand('state', s);
	end
	if s <= mc_reps % use Monte Carlo design
		outsamp = 1;
		model_params = model_params0;
	else 		% sample from parameter space
		outsamp = 0;
		model_params = rand(size(ub,1),1).*(ub-lb) + lb;
	end
	RNGstate = rand('state');
	% break into pieces
	alpha = model_params(1,:);
	beta  = model_params(2,:);
	delta = model_params(3,:);
	gam   = model_params(4,:);
	rho1   = model_params(5,:);
	sigma1 = model_params(6,:);
	rho2   = model_params(7,:);
	sigma2 = model_params(8,:);
	nss   = model_params(9,:);

	% the psi implied by other parameters
	c1 = ((1/beta + delta - 1)/alpha)^(1/(1-alpha));
	kss = nss/c1;
	css = kss * (c1^(1-alpha) - delta);
	c2 = (css)^(-gam/alpha);
	psi = (1-alpha)*((c1/c2)^(-alpha));

	% save parameters and random number state
    if first_time
        save parameterfile1  alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
        save state1 RNGstate;
	dynare CK_Example noclearall;
        first_time = false;
    else
        set_param_value('alppha',alpha)
        set_param_value('betta',beta)
        set_param_value('delta',delta)
        set_param_value('gam',gam)
        set_param_value('rho1',rho1)
        set_param_value('sigma1',sigma1)
        set_param_value('rho2',rho2)
        set_param_value('sigma2',sigma2)
        set_param_value('nss',nss)
        ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
        set_dynare_seed(ss);
        info = stoch_simul(var_list_);
    end
	% get a simulation of length 40 or 160 (10 or 40 years quarterly), and compute aux. statistic
	data = [y c n MPK MPL];
	data = data(101:260,:);   % 40 years
	% data = data(101:140,:); % 10 years
	Z = aux_stat(data);
	contrib = [outsamp model_params' psi Z'];
	if (s==1)
		simdata = zeros(mc_reps+S, columns(contrib));
	end
	simdata(s,:) = contrib;
	if (~mod(s,10))
		fprintf('\n\n##################################################################################\n', s, mc_reps+S);
		fprintf('### Estimation of DSGE model by SBIL: %d simulations of %d total done so far ###\n', s, mc_reps+S);
		fprintf('##################################################################################\n', s, mc_reps+S);
	end	
end

% Now the results
clc;
fprintf('Estimation of DSGE model by SBIL: simulations done, now doing nearest neighbors regression\n');

% find in and out of sample
design = simdata(:,1); % 1 for design, 0 for from param space
insamp = ~design;
design = ~insamp; % make it boolean

% design
s = simdata(design,:);
fprintf('%d full size of design sample: \n', size(s,1));
test = any(isnan(s'));
test = (test ==0);
s = s(test,:);
test = any(isinf(s'));
test = (test ==0);
s = s(test,:);
test = (s(:,12) ~= -1000);  % this is the bad value code
s = s(test,:);
fprintf('%d valid rows of design sample (drops are a problem)\n', size(s,1));

% param space
simdata = simdata(insamp,:);
fprintf('%d full size of param space sample: \n', size(simdata,1));
test = any(isnan(simdata'));
test = (test ==0);
simdata = simdata(test,:);
test = any(isinf(simdata'));
test = (test ==0);
simdata = simdata(test,:);
test = (simdata(:,12) ~= -1000); % this is the bad value code
simdata = simdata(test,:);
fprintf('%d valid rows of param space sample (drops are not a problem)\n', size(simdata,1));
% drop off extremes to reduce NN computations: this is always clean, so don't report reduction is sample
test = max(abs(simdata'))' < 10;
simdata = simdata(test,:);
% put two together again
simdata = [s; simdata];
% update this
design = simdata(:,1); % 1 for design, 0 for from param space
insamp = ~design;
design = ~insamp; % make it boolean

% scale aux stats
Z = simdata(:,12:size(simdata,2));
Z = Z*inv(diag(std(Z)));
Z_in = Z(insamp,:);
Z_out = Z(design,:);

% parameters
theta = simdata(:,2:11);
theta_in = theta(insamp,:);

% get information for psi
psi = theta_in(1:size(theta_in,1),10);
psi_pmean = mean(psi);
psi_psdev = std(psi);

% true parameter values
theta0 = theta(design,:);
theta0 = theta0(1,:);

n = size(Z_in,1);
k = floor(1.0*n^0.25);
n = size(Z_out,1);

%%%%%%%%%%%%%% select statistics for each parameter and find neighbors %%%%%%%%%%%%%%%%%%%%%%%%%
% alpha beta delta nss
selection = [8:12]; % means of the variables
Z_in1 = Z_in(:,selection);
Z_out1 = Z_out(:,selection);
idx = knnsearch(Z_out1, Z_in1, k);
thetahat_abdn = zeros(size(Z_out1,1),size(theta_in,2));
for i = 1:size(Z_out1,1)
	thetahat_abdn(i,:) = mean(theta_in(idx(i,:),:));
end

% gam
selection = [5 14 19:21]; % slope of MPL equation, plus correlations c with other variables
Z_in1 = Z_in(:,selection);
Z_out1 = Z_out(:,selection);
idx = knnsearch(Z_out1, Z_in1, k);
thetahat_gam = zeros(size(Z_out1,1),size(theta_in,2));
for i = 1:size(Z_out1,1)
	thetahat_gam(i,:) = mean(theta_in(idx(i,:),:));
end

% rho1
selection = [2 28:32]; % rho from production equatioon, plus autocovs. of output
Z_in1 = Z_in(:,selection);
Z_out1 = Z_out(:,selection);
idx = knnsearch(Z_out1, Z_in1, k);
thetahat_rho1 = zeros(size(Z_out1,1),size(theta_in,2));
for i = 1:size(Z_out1,1)
	thetahat_rho1(i,:) = mean(theta_in(idx(i,:),:));
end


% sig1
selection = [3]; % sig from production equation
Z_in1 = Z_in(:,selection);
Z_out1 = Z_out(:,selection);
idx = knnsearch(Z_out1, Z_in1, k);
thetahat_sig1 = zeros(size(Z_out1,1),size(theta_in,2));
for i = 1:size(Z_out1,1)
	thetahat_sig1(i,:) = mean(theta_in(idx(i,:),:));
end

% rho2
selection = [6 33:37]; % rho from MPL equation, plus autocovs. of cons.
Z_in1 = Z_in(:,selection);
Z_out1 = Z_out(:,selection);
idx = knnsearch(Z_out1, Z_in1, k);
thetahat_rho2 = zeros(size(Z_out1,1),size(theta_in,2));
for i = 1:size(Z_out1,1)
	thetahat_rho2(i,:) = mean(theta_in(idx(i,:),:));
end

% sig2
selection = [7];  % sig from MPL equation
Z_in1 = Z_in(:,selection);
Z_out1 = Z_out(:,selection);
idx = knnsearch(Z_out1, Z_in1, k);
thetahat_sig2 = zeros(size(Z_out1,1),size(theta_in,2));
for i = 1:size(Z_out1,1)
	thetahat_sig2(i,:) = mean(theta_in(idx(i,:),:));
end

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
clabels = char('true', 'mean', 'pmean', 'sdev.','psdev', 'bias', 'pbias','rmse', 'prmse');
rlabels = char(...
'alpha',...
'beta',...
'delta',...
'gam',...
'rho1',...
'sig1',...
'rho2',...
'sig2',...
'nss',...
'psi'...
);
fprintf('\n\nSBIL estimation results - pmean, etc. correspond to prior mean\n');
prettyprint([theta0; m; priormean; s; priorsdev; b ; priorbias; rmse; priorrmse]', rlabels, clabels);
fprintf('number of Monte Carlo reps: %d\n', size(Z_out,1));
fprintf('\n\n');




