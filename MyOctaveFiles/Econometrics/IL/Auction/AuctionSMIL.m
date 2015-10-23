1;

function d = smil_obj(theta, Z_out, Ztheta_in, k)
	Ztheta_out = [theta' Z_out];
	[idx, d] = nearest_neighbors(Ztheta_out, Ztheta_in, k);
	d = d(k);
endfunction

load simdata.80;
%simdata = simdata(1:50000,:);

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

Ztheta_in = [theta_in Z_in];


% these may look tight, but they are well outside
% the limits of the SBIL estimator. Because SMIL
% is theoretically close to SBIL, no need to go
% further out
lb = [.1; 0.1];
ub = [0.9; 1.1];
% start value
theta = rand(2,1).*(ub-lb) + lb;
nt = 3;
ns = 1;
rt = 0.5; # careful - this is too low for many problems
maxevals = 2000;
neps = 5;
functol = 1e-10;
paramtol = 1e-3;
verbosity = 2; # only final results. Inc
minarg = 1;
control = { lb, ub, nt, ns, rt, maxevals, neps, functol, paramtol, verbosity, 1};

n = rows(Z_in);
k = floor(1.5*n^0.25);
contribs = zeros(rows(Z_out),2);
for i = 1:rows(Z_out)
	Z_out_i = Z_out(i,:);
	% get the pool of possible Zs for this Zn
	[idx, d] = nearest_neighbors(Z_out_i, Z_in, 100*k);
	zz = Ztheta_in(idx(:),:);
	idx = idx(k);
	[thetahat, obj_value, convergence] = samin('smil_obj', {theta, Z_out_i, zz, k}, control);
	contribs(i,:) = [thetahat'];

	theta0 = [0.5 0.5];
	bias = mean(contribs(1:i,:) - theta0);
	rmse = sqrt(mean((contribs(1:i,:) - theta0).^2));
	bias
	rmse
endfor
save AuctionSMIL.out contribs;
