close all;
load thetahat.d1;
thetahat = thetahat;
load priorpsi.d1;
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

lb = lb_param_ub(:,1);
ub = lb_param_ub(:,3);
lb = [lb; min(priorpsi)];
ub = [ub; max(priorpsi)];
theta0 = lb_param_ub(:,2);
% recover estimate of psi
beta = theta0(2,:);
delta = theta0(3,:);
alpha = theta0(1,:);
nss = theta0(9,:);
gam = theta0(4,:);
c1 = ((1/beta + delta - 1)/alpha)^(1/(1-alpha));
kss = nss/c1;
css = kss * (c1^(1-alpha) - delta);
c2 = (css)^(-gam/alpha);
psi = (1-alpha)*((c1/c2)^(-alpha));
theta0 = [theta0; psi];

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

for i = 1:columns(thetahat)
	if i > 1 figure; end
	tmp = thetahat(:,i);
	x = 1:100;
	x = (ub(i)-lb(i))/100*x + lb(i);
	x = x(:);
	d = kernel_density(x, thetahat(:,i),2, '__kernel_epanechnikov', true );
	plot(x, d, "b", 'LineWidth', 4);
	hold on;
	if i==columns(thetahat)
		h = kernel_density(x, priorpsi(1:100000,:),2, '__kernel_epanechnikov', true );
	else
	   	h = ones(size(x))/(ub(i)-lb(i));
	endif
	plot(x, h, "b.", 'LineWidth', 4);
	h = max(d);
	plot([theta0(i);theta0(i)], [0;h], 'b', 'LineWidth', 4);
	of = rlabels(i,:);
	print(of, '-dpng');
end
