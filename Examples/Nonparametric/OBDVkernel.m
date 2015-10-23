# This does kernel regression estimation for OBDV data

# The MEPS data

# Define dep and expl vbls
#
# 	The dep. vbls, in corresponding column
# 	Office based doctor visits	1
# 	Outpatient doctor visits	2
# 	Emergency room visits		3
# 	Inpatient visits		4
# 	Dental visits			5
# 	Prescriptions			6
#

1;


which_dep = 1;
if (which_dep == 1) printf("\nOBDV\n"); endif
if (which_dep == 2) printf("\nOPV");endif
if (which_dep == 3) printf("\nIPV");endif
if (which_dep == 4) printf("\nERV");endif
if (which_dep == 5) printf("\nDV");endif
if (which_dep == 6) printf("\nPRESCR");endif

load meps1996.data;

# data = data(1:400,:); # cut down on size to speed this up (stupid idea in real life!)
y = data(:,which_dep);
x = data(:,7:12);
n = rows(x);

# evaluation points for plot: vbls at means, except AGE goes from 18 - 65
evalp = mean(x);
evalp = repmat(evalp, 48,1);
evalp(:,4) = (18:65)';

ww = "";
fit = kernel_regression(evalp, y, x, ww);
plot(evalp(:,4),fit);
title("Kernel fit, dep. var. versus AGE");
xlabel("Age");
grid("on");
legend("off");
axis([18 65]);
%print("kernelfit.png", "-dpng");
