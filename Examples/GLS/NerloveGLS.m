# Do a GLS correction for the Nerlove Cobb-Douglas model with
# all coefficients common except constant and output, which varies across
# five groups ordered from smallest 29 to largest 29
#
load nerlove.data;

data = data(:,2:6);
data = log(data);
n = rows(data);
y = data(:,1);

# let coefficient on output vary across 5 groups, from smallest to largest
output = data(:,2);
mask = kron(eye(5), ones(29,1));
output = diag(output)*mask;

# assemble the regressor matrix
x = [mask output data(:,3:5)];
k = columns(x);
k = columns(x);

# OLS results, used to consistently estimate the st. dev. of group errors
names = char("constant1", "constant2", "constant3", "constant4", "constant5", "output1" , "output2" , "output3" , "output4" , "output5" ,"labor", "fuel", "capital");
[junk, junk, e] = mc_ols(y, x, names, 1);
e = e .^2;
e = reshape(e, 29,5);
sig = mean(e);

correction = 1 ./ sqrt(sig);
correction = diag(correction)*ones(5,29);
correction = vec(correction');


# the GLS correction: weight by inverse of st. dev.
y = y .* correction;
x = diag(correction)*x;
[b, varb] = mc_ols(y, x, names);

# Plot residuals to check for residual het.
PlotResiduals(y,x,"nerlove_residuals_gls.eps");

# Test HOD1: the input price coefficients sum to 1
R = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1];
r = 1;
printf("\nTesting HOD1\n");
WaldTest(b, varb, R, r);

# this is commented so as not to overwrite the figure used in notes
#PlotResiduals(y, x, 1);
# 