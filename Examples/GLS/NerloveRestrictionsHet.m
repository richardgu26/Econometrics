# Tests CRTS and HOD1 for the Nerlove Cobb-Douglas model
# The constant and output coefficient varies across 5 firm size groups

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

[b, varb] = mc_ols(y,x,"", 1);

# First HOD1: the input price coefficients sum to 1
R = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1];
r = 1;
printf("\nTesting HOD1\n");
WaldTest(b, varb, R, r);

# Now CRTS: need all output coefficients to be 1
R = [zeros(5,5) eye(5) zeros(5,3)];
r = ones(5,1);
printf("\nTesting CRTS\n");
WaldTest(b, varb, R, r);
