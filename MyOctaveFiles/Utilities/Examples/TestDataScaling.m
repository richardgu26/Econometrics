# Verifies the ScaleData and UnscaleParameters functions

# With or without constant - you choose
x = [ones(10,1) randn(10,2)];
#x = randn(10,3);

# Now check that it works
[xs,scalecoefs] = ScaleData(x);
thetas = randn(3,1);
vs = randn(3,3);
[theta, vtheta] = UnscaleParameters(thetas, vs, scalecoefs); 

printf("Check that there's no diff. between scaled and orig.\n");
disp(xs*thetas - x*theta);

# no singularities introduced, we hope?
printf("Check nonsingular - determinant is %d\n", det(xs'*xs));
