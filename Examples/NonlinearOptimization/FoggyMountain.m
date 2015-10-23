


# One BFGS run with poor starting values
theta = [8;-8];
control = {100,2,1,1};
[theta, obj_value, iterations, convergence] = bfgsmin("FoggyMountainObj", {theta}, control);
printf("The result with poor start values\n");
theta'


# Now try simulated annealing
# SA controls
ub = 20*ones(rows(theta),1);
lb = -ub;
nt = 20;
ns = 5;
rt = 0.5; # careful - this is too low for many problems
maxevals = 1e10;
neps = 5;
functol = 1e-10;
paramtol = 1e-3;
verbosity = 1;
minarg = 1;
control = { lb, ub, nt, ns, rt, maxevals, neps, functol, paramtol, verbosity, 1};
# do sa
theta = [8;-8]; # bad start values
[theta, obj_value, convergence] = samin("FoggyMountainObj", {theta}, control);


# Now try 20 random starting values and use a short BFGS on each to get good start values
startvals = 20*rand(2,20) - 10; # 20 start values, random on (-10,10)x(-10,10)
iters = 10;  # number of trial iterations per start value
theta = battery("FoggyMountainObj", {theta}, 1, startvals, iters);
printf("Now try a battery of random start values and \na short BFGS on each, then iterate to convergence\n");
printf("The result using 20 randoms start values\n");
theta'
printf("The true maximizer is near (0.037,0)\n");
