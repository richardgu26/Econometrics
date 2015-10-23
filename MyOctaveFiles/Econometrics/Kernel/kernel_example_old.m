# Copyright (C) 2007 Michael Creel <michael.creel@uab.es>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; If not, see <http://www.gnu.org/licenses/>.

# kernel_example: examples of how to use kernel density and regression functions
# requires the optim and plot packages from Octave Forge
#
# usage: kernel_example;

# sample size (default n = 500) - should get better fit (on average)
# as this increases, supposing that you allow the optimal window width
# to be found, by uncommenting the relevant lines
n = 500;

# set this to greater than 0 to try parallel computations (requires MPITB)
compute_nodes = 1;
nodes = compute_nodes + 1; # count master node

close all;
hold off;

############################################################
# kernel regression example
x = (0:n-1)/n;
x = x';
# generate dependent variable
trueline = sin(2*pi*x) + cos(10*2*pi*x) + .5*randn(n,1);
sig = 0.5;
y = trueline + sig*randn(n,1);
tic;
fit = kernel_regression(x, y, x); # use the default bandwidth
t1 = toc;
printf("\n");
printf("########################################################################\n");
printf("time for kernel regression example using %d data points and %d compute nodes: %f\n", n, nodes, t1);
plot(x, fit, ";fit;", x, trueline,";true;");
grid("on");
title("Example 1: Kernel regression fit");

############################################################
# kernel density example: univariate - fit to Chi^2(3) data
data = sumsq(randn(n,3),2);
# evaluation point are on a grid for plotting
stepsize = 0.1;
grid_x = (-5:stepsize:15)';
bandwidth = 0.55;
# get optimal bandwidth (time consuming, uncomment if you want to try it)
# bandwidth = kernel_optimal_bandwidth(data);
# get the fitted density and do a plot
tic;
dens = kernel_density(grid_x, data, bandwidth, "__kernel_normal", false, false, compute_nodes);
t1 = toc;
printf("\n");
printf("########################################################################\n");
printf("time for univariate kernel density example using %d data points and %d compute nodes: %f\n", n, nodes, t1);
printf("A rough integration under the fitted univariate density is %f\n", sum(dens)*stepsize);
figure();
plot(grid_x, dens, ";fitted density;", grid_x, chi2pdf(grid_x,3), ";true density;");
title("Example 2: Kernel density fit: Univariate Chi^2(3) data");

############################################################
# kernel density example: bivariate
# X ~ N(0,1)
# Y ~ Chi squared(3)
# X, Y are dependent
d = randn(n,3);
data = [d(:,1) sumsq(d,2)];
# evaluation points are on a grid for plotting
stepsize = 0.2;
a = (-5:stepsize:5)'; # for the N(0,1)
b = (-1:stepsize:9)';  # for the Chi squared(3)
gridsize = rows(a);
[grid_x, grid_y] = meshgrid(a, b);
eval_points = [vec(grid_x) vec(grid_y)];
bandwidth = 0.5;
# get optimal bandwidth (time consuming, uncomment if you want to try it)
# bandwidth = kernel_optimal_bandwidth(data);
# get the fitted density and do a plot
tic;
dens = kernel_density(eval_points, data, bandwidth, "__kernel_epanechnikov", true, false, compute_nodes);
t1 = toc;
printf("\n");
printf("########################################################################\n");
printf("time for multivariate kernel density example using %d data points and %d compute nodes: %f\n", n, nodes, t1);
dens = reshape(dens, gridsize, gridsize);
printf("A rough integration under the fitted bivariate density is %f\n", sum(sum(dens))*stepsize^2);
figure();
surf(grid_x, grid_y, dens);
title("Example 3: Kernel density fit: dependent bivariate data");
xlabel("true marginal density is N(0,1)");
ylabel("true marginal density is Chi^2(3)");


# more extensive test of parallel
if compute_nodes > 0 # only try this if parallel is available
	ns =[1000; 2000; 4000];
	printf("\n");
	printf("########################################################################\n");
	printf("kernel regression example with several sample sizes serial/parallel timings\n");
	figure();
	clf;
	title("Compute time versus nodes, kernel regression with different sample sizes");
	xlabel("nodes");
	ylabel("time (sec)");
	hold on;
	ts = zeros(rows(ns), 4);
	for i = 1:rows(ns)
		for nodes = 1:4;
			n = ns(i,:);
			x = 1:n;
			x = x';
			x = 2*x/n;
			# generate dependent variable
			trueline =  x + (x.^2)/2 - 3.1*(x.^3)/3 + 1.2*(x.^4)/4;
			sig = 0.5;
			y = trueline + sig*randn(n,1);
			bandwidth = 0.45;
			tic;
			fit = kernel_regression(x, y, x, bandwidth, "__kernel_normal", false, false, nodes-1);
			t1 = toc;
			ts(i, nodes) = t1;
		endfor
		plot(ts(i,:)');
		
	endfor
	hold off;
	printf("\nThe following table shows compute time by sample size (rows) and MPI ranks (columns)\n");
	printf("The sample sizes are too small to see a good speedup.\nEdit the program to increase them if you have enough memory to allow it\n");
	rlabels = char("1000", "2000", "4000");
	clabels = char("1", "2", "3", "4");
	prettyprint(ts, rlabels, clabels);
endif
