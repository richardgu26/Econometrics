# Copyright (C) 2007, 2009 Michael Creel <michael.creel@uab.es>
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
#
# usage: kernel_example(n, use_mpi, show_plot);
# INPUTS:
# n: sample size (default n = 500)
# use_mpi: (default false) Use MPI or not. To use MPI, you should execute
#	mpirun -np X octave --eval "kernel_example(1000, true);"
#	where X>1 is an integer number of ranks to use.
#	requires open_mpi extensions from Octave Forge
# show_plot (default true) display plot?
function kernel_example(n, use_mpi, show_plot)

	if (nargin < 1) n = 2000; endif
	if (nargin < 2) use_mpi = false; endif
	if (nargin < 3) show_plot = true; endif


	if ((use_mpi) && not(MPI_Initialized)) MPI_Init; endif

	x = (0:n-1)/n;
	x = x';
	# generate dependent variable
	trueline = sin(2*pi*x) + cos(10*2*pi*x);
	sig = 0.1;
	y = trueline + sig*randn(n,1);
	trueline2 = cos(2*pi*x) + sin(10*2*pi*x);
	sig = 0.1;
	y2 = trueline2 + sig*randn(n,1);
	y = [y y2];
	trueline = [trueline trueline2];


	if use_mpi
		CW = MPI_Comm_Load("NEWORLD");
		nodes = MPI_Comm_size(CW);
		if (nodes < 2) MPI_Finalize; error("kernel_example: when using MPI, must specify at least 2 ranks"); endif
		# copy the same data over, so that the plot will be continuous
		# don't time this part
		if (not(MPI_Comm_rank(CW)) && (nodes > 1))
			for i=1:nodes-1
				MPI_Send(y, i, 42, CW);
			endfor
		else
			y = MPI_Recv(0, 42, CW);
		endif
	endif
	tic; # timing only the kernel regression
	fit = kernel_regression(x, y, x, 0.02, "__kernel_normal", true, use_mpi, false, true); # use the default bandwidth
	treg = toc;

		
	############################################################
	# kernel density example: x (from above) is uniform (0,1)
	# evaluation point are on a grid for plotting
	bandwidth = 0.55;
	# get optimal bandwidth (time consuming, uncomment if you want to try it)
	# bandwidth = kernel_optimal_bandwidth(data);
	# get the fitted density and do a plot
	tic;
	dens = kernel_density(x, x, bandwidth, "__kernel_normal", true, use_mpi, false, true);
	tdens = toc;


	
	############################################################
	# kernel density example: bivariate
	data = [x y(:,1)];
	# evaluation points are on a grid for plotting
	a = (0:20-1)/20;
	a = a';
	b = 7*a-3.5;
	gridsize = rows(a);
	[grid_x, grid_y] = meshgrid(a, b);
	eval_points = [vec(grid_x) vec(grid_x)];
	bandwidth = 0.3;
	# get optimal bandwidth (time consuming, uncomment if you want to try it)
	# bandwidth = kernel_optimal_bandwidth(data);
	# get the fitted density and do a plot
	tic;
	dens2 = kernel_density(eval_points, data, bandwidth, "__kernel_epanechnikov", false, use_mpi, false, true);
	tdens2 = toc;

	
	if (use_mpi)
		if not(MPI_Comm_rank(CW))
			printf("time for kernel regression example using %d data points on %d nodes: %f\n", n, nodes, treg);
			printf("time for kernel density example using %d data points on %d nodes: %f\n", n, nodes, tdens);
			printf("time for bivariate kernel density example using %d data points on %d nodes: %f\n", n, nodes, tdens2);
			if show_plot
				close all;
				plot(x, fit, ";fit;", x, trueline,";true;");
				grid("on");
				title("Example 1: Kernel regression fit");
				figure;
				plot(x, dens, ";fitted density;");
				grid("on");
				title("Kernel density fit");
				dens = reshape(dens2, gridsize, gridsize);
				figure();
				surf(grid_x, grid_y, dens);
				title("Bivariate kernel density fit");
			endif
		endif
	else
		printf("time for kernel regression example using %d data points: %f\n", n, treg);
		printf("time for kernel density example using %d data points: %f\n", n, tdens);
		printf("time for bivariate kernel density example using %d data points: %f\n", n, tdens2);
		if show_plot
			close all;
			plot(x, fit, ";fit;", x, trueline,";true;");
			grid("on");
			title("Example 1: Kernel regression fit");
			figure;
			plot(x, dens, ";fitted density;");
			grid("on");
			title("Kernel density fit");
			dens = reshape(dens2, gridsize, gridsize);
			figure();
			surf(grid_x, grid_y, dens);
			title("Bivariate kernel density fit");
		endif
	endif

	if (use_mpi) MPI_Finalize; endif
	 
endfunction
