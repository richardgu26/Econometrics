# Copyright (C) 2010 Michael Creel <michael.creel@uab.es>
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

# kernel_local_linear_example: example of how to use local linear kernel
# estimation. The default window widths are much too large for a good fit. Note
# how the local linear approach fits the edges much better.
#
# usage: kernel_local_linear_example(n, use_mpi, show_plot);
# INPUTS:
# n: sample size (default n = 1000)
# use_mpi: (default false) Use MPI or not. To use MPI, you should execute
#	mpirun -np X octave --eval "kernel_example(1000, true);"
#	where X is an integer number of ranks to use.
#	requires open_mpi extensions from Octave Forge
# show_plot (default true) display plot?
function kernel_local_linear_example(n, use_mpi, show_plot)

	if (nargin < 1) n = 1000; endif
	if (nargin < 2) use_mpi = false; endif
	if (nargin < 3) show_plot = true; endif


	if ((use_mpi) && not(MPI_Initialized)) MPI_Init; endif

	x = (0:n-1)/n;
	x = x';
	# generate dependent variable
	trueline = sin(2*pi*x) + cos(10*2*pi*x);
	sig = 0.5;
	y = trueline + sig*randn(n,1);
	if use_mpi
		CW = MPI_Comm_Load("NEWORLD");
		nodes = MPI_Comm_size(CW);
		# copy the same data over, so that the plot will be continuous
		# don't time this part
		if not(MPI_Comm_rank(CW))
			for i=1:nodes-1
				MPI_Send(y, i, 42, CW);
			endfor
		else
			y = MPI_Recv(0, 42, CW);
		endif
	endif
	tic; # timing only the kernel regression
	fit  = kernel_local_linear(x, y, x, "", "__kernel_epanechnikov", true, use_mpi); # use the default bandwidth
	fit2 = kernel_regression(x, y, x, "", "__kernel_epanechnikov", true, use_mpi); # use the default bandwidth
	t = toc;

	if (use_mpi)
		if not(MPI_Comm_rank(CW))
			printf("time for kernel regression example using %d data points on %d nodes: %f\n", n, nodes, t);
			if show_plot
				plot(x, fit, ";local linear;", x, fit2, ";kernel;", x, trueline,";true;");
				grid("on");
				title("Example 1: Kernel regression fit");
			endif
		endif
	endif

	if not(use_mpi)
		printf("time for kernel regression example using %d data points: %f\n", n, t);
		if show_plot
			plot(x, fit, ";local linear;", x, fit2, ";kernel;", x, trueline,";true;");
			grid("on");
			title("Example 1: Kernel regression fit");
		endif
	endif

	if (use_mpi && show_plot) pause(10); endif
	if (use_mpi) MPI_Finalize; endif

endfunction
