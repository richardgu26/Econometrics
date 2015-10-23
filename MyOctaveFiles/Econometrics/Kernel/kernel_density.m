# Copyright (C) 2010 Michael Creel <michael.creel@uab.es>
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

# kernel_density: kernel density estimator
#
# usage:
# 	dens = kernel_density(eval_points, data, bandwidth)
#
# inputs:
#	eval_points: PxK matrix of points at which to calculate the density
#	data: NxK vector of observations of the G dependent variables
#	bandwidth (optional): positive scalar, the smoothing parameter.
#		Default is N ^ (-1/(4+K))
#	kernel (optional): string. Name of the kernel function. Default is
#		Gaussian kernel.
#	prewhiten bool (optional): default true. If true, rotate data
# 		using Choleski decomposition of inverse of covariance,
#		to approximate independence after the transformation, which
#		makes a product kernel a reasonable choice.
#	use_mpi: (optional) Whether or not to use MPI. Require openmpi_ext. default false.
#		rank 0 does not do computations, it gathers results. So you should always
#		use at least 2 MPI ranks if this option is set to true
#	do_cv: bool (optional). default false. If true, calculate leave-1-out
#		 fit to calculate the cross validation score
#	verbose: bool (optional) default false. monitor receipt of results, if using MPI
# outputs:
#	dens: Px1 vector: the fitted density value at each of the P evaluation points.
#

function dens = kernel_density(eval_points, data, bandwidth, kernel, prewhiten, use_mpi, do_cv, verbose)
	if nargin < 2; error("kernel_density: at least 2 arguments are required"); endif

	[n, k] = size(data);
	nn = rows(eval_points);

	# set defaults for optional args
	if (nargin < 3) bandwidth = (n ^ (-1/(4+k))); endif	# bandwidth - see Li and Racine pg. 66
	if (nargin < 4) kernel = "__kernel_normal"; endif 	# what kernel?
	if (nargin < 5) prewhiten = true; endif 		# automatic prewhitening?
	if (nargin < 6)	use_mpi = false; endif 			# compute nodes to use if parallel
	if (nargin < 7)	do_cv = false; endif 			# ordinary or leave-1-out
	if (nargin < 8)	verbose = false; endif 			# monitor receipt of results? (only for use_mpi=true)


	if (strcmp(bandwidth,"")) bandwidth = (n ^ (-1/(4+k))); endif # allows use of ""
	if (strcmp(kernel,"")); kernel = "__kernel_normal"; endif # allows use of ""

	if prewhiten
		H = bandwidth*chol(cov(data));
 	else
		H = bandwidth;
	endif
	H_inv = pinv(H);

	# weight by inverse bandwidth matrix
	eval_points = eval_points*H_inv;
	data = data*H_inv;


	# ordinary single core fit
	if not(use_mpi) #
		if do_cv
			z = kernel_density_nodes(eval_points, data, kernel, 1);
		else	
			z = kernel_density_nodes(eval_points, data, kernel);
		endif
	else
		# some details
		if not(MPI_Initialized) error("kernel_regression: attempted to use MPI without initializing"); endif
		use_mpi = true;
		CW = MPI_Comm_Load("NEWORLD");
		nodes = MPI_Comm_size(CW);
		node = MPI_Comm_rank(CW);
		mytag = 48;

		if node # compute nodes
			while true
				limits = MPI_Recv(0, mytag, CW);
				startblock = limits(1,:);
				if (startblock < 0)	break; endif # frontend sends this when all done
				endblock = limits(2,:);
				myeval = eval_points(startblock:endblock,:);
				if do_cv
					z = kernel_density_nodes(myeval, data, kernel, startblock);
				else	
					z = kernel_density_nodes(myeval, data, kernel);
				endif	
				MPI_Send(limits, 0, mytag, CW);
				MPI_Send(z, 0, mytag+1, CW);
			endwhile
		else # frontend

		# block size: large enough for efficiency, small enough for load balance 
		blocksize = floor(min(5e7/(n*k), nn/(nodes-1)));
			endblock = 0;
			# send out initial assignments
			for i = 1:nodes-1
				startblock = endblock + 1;
				endblock = endblock + blocksize;
				if (endblock > nn) endblock = nn; endif # that means we're done
				if (startblock > nn) startblock = -1; endif # that will signal nodes they're done
				limits = [startblock; endblock];
				MPI_Send(limits, i, mytag, CW);
			endfor

			done=false;
			while not(done)
				for i = 1:nodes-1
					limits = MPI_Recv(i, mytag, CW);
					contrib = MPI_Recv(i, mytag+1, CW);
					mystartblock = limits(1,:);
					myendblock = limits(2,:);
					if verbose
						if not(node)
							printf("Just received fits for points %d though %d out of %d\n", ...
						 mystartblock, myendblock, nn);
						endif
					endif
					z(mystartblock:myendblock,:) = contrib;
					startblock = endblock + 1;
					endblock = endblock + blocksize;
					if (endblock > nn) endblock = nn; endif # that means we're done
					if (startblock > nn) startblock = -1; endif # that will signal nodes they're done
					limits = [startblock; endblock];
					MPI_Send(limits, i, mytag, CW);
					if (myendblock == nn);
						done=true;
						break;
					endif # that means all have been received
				endfor
			endwhile
		endif
	endif
	dens = z*det(H_inv);
endfunction
