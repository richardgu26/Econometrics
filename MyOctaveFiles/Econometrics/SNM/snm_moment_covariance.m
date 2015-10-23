# Copyright (C) 2008 Michael Creel <michael.creel@uab.es>
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
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

#
# snm_moment_covariance: estimates the covariance matrix of moments conditions
# by Monte Carlo
# usage: omega = snm_moment_covariance(theta, data, momentargs, reps, verbose)
#
# inputs:
#    theta: column vector initial parameters
#    data: data matrix
#    momentargs: (cell) additional inputs needed to compute moments.
#	      May be empty ("")
#    reps: number of replications to use (default 1000)
#    verbose: (boolean) progress monitor printed to screen (false by default)
#
# output:
#    omega: estimated covariance matrix of average moments

function omega = snm_moment_covariance(theta, data, momentargs, reps, verbose)

	if nargin < 3 error("snm_covariances: 3 arguments required"); endif
	if nargin < 4 reps = 1000; endif # default number of reps
	if nargin < 5 verbose = false; endif

	model = momentargs{1};
	modelargs = momentargs{2};
	S = modelargs{2};
	do_cond = momentargs{7};
	do_uncond = momentargs{8};
	realdata = data;
	n = rows(data);


	for rep = 1:reps
		# get bootstrap data (size n)
		modelargs{1} = "";
		modelargs{2} = n;
		modelargs{4} = true; # make instruments too
		[data, L, K, junk] = feval(model, theta, modelargs);
		M = columns(data) - L -K;
	
		# split data into endogs, conditioning variables, and instruments
		endogs = data(:,1:L);
		condvars = data(:,L+1:L+K);
		instruments = data(:,L+K+1:L+K+M);
		[endogs, condvars, P] = snm_dataprep(endogs, condvars);
		data = [endogs, condvars, instruments];

		# set up a random draws of size S
		modelargs{1} = "";
		modelargs{2} = S;
		modelargs{4} = false; # no instruments here
		[junk, junk, junk, snm_draws] = feval(model, theta, modelargs);
		modelargs{1} = snm_draws;
		momentargs{2} = modelargs;
		momentargs{6} = P;

		# errors for this draw
		[ce, ue] = snm_residuals(theta, data, momentargs);

		if do_cond & do_uncond
			# interact instruments with errors
			M = columns(instruments);
			m = zeros(n, L*M);
			for i=1:L
				m(:,i*M-M+1:i*M) = instruments .* repmat(ce(:,i), 1, M); # errors interacted with instruments
			endfor
			if (rep == 1) ms = zeros(reps, L + L*M); endif
			ms(rep,:) = mean([m ue]);
		elseif do_cond & (!do_uncond)
			# interact instruments with errors
			M = columns(instruments);
			m = zeros(n, L*M);
			for i=1:L
				m(:,i*M-M+1:i*M) = instruments .* repmat(ce(:,i), 1, M); # errors interacted with instruments
			endfor
			if (rep == 1) ms = zeros(reps, L*M); endif
			ms(rep,:) = mean(m);
		else
			if (rep == 1) ms = zeros(reps, L); endif
			ms(rep,:) = mean(ue);
		endif
  		if (verbose) printf("snm_moment_covariance: %d of %d replications are done\n", rep, reps); endif
	endfor
	omega = n*cov(ms);

endfunction
