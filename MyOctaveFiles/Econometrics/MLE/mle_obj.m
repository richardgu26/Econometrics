# Copyright (C) 2003,2004,2005  Michael Creel <michael.creel@uab.es>
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

## usage: [obj_value, score] = mle_obj(theta, data, model, modelargs, nslaves)
##
## Returns the average log-likelihood for a specified model
## This is for internal use by mle_estimate

# adapted to use lbfgs 10-May-2010. Still need to parallelize using openmpi_ext


function [obj_value, score] = mle_obj(theta, data, model, modelargs, nslaves)

	n = rows(data);   
	if (nargin < 5) nslaves = 0; endif
	if (nslaves > 0)
		global NSLAVES PARALLEL NEWORLD NSLAVES TAG;

		nn = floor(n/(NSLAVES + 1)); # number of obsns per slave

		# The command that the slave nodes will execute
    		cmd=['contrib = mle_obj_nodes(theta, data, model, modelargs, nn); ',...	
        		'MPI_Send(contrib,0,TAG,NEWORLD);'];	

		# send items to slaves
		NumCmds_Send({"theta", "nn", "cmd"}, {theta, nn, cmd});

		# evaluate last block on master while slaves are busy
  		obj_value = mle_obj_nodes(theta, data, model, modelargs, nn);

		# collect slaves' results
		contrib = 0.0; # must be initialized to use MPI_Recv
  		for i = 1:NSLAVES
			MPI_Recv(contrib,i,TAG,NEWORLD);
			obj_value = obj_value + contrib;
		endfor

		# compute the average
  		obj_value = - obj_value / n;
  		score = "na"; # fix this later to allow analytic score in parallel
		
	else # serial version
		# to use lbfgsb
		[contribs, score] = feval(model, theta, data, modelargs);
		obj_value = - mean(contribs);
		if (isnan(obj_value)||isinf(obj_value)) obj_value = realmax; endif
		if nargout > 1
		  if isnumeric(score) # model may pass analytic score, otherwise it will be "na"
		    score = - mean(score)';
		  else # compute numeric score if we need to
		    score = numgradient(model, {theta, data, modelargs});
		    score = - mean(score)';
		  endif
		 endif
 		 # to use bfgsmin
# 		[contribs, score] = feval(model, theta, data, modelargs);
# 		obj_value = - mean(contribs);
# 		if isnumeric(score) # model may pass analytic score, otherwise it will be "na"
# 		    score = - mean(score)';
# 		endif
		


	endif

	# let's bullet-proof this in case the model goes nuts
	if (((abs(obj_value) == Inf)) || (isnan(obj_value)))
		obj_value = realmax/10;
	endif	    

endfunction
