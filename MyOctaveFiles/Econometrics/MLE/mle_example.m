## Copyright (C) 2003,2004  Michael Creel <michael.creel@uab.es>
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

## Example to show how to use MLE functions



# Generate data
n = 200; # how many observations?

# the explanatory variables: note that they have unequal scales
x = [ones(n,1) -rand(n,1) randn(n,1)];
theta = 1:3; # true coefficients are 1,2,3
theta = theta';

lambda = exp(x*theta);
y = poissrnd(lambda); # generate the dependent variable

####################################
# define arguments for mle_results #
####################################

# starting values
theta = zeros(3,1);
# data
data = [y, x];
# name of model to estimate
model = "poisson";
modelargs = {0}; # if this is zero the function gives analytic score, otherwise not
# parameter names
names = char("beta1", "beta2", "beta3");
mletitle = "Poisson MLE trial"; # title for the run

# controls for bfgsmin: 30 iterations is not always enough for convergence
control = {500,1};

# This displays the results
printf("\n\nanalytic score, unscaled data\n");
tic;
[theta, V, obj_value, infocrit] = mle_results(theta, data, model, modelargs, names, mletitle, 0, control);
toc

# This just calculates the results, no screen display
printf("\n\nanalytic score, unscaled data, no screen display\n");
theta = zeros(3,1);
tic;
[theta, obj_value, convergence] = mle_estimate(theta, data, model, modelargs, control);
printf("obj_value = %f, to confirm it worked\n", obj_value);
toc

# This show how to scale data during estimation, but get results
# for data in original units (recommended to avoid conditioning problems)
# This usually converges faster, depending upon the data
printf("\n\nanalytic score, scaled data\n");
[scaled_x, unscale] = scale_data(x);
data = [y, scaled_x];
theta = zeros(3,1);
tic;
[theta, V, obj_value, infocrit] = mle_results(theta, data, model, modelargs, names, mletitle, unscale, control);
toc

# Example using numeric score
printf("\n\nnumeric score, scaled data\n");
theta = zeros(3,1);
modelargs = {1}; # set the switch for no score
tic;
[theta, V, obj_value, infocrit] = mle_results(theta, data, model, modelargs, names, mletitle, unscale, control);
toc
# Example doing estimation in parallel on a cluster (requires MPITB)
# uncomment the following if you have MPITB installed
#  theta = zeros(3,1);
#  nslaves = 1;
#  title = "MLE estimation done in parallel";
#  [theta, V, obj_value, infocrit] = mle_results(theta, data, model, modelargs, names, mletitle, unscale, control, nslaves);
