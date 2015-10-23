# this is the funtion that appears in the
# FoggySurface and FoggyContour graphics,
# It is the Octave sombrero, tilted a bit
function [z, score] = FoggyMountainObj(theta)
	theta1 = theta(1,:);
	theta2 = theta(2,:);
   	r = sqrt (theta1^2 + theta2^2) + eps;
	z = sin(r) / r;
	z = z + theta1/80 - 0.1*(theta1/16)^2;
	z = -z/10; # switch to minimization
	score = "na";
endfunction	
