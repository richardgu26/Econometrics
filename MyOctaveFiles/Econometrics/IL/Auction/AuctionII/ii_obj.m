# CU-GMM objective function for indirect inference
function [obj_value score] = ii_obj(theta, randdraws1, randdraws2, Z, overid)

	% simulated aux stats
	Zs = aux_stat(theta, randdraws1, randdraws2);

	if !overid
		Zs = Zs(:,1:2);
		Z = Z(:,1:2);
	endif

	% moment conditions
	Zbar = mean(Zs);
	m = Z - Zbar; % row vector

	% weight matrix
	if overid
	  e = Zs - Z;
	  V = e'*e/rows(e);
	  W = inv(V);
	else
	  W = eye(2);
	endif  
	
	obj_value = m*W*m';
	score = "";
endfunction

