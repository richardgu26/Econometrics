# This will work for orders up to and including k_gam=5
function [logdensity, score, prediction] = PoissonSNP(theta, data, otherargs)
	y = data(:,1);
	x = data(:,2:columns(data));
	[n, k] = size(x);
	beta = theta(1:k,:);
	gam = theta(k + 1:rows(theta),:);
	
# automatically scale the polynomial coefficients for (hopefully)
# better numeric stability
	for i = 1 : rows(gam)
		gam(i,:) = gam(i,:)/(10^(i+1));
	endfor
	
	lambda = exp(x*beta);

	k_gam = rows(gam);
	if k_gam > 5 error("Only implemented for k_gam <=5"); endif
	m = zeros(n,2*k_gam+1);

	if (k_gam >= 1) 
		m(:,1) = lambda;
		m(:,2) = lambda.*(1+lambda);
		m(:,3) = lambda.*(1 + 3 .* lambda + lambda.^2);

	
	endif

	if (k_gam >= 2)
		m(:,4) = lambda.*(1 + 7 .* lambda + 6 .* lambda.^2 + lambda.^3);
		m(:,5) = lambda.*(1 + 15 .* lambda + 25 .* lambda.^2 + 10 .* lambda.^3 + lambda.^4);

	endif

	if (k_gam >= 3)
		m(:,6) = lambda.*(1 + 31 .* lambda + 90 .* lambda.^2 + 65 .* lambda.^3 + 15 .* lambda.^4 + lambda.^5);
		m(:,7) = lambda.*(1 + 63 .* lambda + 301 .* lambda.^2 + 350 .* lambda.^3 + 140 .* lambda.^4 + \
		            21 .* lambda.^5 + lambda.^6);

	endif

	if (k_gam >= 4)
		m(:,8) = lambda.*(1 + 127 .* lambda + 966 .* lambda.^2 + 1701 .* lambda.^3 + 1050 .* lambda.^4 + \
		           266 .* lambda.^5 + 28 .* lambda.^6 + lambda.^7);

		m(:,9) = lambda.*(1 + 255 .* lambda + 3025 .* lambda.^2 + 7770 .* lambda.^3 + 6951 .* lambda.^4 + \
		           2646 .* lambda.^5 + 462 .* lambda.^6 + 36 .* lambda.^7 + lambda.^8);
			
	endif
	
	if (k_gam >= 5)
		m(:,10) = lambda.*(1 + 511 .* lambda + 9330 .* lambda.^2 + 34105 .* lambda.^3 + 42525 .* lambda.^4 + \
		           22827 .* lambda.^5 + 5880 .* lambda.^6 + 750 .* lambda.^7 + 45 .* lambda.^8 + lambda.^9);

		m(:,11) = lambda.*(1 + 1023 .* lambda + 28501 .* lambda.^2 + 145750 .* lambda.^3 + 246730 .* lambda.^4 + \
		           179487 .* lambda.^5 + 63987 .* lambda.^6 + 11880 .* lambda.^7 + 1155 .* lambda.^8 + 55 .* lambda.^9 + lambda.^10);
		
	endif	
	
	
	if (k_gam == 1)
		norm_factor = 1 + gam(1,:) .* (2 .* m(:,1) + gam(1,:) .* m(:,2));
	elseif (k_gam == 2)
		norm_factor = 1 + gam(1,:) .^ 2 .* m(:,2) + 2 .* gam(1,:) .* (m(:,1) + gam(2,:) .* m(:,3)) + \
		  		gam(2,:) .* (2 .* m(:,2) + gam(2,:) .* m(:,4));	
	elseif (k_gam == 3)
		norm_factor = 1 + gam(1,:) .^ 2 .* m(:,2) + 2 .* gam(2,:) .* m(:,2) + 2 .* gam(3,:) .* m(:,3) + gam(2,:) .^ 2 .* m(:,4) + \
		  		2 .* gam(1,:) .* (m(:,1) + gam(2,:) .* m(:,3) + gam(3,:) .* m(:,4)) + 2 .* gam(2,:) .* gam(3,:) .* m(:,5) + \
		  		gam(3,:) .^ 2 .* m(:,6);
	elseif (k_gam == 4)
		norm_factor = 1 + gam(1,:) .^ 2 .* m(:,2) + 2 .* gam(2,:) .* m(:,2) + 2 .* gam(3,:) .* m(:,3) + gam(2,:) .^ 2 .* m(:,4) + \
		  			2 .* gam(4,:) .* m(:,4) + 2 .* gam(2,:) .* gam(3,:) .* m(:,5) + \
		  			2 .* gam(1,:) .* (m(:,1) + gam(2,:) .* m(:,3) + gam(3,:) .* m(:,4) + gam(4,:) .* m(:,5)) + \
		  			gam(3,:) .^ 2 .* m(:,6) + 2 .* gam(2,:) .* gam(4,:) .* m(:,6) + 2 .* gam(3,:) .* gam(4,:) .* m(:,7) + \
		  			gam(4,:) .^ 2 .* m(:,8);
	else
		norm_factor = 1 + gam(1,:) .^ 2 .* m(:,2) + 2 .* gam(2,:) .* m(:,2) + 2 .* gam(3,:) .* m(:,3) + gam(2,:) .^ 2 .* m(:,4) + \
					  	2 .* gam(4,:) .* m(:,4) + 2 .* gam(2,:) .* gam(3,:) .* m(:,5) + 2 .* gam(5,:) .* m(:,5) + gam(3,:) .^ 2 .* m(:,6) + \
		  				2 .* gam(2,:) .* gam(4,:) .* m(:,6) + 2 .* gam(1,:) .* \
		   				(m(:,1) + gam(2,:) .* m(:,3) + gam(3,:) .* m(:,4) + gam(4,:) .* m(:,5) + gam(5,:) .* m(:,6)) + \
		  				2 .* gam(3,:) .* gam(4,:) .* m(:,7) + 2 .* gam(2,:) .* gam(5,:) .* m(:,7) + gam(4,:) .^ 2 .* m(:,8) + \
		  				2 .* gam(3,:) .* gam(5,:) .* m(:,8) + 2 .* gam(4,:) .* gam(5,:) .* m(:,9) + gam(5,:) .^ 2 .* m(:,10);
	endif

	norm_factor = eps + norm_factor;

	temp = [1 ; gam];		   	# gam_0 = 1
	prediction = zeros(n,1);
	for i = 0:k_gam # 	// Cameron and Johansson eqn 4 (could modify easily to get variance, etc.
		for j = 0:k_gam
			prediction = prediction + temp(i+1,:)*temp(j+1,:)*m(:,i+j+1);
		endfor
	endfor
	prediction  = prediction  ./ norm_factor;
	
	log_base_density = - lambda	+ y .* (x*beta) - mc_lgamma(y+1);
	poly = ones(n,1);
	for i = 1:k_gam
	poly = poly + (y .^ (i))*gam(i,:);  
	endfor
	polysq = eps + poly .^2; # square the polynomial

	logdensity = log(polysq) + log_base_density - log(norm_factor);
	score = "na";
endfunction





