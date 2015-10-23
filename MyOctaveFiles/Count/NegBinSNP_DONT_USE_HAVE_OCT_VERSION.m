# This will work for orders up to and including k_gam=5
function [logdensity, score, prediction] = NegBinSNP(theta, data, otherargs)
	theta = parameterize(theta, {"NegBinSNP", otherargs});
	y = data[:,1];
	x = data[:,2:columns(data)];
	negbin_type = otherargs{1};
	[n, k] = size(x);
	beta = theta(1:k,:);
	alpha = exp(theta(k+1,:));
	gam = theta(k + 2:rows(theta),:);
	
# automatically scale the polynomial coefficients for (hopefully)
# better numeric stability
	for i = 1 : rows(gam)
		gam(i,:) = gam(i,:)/(10^(i+1));
	endfor
	
	lambda = mc_exp(x*beta);

	if (negbin_type == 1)
		psi = lambda / alpha;
		lampsi = alpha;  # lam/psi ratio
	else
		psi = ones(n,1) / alpha;
		lampsi = lambda*alpha;  # lam/psi ratio
	endif


	k_gam = rows(gam);
	if k_gam > 5 error("Only implemented for k_gam <=4"); endif
	m = zeros(n,11);

		m(:,1) = lambda;

		m(:,2) = (lambda .* (lambda + psi + lambda .* psi)) ./ psi;

		m(:,3) = (lambda .* (psi .^ 2 + 3 .* lambda .* psi .* (1 + psi) + lambda .^ 2 .* (2 + 3 .* psi + psi .^ 2))) ./ \
  		psi .^ 2;

		m(:,4) = (lambda .* (psi .^ 3 + 7 .* lambda .* psi .^ 2 .* (1 + psi) + \
    		6 .* lambda .^ 2 .* psi .* (2 + 3 .* psi + psi .^ 2) + \
			lambda .^ 3 .* (6 + 11 .* psi + 6 .* psi .^ 2 + psi .^ 3))) ./ psi .^ 3;

		m(:,5) = (lambda .* (psi .^ 4 + 15 .* lambda .* psi .^ 3 .* (1 + psi) + \
      		25 .* lambda .^ 2 .* psi .^ 2 .* (2 + 3 .* psi + psi .^ 2) + \
      		10 .* lambda .^ 3 .* psi .* (6 + 11 .* psi + 6 .* psi .^ 2 + psi .^ 3) + \
      		lambda .^ 4 .* (24 + 50 .* psi + 35 .* psi .^ 2 + 10 .* psi .^ 3 + psi .^ 4))) ./ psi .^ 4;
		
		m(:,6) = (lambda .* (psi .^ 5 + 31 .* lambda .* psi .^ 4 .* (1 + psi) + \
      		90 .* lambda .^ 2 .* psi .^ 3 .* (2 + 3 .* psi + psi .^ 2) + \
		    65 .* lambda .^ 3 .* psi .^ 2 .* (6 + 11 .* psi + 6 .* psi .^ 2 + psi .^ 3) + \
      		15 .* lambda .^ 4 .* psi .* (24 + 50 .* psi + 35 .* psi .^ 2 + 10 .* psi .^ 3 + psi .^ 4) + \
      		lambda .^ 5 .* (120 + 274 .* psi + 225 .* psi .^ 2 + 85 .* psi .^ 3 + 15 .* psi .^ 4 + \
         	psi .^ 5))) ./ psi .^ 5;
		
	 	m(:,7) = (lambda .* (psi .^ 6 + 63 .* lambda .* psi .^ 5 .* (1 + psi) + \
		    301 .* lambda .^ 2 .* psi .^ 4 .* (2 + 3 .* psi + psi .^ 2) + \
      		350 .* lambda .^ 3 .* psi .^ 3 .* (6 + 11 .* psi + 6 .* psi .^ 2 + psi .^ 3) + \
      		140 .* lambda .^ 4 .* psi .^ 2 .* (24 + 50 .* psi + 35 .* psi .^ 2 + 10 .* psi .^ 3 + psi .^ 4) + \
      		21 .* lambda .^ 5 .* psi .* (120 + 274 .* psi + 225 .* psi .^ 2 + 85 .* psi .^ 3 + \
         	15 .* psi .^ 4 + psi .^ 5) + \
      		lambda .^ 6 .* (720 + 1764 .* psi + 1624 .* psi .^ 2 + 735 .* psi .^ 3 + 175 .* psi .^ 4 + \
         	21 .* psi .^ 5 + psi .^ 6))) ./ psi .^ 6;
			
		m(:,8) = (lambda .* (psi .^ 7 + 127 .* lambda .* psi .^ 6 .* (1 + psi) + \
      		966 .* lambda .^ 2 .* psi .^ 5 .* (2 + 3 .* psi + psi .^ 2) + \
      		1701 .* lambda .^ 3 .* psi .^ 4 .* (6 + 11 .* psi + 6 .* psi .^ 2 + psi .^ 3) + \
      		1050 .* lambda .^ 4 .* psi .^ 3 .* (24 + 50 .* psi + 35 .* psi .^ 2 + 10 .* psi .^ 3 + \
         	psi .^ 4) + 266 .* lambda .^ 5 .* psi .^ 2 .* \
       		(120 + 274 .* psi + 225 .* psi .^ 2 + 85 .* psi .^ 3 + 15 .* psi .^ 4 + psi .^ 5) + \
      		28 .* lambda .^ 6 .* psi .* (720 + 1764 .* psi + 1624 .* psi .^ 2 + 735 .* psi .^ 3 + \
         	175 .* psi .^ 4 + 21 .* psi .^ 5 + psi .^ 6) + \
      		lambda .^ 7 .* (5040 + 13068 .* psi + 13132 .* psi .^ 2 + 6769 .* psi .^ 3 + \
         	1960 .* psi .^ 4 + 322 .* psi .^ 5 + 28 .* psi .^ 6 + psi .^ 7))) ./ psi .^ 7;
	
		m(:,9) = (lambda .* (psi .^ 8 + 255 .* lambda .* psi .^ 7 .* (1 + psi) + \
      		3025 .* lambda .^ 2 .* psi .^ 6 .* (2 + 3 .* psi + psi .^ 2) + \
      		7770 .* lambda .^ 3 .* psi .^ 5 .* (6 + 11 .* psi + 6 .* psi .^ 2 + psi .^ 3) + \
      		6951 .* lambda .^ 4 .* psi .^ 4 .* (24 + 50 .* psi + 35 .* psi .^ 2 + 10 .* psi .^ 3 + \
         	psi .^ 4) + 2646 .* lambda .^ 5 .* psi .^ 3 .* \
       		(120 + 274 .* psi + 225 .* psi .^ 2 + 85 .* psi .^ 3 + 15 .* psi .^ 4 + psi .^ 5) + \
      		462 .* lambda .^ 6 .* psi .^ 2 .* (720 + 1764 .* psi + 1624 .* psi .^ 2 + 735 .* psi .^ 3 + \
         	175 .* psi .^ 4 + 21 .* psi .^ 5 + psi .^ 6) + \
      		36 .* lambda .^ 7 .* psi .* (5040 + 13068 .* psi + 13132 .* psi .^ 2 + 6769 .* psi .^ 3 + \
         	1960 .* psi .^ 4 + 322 .* psi .^ 5 + 28 .* psi .^ 6 + psi .^ 7) + \
      		lambda .^ 8 .* (40320 + 109584 .* psi + 118124 .* psi .^ 2 + 67284 .* psi .^ 3 + \
         	22449 .* psi .^ 4 + 4536 .* psi .^ 5 + 546 .* psi .^ 6 + 36 .* psi .^ 7 + psi .^ 8))) ./ \
  			psi .^ 8;
		
		m(:,10) = (lambda .* (psi .^ 9 + 511 .* lambda .* psi .^ 8 .* (1 + psi) + \
      		9330 .* lambda .^ 2 .* psi .^ 7 .* (2 + 3 .* psi + psi .^ 2) + \
      		34105 .* lambda .^ 3 .* psi .^ 6 .* (6 + 11 .* psi + 6 .* psi .^ 2 + psi .^ 3) + \
      		42525 .* lambda .^ 4 .* psi .^ 5 .* (24 + 50 .* psi + 35 .* psi .^ 2 + 10 .* psi .^ 3 + \
         	psi .^ 4) + 22827 .* lambda .^ 5 .* psi .^ 4 .* \
       		(120 + 274 .* psi + 225 .* psi .^ 2 + 85 .* psi .^ 3 + 15 .* psi .^ 4 + psi .^ 5) + \
      		5880 .* lambda .^ 6 .* psi .^ 3 .* (720 + 1764 .* psi + 1624 .* psi .^ 2 + 735 .* psi .^ 3 + \
         	175 .* psi .^ 4 + 21 .* psi .^ 5 + psi .^ 6) + \
      		750 .* lambda .^ 7 .* psi .^ 2 .* (5040 + 13068 .* psi + 13132 .* psi .^ 2 + \
         	6769 .* psi .^ 3 + 1960 .* psi .^ 4 + 322 .* psi .^ 5 + 28 .* psi .^ 6 + psi .^ 7) + \
      		45 .* lambda .^ 8 .* psi .* (40320 + 109584 .* psi + 118124 .* psi .^ 2 + \
         	67284 .* psi .^ 3 + 22449 .* psi .^ 4 + 4536 .* psi .^ 5 + 546 .* psi .^ 6 + \
         	36 .* psi .^ 7 + psi .^ 8) + \
      		lambda .^ 9 .* (362880 + 1026576 .* psi + 1172700 .* psi .^ 2 + 723680 .* psi .^ 3 + \
         	269325 .* psi .^ 4 + 63273 .* psi .^ 5 + 9450 .* psi .^ 6 + 870 .* psi .^ 7 + \
         	45 .* psi .^ 8 + psi .^ 9))) ./ psi .^ 9;

		m(:,11) = (lambda .* (psi .^ 10 + 1023 .* lambda .* psi .^ 9 .* (1 + psi) + \
      		28501 .* lambda .^ 2 .* psi .^ 8 .* (2 + 3 .* psi + psi .^ 2) + \
      		145750 .* lambda .^ 3 .* psi .^ 7 .* (6 + 11 .* psi + 6 .* psi .^ 2 + psi .^ 3) + \
      		246730 .* lambda .^ 4 .* psi .^ 6 .* \
       		(24 + 50 .* psi + 35 .* psi .^ 2 + 10 .* psi .^ 3 + psi .^ 4) + \
      		179487 .* lambda .^ 5 .* psi .^ 5 .* \
       		(120 + 274 .* psi + 225 .* psi .^ 2 + 85 .* psi .^ 3 + 15 .* psi .^ 4 + psi .^ 5) + \
      		63987 .* lambda .^ 6 .* psi .^ 4 .* (720 + 1764 .* psi + 1624 .* psi .^ 2 + 735 .* psi .^ 3 + \
         	175 .* psi .^ 4 + 21 .* psi .^ 5 + psi .^ 6) + \
      		11880 .* lambda .^ 7 .* psi .^ 3 .* (5040 + 13068 .* psi + 13132 .* psi .^ 2 + \
         	6769 .* psi .^ 3 + 1960 .* psi .^ 4 + 322 .* psi .^ 5 + 28 .* psi .^ 6 + psi .^ 7) + \
      		1155 .* lambda .^ 8 .* psi .^ 2 .* (40320 + 109584 .* psi + 118124 .* psi .^ 2 + \
         	67284 .* psi .^ 3 + 22449 .* psi .^ 4 + 4536 .* psi .^ 5 + 546 .* psi .^ 6 + \
         	36 .* psi .^ 7 + psi .^ 8) + \
      		55 .* lambda .^ 9 .* psi .* (362880 + 1026576 .* psi + 1172700 .* psi .^ 2 + \
         	723680 .* psi .^ 3 + 269325 .* psi .^ 4 + 63273 .* psi .^ 5 + 9450 .* psi .^ 6 + \
         	870 .* psi .^ 7 + 45 .* psi .^ 8 + psi .^ 9) + \
      		lambda .^ 10 .* (3628800 + 10628640 .* psi + 12753576 .* psi .^ 2 + \
         	8409500 .* psi .^ 3 + 3416930 .* psi .^ 4 + 902055 .* psi .^ 5 + 157773 .* psi .^ 6 + \
         	18150 .* psi .^ 7 + 1320 .* psi .^ 8 + 55 .* psi .^ 9 + psi .^ 10))) ./ psi .^ 10;
	
	
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
	
	log_base_density = mc_lgamma(y + psi) - mc_lgamma(psi) - mc_lgamma(y+1) \
				+ psi .* log(1 ./ (lampsi + 1)) + y .* log(lampsi ./ (lampsi + 1));
	
	poly = ones(n,1);
	for i = 1:k_gam
	poly = poly + (y .^ (i))*gam(i,:);  
	endfor
	polysq = eps + poly .^2; # square the polynomial

	logdensity = log(polysq) + log_base_density - log(norm_factor);
	score = "na";
endfunction





