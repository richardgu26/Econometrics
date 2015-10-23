# digamma   Computes the digamma (also called psi) function (d log gamma).
# 		  Use:  p = digamma(x) where x is a scalar or vector.
# 
#  H2M/cnt Toolbox, Version 2.0
#  Olivier Cappé, 16/08/2001 - 17/08/2001
#  ENST Dpt. TSI / LTCI (CNRS URA 820), Paris
# 
#  Adapted from GAUSS code by Paul L. Fackler available at
#   http://www.american.edu/academic.depts/cas/econ/gaussres/pdf/loggamma.src
# 
#  Formula 6.3.18 with recurrence formula 6.3.5 from Abromowitz and
#  Stegun, Dover (1965)
# Note by Creel - this is GPL'd, and was found as part of http://www.tsi.enst.fr/~cappe/h2m/
function p = mc_digamma(x)

	if (any(x <= 0))
	  error('digamma requires positive arguments.');
	endif

	x = x+6;
	p = 1./(x.^2);
	p = (((0.004166666666667*p-0.003968253986254).*p + 0.008333333333333).*p-0.083333333333333).*p;
	p = p+log(x)-0.5./x-1./(x-1)-1./(x-2)-1./(x-3)-1./(x-4)-1./(x-5)-1./(x-6);
endfunction	
