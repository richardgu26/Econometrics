# This is a smooth replacement for log() that is defined on the
# real line. This is useful to avoid crashes during intermediate
# iterations while doing nonlinear estimation.
# Reference: A.R. Gallant, 1987, Nonlinear Statistical Models, pg. 319.
#
# Michael.Creel@uab.es
# 13/3/2000

function y = slog(z)
    a =  - 299.999999999999886;
    b = 5667638086.9808321;
    c =  - 28288190434904165;
    t1 = z <= 0;
    t2 = (0. < z) & (z <= 1e-7);
    t3 = 1 - t1 - t2;	
    t1 = t1 .* (a + b * z);
    t2 = t2 .* (a + b * z + c * z .^ 2);
    y = t3 .* z + (1 - t3);
    t3 = t3 .* log(y);
    y = t1 + t2 + t3;
endfunction