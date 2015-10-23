## usage:  results2tex (betas, sterrs, file, cname, rname, prec, align, lines)
##
## Saves estimated coefficients and st. errors to LaTeX table file.
##
## betas and sterrs are KxM matrices of estimated coefficients
## and st. errors, where K is the number of coefs. and M is number
## of models. If some models do not contain all coefficient, the
## corresponding value of betas should be set to nan, and the
## standard error can be arbitrary (it won't be printed).
##
## cname and rname are matrices whose rows are the column
## respectively row headings, if c/rname == "" (default), no headings
## are printed.
##
## prec gives the number of digits printed after the
## comma. If prec<0 (default) the number is printed unformatted
##
## align gives the alignment of the data, default is "r" (flush right). If
##
## lines == 0 (default) or == "none", no lines are printed in the
## tabular; if lines == 2 or == "full" there are lines around each cell,
## otherwise there are lines at the border and between the headings and
## the data.
## Author:  Michael Creel <Michael.creel@uab.es>
## Based in part on save2tex.m, by Andreas Wei <Andreas.Weingessel@ci.tuwien.ac.at>
## Description:  Save to a file in a LaTeX tabular environment

# Example of usage that shows how missing elements are treated
# a = rand(3,4);
# b = rand(2,1);
# b = postpad(b,rows(a),nan);
# data = [b, a];
# sterr = rand(size(data));
# rname = char("r1","r2","r3");
# cname = char("c1","c2","c3","c4","c5");
#
# results2tex(data, sterr, "test.tex", cname, rname, 5, "r", 1);

function results2tex(betas, sterrs, file, cname, rname, prec, align, lines)

	if ((nargin == 1) || (nargin > 8))
    		printf("usage: save2tex (X, file, cname, rname, prec, align, lines)\n");
	endif

	if (nargin < 8) lines = 0; endif
	if (nargin < 7) align = "r"; endif  # center align seems to get decimals
	if (nargin < 6) prec = 4; endif
	if (nargin < 5) rname = ""; endif
	if (nargin < 4) cname = ""; endif


	if (ischar (lines))
    		if (strcmp (lines, "full"))
    	  		lines = 2;
	    	elseif (strcmp (lines, "none"))
    		  	lines = 0;
	    	else
    		  	lines = 1;
	    	endif
	endif

	nr = rows(betas);
	nc = columns(betas);
	is_rn = columns(rname);
	is_cn = columns(cname);

	if ((is_rn) && (rows(rname) != nr))
    		error ("results2tex: Numbers of rows and row names do not match.");
	endif

	if ((is_cn) && (rows(cname) != nc))
    		error ("results2tex: Numbers of columns and column names do not match.");
	endif


	### open output file
	FN = fopen (file, "w");
	if (FN < 0)
    		error ("save2tex: Can not open File %s", file);
	endif

	### create format line
	LSTR = "";
	LFSTR = "";
	if (lines)
    		LSTR = "|";
    		if (lines == 2)
    	  		LFSTR = "|";
    		endif
	endif
	STR = ["\\begin{tabular}{", LSTR];
	if (is_rn)
    		STR = [STR, "l", LSTR];
	endif

	for i=1:nc
    		STR = [STR, align, LFSTR];
	endfor

	if (lines == 2)
    		STR = [STR, "}\n"];
	else
    		STR = [STR, LSTR, "}\n"];
	endif

	fprintf (FN, STR);

	if (lines)
    		fprintf (FN, "\\hline\n");
	endif


	### print column headers
	if (is_cn)
    		if (is_rn)
    	  		STR = "  & ";
    		else
    	  		STR = "  ";
    		endif

		for i = 1:nc
    	  		STR = [STR, cname(i,:), " & "];
    		endfor

		les=length(STR);
    		fprintf(FN, [STR(1:les-2), "\\\\\n"]);

		if (lines)
    	  		fprintf (FN, "\\hline\n");
    		endif
	endif

	test = abs(betas ./ sterrs);
	pvalue = 2*(1 - normal_cdf(test));

	### print data
	for i = 1:nr
    		STR = [rname(i,:), "  & "];
		STR2 = "  & "; # no row labels for st. errs or pvalues
		STR3 = STR2;
  		nanform = sprintf(" %%s  & ");
		nanform2 = sprintf("\\small%%s & ");
		nanform3 = sprintf("\\small%%s & ");

  		form = sprintf(" %%.%df  & ", prec);
		form2 = sprintf("\\small(%%.%df) & ", 3);
		form3 = sprintf("\\small[%%.%df] & ", 3);


	    	for j = 1:nc
			if isnan(betas(i,j))
    				STR = [STR, sprintf(nanform,"na")];
     				STR2 = [STR2, sprintf(nanform2," ")];
     				STR3 = [STR3, sprintf(nanform3," ")];
			else
     				STR = [STR, sprintf(form,betas(i,j))];
     				STR2 = [STR2, sprintf(form2,sterrs(i,j))];
     				STR3 = [STR3, sprintf(form3,pvalue(i,j))];
			endif
   		endfor

		les = length(STR);
		les2 = length(STR2);
		les3 = length(STR3);

		if (i < nr)
	   	  	if (lines == 2)
				fprintf(FN, [STR(1:les-2), "\\\\\n\\hline\n"]);
				fprintf(FN, [STR2(1:les2-2), "\\\\\n\\hline\n"]);
				fprintf(FN, [STR3(1:les3-2), "\\\\\n\\hline\n"]);
      			else
	     			fprintf(FN, [STR(1:les-2), "\\\\\n"]);
     				fprintf(FN, [STR2(1:les2-2), "\\\\\n"]);
     				fprintf(FN, [STR3(1:les3-2), "\\\\\n"]);
      			endif
   		else
   	  		if (lines)
				fprintf(FN, [STR(1:les-2),   "\\\\\n"]);
				fprintf(FN, [STR2(1:les2-2),   "\\\\\n"]);
				fprintf(FN, [STR3(1:les3-2), "\\\\\n\\hline\n"]);
   	  		else
     				fprintf(FN, [STR(1:les-2), "\\\\\n"]);
     				fprintf(FN, [STR2(1:les2-2), "\\\\\n"]);
     				fprintf(FN, [STR3(1:les3-2), "\n"]);
   	  		endif
   		endif
	endfor
	### finish output
	fprintf(FN, "\\end{tabular}\n");
	fprintf(FN, "\n \\small ( )  = standard errors; [ ] = p-values\n");

	fclose(FN);
endfunction
