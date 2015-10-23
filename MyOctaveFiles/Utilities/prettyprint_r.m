	# this prints matrices with row labels but no column labels
	function prettyprint_r(mat, rlabels)


		# now print the row labels and rows
		printf("\n");
		k = rows(mat);
		for i = 1:k
			if ischar(rlabels(i,:))
				printf(rlabels(i,:));
			else
				printf("%i", rlabels(i,:));
			endif
			printf("  %8.4f", mat(i,:));
			printf("\n");
		endfor
	endfunction
