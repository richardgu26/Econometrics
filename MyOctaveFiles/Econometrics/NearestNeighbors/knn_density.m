# parallelize the loop
# can kd (a structure, I think) be sent by MPI?
function dens = knn_density(eval_points, data, k, prewhiten, do_cv)
ann;
	if nargin < 2; error("knn_density: at least 2 arguments are required"); endif

	[n, G] = size(data);
	nn = rows(eval_points);

	# set defaults for optional args
	if (nargin < 3) k = floor(1.5*(n^(0.25))); endif	# number of neighbors
	if !isnumeric(k) k = floor(1.5*(n^(0.25))); endif # allow using "" as a placeholder  	
	if (nargin < 4) prewhiten = false; endif 	# automatic prewhitening?
	if (nargin < 5)	do_cv = false; endif

	if prewhiten
		[data scalecoefs] = scale_data([data; eval_points]);
		eval_points = data(n+1:n+nn,:);
		data = data(1:n,:);
	endif

	dens = zeros(nn,1);
	c = (pi^(G/2))/gamma((G+2)/2);
	kd = ANNkd_tree(data);
	for i = 1:nn 
		[nn_idx,dd] = kd.annkPriSearch(eval_points(i,:), k);
		d = dd(k);
		dens(i) = k / (c*n*d^G);
	endfor
endfunction
