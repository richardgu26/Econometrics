% this illustrates the Cholesky decomposition

% create a random p.d. covariance matrix
V = randn(3,3);
V = V'*V; % make it pos. def.
printf("true covariance matrix\n");
disp(V);
printf("press enter to continue\n");
pause;

% take Cholesky decomp.
P = chol(V);
printf("P=chol(V): P'P-V\n");
disp(P'*P - V);
printf("press enter to continue\n");
pause;


% how to sample from N(0,V)
e = randn(1000,3)*P;
printf("sample covariance matrix of 1000 draws from N(0,V)\n");
disp(cov(e));
printf("press enter to continue\n");
pause;


% verify that GLS transformation works
P = chol(inv(V));
printf("P=chol(inv(V)): P*V*P'\n");
disp(P*V*P');
printf("press enter to continue\n");
pause;

e = e*P';
printf("sample covariance matrix of transformed errors\n");
disp(cov(e));



