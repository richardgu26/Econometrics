for i = 1:1000
	a = rand(1000,1000);
	a = inv(a'*a);
	printf("done with iteration %d, press CTRL-C to break\n", i);
	i = i + 1;
endfor