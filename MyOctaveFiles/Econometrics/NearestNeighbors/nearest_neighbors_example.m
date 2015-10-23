a = reshape(1:10,2,5)';
b = rand(100,2)*10;
k = 3;
ind = nearest_neighbors(a, b, 3);

for i = 1:rows(a)
	printf("target %d:\n\n", i);
	disp(a(i,:));
	printf("neighbors\n");
	disp(b(ind(i,:),:));
endfor
