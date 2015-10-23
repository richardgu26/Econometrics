
if not(MPI_Initialized) MPI_Init; end
CW = MPI_Comm_Load("NEWORLD");
node = MPI_Comm_rank(CW);
nodes = MPI_Comm_size(CW);
mytag = 48;

a = reshape(1:10,2,5)';
b = rand(1000000,2)*10;
k = 1;
ind = nearest_neighbors(a, b, 3);
node

for i = 1:rows(a)
	printf("target %d:\n\n", i);
	disp(a(i,:));
	printf("neighbors\n");
	disp(b(ind(i,:),:));
endfor
if !MPI_Finalized MPI_Finalize; endif
