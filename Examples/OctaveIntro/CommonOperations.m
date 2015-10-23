echo on # this causes the commands to be echoed to the output, 
		# so you can see what commands generate what output

# vector of ones
a = ones(10,1); 
a

# vector of zeros
b = zeros(10,1);
b

# element-by-element multiply two vectors
c = a .* b;
c


# create a increasing sequence (this creates a row vector)
d = 1:10;
d

# put vectors next to eachother to define a matrix
a = [a b d']; # square brackets are used to define matrices
a

# stack two matrices on top of one another (note that Octave is case sensitive)
A = ones(2,2);
B = rand(2,2); # these are U(0,1) draws
C = [A; B]; # the semicolon between A and B means B goes under A
			# with a space or a comma it would go beside A
C


# diagonal multiply
# (multiply each element in a row of a matrix by the corresponding
# element of a vector)
# this is useful when forming moment conditions that are errors
# multiplied by instruments
e =  a.*d';
e


# inversion
A = [5,1;1,5]
B = inverse(A)
A*B

# a for loop, and incidentally, exponentiation
for i = 1:10
	10^i
endfor	

# calculate OLS coefficients
y = rand(100,1);
x = [ones(100,1) rand(100,1)];
beta = inverse(x'*x)*x'*y


# test whether some logical statements hold
a = randn(10,1); # N(0,1) random draws
b = a > 0; # element by element comparison of vectors
[a b]

c = C == 1; # test equality (element by element matrix compare to scalar)
c

# reshape a matrix
D = 1:25;
D = reshape(D,5,5);
D

# extract elements out of a matrix
E = D(1:3,2:5);
E
E = diag(D); # useful for calculating t-statistics
E

# lag a time series (missing elements filled with 1's)
E = lag(E,1) # this is another of my utilities
