# create a matrix filled with draws from N(0,1)
x = randn(10,3);

# this is how to save a matrix to disc:
save x;

# To convince you that the load works, erase it
clear x;

# Now load it, and display
load x;
x
