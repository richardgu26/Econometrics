load AuctionII.80new
printf("total rows read in: %d\n", rows(AuctionII));
test = (AuctionII(:,3)==1) & (AuctionII(:,4) ==1);
AuctionII = AuctionII(test,:);
printf("valid rows: %d\n", rows(AuctionII));

printf("number of Monte Carlo reps: %d\n", rows(AuctionII));
est = AuctionII(:,7:10);
m = mean(est);
s = std(est);
t = [0.5 0.5 0.5 0.5];
b = m - t;
rmse = sqrt(b.^2 + s.^2);

m
s
b
rmse
