# Read in the data

load meps1996.data;

#	The variables and their columns
#	VARIABLE						column
# 	Office based doctor visits	    1
# 	Outpatient doctor visits		2
# 	Emergency room visits			3
# 	Inpatient visits				4
# 	Dental visits					5
# 	Prescriptions					6
#  	PUBLIC_INS  					7
# 	PRIVATE_INS						8
# 	SEX								9
# 	AGE								10
# 	EDUC							11
# 	INCOME							12


names = char("OBDV", "OPV","IPV","ERV","DV","RX","PUB","PRIV","SEX","AGE","EDUC","INC");
n = rows(data);
printf("MEPS data, 1996, complete data set statistics\n");
printf("%d observations\n",n);
dstats(data, names);
