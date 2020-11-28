setwd("C:/Users/Nils/Documents/Studium/TUM/09_Masterarbeit/master-thesis-tests")

date = "2020-11-28"
save = TRUE

#################
# Sizes setup 1 #
#################

name = paste("results/", date, "_star-tree_setup-1_n=1000_m=20.pdf", sep="")

# Load data
#s1_grouping = readRDS("results/grouping/sizes/2020-11-21-13-17_star-tree_setup=1_n=1000_m=20.rds")
s1_run_over = readRDS("results/run-over/sizes/2020-11-21-12-29_star-tree_setup=1_n=1000_m=20.rds")
s1_symmetric = readRDS("results/symmetric/sizes/2020-11-21-11-56_star-tree_setup=1_n=1000_m=20.rds")

# Set variables
title = "Emprical test sizes vs. nominal test levels for different test strategies. \n Star tree - setup 1"
subtitle = "500 experiments with m = 20 and n = 1000."
alphas = seq(0.01, 0.99, 0.01)


# create pdf file
if (save){
  pdf(name) 
}

# plot (TODO: add incomplete U-stat)
plot(alphas, s1_grouping, 
     xlab="Nominal level", ylab="Emprical test size", main=title, sub=subtitle,
     type="p", pch=1)
points(alphas, s1_run_over, type="p", pch=16)
points(alphas, s1_symmetric, type="p", pch=18)
abline(coef = c(0,1))
legend("topleft", legend = c("grouping", "run-over", "symmetric"), bty = "n", lwd = 1, 
       cex = 1.2, lty = c(NA, NA, NA), pch = c(1, 16, 18))

# close pdf file
if (save){
  dev.off() 
}


#################
# Sizes setup 2 #
#################

name = paste("results/", date, "_star-tree_setup-2_n=1000_m=20.pdf", sep="")

# Load data
s2_factanal = readRDS("results/factanal/sizes/2020-11-28-14-06_star-tree_setup=2_n=1000_m=20.rds")
s2_run_over = readRDS("results/run-over/sizes/2020-11-21-14-26_star-tree_setup=2_n=1000_m=20.rds")
s2_U_stat= readRDS("results/U-stat/sizes/2020-11-24-09-52_star-tree_setup=2_n=1000_m=20.rds")
s2_symmetric = readRDS("results/symmetric/sizes/2020-11-21-13-29_star-tree_setup=2_n=1000_m=20.rds")

# Set variables
title = "Emprical test sizes vs. nominal test levels for different test strategies. \n Star tree - setup 2"
subtitle = "500 experiments with m = 20 and n = 1000."
alphas = seq(0.01, 0.99, 0.01)


# create pdf file
if (save){
  pdf(name) 
}

# plot (TODO: add incomplete U-stat)
plot(alphas, s2_run_over, 
     xlab="Nominal level", ylab="Emprical test size", main=title, sub=subtitle,
     type="p", pch=16)
points(alphas, s2_factanal, type="p", pch=1)
points(alphas, s2_U_stat, type="p", pch=3)
points(alphas, s2_symmetric, type="p", pch=18)
abline(coef = c(0,1))
legend("bottomright", legend = c("likelihood-ratio", "run-over", "incomplete U-statistic", "grouping"), bty = "n", lwd = 1, 
       cex = 1.2, lty = c(NA, NA, NA, NA), pch = c(1, 16, 3, 18))

# close pdf file
if (save){
  dev.off() 
}


