setwd("C:/Users/Nils/Documents/Studium/TUM/09_Masterarbeit/master-thesis-tests")

save = TRUE

#################
# Sizes setup 1 #
#################

name = paste("results/FINAL-PLOTS/", "sizes_star-tree_setup-1_n=1000_m=20.pdf", sep="")

# Load data
s1_LR = readRDS("results/star_tree/LR/sizes/2020-11-29-00-27_star-tree_setup=1_n=1000_m=20.rds")
s1_grouping = readRDS("results/star_tree/grouping/sizes/2020-11-21-11-56_star-tree_setup=1_n=1000_m=20.rds")
s1_run_over = readRDS("results/star_tree/run-over/sizes/2020-11-21-12-29_star-tree_setup=1_n=1000_m=20.rds")
s1_U_stat= readRDS("results/star_tree/U-stat/sizes/2020-12-17-17-34_star-tree_setup=1_n=1000_m=20.rds")


# Set variables
alphas = seq(0.01, 0.99, 0.01)


# create pdf file
if (save){
  pdf(name, width=8,height=6) 
}


plot(alphas, s1_run_over, 
     xlab="Nominal level", ylab="Emprical test size", #main=title, sub=subtitle,
     type="p", pch=2, cex=0.8,  cex.lab=1.3, cex.axis=1.3, cex.main=1.3, cex.sub=1.3)
points(alphas, s1_LR, type="p", pch=4, cex=0.7)
points(alphas, s1_U_stat, type="p", pch=15, cex=0.7)
points(alphas, s1_grouping, type="p", pch=1, cex=0.9)
abline(coef = c(0,1))
legend("bottomright", legend = c("LR", "Test 1", "Test 2", "Test 3"), bty = "n", lwd = 1, 
       cex = 1.3, lty = c(NA, NA, NA, NA), pch = c(4, 1, 2, 15))

# close pdf file
if (save){
  dev.off() 
}


#################
# Sizes setup 2 #
#################

name = paste("results/FINAL-PLOTS/", "sizes_star-tree_setup-2_n=1000_m=20.pdf", sep="")

# Load data
s2_LR = readRDS("results/star_tree/LR/sizes/2020-11-28-14-06_star-tree_setup=2_n=1000_m=20.rds")
s2_run_over = readRDS("results/star_tree/run-over/sizes/2020-11-21-14-26_star-tree_setup=2_n=1000_m=20.rds")
s2_U_stat= readRDS("results/star_tree/U-stat/sizes/2020-11-24-09-52_star-tree_setup=2_n=1000_m=20.rds")
s2_grouping = readRDS("results/star_tree/grouping/sizes/2020-11-21-13-29_star-tree_setup=2_n=1000_m=20.rds")

# Set variables
alphas = seq(0.01, 0.99, 0.01)


# create pdf file
if (save){
  pdf(name, width=8,height=6) 
}


plot(alphas, s2_run_over, 
     xlab="Nominal level", ylab="Emprical test size", #main=title, sub=subtitle,
     type="p", pch=2, cex=0.8,  cex.lab=1.3, cex.axis=1.3, cex.main=1.3, cex.sub=1.3)
points(alphas, s2_LR, type="p", pch=4, cex=0.7)
points(alphas, s2_U_stat, type="p", pch=15, cex=0.7)
points(alphas, s2_grouping, type="p", pch=1, cex=0.9)
abline(coef = c(0,1))
legend("bottomright", legend = c("LR", "Test 1", "Test 2", "Test 3"), bty = "n", lwd = 1, 
       cex = 1.3, lty = c(NA, NA, NA, NA), pch = c(4, 1, 2, 15))

# close pdf file
if (save){
  dev.off() 
}




##########
# vary B #
##########

# Set variables
alphas = seq(0.01, 0.99, 0.01)

# create pdf file
name = paste("results/FINAL-PLOTS/", "vary-B_1.pdf", sep="")
if (save){
  pdf(name, width=8,height=6) 
}

plot(alphas, B3, 
     xlab="Nominal level", ylab="Emprical test size", #main=title, sub=subtitle,
     type="p", pch=15, cex=0.7,  cex.lab=1.3, cex.axis=1.3, cex.main=1.3, cex.sub=1.3)
points(alphas, B4, type="p", pch=0, cex=0.7)
points(alphas, B5, type="p", pch=4, cex=0.7)
points(alphas, B6, type="p", pch=1, cex=0.9)
points(alphas, B7, type="p", pch=16, cex=0.9)
abline(coef = c(0,1))
legend("bottomright", legend = c("B=3", "B=4", "B=5", "B=6", "B=7"), bty = "n", lwd = 1, 
       cex = 1.3, lty = c(NA, NA, NA, NA, NA), pch = c(15, 0, 4, 1, 16))

# close pdf file
if (save){
  dev.off() 
}


alphas = seq(0.01, 0.2, 0.01)
name = paste("results/FINAL-PLOTS/",  "vary-B_2.pdf", sep="")

# Load data
B3 = readRDS("results/star_tree/run-over/vary-B/2020-11-28-14-51_B=3.rds")[1:20]
B4 = readRDS("results/star_tree/run-over/vary-B/2020-11-28-15-14_B=4.rds")[1:20]
B5 = readRDS("results/star_tree/run-over/vary-B/2020-11-28-15-37_B=5.rds")[1:20]
B6 = readRDS("results/star_tree/run-over/vary-B/2020-11-28-15-58_B=6.rds")[1:20]
B7 = readRDS("results/star_tree/run-over/vary-B/2020-11-28-16-17_B=7.rds")[1:20]

# create pdf file
if (save){
  pdf(name, width=8,height=6) 
}

plot(alphas, B3, 
     xlab="Nominal level", ylab="Emprical test size", #main=title, sub=subtitle,
     type="p", pch=15, cex=1.5*0.7,  cex.lab=1.3, cex.axis=1.3, cex.main=1.3, cex.sub=1.3, xlim = c(0,0.2), ylim = c(0,0.25))
points(alphas, B4, type="p", pch=0, cex=1.5*0.7)
points(alphas, B5, type="p", pch=4, cex=1.5*0.7)
points(alphas, B6, type="p", pch=1, cex=1.5*0.9)
points(alphas, B7, type="p", pch=16, cex=1.5*0.9)
abline(coef = c(0,1))
legend("bottomright", legend = c("B=3", "B=4", "B=5", "B=6", "B=7"), bty = "n", lwd = 1, 
       cex = 1.3, lty = c(NA, NA, NA, NA, NA), pch = c(15, 0, 4, 1, 16))

# close pdf file
if (save){
  dev.off() 
}
