#################
# Set variables #
#################
alphas = seq(0.01, 0.99, 0.01)
save = TRUE
#setwd("/dss/dsshome1/lxc0D/ge73wex3/master-thesis-tests")
setwd("C:/Users/Nils/Documents/git/master-thesis-tests")

width=8
height=6

lr_pch = 4
gr_pch = 1
ro_pch = 2
u_pch = 15
  
lr_cex = 0.7
gr_cex = 0.9
ro_cex = 0.8
u_cex = 0.7

lr_name = "LR"
gr_name = "Test 1"
#ro_name = "Test 2"
u_name = "Test 2"


cex_lab = 1.3
cex_axis = 1.3
cex_main = 1.3
cex_sub = 1.3
cex_legend = 1.3

lwd_legend = 1

###############
# p values LR #
###############

S1M20 = readRDS("results/star_tree/LR/p-values/setup=1_n=500_m=20_exp=5000.rds")
name = paste("results/final_plots/", "pval-LR-setup1-m20.pdf", sep="")
if (save){pdf(name, width=8, height=8)}
par(mar=c(7,4,4,1)+.1)
hist(S1M20, breaks=seq(0,1,0.05), xlab="p-value", ylab="", freq=FALSE, main="", 
     cex.lab=2*cex_lab, cex.axis=2*cex_axis, cex.main=2*cex_main, cex.sub=2*cex_sub, mgp=c(5,2,0)) 
if (save){dev.off()}

S1M200 = readRDS("results/star_tree/LR/p-values/setup=1_n=500_m=200_exp=5000.rds")
name = paste("results/final_plots/", "pval-LR-setup1-m200.pdf", sep="")
if (save){pdf(name, width=8, height=8)}
par(mar=c(7,4,4,1)+.1)
hist(S1M200, breaks=seq(0,1,0.05), xlab="p-value", ylab="", freq=FALSE, main="", 
     cex.lab=2*cex_lab, cex.axis=2*cex_axis, cex.main=2*cex_main, cex.sub=2*cex_sub, mgp=c(5,2,0)) 
if (save){dev.off()}

S2M20 = readRDS("results/star_tree/LR/p-values/setup=2_n=500_m=20_exp=5000.rds")
name = paste("results/final_plots/", "pval-LR-setup2-m20.pdf", sep="")
if (save){pdf(name, width=8, height=8)}
par(mar=c(7,4,4,1)+.1)
hist(S2M20, breaks=seq(0,1,0.05), xlab="p-value", ylab="", freq=FALSE, main="", 
     cex.lab=2*cex_lab, cex.axis=2*cex_axis, cex.main=2*cex_main, cex.sub=2*cex_sub, mgp=c(5,2,0)) 
if (save){dev.off()}




#######################
# Parameter selection #
#######################

############### vary B # 1 ###############
name = paste("results/final_plots/", "vary-B_1.pdf", sep="")

B3 = readRDS("results/star_tree/run-over/vary-B/2020-11-28-14-51_B=3.rds")
B4 = readRDS("results/star_tree/run-over/vary-B/2020-11-28-15-14_B=4.rds")
B5 = readRDS("results/star_tree/run-over/vary-B/2020-11-28-15-37_B=5.rds")
B6 = readRDS("results/star_tree/run-over/vary-B/2020-11-28-15-58_B=6.rds")
B7 = readRDS("results/star_tree/run-over/vary-B/2020-11-28-16-17_B=7.rds")

if (save){pdf(name, width=width,height=height)}

plot(alphas, B3, 
     xlab="Nominal level", ylab="Empirical test size", #main=title, sub=subtitle,
     type="p", pch=15, cex=0.7,  cex.lab=cex_lab, cex.axis=cex_axis, cex.main=cex_main, cex.sub=cex_sub)
lines(alphas, B4, type="p", pch=0, cex=0.7)
points(alphas, B5, type="p", pch=4, cex=0.7)
lines(alphas, B6, type="p", pch=1, cex=0.9)
lines(alphas, B7, type="p", pch=16, cex=0.9)
abline(coef = c(0,1))
legend("bottomright", legend=c("B=3", "B=4", "B=5", "B=6", "B=7"), bty="n", lwd=lwd_legend, 
       cex=cex_legend, lty=c(NA, NA, NA, NA, NA), pch=c(15, 0, 4, 1, 16))

if (save){dev.off() }


############### vary B # 2 ###############
B3 = B3[1:20]
B4 = B4[1:20]
B5 = B5[1:20]
B6 = B6[1:20]
B7 = B7[1:20]

name = paste("results/final_plots/",  "vary-B_2.pdf", sep="")

if (save){pdf(name, width=width,height=height)}

plot(alphas[1:20], B3, 
     xlab="Nominal level", ylab="Empirical test size", #main=title, sub=subtitle,
     type="p", pch=15, cex=1.5*0.7,  cex.lab=cex_lab, cex.axis=cex_axis, cex.main=cex_main, cex.sub=cex_sub, xlim = c(0,0.2), ylim = c(0,0.25))
points(alphas[1:20], B4, type="p", pch=0, cex=1.5*0.7)
points(alphas[1:20], B5, type="p", pch=4, cex=1.5*0.7)
points(alphas[1:20], B6, type="p", pch=1, cex=1.5*0.9)
points(alphas[1:20], B7, type="p", pch=16, cex=1.5*0.9)
abline(coef = c(0,1))
legend("bottomright", legend=c("B=3", "B=4", "B=5", "B=6", "B=7"), bty="n", lwd=lwd_legend, 
       cex=cex_legend, lty=c(NA, NA, NA, NA, NA), pch=c(15, 0, 4, 1, 16))

# close pdf file
if (save){dev.off()}



############### grouping # symmetric ###############
name = paste("results/final_plots/", "grouping_symmetric.pdf", sep="")

sym = readRDS("results/star_tree/grouping/symm vs. nonsym/symmetric_star-tree_setup=2_n=500_m=20.rds")
nonsym = readRDS("results/star_tree/grouping/symm vs. nonsym/nonsymmetric_star-tree_setup=2_n=500_m=20.rds")


if (save){pdf(name, width=width,height=height)}

plot(alphas, sym, 
     xlab="Nominal level", ylab="Empirical test size", #main=title, sub=subtitle,
     type="p", pch=gr_pch, cex=gr_cex,  cex.lab=cex_lab, cex.axis=cex_axis, cex.main=cex_main, cex.sub=cex_sub)
points(alphas, nonsym, type="p", pch=19, cex=gr_cex)
abline(coef = c(0,1))
legend("bottomright", legend=c("h symmetric", "h not symmetric"), bty="n", lwd=lwd_legend, 
       cex=cex_legend, lty=c(NA, NA), pch=c(gr_pch, 19))

if (save){dev.off()}


############### run-over # symmetric ###############
name = paste("results/final_plots/", "run-over_symmetric.pdf", sep="")

sym = readRDS("results/star_tree/run-over/symm vs. nonsymm/symmetric_star-tree_setup=2_n=500_m=20.rds")
nonsym = readRDS("results/star_tree/run-over/symm vs. nonsymm/nonsymmetric_star-tree_setup=2_n=500_m=20.rds")


if (save){pdf(name, width=width,height=height)}

plot(alphas, sym, 
     xlab="Nominal level", ylab="Empirical test size", #main=title, sub=subtitle,
     type="p", pch=ro_pch, cex=ro_cex,  cex.lab=cex_lab, cex.axis=cex_axis, cex.main=cex_main, cex.sub=cex_sub)
points(alphas, nonsym, type="p", pch=17, cex=ro_cex)
abline(coef = c(0,1))
legend("bottomright", legend=c("h symmetric", "h not symmetric"), bty="n", lwd=lwd_legend, 
       cex=cex_legend, lty=c(NA, NA), pch=c(ro_pch, 17))

if (save){dev.off()}



############### U-stat # vary-N ###############
name = paste("results/final_plots/", "vary-N.pdf", sep="")

N1000 = readRDS("results/star_tree/U-stat/vary-N/2020-12-04-23-25_star-tree_setup=2_N=1000.rds")
N2500 = readRDS("results/star_tree/U-stat/vary-N/2020-12-05-00-10_star-tree_setup=2_N=2500.rds")
N5000 = readRDS("results/star_tree/U-stat/vary-N/2020-12-16-22-12_star-tree_setup=2_N=5000.rds")


if (save){pdf(name, width=width,height=height)}

plot(alphas, N1000, 
     xlab="Nominal level", ylab="Empirical test size", #main=title, sub=subtitle,
     type="p", pch=0, cex=u_cex,  cex.lab=cex_lab, cex.axis=cex_axis, cex.main=cex_main, cex.sub=cex_sub)
points(alphas, N2500, type="p", pch=7, cex=u_cex)
points(alphas, N5000, type="p", pch=u_pch, cex=u_cex)
abline(coef = c(0,1))
legend("bottomright", legend=c("n=1000", "n=2500", "n=5000"), bty="n", lwd=lwd_legend, 
       cex=cex_legend, lty=c(NA, NA, NA), pch=c(0, 7, u_pch))

if (save){dev.off()}



############### U-stat # standardization ###############
name = paste("results/final_plots/", "standardization.pdf", sep="")

Upsilon = readRDS("results/star_tree/U-stat/standardization/Upsilon_star-tree_setup=2_n=500_m=20.rds")
Upsilon_h = readRDS("results/star_tree/U-stat/standardization/Upsilon_h_star-tree_setup=2_n=500_m=20.rds")


if (save){pdf(name, width=width,height=height)}

plot(alphas, Upsilon_h, 
     xlab="Nominal level", ylab="Empirical test size", #main=title, sub=subtitle,
     type="p", pch=u_pch, cex=u_cex,  cex.lab=cex_lab, cex.axis=cex_axis, cex.main=cex_main, cex.sub=cex_sub)
points(alphas, Upsilon, type="p", pch=0, cex=u_cex)
abline(coef = c(0,1))
legend("bottomright", legend=c("Upsilon", "Upsilon_h"), bty="n", lwd=lwd_legend, 
       cex=cex_legend, lty=c(NA, NA), pch=c(0, u_pch))

if (save){dev.off()}














###################
# Empirical sizes #
###################

############### sizes # setup 1 ###############
name = paste("results/final_plots/", "sizes_star-tree_setup-1_n=1000_m=20.pdf", sep="")

s1_LR = readRDS("results/star_tree/LR/sizes/2020-11-29-00-27_star-tree_setup=1_n=1000_m=20.rds")
s1_grouping = readRDS("results/star_tree/grouping/sizes/2020-11-21-11-56_star-tree_setup=1_n=1000_m=20.rds")
s1_run_over = readRDS("results/star_tree/run-over/sizes/2020-11-21-12-29_star-tree_setup=1_n=1000_m=20.rds")
s1_U_stat= readRDS("results/star_tree/U-stat/sizes/2020-12-17-17-34_star-tree_setup=1_n=1000_m=20.rds")

if (save){pdf(name, width=width,height=height)}

plot(alphas, s1_run_over, 
     xlab="Nominal level", ylab="Empirical test size", #main=title, sub=subtitle,
     type="p", pch=ro_pch, cex=ro_cex,  cex.lab=cex_lab, cex.axis=cex_axis, cex.main=cex_main, cex.sub=cex_sub)
points(alphas, s1_LR, type="p", pch=lr_pch, cex=lr_cex)
points(alphas, s1_U_stat, type="p", pch=u_pch, cex=u_cex)
points(alphas, s1_grouping, type="p", pch=gr_pch, cex=gr_cex)
abline(coef = c(0,1))
legend("bottomright", legend=c(lr_name, gr_name, ro_name, u_name), bty="n", lwd=lwd_legend, 
       cex=cex_legend, lty=c(NA, NA, NA, NA), pch=c(lr_pch, gr_pch, ro_pch, u_pch))

if (save){dev.off()}


############### sizes # setup 2 ###############
name = paste("results/final_plots/", "sizes_star-tree_setup-2_n=1000_m=20.pdf", sep="")

s2_LR = readRDS("results/star_tree/LR/sizes/2020-11-28-14-06_star-tree_setup=2_n=1000_m=20.rds")
s2_run_over = readRDS("results/star_tree/run-over/sizes/2020-11-21-14-26_star-tree_setup=2_n=1000_m=20.rds")
s2_U_stat= readRDS("results/star_tree/U-stat/sizes/2020-12-17-05-11_star-tree_setup=2_n=1000_m=20.rds")
s2_grouping = readRDS("results/star_tree/grouping/sizes/2020-11-21-13-29_star-tree_setup=2_n=1000_m=20.rds")

if (save){pdf(name, width=width,height=height)}

plot(alphas, s2_run_over, 
     xlab="Nominal level", ylab="Empirical test size", #main=title, sub=subtitle,
     type="p", pch=ro_pch, cex=ro_cex,  cex.lab=cex_lab, cex.axis=cex_axis, cex.main=cex_main, cex.sub=cex_sub)
points(alphas, s2_LR, type="p", pch=lr_pch, cex=lr_cex)
points(alphas, s2_U_stat, type="p", pch=u_pch, cex=u_cex)
points(alphas, s2_grouping, type="p", pch=gr_pch, cex=gr_cex)
abline(coef = c(0,1))
legend("bottomright", legend=c(lr_name, gr_name, ro_name, u_name), bty="n", lwd=lwd_legend, 
       cex=cex_legend, lty=c(NA, NA, NA, NA), pch=c(lr_pch, gr_pch, ro_pch, u_pch))

if (save){dev.off()}


############### sizes # caterpillar n=1000 ###############
name = paste("results/final_plots/", "sizes_cat_binary_n=1000_m=20.pdf", sep="")

LR = readRDS("results/cat_binary/LR/sizes/caterpillar_n=1000.rds")
run_over = readRDS("results/cat_binary/run-over/sizes/2020-12-13-18-02_caterpillar_n=1000.rds")
U_stat= readRDS("results/cat_binary/U-stat/sizes/2020-12-18-12-20_caterpillar_n=1000.rds")
grouping = readRDS("results/cat_binary/grouping/sizes/2020-12-13-17-01_caterpillar_n=1000.rds")

if (save){pdf(name, width=width,height=height)}

plot(alphas, run_over, 
     xlab="Nominal level", ylab="Empirical test size", #main=title, sub=subtitle,
     type="p", pch=ro_pch, cex=ro_cex,  cex.lab=cex_lab, cex.axis=cex_axis, cex.main=cex_main, cex.sub=cex_sub)
points(alphas, LR, type="p", pch=lr_pch, cex=lr_cex)
points(alphas, U_stat, type="p", pch=u_pch, cex=u_cex)
points(alphas, grouping, type="p", pch=gr_pch, cex=gr_cex)
abline(coef = c(0,1))
legend("bottomright", legend=c(lr_name, gr_name, ro_name, u_name), bty="n", lwd=lwd_legend, 
       cex=cex_legend, lty=c(NA, NA, NA, NA), pch=c(lr_pch, gr_pch, ro_pch, u_pch))

if (save){dev.off()}


############### sizes # caterpillar n=500 ###############
name = paste("results/final_plots/", "sizes_cat_binary_n=500_m=20.pdf", sep="")

run_over = readRDS("results/cat_binary/run-over/sizes/2020-12-13-17-41_caterpillar_n=500.rds")
U_stat= readRDS("results/cat_binary/U-stat/sizes/2020-12-18-04-32_caterpillar_n=500.rds")
grouping = readRDS("results/cat_binary/grouping/sizes/2020-12-13-16-56_caterpillar_n=500.rds")

if (save){pdf(name, width=width,height=height)}

plot(alphas, run_over, 
     xlab="Nominal level", ylab="Empirical test size", ylim = c(0,1),
     type="p", pch=ro_pch, cex=ro_cex,  cex.lab=cex_lab, cex.axis=cex_axis, cex.main=cex_main, cex.sub=cex_sub)
points(alphas, U_stat, type="p", pch=u_pch, cex=u_cex)
points(alphas, grouping, type="p", pch=gr_pch, cex=gr_cex)
abline(coef = c(0,1))
legend("bottomright", legend=c(gr_name, ro_name, u_name), bty="n", lwd=lwd_legend, 
       cex=cex_legend, lty=c(NA, NA, NA), pch=c(gr_pch, ro_pch, u_pch))

if (save){dev.off()}















#######################
# High dimensionality #
#######################

############### star-tree # setup 1 ###############
name = paste("results/final_plots/", "sizes_star-tree_setup-1_n=500_m=200.pdf", sep="")

LR = readRDS("results/star_tree/LR/sizes/2020-12-22-14-54_star-tree_setup=1_n=500_m=200.rds")
grouping = readRDS("results/star_tree/grouping/sizes/2020-12-22-14-51_star-tree_setup=1_n=500_m=200.rds")
run_over = readRDS("results/star_tree/run-over/sizes/2020-12-22-16-16_star-tree_setup=1_n=500_m=200.rds")
U_stat= readRDS("results/star_tree/U-stat/sizes/2020-12-22-17-24_star-tree_setup=1_n=500_m=200.rds")

if (save){pdf(name, width=width,height=height)}

plot(alphas, run_over, 
     xlab="Nominal level", ylab="Empirical test size", #main=title, sub=subtitle,
     type="p", pch=ro_pch, cex=ro_cex,  cex.lab=cex_lab, cex.axis=cex_axis, cex.main=cex_main, cex.sub=cex_sub)
points(alphas, LR, type="p", pch=lr_pch, cex=lr_cex)
points(alphas, U_stat, type="p", pch=u_pch, cex=u_cex)
points(alphas, grouping, type="p", pch=gr_pch, cex=gr_cex)
abline(coef = c(0,1))
legend("bottomright", legend=c(lr_name, gr_name, ro_name, u_name), bty="n", lwd=lwd_legend, 
       cex=cex_legend, lty=c(NA, NA, NA, NA), pch=c(lr_pch, gr_pch, ro_pch, u_pch))

if (save){dev.off()}


############### caterpillar tree ###############
name = paste("results/final_plots/", "sizes_cat_binary_n=500_m=200.pdf", sep="")

#LR = readRDS("results/cat_binary/LR/sizes/2020-12-22-14-54_star-tree_setup=1_n=500_m=200.rds")
grouping = readRDS("results/cat_binary/grouping/sizes/2020-12-22-18-06_caterpillar_n=500_m=200.rds")
run_over = readRDS("results/cat_binary/run-over/sizes/2020-12-22-23-52_caterpillar_n=500_m=200.rds")
U_stat= readRDS("results/cat_binary/U-stat/sizes/2020-12-22-18-59_caterpillar_n=500_m=200.rds")

if (save){pdf(name, width=width,height=height)}

plot(alphas, grouping, 
     xlab="Nominal level", ylab="Empirical test size", #main=title, sub=subtitle,
     type="p", pch=gr_pch, cex=gr_cex,  cex.lab=cex_lab, cex.axis=cex_axis, cex.main=cex_main, cex.sub=cex_sub)
#points(alphas, LR, type="p", pch=lr_pch, cex=lr_cex)
points(alphas, run_over, type="p", pch=ro_pch, cex=ro_cex)
points(alphas, U_stat, type="p", pch=u_pch, cex=u_cex)
abline(coef = c(0,1))
legend("bottomright", legend=c(gr_name, ro_name, u_name), bty="n", lwd=lwd_legend, 
       cex=cex_legend, lty=c(NA, NA, NA), pch=c(gr_pch, ro_pch, u_pch))

if (save){dev.off()}














###################
# Only equalities #
###################

############### caterpillar ###############
name = paste("results/final_plots/", "only_equalities_cat_binary_n=500_m=20.pdf", sep="")

run_over = readRDS("results/cat_binary/run-over/sizes/only_equalities_caterpillar_n=500_m=20.rds")
U_stat= readRDS("results/cat_binary/U-stat/sizes/only_equalities_caterpillar_n=500_m=20.rds")
grouping = readRDS("results/cat_binary/grouping/sizes/only_equalities_caterpillar_n=500_m=20.rds")

if (save){pdf(name, width=width,height=height)}

plot(alphas, run_over, 
     xlab="Nominal level", ylab="Empirical test size", #main=title, sub=subtitle,
     type="p", pch=ro_pch, cex=ro_cex,  cex.lab=cex_lab, cex.axis=cex_axis, cex.main=cex_main, cex.sub=cex_sub)
points(alphas, U_stat, type="p", pch=u_pch, cex=u_cex)
points(alphas, grouping, type="p", pch=gr_pch, cex=gr_cex)
abline(coef = c(0,1))
legend("bottomright", legend=c(gr_name, ro_name, u_name), bty="n", lwd=lwd_legend, 
       cex=cex_legend, lty=c(NA, NA, NA), pch=c(gr_pch, ro_pch, u_pch))

if (save){dev.off()}













###################
# Empirical power #
###################

############### fixed-n # star-tree setup 2 ###############
name = paste("results/final_plots/", "power-fixed-n_star-tree_setup-2.pdf", sep="")
H = seq(0.5,10,len=20)

grouping = readRDS("results/star_tree/grouping/power-fixed-n/2020-12-18-19-38_setup=2_n=500_m=20.rds")
#run_over = readRDS("results/star_tree/run-over/power-fixed-n/2020-12-19-03-33_setup=2_n=500_m=20.rds")
U_stat= readRDS("results/star_tree/U-stat/power-fixed-n/2020-12-25-17-16_setup=2_n=500_m=20.rds")


if (save){pdf(name, width=width,height=height)}

plot(H, grouping, 
     xlab="h", ylab="Empirical power", ylim = c(0,1), #main=title, sub=subtitle, 
     type="p", pch=gr_pch, cex=1.5*gr_cex,  cex.lab=cex_lab, cex.axis=cex_axis, cex.main=cex_main, cex.sub=cex_sub)
#points(H, run_over, type="p", pch=ro_pch, cex=1.5*ro_cex)
points(H, U_stat, type="p", pch=u_pch, cex=1.5*u_cex)
legend("bottomright", legend=c(gr_name, u_name), bty="n", lwd=lwd_legend, 
       cex=cex_legend, lty=c(NA, NA), pch=c(gr_pch, u_pch))

if (save){dev.off()}




############### fixed-n # star-tree setup 1 ###############
name = paste("results/final_plots/", "power-fixed-n_star-tree_setup-1.pdf", sep="")
H = seq(1,20,len=20)

LR = readRDS("results/star_tree/LR/power-fixed-n/2021-01-04-18-30_setup=1_n=500_m=20.rds")
grouping = readRDS("results/star_tree/grouping/power-fixed-n/2021-01-04-19-21_setup=1_n=500_m=20.rds")
#run_over = readRDS("results/star_tree/run-over/power-fixed-n/2021-01-05-04-05_setup=1_n=500_m=20.rds")
U_stat= readRDS("results/star_tree/U-stat/power-fixed-n/2021-01-06-14-28_setup=1_n=500_m=20.rds")


if (save){pdf(name, width=width,height=height)}

plot(H, grouping, 
     xlab="h", ylab="Empirical power", ylim = c(0,1), #main=title, sub=subtitle, 
     type="p", pch=gr_pch, cex=1.5*gr_cex,  cex.lab=cex_lab, cex.axis=cex_axis, cex.main=cex_main, cex.sub=cex_sub)
#points(H, run_over, type="p", pch=ro_pch, cex=1.5*ro_cex)
points(H, LR, type="p", pch=lr_pch, cex=1.5*lr_cex)
points(H, U_stat, type="p", pch=u_pch, cex=1.5*u_cex)
legend("topleft", legend=c(lr_name, gr_name, u_name), bty="n", lwd=lwd_legend, 
       cex=cex_legend, lty=c(NA, NA, NA), pch=c(lr_pch, gr_pch, u_pch))

if (save){dev.off()}

############### 2-factor # regular ###############
name = paste("results/final_plots/", "power-2-factor-regular.pdf", sep="")
H = seq(1,20,len=20)

LR = readRDS("results/2-factor/power/LR/vary-alternative_regular_n=500_m=20.rds")
U_stat= readRDS("results/2-factor/power/Ustat/vary-alternative_regular_n=500_m=20.rds")


if (save){pdf(name, width=width,height=height)}

plot(H, LR, 
     xlab="h", ylab="Empirical power", ylim = c(0,1), #main=title, sub=subtitle, 
     type="p", pch=lr_pch, cex=1.5*lr_cex,  cex.lab=cex_lab, cex.axis=cex_axis, cex.main=cex_main, cex.sub=cex_sub)
points(H, U_stat, type="p", pch=u_pch, cex=1.5*u_cex)
legend("topleft", legend=c(lr_name, u_name), bty="n", lwd=lwd_legend, 
       cex=cex_legend, lty=c(NA, NA), pch=c(lr_pch, u_pch))

if (save){dev.off()}

############### 2-factor # singular ###############
name = paste("results/final_plots/", "power-2-factor-singular.pdf", sep="")
H = seq(0.5,10,len=20)


LR = readRDS("results/2-factor/power/LR/vary-alternative_singular_n=500_m=20.rds")
U_stat= readRDS("results/2-factor/power/Ustat/vary-alternative_singular_n=500_m=20.rds")


if (save){pdf(name, width=width,height=height)}

plot(H, LR, 
     xlab="h", ylab="Empirical power", ylim = c(0,1), #main=title, sub=subtitle, 
     type="p", pch=lr_pch, cex=1.5*lr_cex,  cex.lab=cex_lab, cex.axis=cex_axis, cex.main=cex_main, cex.sub=cex_sub)
points(H, U_stat, type="p", pch=u_pch, cex=1.5*u_cex)
legend("topleft", legend=c(lr_name, u_name), bty="n", lwd=lwd_legend, 
       cex=cex_legend, lty=c(NA, NA), pch=c(lr_pch, u_pch))

if (save){dev.off()}




############### fixed-alternative # star-tree setup 2 ###############
name = paste("results/final_plots/", "power-fixed-alternative_star-tree_setup-2.pdf", sep="")
N = seq(250,1200,len=20)

grouping = readRDS("results/star_tree/grouping/power-fixed-alternative/2020-12-21-23-41_setup=2_m=20.rds")
run_over = readRDS("results/star_tree/run-over/power-fixed-alternative/2020-12-22-04-04_setup=2_m=20.rds")
U_stat = readRDS("results/star_tree/U-stat/power-fixed-alternative/2021-01-10-00-50_setup=2_m=20.rds")
U_stat = c(NA,NA,NA,U_stat[1],NA,NA,NA,U_stat[2],NA,NA,NA,U_stat[3],NA,NA,NA,U_stat[4],NA,NA,NA,U_stat[5])

if (save){pdf(name, width=width,height=height)}
plot(N, grouping, 
     xlab="N", ylab="Empirical power", ylim = c(0,1), #main=title, sub=subtitle,
     type="p", pch=gr_pch, cex=1.5*gr_cex,  cex.lab=cex_lab, cex.axis=cex_axis, cex.main=cex_main, cex.sub=cex_sub)
points(N, run_over, type="p", pch=ro_pch, cex=1.5*ro_cex)
points(N, U_stat, type="p", pch=u_pch, cex=1.5*u_cex)
legend("bottomright", legend=c(gr_name, ro_name, u_name), bty="n", lwd=lwd_legend, 
       cex=cex_legend, lty=c(NA, NA, NA), pch=c(gr_pch, ro_pch, u_pch))

if (save){dev.off()}



############### fixed-alternative # star-tree setup 1 ###############
name = paste("results/final_plots/", "power-fixed-alternative_star-tree_setup-1.pdf", sep="")
N = seq(250, 1200,len=20)


LR = readRDS("results/star_tree/LR/power-fixed-alternative/2021-01-04-12-28_setup=1_m=20.rds")
grouping = readRDS("results/star_tree/grouping/power-fixed-alternative/2021-01-04-13-46_setup=1_m=20.rds")
run_over = readRDS("results/star_tree/run-over/power-fixed-alternative/2021-01-05-01-42_setup=1_m=20.rds")
U_stat = readRDS("results/star_tree/U-stat/power-fixed-alternative/2021-01-09-18-45_setup=1_m=20.rds")
U_stat = c(NA,NA,NA,U_stat[1],NA,NA,NA,U_stat[2],NA,NA,NA,U_stat[3],NA,NA,NA,U_stat[4],NA,NA,NA,U_stat[5])

if (save){pdf(name, width=width,height=height)}

plot(N, grouping, 
     xlab="N", ylab="Empirical power", ylim = c(0,1), #main=title, sub=subtitle,
     type="p", pch=gr_pch, cex=1.5*gr_cex,  cex.lab=cex_lab, cex.axis=cex_axis, cex.main=cex_main, cex.sub=cex_sub)
points(N, LR, type="p", pch=lr_pch, cex=1.5*lr_cex)
points(N, run_over, type="p", pch=ro_pch, cex=1.5*ro_cex)
points(N, U_stat, type="p", pch=u_pch, cex=1.5*u_cex)
legend(x=1000, y=0.45, legend=c(lr_name, gr_name, ro_name, u_name), bty="n", lwd=lwd_legend, 
       cex=cex_legend, lty=c(NA, NA, NA, NA), pch=c(lr_pch, gr_pch, ro_pch, u_pch))

if (save){dev.off()}




############### compute-cov # sizes # setup 2 ###############
name = paste("results/final_plots/", "compute-cov_sizes_star-tree_setup-2_n=500_m=10.pdf", sep="")


run_over = readRDS("results/star_tree/run-over/sizes/only_equalities_star-tree_setup=2_n=500_m=10.rds")
run_over_compute_cov = readRDS("results/star_tree/run-over-cov/sizes/2021-01-27-12-40_star-tree_setup=2_n=500_m=10.rds")

if (save){pdf(name, width=width,height=height)}

plot(alphas, run_over, 
     xlab="Nominal level", ylab="Empirical test size", #main=title, sub=subtitle,
     type="p", pch=ro_pch, cex=ro_cex,  cex.lab=cex_lab, cex.axis=cex_axis, cex.main=cex_main, cex.sub=cex_sub)

points(alphas, run_over_compute_cov, type="p", pch=17, cex=ro_cex)
abline(coef = c(0,1))
legend("bottomright", legend=c("batched-mean estimator", "g(S)"), bty="n", lwd=lwd_legend, 
       cex=cex_legend, lty=c(NA, NA), pch=c(ro_pch, 17))

if (save){dev.off()}



############### compute-cov # power # setup 2 ###############
name = paste("results/final_plots/", "compute-cov_power_star-tree_setup-2_n=500_m=10.pdf", sep="")
H = seq(0.5,10,len=20)


run_over = readRDS("results/star_tree/run-over/power-fixed-n/2021-01-27-15-21_setup=2_n=500_m=10.rds")
run_over_compute_cov = readRDS("results/star_tree/run-over-cov/power-fixed-n/2021-01-27-14-36_setup=2_n=500_m=10.rds")

if (save){pdf(name, width=width,height=height)}

plot(H, run_over, 
     xlab="Nominal level", ylab="Empirical power", #main=title, sub=subtitle,
     type="p", pch=ro_pch, cex=1.5*ro_cex,  cex.lab=cex_lab, cex.axis=cex_axis, cex.main=cex_main, cex.sub=cex_sub)

points(H, run_over_compute_cov, type="p", pch=17, cex=1.5*ro_cex)
legend("bottomright", legend=c("batched-mean estimator", "g(S)"), bty="n", lwd=lwd_legend, 
       cex=cex_legend, lty=c(NA, NA), pch=c(ro_pch, 17))

if (save){dev.off()}








