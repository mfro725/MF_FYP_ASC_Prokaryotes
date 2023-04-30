
library(stringr)
library(stats)

file_list <- list.files(path=".",pattern = "\\_rates_1.csv$" )

pdf("21.1_fig_1.pdf")
par(mfrow=c(2,2))

for (fl in file_list) {
	s.fl <- str_replace(fl, '_rates_1.csv', '')
	stop_codon <- gsub('.*_', '',s.fl)
	
	data1 <- read.csv(fl,header =TRUE)
	z <- unlist(data1$Z_score)

	hist(z, xlab = 'Z score', col = 'deepskyblue2', main = stop_codon, xlim = c(-15, 15))

positive_success <- length(z[which(z>0)])
num_trials <- length(z)
stats1 <- binom.test(positive_success,num_trials,p=0.5, alternative = c("two.sided"))
print(stop_codon)
print(stats1)

}



data_file <- read.csv('13.1_merged_simulation_rates.csv',header = TRUE)

observed_tag <- unlist(data_file$observed_tag)
observed_tga <- unlist(data_file$observed_tga)
observed_taa <- unlist(data_file$observed_taa)

total_observed_ASC <- observed_tag + observed_tga +  observed_taa

mean_tag <- unlist(data_file$mean_x)
mean_tga <- unlist(data_file$mean_y)
mean_taa <- unlist(data_file$mean)

mean_expected_ASC <- mean_tag + mean_tga + mean_taa 

st_tag <- unlist(data_file$st_dev_x)
st_tga <- unlist(data_file$st_dev_y)
st_taa <- unlist(data_file$st_dev)

st_dev_ASC <- st_tag + st_tga + st_taa 

z_all <- ((total_observed_ASC - mean_expected_ASC )/ st_dev_ASC)

hist(z_all, xlab = 'Z score', col = 'dodgerblue3', main = 'all', xlim = c(-10, 10), breaks = 10)

positive_success_all <- length(z_all[which(z_all>0)])
num_trials_all <- length(z_all)
stats2 <- binom.test(positive_success_all,num_trials_all,p=0.5, alternative = c("two.sided"))

positive_success_all_2 <- length(z_all[which(z_all>1.96)])
num_trials_all <- length(z_all)
stats3 <- binom.test(positive_success_all_2,num_trials_all,p=0.5, alternative = c("two.sided"))
print(paste("binom test assumes 50:50"))
print(stats2)
print(paste("binom test assumes 80% not deviated"))
print(stats3)

dev.off()

