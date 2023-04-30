#histogram of out of frame - out of frame frequencies for tag, tga, taa and all

library(stringr)

pdf("17.1_fig_4.pdf")
par(mfrow=c(2,2))

data1 <- read.csv('13.3_merged_OOF_CDS_rates.csv',header = TRUE)

observed_taa_frequency <- unlist(data1$frequency_taa_coding)
expected_taa_frequency <- unlist(data1$frequency_taa_OOF) 
observed_vs_expected_taa <- observed_taa_frequency - expected_taa_frequency
hist(observed_vs_expected_taa, breaks = 35, xlab = 'observed - expected', col = 'deepskyblue2', main = 'taa',xlim = c(-0.010, 0.005))


observed_tag_frequency <- unlist(data1$frequency_tag_coding)
expected_tag_frequency <- unlist(data1$frequency_tag_OOF) 
observed_vs_expected_tag <- observed_tag_frequency - expected_tag_frequency
hist(observed_vs_expected_tag, breaks = 30,xlab = 'observed - expected', col = 'deepskyblue2', main = 'tag',xlim = c(-0.010, 0.005))

observed_tga_frequency <- unlist(data1$frequency_tga_coding)
expected_tga_frequency <- unlist(data1$frequency_tga_OOF) 
observed_vs_expected_tga <- observed_tga_frequency - expected_tga_frequency
hist(observed_vs_expected_tga, breaks = 30, xlab = 'observed - expected', col = 'deepskyblue2', main = 'tga',xlim = c(-0.010, 0.005))



observed_all_frequency <- unlist(data1$frequency_ASC_coding)
expected_all_frequency <- unlist(data1$frequency_ASC_OOF) 
observed_vs_expected_all <- observed_all_frequency - expected_all_frequency
hist(observed_vs_expected_all, breaks = 30, xlab = 'observed - expected', col = 'dodgerblue3', main = 'all',xlim = c(-0.02, 0.02))


stats1 <- t.test(observed_taa_frequency,expected_taa_frequency, paired = TRUE)
print(paste('taa'))
print(stats1)

stats2 <- t.test(observed_tag_frequency,expected_tag_frequency, paired = TRUE)
print(paste('tag'))
print(stats2)


stats3 <- t.test(observed_tga_frequency,expected_tga_frequency, paired = TRUE)
print(paste('tga'))
print(stats3)


stats4 <- t.test(observed_all_frequency,expected_all_frequency, paired = TRUE)
print(paste('all'))
print(stats4)


dev.off()