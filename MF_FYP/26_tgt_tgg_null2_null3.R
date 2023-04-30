
#histogram of coding - non-coding frequencies for tgg

library(stringr)
library(stats)


pdf("26.1_tgt_tgg_appendix2.pdf")
par(mfrow=c(2,2))

#TGG OOF null model (3)

data <- read.csv('25.1_tgt_tgg_data.csv',header = TRUE)

observed_tgg_frequency_1 <- unlist(data$frequency_tgg_coding)
expected_tgg_frequency_1 <- unlist(data$frequency_tgg_OOF)

stats_tgg_1 <- t.test(observed_tgg_frequency_1,expected_tgg_frequency_1, paired = TRUE)

observed_vs_expected_tgg_1 <- observed_tgg_frequency_1 - expected_tgg_frequency_1
hist(observed_vs_expected_tgg_1, breaks = 30, xlab = 'observed - expected', col = 'deepskyblue2', main = 'tgg: OOF null model',xlim = c(-0.01, 0.01))

print(paste('tgg: OOF null model'))
print(stats_tgg_1)


#TGG NCDS null model (2)


observed_tgg_frequency_2 <- unlist(data$frequency_tgg_coding)
expected_tgg_frequency_2 <- unlist(data$frequency_tgg_non_coding)


stats_tgg_2 <- t.test(observed_tgg_frequency_2,expected_tgg_frequency_2, paired = TRUE)

observed_vs_expected_tgg_2 <- observed_tgg_frequency_2 - expected_tgg_frequency_2

hist(observed_vs_expected_tgg_2, breaks = 60, xlab = 'observed - expected', col = 'deepskyblue2', main = 'tgg: NCDS null model',xlim = c(-0.05, 0.05))

print(paste('tgg: NCDS null model'))
print(stats_tgg_2)


#TGT OOF null model (3)


observed_tgt_frequency_1 <- unlist(data$frequency_tgt_coding)
expected_tgt_frequency_1 <- unlist(data$frequency_tgt_OOF)

stats_tgt_1 <- t.test(observed_tgt_frequency_1,expected_tgt_frequency_1, paired = TRUE)

observed_vs_expected_tgt_1 <- observed_tgt_frequency_1 - expected_tgt_frequency_1
hist(observed_vs_expected_tgt_1, breaks = 30, xlab = 'observed - expected', col = 'deepskyblue2', main = 'tgt: OOF null model',xlim = c(-0.01, 0.01))

print(paste('tgt: OOF null model'))
print(stats_tgt_1)

                
#TGT NCDS null model (2)


observed_tgt_frequency_2 <- unlist(data$frequency_tgt_coding)
expected_tgt_frequency_2 <- unlist(data$frequency_tgt_non_coding)


stats_tgt_2 <- t.test(observed_tgt_frequency_2,expected_tgt_frequency_2, paired = TRUE)

observed_vs_expected_tgt_2 <- observed_tgt_frequency_2 - expected_tgt_frequency_2
hist(observed_vs_expected_tgt_2, breaks = 30, xlab = 'observed - expected', col = 'deepskyblue2', main = 'tgt: NCDS null model',xlim = c(-0.05, 0.05))

print(paste('tgt: NCDS null model'))
print(stats_tgt_2)




dev.off()