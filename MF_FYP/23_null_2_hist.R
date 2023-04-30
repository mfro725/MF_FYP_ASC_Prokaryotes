
#histogram of coding - non-coding frequencies for tag, tga, taa and all

library(stringr)
library(stats)


pdf("23.1_fig_3.pdf")
par(mfrow=c(2,2))


data2 <- read.csv('13.2_merged_NCDS_CDS_rates.csv',header = TRUE)

observed_all_frequency <- unlist(data2$frequency_tag_coding)
expected_all_frequency <- unlist(data2$frequency_tag_non_coding) 

x = observed_all_frequency - expected_all_frequency

stats_tag <- t.test(observed_all_frequency,expected_all_frequency, paired = TRUE)


observed_all_frequency <- unlist(data2$frequency_tga_coding)
expected_all_frequency <- unlist(data2$frequency_tga_non_coding) 
x = observed_all_frequency - expected_all_frequency
stats_tga <- t.test(observed_all_frequency,expected_all_frequency, paired = TRUE)


observed_all_frequency <- unlist(data2$frequency_taa_coding)
expected_all_frequency <- unlist(data2$frequency_taa_non_coding) 
x = observed_all_frequency - expected_all_frequency
stats_taa <- t.test(observed_all_frequency,expected_all_frequency, paired = TRUE)


observed_all_frequency <- unlist(data2$frequency_ASC_coding)
expected_all_frequency <- unlist(data2$frequency_ASC_non_coding) 
x = observed_all_frequency - expected_all_frequency
stats <- t.test(observed_all_frequency,expected_all_frequency, paired = TRUE)
p_value = round((stats$p.value),3)
p = as.character(p_value)
label = paste('p =', p)


data1 <- read.csv('13.2_merged_NCDS_CDS_rates.csv',header = TRUE)

observed_taa_frequency <- unlist(data1$frequency_taa_coding)
expected_taa_frequency <- unlist(data1$frequency_taa_non_coding) 
observed_vs_expected_taa <- observed_taa_frequency - expected_taa_frequency
hist(observed_vs_expected_taa, breaks = 30, xlab = 'observed - expected', col = 'deepskyblue2', main = 'taa',xlim = c(-0.05, 0.05))

observed_tag_frequency <- unlist(data1$frequency_tag_coding)
expected_tag_frequency <- unlist(data1$frequency_tag_non_coding) 
observed_vs_expected_tag <- observed_tag_frequency - expected_tag_frequency
hist(observed_vs_expected_tag, breaks = 10,xlab = 'observed - expected', col = 'deepskyblue2', main = 'tag',xlim = c(-0.05, 0.05))


observed_tga_frequency <- unlist(data1$frequency_tga_coding)
expected_tga_frequency <- unlist(data1$frequency_tga_non_coding) 
observed_vs_expected_tga <- observed_tga_frequency - expected_tga_frequency
hist(observed_vs_expected_tga, breaks = 30, xlab = 'observed - expected', col = 'deepskyblue2', main = 'tga',xlim = c(-0.05, 0.05))



observed_all_frequency <- unlist(data1$frequency_ASC_coding)
expected_all_frequency <- unlist(data1$frequency_ASC_non_coding) 
observed_vs_expected_all <- observed_all_frequency - expected_all_frequency
hist(observed_vs_expected_all, breaks = 30, xlab = 'observed - expected', col = 'dodgerblue3', main = 'all',xlim = c(-0.10, 0.05))

print(paste('tag'))
print(stats_tag)
print(paste('tga'))
print(stats_tga)
print(paste('taa'))
print(stats_taa)
print(paste('all'))
print(stats)
                    
mtext(label,          
      side = 4,
      cex = 0,5)


dev.off()