#scatter graph of Z scores (1000 sim) vs GC3 content of genomes
library(stringr)

pdf("15.1_fig_2.pdf")


data_file <- read.csv('13.1_merged_simulation_rates.csv',header = TRUE)

observed_tag <- unlist(data_file$observed_tag)
observed_tga <- unlist(data_file$observed_tga)
observed_taa <- unlist(data_file$observed_taa)

total_observed_ASC <- observed_tag + observed_tga +  observed_taa

mean_tag <- unlist(data_file$mean_x)
mean_tga <- unlist(data_file$mean_y)
mean_taa <- unlist(data_file$mean)

 mean_expected_ASC <- mean_tag + mean_tga + mean_taa

st_dev_tag <- unlist(data_file$st_dev_x)
st_dev_tga <- unlist(data_file$st_dev_y)
st_dev_taa <- unlist(data_file$st_dev)

st_dev_ASC <- st_dev_tag + st_dev_tga + st_dev_taa 

z_all <- ((total_observed_ASC - mean_expected_ASC )/ st_dev_ASC)
x <- unlist(data_file$GC3)
y <- unlist(z_all)

plot(x, y,col = 'deepskyblue1',  xlab = "GC3 Content", ylab = "Z Score", pch = 19, frame = FALSE,ylim = c(-10, 5))
    
# Add regression line

abline(lm(y ~ x, data = data_file), col = "black", lwd = 3)

spearmans_test <- cor.test(x, y,source('~/Docs/Python_Script/Maya_Frost_FYP/15_GC3_vs_Z_scatter.R'), chdir = TRUE,
         method = "spearman")

print(spearmans_test)

dev.off()
