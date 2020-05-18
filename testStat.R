# LOAD DATASETS PACKAGES ###################################

library(datasets)  # Load/unload base packages manually
library(readr)


# LOAD DATA ################################################

my_data <- read.table("DATA/datFft_50_5000.txt", header=T,sep="\t")


my_data


# PLOTS ####################################################

# Good to first check univariate distributions
hist(my_data$img)
hist(my_data$mean)

# Basic X-Y plot for two quantitative variables
with(my_data[my_data$img==1,], plot(my_data$PSNR, my_data$SSIM))

# Add some options
with(my_data[my_data$img==1,],plot(my_data$alpha_0, my_data$PSNR,
     pch = 0.1,         # Solid circle
     cex = 0.1,        # Make 150% size
     col = "#cc0000"))  # Red
     #main = "alpha_1 vs SSIM",
     #xlab = "alpha_1 value",
     #ylab = "SSIM"

plot(my_data)
# CLEAN UP #################################################

# Clear packages
detach("package:datasets", unload = TRUE)  # For base

# Clear plots
dev.off()  # But only if there IS a plot

# Clear console
cat("\014")  # ctrl+L

# Clear mind :)