# LOAD DATASETS PACKAGES ###################################

library(datasets)  # Load/unload base packages manually
library(readr)


# LOAD DATA ################################################

my_data <- read.table("DATA/datCondat_10_10.txt", header=T,sep=",")


my_data


# PLOTS ####################################################

# Good to first check univariate distributions
hist(my_data$img)
hist(my_data$mean)

# Basic X-Y plot for two quantitative variables
plot(my_data$alpha_1, my_data$alpha_0)

# Add some options
plot(my_data$alpha_0, my_data$PSNR,
     pch = 19,         # Solid circle
     cex = 1.5,        # Make 150% size
     col = "#cc0000",  # Red
     main = "\alpha_0 vs PSNR",
     xlab = "\alpha_0 value",
     ylab = "PSNR")

plot(my_data)
# CLEAN UP #################################################

# Clear packages
detach("package:datasets", unload = TRUE)  # For base

# Clear plots
dev.off()  # But only if there IS a plot

# Clear console
cat("\014")  # ctrl+L

# Clear mind :)