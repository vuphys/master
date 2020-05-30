# LOAD DATASETS PACKAGES ###################################

library(datasets)  # Load/unload base packages manually
library(readr)


# LOAD DATA ################################################

my_data <- read.table("DATA/da_ana/datFft_10_1000_2020-05-27T060716.txt", header=T,sep="\t")


my_data


# PLOTS ####################################################

# Good to first check univariate distributions
hist(my_data$img)
hist(my_data$mean)

# Basic X-Y plot for two quantitative variables
with(my_data[my_data$img==1,], plot(my_data$PSNR, my_data$SSIM))

# Add some options
with(my_data[my_data$img==1,],plot(my_data$alpha_0, my_data$VIF,
     pch = 0.1,         # Solid circle
     cex = 0.1,        # Make 150% size
     col = "#cc0000"))  # Red
     #main = "alpha_1 vs SSIM",
     #xlab = "alpha_1 value",
     #ylab = "SSIM"

# withouth constrain
plot(my_data$alpha_1, my_data$FSIM,
            pch = 0.05,         # Solid circle
            cex = 0.05,
          col = "#cc0000") # Make 150% size
            #main = "alpha_1 vs SSIM",
#xlab = "alpha_1 value",
#ylab = "SSIM"
install.packages('latex2exp')

library(ggplot2) #plot in R
library(latex2exp) #Latex in R
library(plotly) #plot 3d


# Test code:
ggplot(my_data, aes(x=alpha_1, y=SSIM))
ggplot(my_data, aes(x=alpha_1, y=SSIM, color=img)) + geom_point() + geom_smooth()

ggplot(my_data) + geom_point(aes(x=alpha_1, y=SSIM, color=img)) + geom_smooth(aes(x=alpha_1, y=SSIM, color=img)) 

ggplot(my_data, aes(x=alpha_1, y=SSIM, color=img)) + geom_point()

# CODE USE FOR PLOT!!!!
gg <- ggplot(my_data, aes(x=alpha_0, y=VIF,
        color=img)) + geom_point() + labs(title=TeX("$\\alpha_0$ vs VIF"), x=TeX("$\\alpha_0$"), y="VIF")
print(gg)
gg + facet_wrap(std~img) #divide in several plots

# 3d plot
plot_ly(x=my_data$alpha_1, y=my_data$alpha_0, z=my_data$MS_SSIM,
        type="scatter3d",
        mode="markers", 
        color=my_data$img,
        marker = list(size = 2))

plot(my_data)
# CLEAN UP #################################################

# Clear packages
detach("package:datasets", unload = TRUE)  # For base

# Clear plots
dev.off()  # But only if there IS a plot

# Clear console
cat("\014")  # ctrl+L

# Clear mind :)