MergeShifts <- function (target) {
        library(tidyverse)
        library(readr)
        
        # get working directory to reset at end of function
        directory <- getwd()
        
        # generate list of rings to analyze (gets all *.csv files)
        ringList <- list.files(directory, pattern = ".csv", recursive = FALSE)
        
        # add data to data frame corresponding for each ring in rings
        df <- data.frame()
        for (i in ringList) {
                tmp <- read_csv(i)
                df <- rbind(df, tmp)
        }

        # saves aggregated data with name_allRings.csv
        df <- filter(df, Target != "thermal")
        df <- filter(df, Target != "0.0")
        write_csv(df, "IdoNOTevenGiveaShit.csv")
        
        # returns working directory to top level
        setwd(directory)
}

library(drc)
dat <- read_csv("IdoNOTevenGiveaShit.csv")
dat <- filter(dat, Target != 0.0)
dat <- filter(dat, Target != 1.6)

fit <- drm(formula = dat$`Net Shift`~dat$Target, fct = L.4())
fit.loess <- glm(dat$`Net Shift`~dat$Target)
plot(fit)
plot(fit.glm)
summary(fit.loess)

plot <- ggplot(dat, aes(y = `Net Shift`, x = Target)) +
        geom_point() + 
        scale_x_log10(breaks = c(1, 10, 100, 1000, 10000)) +
        geom_smooth(method = "drm", fct = L.4(), se = FALSE) +
        theme_few(base_size = 16) +
        labs(x = "Target Concentration (pg/mL)",
             y = expression(paste("Relative Shift ( ",Delta,"pm)"))) +
        ggtitle("CCL8 Calbration")


plot
ggsave(plot, filename = "CCL8.png", width = 9, height = 6)

plot + geom_smooth()
plot + geom_smooth(method = glm, method.args = list(family = "binomial"))
plot



FittingPlot <- function () {
        library(tidyverse)
        library(RColorBrewer)
        library(ggthemes)
        
        #Read Data
        dat_shifts <- read.csv("IdoNOTevenGiveaShit.csv", header=TRUE)
        
        #configure plot and legend
        A1_start <- -300
        A2_start <- 10000
        x0_start <- 500
        p_start <- 0.5
        plots <- ggplot(dat_shifts, aes(x = Target, y = NetShift)) +
                theme_few() +
                ylab(expression(paste("Relative Shift (",Delta,"pm)"))) +
                xlab("Concentration (pg/mL)")+
                scale_x_log10() +
                geom_smooth(method = "nls", se = FALSE,size = 1, 
                            fullrange = TRUE, 
                            mapping = aes(x = dat_shifts$Target, 
                                          y = dat_shifts$NetShift))
        
        plots
        #save
        ggsave(file = "NetShifts.png",plots, width = 10, height = 6)
}

GetValues <- function(){
        
        # loads relevant libraries
        library(dplyr)
        library(tidyr)
        
        
        #Read Data 1
        dat_shifts <- read.csv("IdoNOTevenGiveaShit.csv", header=TRUE)
        x <- dat_shifts$Target
        y <- dat_shifts$NetShift
        
        #Fit Eqn
        A1_start <- 50
        A2_start <- 10000
        x0_start <- 3000
        p_start <- 0.5
        fit <- nls(formula = y ~ (((A1 - A2) / (1 + (x / x0) ^ p)) + A2),
                   start = list(A1 = A1_start, 
                                A2 = A2_start, 
                                x0 = x0_start, 
                                p = p_start),
                   algorithm = "port",
                   lower=c(-400, 0, 0, 0),
                   upper=c(100, 30000, 4000, 1))
        
        #Output
        data.frame(coef(summary(fit)))
        
}


#mydata <- read.csv("http://www.ats.ucla.edu/stat/data/binary.csv")
#head(mydata)
#mydata$rank <- factor(mydata$rank)
#mylogit <- glm(admit ~ gre + gpa + rank, data = mydata, family = "binomial")
#summary(mylogit)
