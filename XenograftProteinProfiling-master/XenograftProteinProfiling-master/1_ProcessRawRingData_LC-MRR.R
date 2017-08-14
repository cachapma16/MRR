GetName <- function(){
        # get the filename from the current working directory
        directory <- basename(getwd())
        
        # directory naming from MRR: "CHIPNAME_gaskGASkETTYPE_DATE"
        # extracts and returns GASKETTYPE from directory name
        name <- unlist(strsplit(directory, split = "_"))
        name <- name[2]
        
        # define name as global variable for use in other functions
        name <<- gsub('gask','',name) # removes "gask" from name
}

AggData <- function(loc = 'plots') {
        library(tidyverse)
        
        # get working directory to reset at end of function
        directory <- getwd()
        
        # change this file name to use alternative ring or group labels
        filename <- "groupNames_unmodified.csv"
        
        # get information of chip layout from github repository
        if (!file.exists("groupNames_unmodified.csv")){
                url <- "https://raw.githubusercontent.com/JamesHWade/XenograftProteinProfiling/master/groupNames_unmodified.csv"
                filename <- basename(url)
                download.file(url, filename)
        }
        
        # read in recipe/chip layout
        recipe <- read_csv(filename)
        colnames(recipe)[1] <- "Target" # rename column to remove byte order mark
        targets <- recipe$Target
        
        # generate list of rings to analyze (gets all *.csv files)
        rings <- list.files(directory, pattern = ".csv", recursive = FALSE)
        idfile <- grepl("group", rings)
        removeFiles <- c("comments.csv", rings[idfile])
        rings <- rings[!rings %in% removeFiles]
        
        # add data to data frame corresponding for each ring in rings
        df <- data.frame()
        for (i in rings) {
                dat <- read_csv(i, col_names = FALSE)
                ringNum <- as.numeric(strsplit(i, "\\.")[[1]][1])
                recipe.col <- which(recipe$Ring == ringNum)
                groupNum <- recipe$Group[recipe.col]
                groupName <- as.character(recipe$Target[[recipe.col]])
                channel <- recipe$Channel[[recipe.col]]
                tmp <- dat[,c(1,2)]
                names(tmp) <- c("Time", "Shift")
                tmp$ring <- ringNum
                tmp$group <- groupNum
                tmp$groupName <- groupName
                tmp$channel <- channel
                tmp$run <- name
                tmp$time_point <- seq(1:nrow(dat))
                df <- rbind(df, tmp)
        }
        
        # renames columns in df
        names(df) <- c("Time", "Shift", "Ring", "Group", "Target", "Channel",
                       "Experiment", "Time Point")
        
        # creates "plots" directory if one does not exist
        if (!file.exists(loc)){dir.create(loc)}
        
        # saves aggregated data with name_allRings.csv
        write_csv(df, paste(loc, '/', name, "_allRings.csv", sep=""))
        
        # returns working directory to top level
        setwd(directory)
}

PlotIndyRings <- function() {
        library(tidyverse)
        library(RColorBrewer)
        library(ggthemes)
        
        # get working directory to reset at end of function
        directory <- getwd()
        
        # use thermally controlled data if desired
        dat <- read_csv(paste0("plots/", name, "_allRings.csv"))
        dat <- filter(dat, Target != "thermal")
        dat$Ring <- factor(dat$Ring)
        
        #separates into subsets corrdinating to each PS injection
        startTime <- 33
        df0 <- subset(dat, Time > startTime & Time < startTime + 30)
        df1 <- subset(dat, Time > startTime + 30 & Time < startTime + 61)
        df2 <- subset(dat, Time > startTime + 61 & Time < startTime + 92)
        df3 <- subset(dat, Time > startTime + 92 & Time < startTime + 123)
        # df4 <- subset(dat, Time > startTime + 123 & Time < startTime + 154)
        # df5 <- subset(dat, Time > startTime + 154 & Time < startTime + 185)
        # df6 <- subset(dat, Time > startTime + 185 & Time < startTime + 216)
        # df7 <- subset(dat, Time > startTime + 216 & Time < startTime + 247)
        
        #add column with standard name 
        df0$Standard <- "Ethyl Acetate"
        df1$Standard <- "1.3 kDa"
        df2$Standard <- "3.5 kDa"
        df3$Standard <- "8.7 kDa"
        # df4$Standard <- "17.6 kDa"
        # df5$Standard <- "35 kDa"
        # df6$Standard <- "130 kDa"
        # df7$Standard <- "304 kDa"
        
        df <- rbind(df0, df1, df2, df3) #, df4, df5, df6, df7)
        
        #normalize time and shift for each injection
        standardList <- unique(df$Standard)
        ringList <- unique(df$Ring)
        
        dat.all <- data.frame()
        for (i in standardList){
                dat.std <- filter(df, Standard == i)
                for (i in ringList){
                        dat.ring <- filter(dat.std, Ring == i)
                        dat.ring$Time <- dat.ring$Time - dat.ring$Time[1]
                        dat.ring$Shift <- dat.ring$Shift - dat.ring$Shift[1]
                        tmp <- dat.ring
                        dat.all <- rbind(dat.all, tmp)
                }
        }
        
        dat.all$Standard <- factor(dat.all$Standard,
                                   levels = c("Ethyl Acetate","1.3 kDa", 
                                              "3.5 kDa", "8.7 kDa", "17.6 kDa", 
                                              "35 kDa", "130 kDa", "304 kDa"))
        
        plot <- ggplot(dat.all, aes(y = Shift, x = Time, color = Ring)) +
                geom_line() + facet_grid(Ring~Standard) + 
                xlab("Time (min)") +
                ylab(expression(paste("Relative Shift (",Delta,"pm)"))) +
                theme_few()
        
        ggsave(plot, filename = paste0("plots/", name, "_IndyRings.png"), 
               width = 8, height = 6)
}

DataSplitting <- function(){
        library(tidyverse)
        library(reshape2)
        library(pracma)
        library(baseline)
        
        dat <- read_csv(paste0("plots/", name, "_allRings.csv"))
        dat <- filter(dat, Target != "thermal")
        
        mean <- function(x) {mean(x, trim = 0.5)}
        
        dat.avg <- dat %>% group_by(Target, `Time Point`) %>%
                summarise_each(funs(mean, sd), c(Time, Shift))
        
        #separates into subsets corrdinating to each PS injection
        startTime <- 3
        df0 <- subset(dat.avg, Time_mean > startTime & 
                              Time_mean < startTime + 30)
        df1 <- subset(dat.avg, Time_mean > startTime + 30 & 
                              Time_mean < startTime + 61)
        df2 <- subset(dat.avg, Time_mean > startTime + 61 & 
                              Time_mean < startTime + 92)
        df3 <- subset(dat.avg, Time_mean > startTime + 92 & 
                              Time_mean < startTime + 123)
        # df4 <- subset(dat.avg, Time_mean > startTime + 123 & 
        #                       Time_mean < startTime + 154)
        # df5 <- subset(dat.avg, Time_mean > startTime + 154 & 
        #                       Time_mean < startTime + 185)
        # df6 <- subset(dat.avg, Time_mean > startTime + 185 & 
        #                       Time_mean < startTime + 216)
        # df7 <- subset(dat.avg, Time_mean > startTime + 216 & 
        #                       Time_mean < startTime + 247)
        
        #add column with standard name 
        df0$Standard <- "Ethyl Acetate"
        df1$Standard <- "1.3 kDa"
        df2$Standard <- "3.5 kDa"
        df3$Standard <- "8.7 kDa"
        # df4$Standard <- "17.6 kDa"
        # df5$Standard <- "35 kDa"
        # df6$Standard <- "130 kDa"
        # df7$Standard <- "304 kDa"
        
        runsList <- list(df0, df1, df2, df3) #, df4, df5, df6, df7)
        
        #smoothing data with whitaker smoothing
        runsListSmooth <- lapply(runsList, function(i){
                i$Shift_mean <- i$Shift_mean - i$Shift_mean[1]
                i$Whit <- whittaker(i$Shift_mean, lambda = 10000, d = 3)
                peakDetectionShift <- baseline.irls(t(i$Whit))
                i$Baseline <- t(peakDetectionShift$corrected)
                i$Time_mean <- i$Time_mean - i$Time_mean[1]
                i$ElutionVolume <- i$Time_mean*0.3
                return(i)
        })
        
        
        df <- bind_rows(runsListSmooth)
        #df.melt <- melt(data = df, id.vars = "Time_mean", measure.vars = c("Baseline", "Whit"))
        #ggplot(df.melt) + geom_line(aes(Time_mean, value, color = variable))
        
        df$Standard <- factor(df$Standard,
                              levels = c("Ethyl Acetate","1.3 kDa", 
                                         "3.5 kDa", "8.7 kDa", "17.6 kDa", 
                                         "35 kDa", "130 kDa", "304 kDa"))
        
        write_csv(x = df, 
                  path = paste0("plots/", name, "inj_combined.csv"))
}

PlotAvgData <- function(loc = 'plots'){
        # loads relevant libraries and plot theme
        library(tidyverse)
        library(RColorBrewer)
        library(ggthemes)
        
        # get working directory to reset at end of function
        directory <- getwd()
        
        # use thermally controlled data if desired
        dat <- read_csv(paste0(loc, "/", name, "inj_combined.csv"))
        dat <- filter(dat, Target != "thermal")
        setwd(loc)
        
        dat$Standard <- factor(dat$Standard,
                               levels = c("Ethyl Acetate", "1.3 kDa", "3.5 kDa",
                                          "8.7 kDa", "17.6 kDa", "35 kDa", 
                                          "130 kDa", "304 kDa") )
        
        # alternative plots with averaged clusters
        plot_theme <- theme_few(base_size = 16)
        
        
        plots.avg <- ggplot(dat, aes(x = ElutionVolume, y = Shift_mean, 
                                     color = Standard)) + 
                geom_line() + 
                xlab("Elution Volume (mL)") +
                ylab(expression(paste("Relative Shift (",Delta,"pm)"))) + 
                geom_ribbon(aes(ymin =Shift_mean - Shift_sd,
                                ymax =Shift_mean + Shift_sd, linetype = NA),
                            fill = "slategrey", alpha = 0) +
                plot_theme
        
        plots.avg.wrap <- plots.avg + facet_wrap(~Standard)
        
        ggsave(plots.avg, filename = paste0(name, "avg_overlapped.png"), 
               width = 8, height = 6)
        ggsave(plots.avg.wrap, filename = 
                       paste0(name, "avg_overlapped_wrapped.png"), 
               width = 8, height = 6)
        
        # smoothing data with rolling (moving) average
        smooth <- ggplot(dat, 
                         aes(y = Baseline, x = ElutionVolume, color = Standard)) +
                geom_line(size = 0.5) + xlab("Elution Volume (mL)") +
                ylab(expression(paste("Relative Shift (",Delta,"pm)"))) + 
                plot_theme
        
        smooth.comp <- smooth + 
                geom_line(aes(y =Shift_mean, x =ElutionVolume), alpha = 1/2)
        
        smooth.wrap <- smooth + facet_wrap(~Standard)
        
        ggsave(smooth, 
               filename = paste0(name, " Smooth.png"), 
               width = 8, height = 6)
        ggsave(smooth.comp, 
               filename = paste0(name, " Smoothed and Averaged.png"), 
               width = 8, height = 6)
        ggsave(smooth.wrap, 
               filename = paste0(name, " Smooth Wrap.png"), 
               width = 8, height = 6)
        
        setwd(directory)
}

AnalyzeData <- function(){
        GetName()
        AggData()
        DataSplitting()
        PlotIndyRings()
        PlotAvgData()
}

BatchAnalyze <- function(){
        dataList <- list.dirs(recursive = FALSE)
        directory <- getwd()
        for (i in dataList){
                setwd(i)
                AnalyzeData()
                setwd(directory)
        }
}