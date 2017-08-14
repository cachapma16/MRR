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
        # load relevant libraries
        library(tidyverse)
        
        # get working directory to reset at end of function
        directory <- getwd()
        
        # change this file name to use alternative ring or group labels
        filename <- "groupNames_unmodified.csv"
        
        # get information of chip layout from github repository
        if (!file.exists(filename)){
                url <- "https://raw.githubusercontent.com/JamesHWade/XenograftProteinProfiling/master/groupNames_unmodified.csv"
                filename <- basename(url)
                download.file(url, filename)
        }
        
        # define recipe as global variable for use in other functions
        recipe <<- read_csv(filename)
        colnames(recipe)[1] <- "Target" # rename column & remove byte order mark
        targets <- recipe$Target
        
        # generate list of rings to analyze (gets all *.csv files)
        rings <- list.files(directory, pattern = ".csv", recursive = FALSE)
        idfile <- grepl("group", rings)
        removeFiles <- c("comments.csv", rings[idfile])
        rings <- rings[!rings %in% removeFiles]
        
        # create empty data frame to store data
        df <- data.frame()
        
        # add data to data frame corresponding for each ring in rings
        for (i in rings) {
                ring <- as.vector(i)
                dat <- read_csv(ring, col_types = cols(), col_names = FALSE)
                time_shift <- dat[ ,1]
                shift <- dat[ ,2]
                ringStr <- strsplit(i, "\\.")[[1]]
                ringNum <- as.numeric(ringStr[1])
                recipe.col <- which(recipe$Ring == ringNum)
                groupNum <- recipe$Group[recipe.col]
                ring <- rep(ringNum, nrow(dat))
                group <- rep(groupNum, nrow(dat))
                groupName <- as.character(recipe$Target[[recipe.col]])
                groupName <- rep(groupName, nrow(dat))
                channel <- recipe$Channel[[recipe.col]]
                channel <- rep(channel, nrow(dat))
                run <- rep(name, nrow(dat))
                time_point <- seq(1:nrow(dat))
                tmp <- data.frame(ring, group, time_shift, shift, groupName, 
                                  channel, run, time_point)
                df <- rbind(df, tmp)
        }
        
        # renames columns in df
        names(df) <- c("Ring", "Group", "Time", "Shift", "Target", "Channel",
                       "Experiment", "Time Point")
        
        # creates "plots" directory if one does not exist
        if (!file.exists(loc)){dir.create(loc)}
        
        # saves aggregated data with name_allRings.csv
        write_csv(df, paste(loc, '/', name, "_allRings.csv", sep=""))
        
        # returns working directory to top level
        setwd(directory)
}

SubtractControl <- function(loc = 'plots', 
                            ch, 
                            cntl){
        #load relevant libraries
        library(tidyverse)
        
        # get working directory to reset at end of function
        directory = getwd()
        
        # get ring data and filter by channel
        dat <- read_csv(paste0(loc, "/", name, "_", "allRings.csv"))
        if (ch != "U"){
                dat <- filter(dat, Channel == ch)
        }
        dat <- filter(dat, Target != "Ignore")
        
        # get thermal control averages
        controls <- filter(dat, Target == cntl)
        ringList <- unique(controls$Ring)
        
        # gets times from first thermal control
        times <- filter(controls, Ring == ringList[1]) %>% select(Time)
        df.controls <- data.frame(times)
        
        # create dataframe with all controls
        for (i in ringList){
                ringShift <- filter(controls, Ring == i) %>% select(Shift)
                names(ringShift) <- paste('Ring', i, sep='')
                df.controls <- cbind(df.controls, ringShift)
        }
        
        # averages thermal controls
        cols <- ncol(df.controls)
        if (length(unique(controls$Ring)) != 1) {
                df.controls$avgControls <- rowMeans(df.controls[,c(2:cols)])
        } else {
                df.controls$avgControls <- df.controls[,c(2:cols)]
        }
        avgControls <- as.vector(df.controls$avgControls)
        
        #subtracts thermal controls from each ring
        ringNames <- unique(dat$Ring)
        for(i in ringNames){
                ringDat <- filter(dat, Ring == i) %>% select(Shift)
                ringTC <- ringDat - avgControls
                dat[dat$Ring == i, 3] <- ringTC
        }
        
        dat <- filter(dat, Ring != 3)
        
        write_csv(dat, paste(loc,"/", name, "_", cntl, "Control", "_ch", ch, 
                             ".csv", sep = ''))   
}

PlotIndyRings <- function(loc = 'plots', 
                          delay = 30) {
        library(tidyverse)
        library(RColorBrewer)
        library(ggthemes)
        
        # get working directory to reset at end of function
        directory <- getwd()
        
        # use thermally controlled data if desired
        dat <- read_csv(paste0(loc, "/", name, "_allRings.csv"))
        dat <- filter(dat, Target != "thermal")
        dat$Ring <- factor(dat$Ring)
        
        #separates into subsets corrdinating to each PS injection
        startTime <- delay
        df0 <- subset(dat, Time > startTime & Time < startTime + 30)
        df1 <- subset(dat, Time > startTime + 30 & Time < startTime + 61)
        df2 <- subset(dat, Time > startTime + 61 & Time < startTime + 92)
        df3 <- subset(dat, Time > startTime + 92 & Time < startTime + 123)
        df4 <- subset(dat, Time > startTime + 123 & Time < startTime + 154)
        df5 <- subset(dat, Time > startTime + 154 & Time < startTime + 185)
        df6 <- subset(dat, Time > startTime + 185 & Time < startTime + 216)
        df7 <- subset(dat, Time > startTime + 216 & Time < startTime + 247)
        
        #add column with standard name 
        df0$Standard <- "Ethyl Acetate"
        df1$Standard <- "1.3 kDa"
        df2$Standard <- "3.5 kDa"
        df3$Standard <- "8.7 kDa"
        df4$Standard <- "17.6 kDa"
        df5$Standard <- "35 kDa"
        df6$Standard <- "130 kDa"
        df7$Standard <- "304 kDa"
        
        df <- rbind(df0, df1, df2, df3, df4, df5, df6, df7)
        
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
        
        ggsave(plot, filename = paste0(loc, "/", name, "_IndyRings.png"), 
               width = 8, height = 6)
        
        write_csv(df, paste0(loc, '/', name, "_allRings_byRing.csv"))
        
}

CheckRingQuality <- function(loc = 'plots', 
                             varLevel = 10) {
        library(tidyverse)
        
        dat <- read_csv(paste0(loc,"/", name, "_allRings_byRing.csv"))
        
        dat.avg <- dat %>% group_by(Ring) %>%
                summarise_each(funs(var), c(Shift))
        
        ringWinners <- filter(dat.avg, Shift < varLevel) %>% select(Ring)
        write_csv(ringWinners, paste0(loc, '/', name, "_ringWinners.csv"))
}

DataSplitting <- function(winner, 
                          loc = 'plots', 
                          delay = 30){
        library(tidyverse)
        library(zoo)
        library(pracma)
        
        dat <- read_csv(paste0(loc, '/', name, "_allRings.csv"))
        
        if (winner == TRUE) {
                winners <- read_csv(paste0(loc, '/', name, "_ringWinners.csv"))
                rings <- as.vector(unlist(winners))
                dat <- filter(dat, grepl(paste(rings, collapse = "|"), Ring))
        }
        
        kSmooth <- 7
        flSavgol <- 31
        fSavgol <- 4
        dSavgol <- 0
        mean <- function(x) {mean(x, trim = 0.5)}
        
        dat.avg <- dat %>% group_by(Target, `Time Point`) %>%
                summarise_each(funs(mean, sd), c(Time, Shift))
        dat.avg$Smooth <- rollmean(dat.avg$Shift_mean, 
                                   k = kSmooth, 
                                   fill = "extend")
        dat.avg$Savgol <- savgol(dat.avg$Shift_mean, 
                                 fl = flSavgol, 
                                 forder = fSavgol, 
                                 dorder = dSavgol)
        
        #separates into subsets corrdinating to each PS injection
        startTime <- delay
        df0 <- subset(dat.avg, Time_mean > startTime & 
                              Time_mean < startTime + 30)
        df1 <- subset(dat.avg, Time_mean > startTime + 30 & 
                              Time_mean < startTime + 61)
        df2 <- subset(dat.avg, Time_mean > startTime + 61 & 
                              Time_mean < startTime + 92)
        df3 <- subset(dat.avg, Time_mean > startTime + 92 & 
                              Time_mean < startTime + 123)
        df4 <- subset(dat.avg, Time_mean > startTime + 123 & 
                              Time_mean < startTime + 154)
        df5 <- subset(dat.avg, Time_mean > startTime + 154 & 
                              Time_mean < startTime + 185)
        df6 <- subset(dat.avg, Time_mean > startTime + 185 & 
                              Time_mean < startTime + 216)
        df7 <- subset(dat.avg, Time_mean > startTime + 216 & 
                              Time_mean < startTime + 247)
        
        #add column with standard name 
        df0$Standard <- "Ethyl Acetate"
        df1$Standard <- "1.3 kDa"
        df2$Standard <- "3.5 kDa"
        df3$Standard <- "8.7 kDa"
        df4$Standard <- "17.6 kDa"
        df5$Standard <- "35 kDa"
        df6$Standard <- "130 kDa"
        df7$Standard <- "304 kDa"
        
        #smoothing data with rolling (moving) average
        df0$Time_mean <- df0$Time_mean - df0$Time_mean[1]
        df0$Shift_mean <- df0$Shift_mean - df0$Shift_mean[1]
        df0$Smooth <- df0$Smooth - df0$Smooth[1]
        df0$Savgol <- df0$Savgol - df0$Savgol[1]
        
        df1$Time_mean <- df1$Time_mean - df1$Time_mean[1]
        df1$Shift_mean <- df1$Shift_mean - df1$Shift_mean[1]
        df1$Smooth <- df1$Smooth - df1$Smooth[1]
        df1$Savgol <- df1$Savgol - df1$Savgol[1]
        
        df2$Time_mean <- df2$Time_mean - df2$Time_mean[1]
        df2$Shift_mean <- df2$Shift_mean - df2$Shift_mean[1]
        df2$Smooth <- df2$Smooth - df2$Smooth[1]
        df2$Savgol <- df2$Savgol - df2$Savgol[1]
        
        df3$Time_mean <- df3$Time_mean - df3$Time_mean[1]
        df3$Shift_mean <- df3$Shift_mean - df3$Shift_mean[1]
        df3$Smooth <- df3$Smooth - df3$Smooth[1]
        df3$Savgol <- df3$Savgol - df3$Savgol[1]
        
        df4$Time_mean <- df4$Time_mean - df4$Time_mean[1]
        df4$Shift_mean <- df4$Shift_mean - df4$Shift_mean[1]
        df4$Smooth <- df4$Smooth - df4$Smooth[1]
        df4$Savgol <- df4$Savgol - df4$Savgol[1]
        
        df5$Time_mean <- df5$Time_mean - df5$Time_mean[1]
        df5$Shift_mean <- df5$Shift_mean - df5$Shift_mean[1]
        df5$Smooth <- df5$Smooth - df5$Smooth[1]
        df5$Savgol <- df5$Savgol - df5$Savgol[1]
        
        df6$Time_mean <- df6$Time_mean - df6$Time_mean[1]
        df6$Shift_mean <- df6$Shift_mean - df6$Shift_mean[1]
        df6$Smooth <- df6$Smooth - df6$Smooth[1]
        df6$Savgol <- df6$Savgol - df6$Savgol[1]
        
        df7$Time_mean <- df7$Time_mean - df7$Time_mean[1]
        df7$Shift_mean <- df7$Shift_mean - df7$Shift_mean[1]
        df7$Smooth <- df7$Smooth - df7$Smooth[1]
        df7$Savgol <- df7$Savgol - df7$Savgol[1]
        
        df <- rbind(df0, df1, df2, df3, df4, df5, df6, df7)
        
        
        df$Standard <- factor(df$Standard,
                              levels = c("Ethyl Acetate","1.3 kDa", 
                                         "3.5 kDa", "8.7 kDa", "17.6 kDa", 
                                         "35 kDa", "130 kDa", "304 kDa"))
        
        write_csv(x = df, 
                  path = paste0(loc, "/", name, "inj_combined.csv"))
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
        
        
        plots.avg <- ggplot(dat, aes(x = Time_mean, y = Shift_mean, 
                                     color = Standard)) + 
                geom_line() + 
                xlab("Time (min)") +
                ylab(expression(paste("Relative Shift (",Delta,"pm)"))) + 
                geom_ribbon(aes(ymin =Shift_mean - Shift_sd,
                                ymax =Shift_mean + Shift_sd, linetype = NA),
                            fill = "slategrey", alpha = 1/4) +
                plot_theme
        
        plots.avg.wrap <- plots.avg + facet_wrap(~Standard)
        
        plots.avg.wrap
        
        ggsave(plots.avg, filename = paste0(name, "avg_overlapped.png"), 
               width = 8, height = 6)
        ggsave(plots.avg.wrap, filename = 
                       paste0(name, "avg_overlapped_wrapped.png"), 
               width = 8, height = 6)
        
        # smoothing data with rolling (moving) average
        smooth <- ggplot(dat, 
                         aes(y = Smooth, x =Time_mean, color = Standard)) +
                geom_line(size = 0.5) + plot_theme
        
        smooth.comp <- smooth + 
                geom_line(aes(y =Shift_mean, x =Time_mean), alpha = 1/2)
        
        smooth.wrap <- smooth + facet_wrap(~Standard)
        
        savgol <- ggplot(dat, 
                         aes(y = Savgol, x =Time_mean, color = Standard)) +
                geom_line(size = 0.5) + plot_theme
        
        savgol + facet_wrap(~Standard)
        
        savgol.comp <- savgol + 
                geom_line(aes(y =Shift_mean, x =Time_mean), alpha = 1/2)
        
        savgol.none <- ggplot(dat) +
                geom_line(aes(y =Shift_mean, x =Time_mean,
                              color = Standard), alpha = 1/2) +
                plot_theme
        
        savgol.wrap <- savgol + facet_wrap(~Standard)
        
        comp <- ggplot(dat) +
                geom_line(aes(y = Savgol, x =Time_mean, group = Standard), 
                          color = "black") +
                geom_line(aes(y = Smooth, x =Time_mean, group = Standard),
                          color = "red") +
                plot_theme
        
        ggsave(smooth, 
               filename = paste0(name, " RM-Smooth.png"), 
               width = 8, height = 6)
        ggsave(smooth.comp, 
               filename = paste0(name, " RM-Smoothing and Averaged.png"), 
               width = 8, height = 6)
        ggsave(smooth.wrap, 
               filename = paste0(name, " RM-Smooth Wrap.png"), 
               width = 8, height = 6)
        ggsave(savgol, 
               filename = paste0(name, " SG-Smooth.png"), 
               width = 8, height = 6)
        ggsave(savgol.comp, 
               filename = paste0(name, " SG-Smoothing and Averaged.png"), 
               width = 8, height = 6)
        ggsave(savgol.none, 
               filename = paste0(name, " No Smooth.png"), 
               width = 8, height = 6)
        ggsave(savgol.wrap, 
               filename = paste0(name, " SG-Smooth Wrap.png"), 
               width = 8, height = 6)
        ggsave(comp, 
               filename = paste0(name, " Smoothing Comparison.png"), 
               width = 8, height = 6)
        
        setwd(directory)
}

AnalyzeData <- function(cntl = FALSE, 
                        loc = 'plots', 
                        winner = FALSE, 
                        delay = 30){
        #setwd(choose.dir())
        GetName()
        AggData(loc)
        if (cntl == TRUE){
                SubtractControl(loc, ch = 'U', cntl = "thermal")
        }
        #PlotIndyRings(loc, delay)
        CheckRingQuality(loc, varLevel = 10)
        #DataSplitting(winner, loc, delay)
        #PlotAvgData(loc)
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
