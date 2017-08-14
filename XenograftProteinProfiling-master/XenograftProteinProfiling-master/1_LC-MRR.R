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
        filename <- "groupNames_allClusters.csv"
        
        # get information of chip layout from github repository
        if (!file.exists("groupNames.csv")){
        url <- "https://raw.githubusercontent.com/JamesHWade/XenograftProteinProfiling/master/groupNames_allClusters.csv"
        filename <- basename(url)
        download.file(url, filename)
        }
        
        # define recipe as global variable for use in other functions
        recipe <<- read_csv(filename)
        colnames(recipe)[1] <- "Target" # rename column to remove byte order mark
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
                dat <- read_csv(ring, col_names = FALSE)
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

SubtractControl <- function(loc = 'plots', ch, cntl){
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
        dat[dat$Ring == i, 4] <- ringTC
        }
        
        write_csv(dat, paste(loc,"/", name, "_", cntl, "Control", "_ch", ch, 
                       ".csv", sep = ''))   
}

PlotRingData <- function(cntl, ch, loc = 'plots'){
        # loads relevant libraries and plot theme
        library(tidyverse)
        library(RColorBrewer)
        
        plot_theme <- theme_bw() + 
                      theme(text = element_text(size = 22),
                      axis.line = element_line(colour = "black"),
                      panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank())
        
        # get working directory to reset at end of function
        directory <- getwd()
        
        # use thermally controlled data if desired
        if (cntl != "raw"){
        dat <- read_csv(paste(loc, "/", name, "_", cntl, "Control", 
                          "_ch", ch,".csv", sep=''), col_types = cols())
        } else if (cntl == "raw") {
        dat <- read_csv(paste(loc, "/", name, "_allRings.csv", sep=''), 
                    col_types = cols())
        }
        
        dat.plot <- filter(dat, Ring != 4)
        
        #configure plot and legend
        plots <- ggplot(dat.plot, aes(x = Time, y = Shift, color = factor(Ring))) + 
                geom_line() +
                labs(color = "Rings", x = "Time (min)", 
                     y = expression(paste("Relative Shift (",Delta,"pm)"))) +
                plot_theme #+ facet_wrap(~Ring)
        
        
        run.1 <- plots + facet_wrap(~Ring)
        run.1
        
        ggsave(plots, filename = "AllRings_FullRun.png", width = 8, height = 6)
        ggsave(run.1, filename = "AllRings_FullRun_Facet.png", width = 8, height = 6)
        
        
        run.2 <- plots + xlim(20, 110)
        run.2
        
        run.3 <- plots + xlim(95, 125) + ylim(-80, -40)
        run.3
        
        ggsave(run.3, filename = "AllRings_Subset.png", width = 8, height = 6)
        
        
        if (cntl == "raw"){
                plots <- plots + geom_point(size = 1) + facet_grid(.~ Channel)
        }
        

        # alternative plots with averaged clusters

        dat.2 <- filter(dat.plot, Ring == c(10, 12))
        dat.2 <- dat.plot %>% group_by(Target, `Time Point`) %>% 
                summarise_each(funs(mean, sd), c(Time, Shift))

        plots.avg <- ggplot(dat.2, aes(x = Time_mean, y = Shift_mean, color = Target)) + 
                geom_line(size = 1) + 
                plot_theme + xlab("Time (min)") +
                ylab(expression(paste("Relative Shift (",Delta,"pm)"))) + 
                geom_ribbon(aes(ymin = Shift_mean - Shift_sd,
                                    ymax = Shift_mean + Shift_sd, linetype = NA),
                                    fill = "slategrey", alpha = 1/4) +
                theme(legend.position = "none")
        
        ggsave(plots.avg, filename = "Average_FullRun.png", width = 8, height = 6)
        
        
        
        
        
        # smoothing data with rolling (moving) average
        library(zoo)
        library(reshape2)
        dat.plot$Smooth <- rollmean(dat.plot$Shift, k = 51, fill = "extend")
        dat.3 <- dat.plot[,-c(1,2,5,6,7,8)]
        dat.melt <- melt(dat.3, id.vars = "Time")
        
        
        ggplot(dat.plot, aes(x = Time, y = Shift)) + geom_line() + plot_theme +
                xlim(80,110)
        ggplot(dat.melt, aes(x = Time, y = value, color = variable)) + geom_line(size = 1) +
                plot_theme + xlim(30, 60)
        
        
        # smooth with subset
        dat.subset <- filter(dat.plot, between(Time, 60, 90))
        dat.subset$Shift <- dat.subset$Shift - dat.subset$Shift[1]
        dat.subset$Time <- dat.subset$Time - dat.subset$Time[1]
        dat.subset$Smooth <- rollmean(dat.subset$Shift, k = 21, fill = "extend")
        dat.smoothsub <- dat.subset[,-c(1,2,5,6,7,8)]
        dat.submelt <- melt(dat.smoothsub, id.vars = "Time")
        
        dat.a <- filter(dat.smoothsub, between(Time, 20, 21))
        sd(dat.a$Shift)/sd(dat.a$Smooth)
        
        3*sd(dat.a$Smooth)
        3*sd(dat.a$Shift)
        
        subplot <- ggplot(dat.submelt, aes(x = Time, y = value, color = variable)) + 
                geom_line(size = 1) + plot_theme + xlab("Time (min)") +
                ylab(expression(paste("Relative Shift (",Delta,"pm)")))
        
        subplot
        
        subplot.f <- subplot + facet_grid(variable~.)
        subplot.f
        
        
        ggsave(subplot, filename = "SingleRing_Smoothed_Subset.png", width = 8, height = 6)
        ggsave(subplot.f, filename = "singleRing_Smoothed_Subset_Facet.png", width = 8, height = 6)
        
        
        
        dat.smooth <- filter(dat, Target != "thermal")
        plots.smooth <- ggplot(dat.smooth, aes(x = Time, y = Shift)) +
                geom_smooth(se = TRUE, method = ) + plot_theme
        
        plots.smooth + xlim(30,60)
        
        plots.avg + xlim(0, 30)
        plots.avg + xlim(30, 60)
        plots.avg + xlim(60, 90)
        plots.avg + xlim(95, 125) + ylim(-80, -40)
        plots.avg + xlim(120, 150)
        plots.avg + xlim(150, 180)
        plots.avg + xlim(180, 210)
        
        
        #save plot, uncomment to save
        filename <- paste0(name, "_", cntl, "Control", "_ch", ch)
        #filename2 <- paste0(name, "_", cntl, "Control", "_ch", ch, "_avg")
        setwd(loc)
        ggsave(plots, file = paste0(filename, ".png"), width = 8, height = 6)
        ggsave(plots, file = paste0(filename, ".pdf"), width = 8, height = 6)
        setwd(directory)
}

GetName()
AggData()
SubtractControl(ch = 2, cntl = "thermal")
PlotRingData(ch = 2, cntl = "thermal") 
